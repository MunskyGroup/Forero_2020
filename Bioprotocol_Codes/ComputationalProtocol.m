%% Data Normalization Pipeline and Model Fitting for Forero2020

% This code goes with the bio-protocol paper:

% Visualization, quantification, and prediction of endogenous RNA 
% polymerase II phosphorylation at a single-copy gene in living cells
% -------------------------------------------------------------------------
% Linda S. Forero-Quintero, William Raymond, 
% Brian Munsky & Timothy J. Stasevich 
% -------------------------------------------------------------------------

% adapted from:

% Live-cell imaging reveals the spatiotemporal organization of 
% endogenous RNA polymerase II phosphorylation at a single gene
% -------------------------------------------------------------------------
% Linda S. Forero-Quintero, William Raymond, Tetsuya Handa, 
% Matthew N. Saxton, Tatsuya Morisaki,
% Hiroshi Kimura, Edouard Bertrand, Brian Munsky & Timothy J. Stasevich 
% Nature Communications volume 12, Article number: 3158 (2021) 
% https://doi.org/10.1038/s41467-021-23417-0
% -------------------------------------------------------------------------

clc
clear 
close all

%% Section 1 - Load an Normalize Trajectory Data, and Calculate Normalized Correlation Functions
%  Bio-Protocol Section E

% Returns 2 structures: IntensityData, NormalizedCorrelationData

% IntensityData.Name.normalized_data - normalized intensity trajectories
% IntensityData.Name.raw_data - raw intensity trajectories

% NormalizedCorrelationData.Name_Name2.normalized_corrs - all correlations across the two channels
% NormalizedCorrelationData.Name_Name2.mean_corr  - mean correlation across the two channels
% NormalizedCorrelationData.Name_Name2.sem_corr - The standard error of mean
% NormalizedCorrelationData.Name_Name2.G0_mean - The G0 used for normalizations

% User Defined Parameters
FileName = 'Raw_Intensities Analysis-BLC_W_BG_WO_RunAve_May18';  % Intensity Trajectories (see Section XXX in Protocol).
DataCells = {'C2:C4001','D2:D4001','E2:E4001'};  % Cells or area containing intensity trajectories to normalize 
NChannels = length({'C2:C4001','D2:D4001','E2:E4001'});
ChannelNames = {'CTD','Ser5ph','mRNA'};
NCells = 20;                   % number of cells in the data
LTrajectory = 200;             % number of time points in each trajectory
NormalizationQuantile = 0.95;  % the percentile to normalize as 1 (0-1)
MaxValue = 1.5;                % the maximum value after normalization 
RezeroAutoCorrelation = true; % option to rezero autocorrelation tails
G0_type = 'G0_intp';           % G0 option type: G0 (raw g0), none (divide by 1), G0_intp (Interpolate G0 from G1,G2,G3), G1
CC_Norm = 'ACC_SQRT';           % Normalization type for CC, ACC_SQRT (use the sqrt(g0_signal1)*sqrt(g0_signal2))
CC_G0 = 'G0';                 % G0 option for cross correlations: G0 (use the center delay), max (use the max value), none (use none)
CC_delays = 10;                % How many time points to calculate delay out too for cross correlations
ACC_delays = 31;               % How many time points to calculate delay out too for auto correlations
RezeroPoints = 9;              % how many final points to rezero by for autocorrelations


% FISHQUANT / Quantification of mRNA amount
MeanSpotCount = 15.5;
VarianceSpotCount = .93;

[IntensityData, NormalizedCorrelationData] = get_normalized_ccs_from_xls(FileName, DataCells, NChannels, ...
                                                        ChannelNames, NCells, LTrajectory, NormalizationQuantile, ...
                                                        MaxValue, RezeroAutoCorrelation, RezeroPoints, G0_type, CC_G0, CC_Norm, CC_delays, ACC_delays);

%% Section 2 - Model Specification
% Specify the model here, see Forero 2020 methods section.

% selected model from forero2020, see figure 3c for model diagram
Nstates = 3;
parameters = [  0.433596648003337, ... % kon 1/min
                1, ...              % koff 1/min
                0.666624850282066, ... % kesc 1/min
                0.198727194675144, ... % kproc 1/min
                15.4054402060124, ... % kin 1/min
                0.777851806420281, ... % kout 1/min
                1, ...                 % frac  percentage
                1.98448561566709, ...  % eta_rnap
                1.42073542241570, ...  % eta_ser5
                0.408137598525255, ... % eta_ts
                0.924805163389989, ... % scale_rnap (scale of shot noise)
                1.22915843769510, ...  % scale_ser5 (scale of shot noise)
                1.14124029524139  ];   % scale_mrna (scale of shot noise)

kon = parameters(1);
koff = parameters(2);
kesc = parameters(3);
kproc = parameters(4);
kin =  parameters(5);
kout = parameters(6);
frac = parameters(7);

eta_rnap = parameters(8);
eta_ser5 = parameters(9);
eta_ts = parameters(10);

sc1 = parameters(11);
sc2 = parameters(12);
sc3 = parameters(13);

% S*W1(x) + W0 = dx/dt (system of linear ODEs)
% stoichiometry matrix
S =  [[1    -1     0     0     0     0];
      [0     0     1    -1    -1     0];
      [0     0     0     0     1    -1]];   

% W1 (propensities)

W1 = @(p)[[-p(1)           0         0];
          [p(2)*1000           0         0];
          [p(5)*1000           0         0];
          [    0         p(6)       0];
          [    0       p(7)*p(3)    0];
          [    0          0       p(4)];];   

% independent of x
W0 = @(p)[[p(1); 0; 0; 0; 0; 0;]];
b =  @(p)[p(1); 0; 0;];

% signal transformation matrix c, this matrix converts states x to intensity channels
% e.g. state X1 adds 1 to both R channel and G channel since the RNAP2 has ser5ph and ctd

c = @(p)[[0,   1,  1]; 
         [0, p(7), 1];
         [0,   0,  1]];

noise_parameters = @(p)[p(8), p(9), p(10), p(11), p(12), p(13)];


%% Section 3 - Function to solve model [a,b,c] = solve_model(S,w,pars)
% Function that solves the defined model for Correlations, Means, and Variances


time_vec = [0:30]; % Time span to solve over (minutes)
[ModelCorrelations, ModelMeans, ModelVariances] = solve_model_from_mats(S, W1, W0, c, b, ...
                                                                        noise_parameters, time_vec, parameters);


%% Section 4 - Function to compare model to data

% User Defined Parameters
cc_range = 15; % cross correlation delays to include in the log likelihood calculation
acc_range = 20; % auto correlation delays to include in the log likelihood calculation
channels_quantified = [3]; % which channel was quantified for spot counts? in our case 3 for mRNA
DataMeans = [0,0,15.5]; % Data means for spot counts, 0 if not quantified
DataSEMs = [0,0,.93];   % SEM of the spots quantified, 0 if not quantified
Nspots = [0,0,130];     % Number of spots quantified, 0 if not quantified


par_fixed = parameters; % Which parameters are fixed
par_changed = [1,3:6,8:13]; %which parameters are free to change

[ModelLikelihood] = get_log_likelihood(parameters(par_changed), NormalizedCorrelationData, DataMeans, ...
                                       DataSEMs, S, W1, W0, c, b, ...
                                       noise_parameters, time_vec, cc_range, acc_range, ...
                                       channels_quantified, Nspots, ChannelNames, par_fixed, par_changed);



%% Section 5 - Maximize Likelihood

% Returns 1xN best parameters and MLE value

% User Defined Parameters
search_chains = 20; %Number of iterations to run the GA/fminsearch combo
GA_pop = 200; %number for GA algorithm populations
save_file_name = 'tmp_best_par'; %where to save the best parameter set file

Constraints.LB=-5*ones(size(par_changed)); %constraints of parameters Upper and Lower
Constraints.UB=5*ones(size(par_changed));

% Get loglikelihood function only varying by parameters
get_LL = @(parameters) sum(get_log_likelihood(10.^parameters, NormalizedCorrelationData, DataMeans, ...
                                       DataSEMs, S, W1, W0, c, b, ...
                                       noise_parameters, time_vec, cc_range, acc_range, ...
                                       channels_quantified, Nspots, ChannelNames, par_fixed, par_changed));


[ParBest,MLE] = get_MLE(get_LL, parameters, par_changed, search_chains, Constraints, GA_pop, save_file_name);

%% Section 6 - Compute BIC/AIC

% Returns BIC and AIC.

% User Defined Parameters
k = 5; % Number of parameters
n = 8; % Observed datapoints or independent data sets (6 correlations + 2 moments of spots)

BIC = k*log(n) - 2*-MLE; %Bayes information criterion
AIC = k - 2*-MLE; % Akaike information criterion

%% Section 7 - Run Metropolis Hastings for Parameter Uncertainty Quantification

% User Defined Parameters
delta = 0.03*ParBest(par_changed); % percentage of change of each parameter
proprnd = @(x)(x+delta.*randn(size(x)).*(randn(size(x))>0.5)); % Proposal function
parnames = {'kon','na','kesc','kproc','beta','kout','na','eta_ctd','eta_ser5','eta_mrna',...
    'sc_ctd','sc_ser5','sc_mrna'}; %Parameter names

nsamples = 5000; % n samples per chain
nchains = 1; % n chains total
nsegments = 20; % nsegments per chain
thin = 20; % mh thinning rate
save_file_name = 'met_hast_pars_example_'; %names of the mh chain files
pars_to_plot = [1,3,4,5,6];
fignum = 5;

%if parfor not available replace with for
par_opt = log10(ParBest(par_changed));
parfor i=1:nchains    
    sv_file = [save_file_name,num2str(i)];
    run_chain(i, par_opt, sv_file, get_LL, proprnd, thin, nsamples, nsegments) 
end

plot_mh_from_files(parnames, save_file_name, nchains, nsegments, pars_to_plot, fignum)


%% Section 8 - Make Plots of Correlations
% Compare the model correlations with the new parameters to the data 
% correlations.

%recalculate correlations with newly fit parameters.
time_vec = [0:30]; 
[ModelCorrelations, ModelMeans, ModelVariances] = solve_model_from_mats(S, W1, W0, c, b, ...
                                                                        noise_parameters, time_vec, ParBest);

plot_correlations(ModelCorrelations, NormalizedCorrelationData, ChannelNames)




%% Section 9 - Sample Model Intensities 
% Function that generates sample intensity trajectories from the defined model

% Returns: IntensityData structure
% IntensityData.Name.normalized_data 
% IntensityData.Name.raw_data 

% User Defined parameters
number_trajectories = 20; %number of samples
time_vec = [0:200]; % time vector to sample over
zero_values = [30/200, 23/200,5/200]; % approximation of zero from real data for each signal type (how many n points out of 200 were under 0?)

[MoleculeTrajectories] = sample_model_from_mats(S, W1, W0, c, time_vec, parameters, ...
                      number_trajectories);

[IntensityTrajectories] = add_shot_noise(MoleculeTrajectories, parameters, noise_parameters, zero_values);
[IntensityTrajectories] = normalize_simulated_trajectories(IntensityTrajectories, NChannels, ...
                                                        ChannelNames, NormalizationQuantile, ...
                                                        MaxValue);

figure()
plot(NormalizedIntensityTrajectories.CTD.normalized_data(:,1),'r'); hold on;
plot(NormalizedIntensityTrajectories.Ser5ph.normalized_data(:,1),'g');
plot(NormalizedIntensityTrajectories.mRNA.normalized_data(:,1),'b');
legend('CTD','Ser5ph','mRNA')
title('Example Simulated Intensity Trace')

%% Section 10 - Compare Model Intensities with Data Intensities
% Generate a plot comparing normalized model and data intensities with histograms

% User Defined parameters
nbins = 30; %n number of bins for the histograms

plot_histograms(IntensityTrajectories, IntensityData, ChannelNames, nbins)

%% Section 11 - Predict RNAP2 Locations (e.g, ChIP)
% Simulate a ChIP experiment for the defined model for Nsites + 2 
% (front without escaping / back TS processing). Internal sites are assumed
% to have a uniform distribution of RNAP. Residence time is defined as
% 1/processing rate for the time RNAP are actively elongating or processing.

% Returns: 1xN mean and 1xN std for RNAP at each ChIP site.

% User Defined parameters
time_vec = [0:600];
number_samples = 100;  % number of times to run the simulated ChIP
number_trajectories = 100; % number of TS sites per simulated ChIP
n_sites = 6;           % Number of Internal ChIP sites 
gene_size = 5.2;       % Kilobases (kb)
elongation_rate = 4.1; % kb/min
total_residence_time = 1/kproc; % total time rnaps are residing on the gene with TS


[SimulatedChIP, SimulatedChIP_std] = simulate_ChIP(S, W1, W0, c, ...
                      time_vec, parameters, gene_size, elongation_rate, total_residence_time, ...
                      number_trajectories, number_samples, n_sites);
figure();
bar([1:n_sites+2], SimulatedChIP); hold on;
er = errorbar([1:n_sites+2],SimulatedChIP, SimulatedChIP-SimulatedChIP_std, SimulatedChIP+SimulatedChIP_std);
er.Color = [0,0,0];
er.LineStyle = 'none';
title('Simulated ChIP Experiment')

%% Section 12 - Predict Peturbation Assay Dynamics
% Code that simulates perturbations to the transcription model by changing
% parameters by a percentage at a set time.

% User Defined parameters
number_trajectories = 100;
perturbation_time = 700;
time_vec = [0:1000];
MaxValue = 2;
zero_values = [30/200, 23/200,5/200];

% vector to multiply the parameters by after perturbation time
perturbation_vector = ones(size(parameters));
% add changes to parameters here:
perturbation_vector(4) = .1; %lower k_processing by 90%

% Simulated perturbation with SSA
[PerturbedMoleculeTrajectories] = simulated_perturbation(S, W1, W0, c,  time_vec, parameters, perturbation_vector, ...
                      perturbation_time, number_trajectories);
[PerturbedIntensityTrajectories] = add_shot_noise(PerturbedMoleculeTrajectories, parameters, noise_parameters, zero_values);

[PerturbedIntensityTrajectories] = normalize_simulated_trajectories_perturbation(PerturbedIntensityTrajectories, NChannels, ...
                                                        ChannelNames, NormalizationQuantile, ...
                                                        MaxValue, perturbation_time);


figure();
plot(PerturbedIntensityTrajectories.CTD.normalized_data(:,1),'r'); hold on;
plot(PerturbedIntensityTrajectories.Ser5ph.normalized_data(:,1),'g');
plot(PerturbedIntensityTrajectories.mRNA.normalized_data(:,1),'b');
plot([perturbation_time,perturbation_time],[-1,MaxValue],'k--',LineWidth=2)
legend('CTD','Ser5ph','mRNA')
title('Example Simulated Intensity Trace w/ Perturbations')

%analytical perturbation of moments with ODE solver

[PerturbedAnalyticalSolution] = solve_perturbed_analytical_solution(S, W1, W0, ...
                                                                     c, b, parameters, perturbation_vector);
figure()
plot([0:.1:33], PerturbedAnalyticalSolution(1,:),'r','LineWidth',2); hold on;
plot([0:.1:33],PerturbedAnalyticalSolution(2,:),'g--','LineWidth',2);
plot([0:.1:33],PerturbedAnalyticalSolution(3,:),'b','LineWidth',2);
plot([2.9,2.9],[0,25],'k--')
legend('CTD','Ser5ph','mRNA')
title('Perturbation Analytical Solution')
xlabel('Minutes')
ylabel('Mean Molecule Count')
