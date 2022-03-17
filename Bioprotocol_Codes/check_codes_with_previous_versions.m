clc
clear 
close all

addpath('./Functions','-end')
addpath('../Functions','-end')
addpath('../Data_files','-end')
%%

[mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap] = load_normalization_variance_gui(1,'G0_intp');

sig_dat = [rnap.mn_ac(1),ser5_rnap.mn_cc(11),mrna_rnap.mn_cc(11);...
    ser5_rnap.mn_cc(11),ser5.mn_ac(1),mrna_ser5.mn_cc(11);...
    mrna_rnap.mn_cc(11),mrna_ser5.mn_cc(11),mrna.mn_ac(1)];
sig_dat_sem = [rnap.sem_ac(1),ser5_rnap.sem_cc(11),mrna_rnap.sem_cc(11);...
    ser5_rnap.sem_cc(11),ser5.sem_ac(1),mrna_ser5.sem_cc(11);...
    mrna_rnap.sem_cc(11),mrna_ser5.sem_cc(11),mrna.sem_ac(1)];

mrna_rnap.sem_cc = mrna_rnap.sem_cc/mrna_rnap.mn_cc(11);
mrna_rnap.mn_cc = mrna_rnap.mn_cc/mrna_rnap.mn_cc(11);

mrna_ser5.sem_cc = mrna_ser5.sem_cc/mrna_ser5.mn_cc(11);
mrna_ser5.mn_cc = mrna_ser5.mn_cc/mrna_ser5.mn_cc(11);

ser5_rnap.sem_cc = ser5_rnap.sem_cc/ser5_rnap.mn_cc(11);
ser5_rnap.mn_cc = ser5_rnap.mn_cc/ser5_rnap.mn_cc(11);

mrna.mn_ac

load best_simple_pars

parameters(2) = 1;
parameters(7) = 1;
par_fixed = parameters;
par_changed = [1,3:13];
par_opt = log10(par_fixed(par_changed));
get_err = @(pars)sum(get_log_l_simplified(10.^pars,mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap,par_fixed,par_changed));

old_parameters = parameters;
old_LL = get_log_l_simplified(10.^par_opt,mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap,par_fixed,par_changed);

%%


% User Defined Parameters
FileName = 'Raw_Intensities Analysis-BLC_W_BG_WO_RunAve_May18';  % Intensity Trajectories (see Section XXX in Protocol).
DataCells = {'C2:C4001','D2:D4001','E2:E4001'};  % Cells or area containing intensity trajectories to normalize 
NChannels = length(DataCells); % Number of channels (3 in our case)
ChannelNames = {'CTD','Ser5ph','mRNA'}; 
NCells = 20;                   % number of cells in the data
LTrajectory = 200;             % number of time points in each trajectory
NormalizationQuantile = 0.95;  % the percentile to normalize as 1 (0-1)
MaxValue = 1.5;                % the maximum value after normalization 
RezeroAutoCorrelation = true;  % option to rezero autocorrelation tails
G0_type = 'G0_intp';           % G0 option type: G0 (raw g0), none (divide by 1), G0_intp (Interpolate G0 from G1,G2,G3), G1
CC_Norm = 'ACC_SQRT';          % Normalization type for CC, ACC_SQRT (use the sqrt(g0_signal1)*sqrt(g0_signal2))
CC_G0 = 'G0';                  % G0 option for cross correlations: G0 (use the center delay), max (use the max value), none (use none)
CC_delays = 10;                % How many time points to calculate delay out too for cross correlations
ACC_delays = 31;               % How many time points to calculate delay out too for auto correlations
RezeroPoints = 9;              % how many final points to rezero by for autocorrelations


% FISHQUANT / Quantification of mRNA amount
MeanSpotCount = 15.5;          % Average number of nascent mRNA per transcription site (measured quantity)
VarianceSpotCount = .93;       % Variance in number of nascent mRNA per transcription site (measured quantity)

[IntensityData, NormalizedCorrelationData] = get_normalized_ccs_from_xls(FileName, DataCells, NChannels, ...
                                                        ChannelNames, NCells, LTrajectory, NormalizationQuantile, ...
                                                        MaxValue, RezeroAutoCorrelation, RezeroPoints, G0_type, CC_G0, CC_Norm, CC_delays, ACC_delays);

%% Section 2 - Model Specification
% Specify the model here, see Forero 2020 methods section.

% selected model from forero2020, see figure 3c for model diagram
Nstates = 3;  % Number of species in the model

% Identified parameters for the final model.
parameters = [  0.433596648003337, ... % kon 1/min = omega
                1000, ...              % koff 1/min = fixed at 1000 for instantaneous bursts
                0.666624850282066, ... % kesc 1/min
                0.198727194675144, ... % kproc 1/min
                15.4054402060124, ...  % burst size, b
                0.777851806420281, ... % kout 1/min
                1, ...                 % (frac)tion of RNAP2 in State 2 that are Ser5ph labeled. Fixed at 1.
                1.98448561566709, ...  % eta_rnap
                1.42073542241570, ...  % eta_ser5
                0.408137598525255, ... % eta_ts
                0.924805163389989, ... % scale_rnap (scale of shot noise)
                1.22915843769510, ...  % scale_ser5 (scale of shot noise)
                1.14124029524139  ];   % scale_mrna (scale of shot noise)


% the system of ODES is:
% dx/dt = S*(W1(x) + W0) 
% where:
% S = stoichiometry matrix
S =  [[1    -1     0     0     0     0];
      [0     0     1    -1    -1     0];
      [0     0     0     0     1    -1]];   

% propensity functions have the affine linear form W = W1*x + W0
W1 = @(p)[[-p(1)      0         0];
          [p(2)       0         0];
          [p(5)*p(2)  0         0];
          [    0      p(6)      0];
          [    0      p(7)*p(3) 0];
          [    0      0         p(4)];];   

% independent of x
W0 = @(p)[[p(1); 0; 0; 0; 0; 0;]];

% signal transformation matrix c, this matrix converts states x to intensity channels
% e.g. state X1 adds 1 to both R channel and G channel since the RNAP2 has ser5ph and ctd

c = @(p)[[0,   1,  1]; 
         [0, p(7), 1];
         [0,   0,  1]];

noise_parameters = @(p)[p(8), p(9), p(10), p(11), p(12), p(13)];



% Function that solves the defined model for Correlations, Means, and Variances

time_vec = [0:30]; % Time span to solve over (minutes)
[ModelCorrelations, ModelMeans, ModelVariances] = solve_model_from_mats(S, W1, W0, c, ...
                                                                        noise_parameters, time_vec, parameters);



% User Defined Parameters
cc_range = 15; % cross correlation delay times to include in the log likelihood calculation
acc_range = 20; % auto correlation delay times to include in the log likelihood calculation
channels_quantified = [3]; % which channel was quantified for spot counts? in our case 3 for mRNA
DataMeans = [0,0,15.5]; % Data means for spot counts, 0 if not quantified
DataSEMs = [0,0,.93];   % SEM of the spots quantified, 0 if not quantified
Nspots = [0,0,130];     % Number of spots quantified, 0 if not quantified


par_fixed = parameters; % Which parameters are fixed
par_changed = [1,3:6,8:13]; %which parameters are free to change

[ModelLikelihood] = get_log_likelihood(parameters(par_changed), NormalizedCorrelationData, DataMeans, ...
                                       DataSEMs, S, W1, W0, c, ...
                                       noise_parameters, time_vec, cc_range, acc_range, ...
                                       channels_quantified, Nspots, ChannelNames, par_fixed, par_changed);
        
ModelLikelihoodTotal = sum(ModelLikelihood)
% After running this with the default parameters, the total negative log likelihood of
% the data given the model should be 14.558.



%%

ctd_mn_ac_err = sum(mean(NormalizedCorrelationData.CTD_CTD.mean_corr,1)  - rnap.mn_ac');
ser5_mn_ac_err = sum(mean(NormalizedCorrelationData.Ser5ph_Ser5ph.mean_corr,1)  - ser5.mn_ac');
mrna_mn_ac_err = sum(mean(NormalizedCorrelationData.mRNA_mRNA.mean_corr,1)  - mrna.mn_ac');

mrnactd_mn_ac_err = sum(mean(NormalizedCorrelationData.mRNA_CTD.mean_corr,1)  - mrna_rnap.mn_cc');
mrnaser5_mn_ac_err = sum(mean(NormalizedCorrelationData.mRNA_Ser5ph.mean_corr,1)  - mrna_ser5.mn_cc');
ser5ctd_mn_ac_err = sum(mean(NormalizedCorrelationData.Ser5ph_CTD.mean_corr,1)  - ser5_rnap.mn_cc');

sum([ctd_mn_ac_err, ser5_mn_ac_err, mrna_mn_ac_err, mrnactd_mn_ac_err, mrnaser5_mn_ac_err, ser5ctd_mn_ac_err])
%%

ctd_sem_ac_err = sum(mean(NormalizedCorrelationData.CTD_CTD.sem_corr,1)  - rnap.sem_ac');
ser5_sem_ac_err = sum(mean(NormalizedCorrelationData.Ser5ph_Ser5ph.sem_corr,1)  - ser5.sem_ac');
mrna_sem_ac_err = sum(mean(NormalizedCorrelationData.mRNA_mRNA.sem_corr,1)  - mrna.sem_ac');

mrnactd_sem_ac_err = sum(mean(NormalizedCorrelationData.mRNA_CTD.sem_corr,1)  - mrna_rnap.sem_cc');
mrnaser5_sem_ac_err = sum(mean(NormalizedCorrelationData.mRNA_Ser5ph.sem_corr,1)  - mrna_ser5.sem_cc');
ser5ctd_sem_ac_err = sum(mean(NormalizedCorrelationData.Ser5ph_CTD.sem_corr,1)  - ser5_rnap.sem_cc');

sum([ctd_sem_ac_err, ser5_sem_ac_err, mrna_sem_ac_err, mrnactd_sem_ac_err, mrnaser5_sem_ac_err, ser5ctd_sem_ac_err])


%%

ModelLikelihood - old_LL(1:8)



