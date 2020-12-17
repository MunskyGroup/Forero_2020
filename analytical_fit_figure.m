
%% Plot fit with data for activation/activation transcription during fluct.
%figure('visible', 'off');

%Final figure code that generates all computational model figures used for
%Forero2020 

%Set up figure sizes
clc; close all; clear all;
X_SIZE = 13; Y_SIZE = 15;
figure(1);clf;
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, X_SIZE, Y_SIZE]; % x,y, width, height
global_color = 'k';

subplot(3,1,1)
norm_intensity = 0;
minYlim = -0.1;
maxYlim = 5;
yyaxis left

rng(0)

%% Load in data from excel files and normalize them

norm_signal = 1; %normalized signal is true
seednum = floor(rand*1000);  %get random seed

%Load Cross Correlations with G0 normalization
[~,~,~,mrna_rnap,mrna_ser5,ser5_rnap] = load_normalization_variance_norm_cc(0);  %load normalization 21 pts for cross correlations
%Load Autocorrelations with G0_intp normalization and no rezeroing 
[mrna,ser5,rnap,~,~,~] = load_normalization_variance(1,'G0_intp','none',10);





% 
% model_options.trypt = [0,0,0,0,0,0,0];  
% model_options.burst = 1;
% model_options.processing = 1;
% model_options.elongation_randomness = 0;
% model_options.shot_noise = 1;
% model_options.photobleach=0;
% model_options.freepars = [1:10,12,13,14];
% model_options.calibration = 0;
% 
% model_options.trypt_mu = 0;
% model_options.trypt_var = 0;
% model_options.trypt = [0,0,0,0,0,0,0];
% model_options.trypt_time = 0;
% model_options.trypt_scale =0;

load('best_simple_pars.mat')  %load parameter file
real_valscl = parameters;

%% Generate analytical correlation signals

[sigred,TT,means,sig,TT2,mins] = get_ac_and_cc_mod_simplified(parameters,[0:.1:30]);

parnames = {'kon','kesc','kproc','beta','kout'};
par_changed = [1:5];

mh_pars = [];   %load in the methaste sampling data for plots
mh_vals = [];
ikeep = 1;
i=1;         
while ikeep==1
    try
        fn = ['met_hast_pars_2x_',num2str(i),'.mat'];
        load(fn)
        mh_pars = [mh_pars;mh_smpl];
        mh_vals = [mh_vals;mh_value];
        i=i+1;
    catch
        ikeep=0;
    end    
end


par_fixed = parameters;
par_changed = [1,3:6];
parameters_samp = parameters;

lags1 = [0:31];   

fntsize = 18;   %set up figure 3
close all
figure(1)
subplot(4,3,2)
fig1= gcf;
%set(gcf,'color','k');
fig1.PaperUnits = 'centimeters';
xh = 2*3*6;
yh = 2*14.3/3;
fig1.PaperPosition = [0, 0, 2/3*xh, 2*yh]; % x,y, width, height
%sem_rnap = std(Xdata_sem(:,1:30),1)/(sqrt(20));




plot(TT, sigred(1,:),'k-','linewidth',2)  %plot model rnap_acc
hold on;
%plot the data rnap_acc
b1=errorbar(lags1(1:25), rnap.mn_ac(1:25), rnap.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[256, 0, 0]./256,'Color',[256, 0, 0]./256);hold on;

par_fixed = parameters;
par_changed = [1,3:6];
parameters_samp = parameters;

minvec = sigred(1,:);
maxvec = sigred(1,:);   %randomly sample 50 mh par lines to generate gray lines
for i=1:50
    n = ceil(rand*size(mh_pars,1));
    

    par_samp = mh_pars(n,:);
     parameters_samp(par_changed) = par_samp;
    [sigred_samp,TT_samp,means_samp,sig_samp,TT2_samp,mins_samp] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
    %plot(TT_samp, sigred_samp(1,:),'-','Color',[.2,.2,.2],'linewidth',2)

    for j = 1:length(minvec)               %keep the minimum and maximum points
        if sigred_samp(1,j) > maxvec(1,j)
            maxvec(1,j) = sigred_samp(1,j);
            
        end    
         if sigred_samp(1,j) < minvec(1,j)
            minvec(1,j) = sigred_samp(1,j);
            
        end           
    end
 
end


%fill in the points between min and max MH sampled signals with gray
fill([TT',fliplr(TT')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);
b1=errorbar(lags1(1:25), rnap.mn_ac(1:25), rnap.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[256, 0, 0]./256,'Color',[256, 0, 0]./256);hold on;

plot(TT, sigred(1,:),'k-','linewidth',2) %replot so model acc signal is on top

% Setup labels and limits 
xlim([-.5,25])
plot([0,25],[0,0],'k--')
ylim([-.1,1.3])
set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color','Ycolor',global_color);
title({'\color{black}Autocorrelations'},'FontSize',fntsize,'FontWeight','bold')
set (gca ,'TickLength',[.01,.3],'LineWidth',1);
title({'CTD'},'FontSize',fntsize,'FontWeight','bold')


%%


subplot(4,3,5)      %SER5 ACC figure
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

rnap_color = [138, 25, 38]./256;
ser5_color = [20, 122, 47]./256;

plot(TT,sigred(5,:),'-k','linewidth',2)  %plot model ser5 acc
hold on;
% plot data ser5 acc
b1=errorbar(lags1(1:25), ser5.mn_ac(1:25), ser5.sem_ac(1:25),'s','MarkerSize',5,'MarkerFaceColor',[0, 256, 0]./256,'Color',[0, 256, 0]./256);hold on;


par_fixed = parameters;
par_changed = [1,3:6];
parameters_samp = parameters;

pt = 5;
minvec = sigred(pt,:);
maxvec = sigred(pt,:);

for i=1:50   %sample and get min / max from 50 random mh signals
    n = ceil(rand*size(mh_pars,1));
    
    parameters_samp(par_changed) = par_samp;
    par_samp = mh_pars(n,:);
     parameters_samp(par_changed) = par_samp;
    [sigred_samp,TT_samp,means_samp,sig_samp,TT2_samp,mins_samp] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
    %plot(TT_samp, sigred_samp(5,:),'-','Color',[.8,.8,.8],'linewidth',2)
    
    for j = 1:length(sigred_samp)
        if sigred_samp(pt,j) > maxvec(1,j)
            maxvec(1,j) = sigred_samp(pt,j);
            
        end    
         if sigred_samp(pt,j) < minvec(1,j)
            minvec(1,j) = sigred_samp(pt,j);
            
        end           
    end
    
end
%fill in the gray area and replot other signals on top
fill([TT',fliplr(TT')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);
b1=errorbar(lags1(1:25), ser5.mn_ac(1:25), ser5.sem_ac(1:25),'s','MarkerSize',5,'MarkerFaceColor',[0, 256, 0]./256,'Color',[0, 256, 0]./256);hold on;
plot(TT,sigred(5,:),'-k','linewidth',2)

xlim([-.5,25])   %figure labels and limits
plot([0,25],[0,0],'k--')
ylim([-.1,1.3])
title({'Ser5ph'},'FontSize',fntsize,'FontWeight','bold')
set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
set (gca ,'TickLength',[.01,.3],'LineWidth',1);

%%
subplot(4,3,8)   %mRNA ACC figure
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

rnap_color = [138, 25, 38]./256;
ser5_color = [20, 122, 47]./256;
mrna_color = [22, 68, 148]./256;

plot(TT,sigred(9,:),'-k','linewidth',2)  %model mrna ACC
hold on;
%data mrna ACC
b1=errorbar(lags1(1:25), mrna.mn_ac(1:25), mrna.sem_ac(1:25),'d','MarkerSize',5,'MarkerFaceColor',[0, 0, 256]./256,'Color',[0, 0, 256]./256);hold on;

par_fixed = parameters;
par_changed = [1,3:6];
parameters_samp = parameters;

pt = 9;
minvec = sigred(pt,:);
maxvec = sigred(pt,:);

for i=1:50   %get mh min and max for gray fill
    n = ceil(rand*size(mh_pars,1));
    
    parameters_samp(par_changed) = par_samp;
    par_samp = mh_pars(n,:);
     parameters_samp(par_changed) = par_samp;
    [sigred_samp,TT_samp,means_samp,sig_samp,TT2_samp,mins_samp] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
    %plot(TT_samp, sigred_samp(9,:),'-','Color',[.8,.8,.8],'linewidth',2)
    for j = 1:length(sigred_samp)
        if sigred_samp(pt,j) > maxvec(1,j)
            maxvec(1,j) = sigred_samp(pt,j);
            
        end    
         if sigred_samp(pt,j) < minvec(1,j)
            minvec(1,j) = sigred_samp(pt,j);
            
        end           
    end
    
end

%fill the min max mh pars, and replot signals on top
fill([TT',fliplr(TT')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);
b1=errorbar(lags1(1:25), mrna.mn_ac(1:25), mrna.sem_ac(1:25),'d','MarkerSize',5,'MarkerFaceColor',[0, 0, 256]./256,'Color',[0, 0, 256]./256);hold on;
plot(TT,sigred(9,:),'-k','linewidth',2)

xlim([-.5,25])  %Limits and labels
plot([0,25],[0,0],'k--')
ylim([-.1,1.3])
title({'mRNA'},'FontSize',fntsize,'FontWeight','bold')
set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
set (gca ,'TickLength',[.01,.3],'LineWidth',1);

saveas(gca,'Accs.epsc')  %save the autocorrelations

%%
lags2 = [-10:10]; %time vector for CCs

subplot(4,3,3)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

% Ser5 RNAP cross correlation
plot([-TT(end:-1:2);TT],[sigred(4,end:-1:2),sigred(2,:)],'k-','linewidth',2) %model
hold on;
%data 
b1=errorbar(lags2, ser5_rnap.mn_cc, ser5_rnap.sem_cc,'s','MarkerSize',5,'MarkerFaceColor',[0, 0.67, 1],'Color',[0, 0.67, 1]);hold on;

par_fixed = parameters;
par_changed = [1,3:6];
parameters_samp = parameters;


minvec = [sigred(4,end:-1:2),sigred(2,:)];  %pull the mh min max for fill of gray 
maxvec = [sigred(4,end:-1:2),sigred(2,:)];

for i=1:50
    n = ceil(rand*size(mh_pars,1));
    
    parameters_samp(par_changed) = par_samp;
    par_samp = mh_pars(n,:);
     parameters_samp(par_changed) = par_samp;
    [sigred_samp,TT_samp,means_samp,sig_samp,TT2_samp,mins_samp] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
    %plot([-TT_samp(end:-1:2);TT_samp],[sigred_samp(4,end:-1:2),sigred_samp(2,:)],'-','Color',[.8,.8,.8],'linewidth',2)
    
    svec = [sigred_samp(4,end:-1:2),sigred_samp(2,:)];
    for j = 1:length(svec)
        if svec(1,j) > maxvec(1,j)
            maxvec(1,j) = svec(1,j);
            
        end    
         if svec(1,j) < minvec(1,j)
            minvec(1,j) = svec(1,j);
            
        end           
    end
end

fill([[-TT(end:-1:2);TT]',fliplr([-TT(end:-1:2);TT]')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);

b1=errorbar(lags2, ser5_rnap.mn_cc, ser5_rnap.sem_cc,'s','MarkerSize',5,'MarkerFaceColor',[0, 0.67, 1],'Color',[0, 0.67, 1]);hold on;

plot([-TT(end:-1:2);TT],[sigred(4,end:-1:2),sigred(2,:)],'k-','linewidth',2) %replot signals
plot([0,0],[0,1.3],'--','Color',global_color)

%plotting labels and legends and limits
set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
title({'\color{black}Cross-correlations'},'FontSize',20,'FontWeight','bold')
title({'Ser5ph-CTD'},'FontSize',fntsize,'FontWeight','bold')
xlim([-10,10])
ylim([-.1,1.3]) 
set (gca ,'TickLength',[.01,.3],'LineWidth',1);

%%
subplot(4,3,6)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

plot([-TT(end:-1:2);TT],[sigred(7,end:-1:2),sigred(3,:)],'k-','linewidth',2) % mrna_rnap CC model
hold on;  %mrna_rnap CC data
b1=errorbar(lags2, mrna_rnap.mn_cc, mrna_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[1, .5, 0],'Color',[1, .5, 0]);hold on;

par_fixed = parameters;
par_changed = [1,3:6];
parameters_samp = parameters;

minvec = [sigred(7,end:-1:2),sigred(3,:)];
maxvec = [sigred(7,end:-1:2),sigred(3,:)];


for i=1:50  %load min and max mh_pars for gray bars
    n = ceil(rand*size(mh_pars,1));
    
    parameters_samp(par_changed) = par_samp;
    par_samp = mh_pars(n,:);
     parameters_samp(par_changed) = par_samp;
    [sigred_samp,TT_samp,means_samp,sig_samp,TT2,mins] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
   % plot([-TT_samp(end:-1:2);TT_samp],[sigred_samp(7,end:-1:2),sigred_samp(3,:)],'-','Color',[.8,.8,.8],'linewidth',2)
    
    svec = [sigred_samp(7,end:-1:2),sigred_samp(3,:)];
    for j = 1:length(svec)
        if svec(1,j) > maxvec(1,j)
            maxvec(1,j) = svec(1,j);
            
        end    
         if svec(1,j) < minvec(1,j)
            minvec(1,j) = svec(1,j);
            
        end           
    end
end

fill([[-TT(end:-1:2);TT]',fliplr([-TT(end:-1:2);TT]')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);
b1=errorbar(lags2, mrna_rnap.mn_cc, mrna_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[1, .5, 0],'Color',[1, .5, 0]);hold on;
plot([-TT(end:-1:2);TT],[sigred(7,end:-1:2),sigred(3,:)],'k-','linewidth',2)
plot([0,0],[0,1.3],'--','Color',global_color)

set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
set (gca ,'TickLength',[.01,.3],'LineWidth',1);
title({'mRNA-CTD'},'FontSize',fntsize,'FontWeight','bold')
xlim([-10,10])
ylim([-.1,1.3])

%%

subplot(4,3,9) %mrna_ser5 signal
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

plot([-TT(end:-1:2);TT],[sigred(8,end:-1:2),sigred(6,:)],'k-','linewidth',2)
hold on;
b1=errorbar(lags2, mrna_ser5.mn_cc, mrna_ser5.sem_cc,'d','MarkerSize',5,'MarkerFaceColor',[.5, 0, 1],'Color',[.5, 0, 1]);hold on;



par_fixed = parameters;
par_changed = [1,3:6];
parameters_samp = parameters;

minvec = [sigred(8,end:-1:2),sigred(6,:)];
maxvec = [sigred(8,end:-1:2),sigred(6,:)];

for i=1:50
    n = ceil(rand*size(mh_pars,1));
    
    parameters_samp(par_changed) = par_samp;
    par_samp = mh_pars(n,:);
     parameters_samp(par_changed) = par_samp;
    [sigred_samp,TT_samp,means_samp,sig_samp,TT2_samp,mins_samp] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
    %plot([-TT_samp(end:-1:2);TT_samp],[sigred_samp(8,end:-1:2),sigred_samp(6,:)],'-','Color',[.8,.8,.8],'linewidth',2)
    svec = [sigred_samp(8,end:-1:2),sigred_samp(6,:)];
    for j = 1:length(svec)
        if svec(1,j) > maxvec(1,j)
            maxvec(1,j) = svec(1,j);
            
        end    
         if svec(1,j) < minvec(1,j)
            minvec(1,j) = svec(1,j);
            
        end           
    end
end

fill([[-TT(end:-1:2);TT]',fliplr([-TT(end:-1:2);TT]')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);
b1=errorbar(lags2, mrna_ser5.mn_cc, mrna_ser5.sem_cc,'d','MarkerSize',5,'MarkerFaceColor',[.5, 0, 1],'Color',[.5, 0, 1]);hold on;
plot([-TT(end:-1:2);TT],[sigred(8,end:-1:2),sigred(6,:)],'k-','linewidth',2)
plot([0,0],[0,1.3],'--','Color',global_color)

set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
set (gca ,'TickLength',[.01,.3],'LineWidth',1);
title({'mRNA-Ser5ph'},'FontSize',fntsize,'FontWeight','bold')

xlim([-10,10])
ylim([-.1,1.3])

saveas(gca,'ccs_nodist.epsc')

%% SSA trajectory for bottom of figure 3
kon = parameters(1);
koff = 1000;%parameters(2);
kesc = parameters(3);
kproc = parameters(4);
kin =  1000*parameters(5);
kout = parameters(6);
frac = parameters(7);
eta_rnap = parameters(8);
eta_ser5 = parameters(9);
eta_ts = parameters(10);


Nstates = 3;  %set up propensity and stoich
b = zeros(Nstates,1);
b(1) = kon;

c = zeros(3,3);
c(1,1:3)=[0,1,1];
c(2,1:3)=[0,frac,1];
c(3,3)=1;

S = zeros(3,6);
W1 = zeros(6,3);
W0 = zeros(6,1);
S(1,1) = 1;  W1(1,1) = -kon; W0(1,1) = kon;
S(1,2) = -1; W1(2,1) = koff; 

S(2,3) = 1; W1(3,1) = kin; 
S(2,4) = -1; W1(4,2) = kout; 

S(2:3,5) = [-1;1]; W1(5,2) = kesc*frac; 
S(3,6) = -1; W1(6,3) = kproc; 


x0 = [0,0,0]';
T_array = [0:1:1000];
time_var = 0;
signal_update_rate = 0;

W = @(x) W1*x + W0;
rng(45)
%Solve a singal trajectory 
sol = run_single_SSA_linda(x0,S,W,T_array,time_var,signal_update_rate);  

pol2_ssa = sol(2,:)';
ser5_ssa = sol(2,:)';
ts_ssa = sol(3,:)';

[pol2_ssa,ser5_ssa,ts_ssa,~] = get_model_intesities(sol,eta_rnap,eta_ser5,eta_ts); %convert molecules to signal
[pol2norm,ser5norm,tsnorm] = Normalize_simulated_intensities(.95,pol2_ssa,ser5_ssa,ts_ssa);

X_SIZE = 13; Y_SIZE = 15;
figure(35);clf;
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, X_SIZE, Y_SIZE]; % x,y, width, height
global_color = 'k';

cla;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 1.1*yh]; % x,y, width, height


avpol2 =movmean(pol2norm(end-200:end),3);
avser5 = movmean(ser5norm(end-200:end),3);
avts = movmean(tsnorm(end-200:end),3);
plot(avpol2(1:201),'r','linewidth',2); hold on; plot(avser5(1:201),'g','linewidth',2); plot(avts(1:201),'b','linewidth',2 )
legend('POL2','SER5ph','mRNA')
title('Representative trace','FontSize',fntsize,'FontWeight','bold')
ylabel('Normalized intensity','FontSize',fntsize,'FontWeight','bold')
xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
xlim([0,200])
ylim([-.5,1.5])





saveas(gca, 'trace.epsc')  %save the whole figure



%% Acov taus
T_array = [0:1:1000];

mrna_bootstrapped_tau = [];
ser5_bootstrapped_tau = [];
rnap_bootstrapped_tau = [];

mrna_sim_acov = zeros(400,32);
ser5_sim_acov = zeros(400,32);
rnap_sim_acov = zeros(400,32);

k = 0;
% while k < 400
%     t = [0:1:199];
% 
% 
%     simulated_pol2 = zeros(20,200);
%     simulated_ser5 = zeros(20,200);
%     simulated_mrna = zeros(20,200);
% 
%     for i = 1:20
%         sol = run_single_SSA_linda(x0,S,W,T_array,time_var,signal_update_rate);  
%         [pol2_ssa,ser5_ssa,ts_ssa,~] = get_model_intesities(sol,0,0,0); %convert molecules to signal
%         [pol2norm,ser5norm,tsnorm] = Normalize_simulated_intensities(.95,pol2_ssa,ser5_ssa,ts_ssa);
%         simulated_pol2(i,:) = pol2norm(end-199:end);
%         simulated_ser5(i,:) = ser5norm(end-199:end);
%         simulated_mrna(i,:) = tsnorm(end-199:end);
%     end
% 
%     [mrna_sim,ser5_sim,rnap_sim,~,~,~,G1_rnap_sim,G1_ser5_sim,G1_mrna_sim] = get_simulated_cc([0:1:199],simulated_pol2', simulated_ser5', simulated_mrna',0);
% 
% 
%     decorr_mrna = t(mrna_sim.mn_ac < .01);
%     if length(decorr_mrna) == 0
%         continue
%     end
%     
%     
%     decorr_ser5 = t(ser5_sim.mn_ac < .01);
% 
%     if length(decorr_ser5) == 0
%         continue
%     end
%     decorr_rnap = t(rnap_sim.mn_ac < .01);
% 
%     if length(decorr_rnap) == 0
%         continue
%     end
%     
%     k = k + 1
%     mrna_bootstrapped_tau = [mrna_bootstrapped_tau, decorr_mrna(1)];
%     ser5_bootstrapped_tau = [ser5_bootstrapped_tau, decorr_ser5(1)];
%     rnap_bootstrapped_tau = [rnap_bootstrapped_tau, decorr_rnap(1)];
%     
%     mrna_sim_acov(k,:) = mrna_sim.mn_ac;
%     ser5_sim_acov(k,:) = ser5_sim.mn_ac;
%     rnap_sim_acov(k,:) = rnap_sim.mn_ac;
% 
% end

load('rnap_sim_acov_nonoise.mat')
load('ser5_sim_acov_nonoise.mat')
load('mrna_sim_acov_nonoise.mat')


%%

thresh = .2;
t = [0:1:199];
mrna_bootstrapped_tau = [];
ser5_bootstrapped_tau = [];
rnap_bootstrapped_tau = [];
for i = 1:400
   tmpacov = t(mrna_sim_acov(i,:) < thresh);
 
   mrna_bootstrapped_tau = [mrna_bootstrapped_tau, tmpacov(1)];

   tmpacov = t(ser5_sim_acov(i,:) < thresh);
   ser5_bootstrapped_tau = [ser5_bootstrapped_tau, tmpacov(1)];

   tmpacov = t(rnap_sim_acov(i,:) < thresh);
   rnap_bootstrapped_tau = [rnap_bootstrapped_tau, tmpacov(1)];
 
end

% figure(32)
% histogram(mrna_bootstrapped_tau,'FaceColor','b')
% hold on;
% histogram(ser5_bootstrapped_tau,'FaceColor','g')
% histogram(rnap_bootstrapped_tau,'FaceColor','r')



%%
figure(1)
% subplot(4,3,2)
% plot( [mean(rnap_bootstrapped_tau),mean(rnap_bootstrapped_tau)], [-1,6],'r', 'LineWidth',2)
% x1 = mean(rnap_bootstrapped_tau)-std(rnap_bootstrapped_tau);
% x2 = mean(rnap_bootstrapped_tau)+std(rnap_bootstrapped_tau);
% fill( [x1,x2, x2,x1], [-1,-1,6,6],'r', 'FaceAlpha',.2,'LineStyle','none')
% 
% subplot(4,3,5)
% plot( [mean(ser5_bootstrapped_tau),mean(ser5_bootstrapped_tau)], [-1,6],'g', 'LineWidth',2)
% x1 = mean(ser5_bootstrapped_tau)-std(ser5_bootstrapped_tau);
% x2 = mean(ser5_bootstrapped_tau)+std(ser5_bootstrapped_tau);
% fill( [x1,x2, x2,x1], [-1,-1,6,6],'g', 'FaceAlpha',.2,'LineStyle','none')
% 
% subplot(4,3,8)
% plot( [mean(mrna_bootstrapped_tau),mean(mrna_bootstrapped_tau)], [-1,6],'b', 'LineWidth',2)
% 
% 
% x1 = mean(mrna_bootstrapped_tau)-std(mrna_bootstrapped_tau);
% x2 = mean(mrna_bootstrapped_tau)+std(mrna_bootstrapped_tau);
% fill( [x1,x2, x2,x1], [-1,-1,6,6],'b', 'FaceAlpha',.2,'LineStyle','none')
% 
% 





load('rnap_sim_acov.mat')
load('ser5_sim_acov.mat')
load('mrna_sim_acov.mat')

t = [0:1:199];
mrna_bootstrapped_tau = [];
ser5_bootstrapped_tau = [];
rnap_bootstrapped_tau = [];
for i = 1:400
   tmpacov = t(mrna_sim_acov(i,:) < thresh);
 
   mrna_bootstrapped_tau = [mrna_bootstrapped_tau, tmpacov(1)];

   tmpacov = t(ser5_sim_acov(i,:) < thresh);
   ser5_bootstrapped_tau = [ser5_bootstrapped_tau, tmpacov(1)];

   tmpacov = t(rnap_sim_acov(i,:) < thresh);
   rnap_bootstrapped_tau = [rnap_bootstrapped_tau, tmpacov(1)];
 
end


figure(1)
subplot(4,3,2)
plot( [mean(rnap_bootstrapped_tau),mean(rnap_bootstrapped_tau)], [-1,6],'r', 'LineWidth',2)
x1 = mean(rnap_bootstrapped_tau)-std(rnap_bootstrapped_tau);
x2 = mean(rnap_bootstrapped_tau)+std(rnap_bootstrapped_tau);
fill( [x1,x2, x2,x1], [-1,-1,6,6],'r', 'FaceAlpha',.2,'LineStyle','--','edgecolor','r')

subplot(4,3,5)
plot( [mean(ser5_bootstrapped_tau),mean(ser5_bootstrapped_tau)], [-1,6],'g', 'LineWidth',2)
x1 = mean(ser5_bootstrapped_tau)-std(ser5_bootstrapped_tau);
x2 = mean(ser5_bootstrapped_tau)+std(ser5_bootstrapped_tau);
fill( [x1,x2, x2,x1], [-1,-1,6,6],'g', 'FaceAlpha',.2,'LineStyle','--','edgecolor','g')

subplot(4,3,8)
plot( [mean(mrna_bootstrapped_tau),mean(mrna_bootstrapped_tau)], [-1,6],'b', 'LineWidth',2)


x1 = mean(mrna_bootstrapped_tau)-std(mrna_bootstrapped_tau);
x2 = mean(mrna_bootstrapped_tau)+std(mrna_bootstrapped_tau);
fill( [x1,x2, x2,x1], [-1,-1,6,6],'b', 'FaceAlpha',.2,'LineStyle','--','edgecolor','b')




saveas(gca, 'corrs_with_dwell.epsc')  %save the whole figure
return

%% Supplemental Figure
% Plot the molecule signals with bursting highlighted and filled in 

figure(17)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, yh]; % x,y, width, height
T_array = [0:.1:200];
sol = run_single_SSA_linda(x0,S,W,T_array,time_var,signal_update_rate);

pol2_ssa = (sol(2,:) + sol(3,:))';
ser5_ssa = (sol(2,:)+ sol(3,:))';
ts_ssa = sol(3,:)';


bubble_top = pol2_ssa(end-500:end);
bubble_bottom = ts_ssa(end-500:end);

p2 = pol2_ssa(end-500:end)-ts_ssa(end-500:end);
on = (p2 > 5);
fill([T_array(1:501),fliplr(T_array(1:501))], [(on*500)'-500, fliplr((on*500)')],[220,220,220]./256,'linestyle','None','HandleVisibility','off' )
hold on;
fill([T_array(1:501),fliplr(T_array(1:501))], [bubble_top', fliplr(bubble_bottom')], [255, 212, 212]./256, 'linestyle','None')


plot(T_array(1:501),pol2_ssa(end-500:end),'r','linewidth',2); hold on;  plot(T_array(1:501),ts_ssa(end-500:end),'b','linewidth',2 )

ylim([0,150])
xlabel({'Time (min)'},'FontSize',fntsize,'FontWeight','bold')
ylabel({'Simulated molecule counts'},'FontSize',fntsize,'FontWeight','bold')
a = legend({'Burst','CTD','mRNA'},'FontSize',fntsize,'FontWeight','bold');
%set (gca ,'TickLength',[.01,.3],'LineWidth',2);
set(a, 'Box', 'off');
saveas(gca,'Burst_traj.epsc')




%% histograms

data = xlsread('Nasc_TranscriptsCount.xlsx');
mrna_hist = data(:,1);
%pol2_hist = data(:,1);

mrna_hist = mrna_hist(1:end-5);
%pol2_hist = pol2_hist(1:end-5);



pol2_ssa_traj = [];
ser5_ssa_traj = [];
ts_ssa_traj = [];
mrnai_traj = [];
rnap_c = [];
ser5_c = [];
mrna_c = [];
for i = 1:1000
    sol = run_single_SSA_linda(x0,S,W,[0:1:2000],time_var,signal_update_rate);
    pol2_ssa = sol(2,:)';
    ser5_ssa = sol(2,:)';
    ts_ssa = sol(3,:)';
    mrnai = ts_ssa;
    
 [pol2_ssa,ser5_ssa,ts_ssa,~] = get_model_intesities(sol,eta_rnap,eta_ser5,eta_ts);   

    %mrnai = ts_ssa;



    [pol2norm,ser5norm,tsnorm] = Normalize_simulated_intensities(.95,pol2_ssa,ser5_ssa,ts_ssa);

    
    pol2_ssa_traj = [pol2_ssa_traj, pol2norm(1:200:end)];
    ser5_ssa_traj = [ser5_ssa_traj, ser5norm(1:200:end)];
    ts_ssa_traj = [ts_ssa_traj, ts_ssa(1:200:end)];
    mrnai_traj = [mrnai_traj, mrnai(1:200:end)];
    
  

end


pol2_traj = pol2_ssa_traj';
ser5_traj = ser5_ssa_traj';
ts_traj = ts_ssa_traj';

[rnap_I, ser5_I, mrna_I] = Normalize_raw_intensities(.95);   %get the normalized intensities from the data 95 percentile
rnap_data = reshape(rnap_I(1:20:end,:),1,200);  %rnap_data sorted by 20 time incrementes (decorrelated) this was for histograms
ser5_data = reshape(ser5_I(1:20:end,:),1,200);
mrna_data = reshape(mrna_I(1:20:end,:),1,200);

%xh = 6.3;
%yh = 6.4;

%%
fntsize = 18;


figure(11)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0,xh/3, 3*yh/3]; % x,y, width, height



[x,n] = hist(reshape(pol2_traj.',1,[]),30);
[11, 252, 3]./256
[218, 51, 255]./256
x1 = histogram(reshape(pol2_traj.',1,[]),n,'Normalization','probability','FaceColor',[11, 252, 3]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[11, 252, 3]./256)%,'DisplayStyle','stairs')
hold on;
x2 = histogram(reshape(pol2_traj.',1,[]),n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256,'DisplayStyle','stairs')

x3 = histogram(rnap_data ,n,'Normalization','probability','FaceColor',[218, 51, 255]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[218, 51, 255]./256)%,'DisplayStyle','stairs')

x4 = histogram(rnap_data ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')

a = legend(gca, [x1,x3],{'Model','Data'},'Location','Best')
set(a, 'Box', 'off');

set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
set(a,'TextColor',global_color);

ylim([0,.15])
xlim([-1.3,1.5])
set(gca,'linewidth',2)
saveas(gca,'pol2_hist.epsc')
set(gca,'YTickLabel',[],'XtickLabel',[])

saveas(gca,'pol2_hist_nolabels.epsc')
%%
% 
% title({'CTD Normalized Intensity'},'FontSize',fntsize,'FontWeight','bold','Color',global_color)
% xlabel({'Intensity (Norm)'},'FontSize',fntsize,'FontWeight','bold')
% ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
% 
% ylim([minYlim maxYlim])
% set (gca ,'TickLength',[.01,.3],'LineWidth',1);
% set (gca ,'FontSize',fntsize,'FontName', 'Arial');
% 
% %set(gca,'Color','k')
% %set(gcf,'Color','k')
% set(gca,'YTickLabel',[],'XtickLabel',[])

% saveas(gca, 'pol2hist.epsc')

figure(10)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh/3, 3*yh/3]; % x,y, width, height
[x,n] = hist(reshape(ser5_traj.',1,[]),30);
x1 = histogram(reshape(ser5_traj.',1,[]),n,'Normalization','probability','FaceColor',[11, 252, 3]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[11, 252, 3]./256)%,'DisplayStyle','stairs')
hold on;
x2 = histogram(reshape(ser5_traj.',1,[]),n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256,'DisplayStyle','stairs')

x3 = histogram(ser5_data ,n,'Normalization','probability','FaceColor',[218, 51, 255]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[218, 51, 255]./256)%,'DisplayStyle','stairs')

x4 = histogram(ser5_data ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')
%a = legend(gca,[x1,x3], {'Model','Data'},'Location','Best')
%set(a, 'Box', 'off');
ylim([0,.15])
xlim([-1.3,1.5])
set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
set(a,'TextColor',global_color);
set(gca,'linewidth',2)
saveas(gca, 'ser5hist.epsc')
set(gca,'YTickLabel',[],'XtickLabel',[])
saveas(gca,'SER5_dist_nolabels.epsc') 

%%

% %legend('Model','Data')
% title({'SER5 Normalized Intensity'},'FontSize',fntsize,'FontWeight','bold','Color',global_color)
% xlabel({'Intensity (Norm)'},'FontSize',fntsize,'FontWeight','bold')
% ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
% %ylim([minYlim maxYlim])
% %ylim([0,.2])
% set (gca ,'TickLength',[.01,.3],'LineWidth',1);
% set (gca ,'FontSize',fntsize,'FontName', 'Arial');


%set(gca,'Color','k')
%set(gcf,'Color','k')

figure(12)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh/3, 3*yh/3]; % x,y, width, height

[x,n] = hist(mrnai_traj,30)

x1 = histogram(mrnai_traj,n,'Normalization','probability','FaceColor',[11, 252, 3]./256,'FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256)%'DisplayStyle','stairs')
hold on;
%histogram(mrna_hist ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')
x2 = histogram(mrnai_traj,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256,'DisplayStyle','stairs')

x3 = histogram(mrna_hist ,n,'Normalization','probability','FaceColor',[218, 51, 255]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[218, 51, 255]./256)%,'DisplayStyle','stairs')

x4 = histogram(mrna_hist ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')

%a = legend(gca,[x1,x3], {'Model','Data'},'Location','Best')
%set(a, 'Box', 'off');

set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
set(a,'TextColor',global_color);
set(gca,'linewidth',2)
saveas(gca, 'tshist_labels.epsc') 
set(gca,'YTickLabel',[],'XtickLabel',[])
saveas(gca, 'tshist_nolabels.epsc') 

% title({'mRNA count'},'FontSize',fntsize,'FontWeight','bold','Color',global_color)
% xlabel({'N molecules'},'FontSize',fntsize,'FontWeight','bold')
% ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
% %ylim([minYlim maxYlim])
% set (gca ,'TickLength',[.01,.3],'LineWidth',1);
% set (gca ,'FontSize',fntsize,'FontName', 'Arial');
%set(gca,'Color','k')



%% MH plots

addpath ../Data_files/
addpath ../

parnames = {'kon','kesc','kproc','beta','kout'};
parnames = {'beta','omega','k out','k esc','k complete'};

par_changed = [1:5];

mh_pars = [];
mh_vals = [];
ikeep = 1;
i=1;
while ikeep==1
    try
        fn = ['met_hast_pars_2x_',num2str(i),'.mat'];
        load(fn)
        mh_pars = [mh_pars;mh_smpl];
        mh_vals = [mh_vals;mh_value];
        i=i+1;
    catch
        ikeep=0;
    end    
end

Np = 5;

mh_pars = mh_pars(:,[4,1,5,2,3]);

sz = 4*(1+max(mh_vals)-mh_vals);
figure(6)
xh = 2*3*6;
yh = 2*14.3/3;
fntsize = 18;
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height
for i=1:Np
    subplot(Np,Np,(i-1)*Np+i)

    H = histogram(mh_pars(:,i),20,'Normalization','Probability'); hold on;
    %H = histogram(mh_pars(:,i),20,'Normalization','Probability'); hold on;

    X =  sort(mh_pars(:,i));
    low = X(floor(length(X)/10));
    high = X(ceil(9*length(X)/10));
    
    plot(mh_pars(1,i)*[1,1],get(gca,'ylim'),'r--','linewidth',3)

    title([parnames{par_changed(i)}])

    
    if i==1
        ylabel({parnames{par_changed(i)}},'FontSize',fntsize,'FontWeight','bold');
    elseif i==Np
        xlabel({parnames{par_changed(i)}},'FontSize',fntsize,'FontWeight','bold');
    end
    %     hold on
    %     pr = max(0,log(H.BinEdges/100));
    %     switch parnames{par_changed(i)}
    %         case 'kproc'
    %             pr = pr + log(H.BinEdges/(1/(103/60)));
    %         case 'ke'
    %             pr = pr + log(H.BinEdges/4.1);
    %     end
    
    %     [qw,qe] =max(H.Values);
    
    %     pr = max(H.Values)*exp(-pr)/exp(-pr(qe));
    %     plot(H.BinEdges,pr,'r','linewidth',2)

        set (gca ,'TickLength',[.01,.3],'LineWidth',2);
    set (gca ,'FontSize',fntsize,'FontName', 'Arial');
    for j=i+1:Np
        subplot(Np,Np,(j-1)*Np+i)

        scatter(mh_pars(1:100:end,i),mh_pars(1:100:end,j),sz(1:100:end),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerFaceColor',[50, 139, 191]./256,'MarkerEdgeColor',[50, 139, 191]./256); hold on
        scatter(mh_pars(1:100:end,i),mh_pars(1:100:end,j),sz(1:100:end),'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'MarkerFaceColor',[1,.5,0]); hold on

        plot(mh_pars(1,i),mh_pars(1,j),'ko','markersize',8,'markerfacecolor','k')
        if i==1
            ylabel({parnames{par_changed(j)}},'FontSize',fntsize,'FontWeight','bold');
        end
        if j==Np
            xlabel({parnames{par_changed(i)}} ,'FontSize',fntsize,'FontWeight','bold');
        end
        set (gca ,'TickLength',[.01,.3],'LineWidth',2);
        set (gca ,'FontSize',fntsize,'FontName', 'Arial');
    end
    
    
end

saveas(gca,'sensitivity.epsc')


kon = mh_pars(:,1);
koff = 1000;
kesc = mh_pars(:,2);
kproc = mh_pars(:,3);
kin = mh_pars(:,4)*koff;
kout = mh_pars(:,5);

mu_ctd = kon./(kon+koff).*kin./(kout+kesc);
mu_ctd(1)
% mu_ctd_on = kin./(kout+kesc);

mu_rna = kon./(kon+koff).*kin./(kout+kesc).*kesc./kproc;
mu_rna(1)
% mu_rna_on = kin./(kout+kesc).*kesc./kproc;

rate_prod = mu_ctd.*kesc;
% rate_prod_on = mu_ctd_on.*kesc;

burst_refractory_period = 1./kon;
burst_duration = 1./koff;
pol2_burst_size = kin./koff;

burst_efficiency = 100*kesc./(kesc+kout);
mrna_burst_size = pol2_burst_size.*burst_efficiency/100;

figure(7);


fntsize = 18;
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height
sbplots = {'mu_ctd','mu_rna','rate_prod',...
    'burst_refractory_period',...
    'pol2_burst_size','burst_efficiency','mrna_burst_size'};
sbplots_leg = {'\mu ctd (mol)','\mu rna (mol)',...
    'rate prod (mol/min)',...
    'burst refractory period (min)',...
    'pol2 burst size (mol)','burst efficiency (%)','mrna burst size (mol)'}
for i=1:7
    subplot(3,3,i); hold off
    histogram(eval(sbplots{i}),'Normalization','probability','FaceColor',[11, 252, 3]./256,'FaceAlpha',.5,'linewidth',.5,'edgecolor',[11, 252, 3]./256 ); hold on
    plot(eval([sbplots{i},'(1)'])*[1,1],get(gca,'ylim'),'k--','linewidth',3)
    
    X =  sort(eval(sbplots{i}));
    low = X(floor(length(X)/10));
    high = X(ceil(9*length(X)/10));
    
    title([sbplots_leg{i},'; ',num2str(eval([sbplots{i},'(1)']),3),...
        ' (',num2str(low,3),', ',num2str(high,3),')'],'FontSize',fntsize-5,'FontWeight','bold')
    if mod(i,3) == 1 
        ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
    end
    if ismember(i,[7,6,5])
        xlabel({'Value'},'FontSize',fntsize,'FontWeight','bold')
    end
    
end

saveas(gca,'par.epsc')
%% new inhib figure

inhibs = ones(size(parameters));
n = n+1;

inhibs(1) = .001; % inhibit kon
pol2_traj = [];
ser5_traj = [];
ts_traj = [];

pol2_traji = [];
ser5_traji = [];
ts_traji = [];

for j = 1:40
    sol = run_single_SSA_linda_inhib(x0,S,W,[0:1:1200],time_var,signal_update_rate,parameters,inhibs,1110);
    
        
    [pol2_ssa,ser5_ssa,ts_ssa,~] = get_model_intesities(sol,eta_rnap,eta_ser5,eta_ts); %convert molecules to signal
    
%     pol2_ssa = sol(2,:)';
%     ser5_ssa = sol(2,:)';
%     ts_ssa = sol(3,:)';
    
    
%  pol2_ssa = sol(2,:)' + sol(3,:)';
% ser5_ssa = pol2_ssa;
% ts_ssa = sol(3,:)';



%pol2_ssa = add_shot_std(pol2_ssa,eta_rnap,30/200, std(pol2_ssa(1:1100)));
%ser5_ssa = add_shot_std(ser5_ssa,eta_ser5,23/200, std(ser5_ssa(1:1100)));
%ts_ssa = add_shot_std(ts_ssa,eta_ts,5/200, std(ts_ssa(1:1100))  );

   % max_pol2 = mean(pol2_ssa(1110-2:1110));
    %max_ser5 = mean(ser5_ssa(1110-2:1110));
    %max_ts = mean(ts_ssa(1110-2:1110));
    
    %min_pol2 = mean(pol2_ssa(end-2:end));
    %min_ser5 = mean(ser5_ssa(end-2:end));
    %min_ts = mean(ts_ssa(end-2:end));
    
    
    %pol2_ssa2 = (pol2_ssa - min_pol2)./(max_pol2-min_pol2);
    %ser5_ssa2 = (ser5_ssa - min_ser5)./(max_ser5-min_ser5);
    %ts_ssa2 = (ts_ssa - min_ts)./(max_ts-min_ts);
    
%     pol2_ssa2 = min_pol2+(max_pol2-min_pol2)*(pol2_ssa-min(pol2_ssa))/( max(pol2_ssa)-min(pol2_ssa));
%     ser5_ssa2 = min_ser5+(max_ser5-min_ser5)*(ser5_ssa-min(ser5_ssa))/( max(ser5_ssa)-min(ser5_ssa));
%     ts_ssa2 = min_ts+(max_ts-min_ts)*(ts_ssa-min(ts_ssa))/( max(ts_ssa)-min(ts_ssa));
    top_pol2 = quantile(pol2_ssa(1:1*1100),.95); 
    top_ser5 = quantile(ser5_ssa(1:1*1100),.95); 
    top_ts = quantile(ts_ssa(1:1*1100),.95); 

    pol2_ssa = pol2_ssa./top_pol2;
    ser5_ssa = ser5_ssa./top_ser5;
    ts_ssa = ts_ssa./top_ts;
    
     pol2_ssa = min(pol2_ssa,1.5);
     ser5_ssa = min(ser5_ssa,1.5);
     ts_ssa = min(ts_ssa,1.5);

    pol2_traj = [pol2_traj; pol2_ssa'];
    ser5_traj = [ser5_traj; ser5_ssa'];
    ts_traj = [ts_traj; ts_ssa'];
    


end
avpol2 = mean(pol2_traj,1);
avser5 = mean(ser5_traj,1);
avts = mean(ts_traj,1);

bef_pol2 = mean((avpol2(1*1100:1110)));
bef_ser5 = mean((avser5(1*1100:1110)));
bef_ts = mean((avts(1*1100:1110)));
avpol2 = avpol2./mean((avpol2(1*1100:1110)));
avser5 = avser5./mean((avser5(1*1100:1110)));
avts = avts./mean((avts(1*1100:1110)));
per_len = 40;
t = [-10:1:per_len-10];


figure
hold on;
avpol2 = mean(pol2_traj,1);
avser5 = mean(ser5_traj,1);
avts = mean(ts_traj,1);

bef_pol2 = mean((avpol2(1*1100-2:1110)));
bef_ser5 = mean((avser5(1*1100-2:1110)));
bef_ts = mean((avts(1*1100-2:1110)));

after_pol2 = mean((avpol2(end-2:end))); 
after_ser5 = mean((avser5(end-2:end))); 
after_ts = mean((avts(end-2:end))); 

avpol2 = (avpol2- after_pol2)./(bef_pol2-after_pol2);
avser5 = (avser5- after_ser5)./(bef_ser5-after_ser5);
avts = (avts- after_ts)./(bef_ts-after_ts);
%plot(t, avpol2(1*1100:1*1100+1*per_len),'r','LineWidth',1);  plot(t, avser5(1*1100:1*1100+1*per_len),'g','LineWidth',1); plot(t, avts(1*1100:1*1100+1*per_len),'b','LineWidth',1)

plot(t, avpol2(1*1100:1*1100+1*per_len),'Color','r','LineWidth',1);  plot(t, avser5(1*1100:1*1100+1*per_len),'Color','g','LineWidth',1); plot(t, avts(1*1100:1*1100+1*per_len),'b','LineWidth',1)
plot([0,0],[-.1,1.1],'k--','LineWidth',2);

%%
p2_x = pol2_traj(:,1107:1140);
s5_x =  ser5_traj(:,1107:1140);
m_x = mean(ts_traj(:,1107:1140),1);
xx = 0:1:33;
modelFun = @(p,x) atan(x.*p(1)./p(2)./(p(1)^2-x.^2));
startingVals = [1 1]; % try with any other starting values. 
coefEsts = nlinfit(xx,m_x, modelFun, startingVals);

figure
plot(m_x);
hold on;
plot(modelFun(coefEsts,xx));


%% Inhibs 
n = 0;
figure(9)
titles = {'Blocking k_{on}','Blocking \beta','Blocking k_{escape}','Blocking k_{esc} and 30% k_{proc}','Blocking k_{on} and 30% k_{proc} ','Blocking \beta and 30% k_{proc}'};
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height
for i = 1:6
subplot(2,3,i)

inhibs = ones(size(parameters));
n = n+1;
switch n
    case 1  %inhibit kon
        inhibs(1) = .01;
    case 2  %inhibit beta
        
        inhibs(6) = 1000;
    case 3 %inhibit kescape
        inhibs(3) = .01;
    case 4 %reduce esc and proc
        inhibs(4) = .3;
        inhibs(3) = .01;
        
    case 5%reduce on and kproc
        inhibs(1) = .01;
        inhibs(4) = .3;
        
    case 6 %reduce kproc and beta
        inhibs(4) = .3;
        inhibs(6) = 1000;
end
inhibs;
pol2_traj = [];
ser5_traj = [];
ts_traj = [];

pol2_traji = [];
ser5_traji = [];
ts_traji = [];

for j = 1:200
    sol = run_single_SSA_linda_inhib(x0,S,W,[0:1:1200],time_var,signal_update_rate,parameters,inhibs,1110);
    
        
    
    pol2_ssa = sol(2,:)';
    ser5_ssa = sol(2,:)';
    ts_ssa = sol(3,:)';
    
    
 pol2_ssa = sol(2,:)' + sol(3,:)';
ser5_ssa = pol2_ssa;
ts_ssa = sol(3,:)';

pol2_ssa = add_shot_std(pol2_ssa,eta_rnap,30/200, std(pol2_ssa(1:1100)));
ser5_ssa = add_shot_std(ser5_ssa,eta_ser5,23/200, std(ser5_ssa(1:1100)));
ts_ssa = add_shot_std(ts_ssa,eta_ts,5/200, std(ts_ssa(1:1100))  );



    %[pol2_ssa,ser5_ssa,ts_ssa,~] = get_model_intesities(sol,eta_rnap,eta_ser5,eta_ts); 
    
    top_pol2 = quantile(pol2_ssa(1:1*1100),.95); 
    top_ser5 = quantile(ser5_ssa(1:1*1100),.95); 
    top_ts = quantile(ts_ssa(1:1*1100),.95); 

    pol2_ssa = pol2_ssa./top_pol2;
    ser5_ssa = ser5_ssa./top_ser5;
    ts_ssa = ts_ssa./top_ts;
    
     pol2_ssa = min(pol2_ssa,1.5);
     ser5_ssa = min(ser5_ssa,1.5);
     ts_ssa = min(ts_ssa,1.5);

    pol2_traj = [pol2_traj; pol2_ssa'];
    ser5_traj = [ser5_traj; ser5_ssa'];
    ts_traj = [ts_traj; ts_ssa'];
    


end
avpol2 = mean(pol2_traj,1);
avser5 = mean(ser5_traj,1);
avts = mean(ts_traj,1);

bef_pol2 = mean((avpol2(1*1100:1110)));
bef_ser5 = mean((avser5(1*1100:1110)));
bef_ts = mean((avts(1*1100:1110)));
avpol2 = avpol2./mean((avpol2(1*1100:1110)));
avser5 = avser5./mean((avser5(1*1100:1110)));
avts = avts./mean((avts(1*1100:1110)));
per_len = 40;
t = [-10:1:per_len-10];


plot(t, pol2_traj(1,1*1100:1*1100+1*per_len)./bef_pol2,'r','LineWidth',1); hold on; plot(t, ser5_traj(1,1*1100:1*1100+1*per_len)./bef_ser5,'g','LineWidth',1); plot(t, ts_traj(1,1*1100:1*1100+1*per_len)./bef_ts,'b','LineWidth',1)
if i == 1
    %legend({'CTD','Ser5ph','mRNA'})
end
plot(t, pol2_traj(2:5,1*1100:1*1100+1*per_len)./bef_pol2,'r','LineWidth',1); plot(t, ser5_traj(2:5,1*1100:1*1100+1*per_len)./bef_ser5,'g','LineWidth',1); plot(t, ts_traj(2:5,1*1100:1*1100+1*per_len)./bef_ts,'b','LineWidth',1)

plot(t, avpol2(1*1100:1*1100+1*per_len),'r','LineWidth',3);  plot(t, avser5(1*1100:1*1100+1*per_len),'g','LineWidth',3); plot(t, avts(1*1100:1*1100+1*per_len),'b','LineWidth',3)
plot([0,0],[-2,2],'k--','LineWidth',2);

ylim([-1,2])
title(titles{i},'FontSize',fntsize,'FontWeight','bold')
if i == 5
    xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
end
if i == 6
    xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
end

if i == 4
    xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
end

if ismember(i,[1,4])
    ylabel('Normalized Intensity','FontSize',fntsize,'FontWeight','bold')
end




end
saveas(gca,'perturb.epsc')


return 

%%
n = 0;
figure(15)
titles = {'Blocking k_{on}','Blocking k_{escape}','Blocking k_{esc} and 30% k_{proc}','Blocking k_{on} and 30% k_{proc} '};
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

k = 0;
per_len = 40;
for i = 1:4


inhibs = ones(size(parameters));
k = k+1;
switch k
    case 1  %inhibit kon
        inhibs(1) = .000001;

    case 2 %inhibit kescape
        inhibs(3) = .01;
    case 3 %reduce esc and proc
        inhibs(4) = .3;
        inhibs(3) = .01;
        
    case 4%reduce on and kproc
        inhibs(1) = .01;
        inhibs(4) = .3;
        

end
inhibs;
pol2_traj = [];
ser5_traj = [];
ts_traj = [];

pol2_traji = [];
ser5_traji = [];
ts_traji = [];
n = 0;
while n < 7
    
    sol = run_single_SSA_linda_inhib(x0,S,W,[0:1:1200],time_var,signal_update_rate,parameters,inhibs,1110);
    sol(2,1110);
    if sol(3,1110) > means(3)
        n = n+1;
        pol2_ssa = sol(2,:)' + sol(3,:)';
        ser5_ssa = pol2_ssa;
        ts_ssa = sol(3,:)';        

        pol2_ssa = sol(2,:)' + sol(3,:)';
        ser5_ssa = pol2_ssa;
        ts_ssa = sol(3,:)';       
        
        
        pol2_traj = [pol2_traj; pol2_ssa'];
        ser5_traj = [ser5_traj; ser5_ssa'];
        ts_traj = [ts_traj; ts_ssa'];

    end
    



end


newpars = parameters;
newpars = newpars.*inhibs;

inhibs 
kon = newpars(1);
koff = 1000;%parameters(2);
kesc = newpars(3);
kproc = newpars(4);
kin =  1000*newpars(5);
kout = newpars(6);
frac = newpars(7);
eta_rnap = newpars(8);
eta_ser5 = newpars(9);
eta_ts = newpars(10);

Nstates = 3;
b = zeros(Nstates,1);
b(1) = kon;
b(2) = 0;
b(3) = 0;

c = zeros(3,3);
c(1,1:3)=[0,1,1];
c(2,1:3)=[0,frac,1];
c(3,3)=1;

S = zeros(3,6);
W1 = zeros(6,3);
W0 = zeros(6,1);
S(1,1) = 1;  W1(1,1) = -kon; W0(1,1) = kon;
S(1,2) = -1; W1(2,1) = koff; 

S(2,3) = 1; W1(3,1) = kin; 
S(2,4) = -1; W1(4,2) = kout; 

S(2:3,5) = [-1;1]; W1(5,2) = kesc*frac; 
S(3,6) = -1; W1(6,3) = kproc; 

A=S*W1;
tode = 0:.1:per_len-10;
on = [];
ctd_ode = [];
ts_ode = [];
for j = 1:length(tode)
    
    %P = expm(A*tode(j)+ W0)*b;
%      if j == 1
%          b(2) = means(2);
%          b(3) = means(3);
%      else
%          b(2) = 0;
%          b(3) = 0;         
%      end

     EX = -A\b;
     means_new = c*EX;
     %P =  A^-1*(-eye(3) + expm(A*tode(j)))*b - x0;
     P= A^-1*(-eye(3) + expm(A*tode(j)))*b  + expm(A*tode(j))*[0;4.5649;15.2925];
    on = [on, P(1)];
    ctd_ode = [ctd_ode,P(2)];
    ts_ode = [ts_ode,P(3)];
end

t = [-10:1:per_len-10];

subplot(3,4,i)


plot(t, pol2_traj(1:n,1*1100:1*1100+1*per_len),'r','LineWidth',1); 
hold on;
plot(tode,ctd_ode+ts_ode,'k-','LineWidth',2)
% 
% plot(t, avpol2(1*1100:1*1100+1*per_len),'r','LineWidth',3); 
% 
% fill([t, fliplr(t)], [avpol2(1*1100:1*1100+1*per_len)+stdpol2(1*1100:1*1100+1*per_len),fliplr(avpol2(1*1100:1*1100+1*per_len)-stdpol2(1*1100:1*1100+1*per_len))],'r', 'FaceAlpha',.5)
% 
% plot(t, avpol2(1*1100:1*1100+1*per_len)-stdpol2(1*1100:1*1100+1*per_len),'r','LineWidth',1); 
% plot(t, avpol2(1*1100:1*1100+1*per_len)+stdpol2(1*1100:1*1100+1*per_len),'r','LineWidth',1); 
% 
% plot(t, movmean(pol2_traj(1,1*1100:1*1100+1*per_len)./bef_pol2,3),'k','LineWidth',2); 

plot([0,0],[-10,150],'k--','LineWidth',1);
plot([-10,per_len-10],[0,0],'k--','LineWidth',1);
plot([-5,0],[19.8575,19.8575],'k-','LineWidth',2)
ylim([-10,110])
xlim([-5,per_len-10])


%ylim([-1,2])
title(titles{i},'FontSize',fntsize,'FontWeight','bold')


if i == 1
    
   ylabel('CTD Count','FontSize',fntsize,'FontWeight','bold') 
end

subplot(3,4,i+4)


 plot(t,ser5_traj(1:n,1*1100:1*1100+1*per_len),'g','LineWidth',1); 
 hold on;
 plot(tode,ctd_ode+ts_ode,'k-','LineWidth',2)
 
%  plot(t, avser5(1*1100:1*1100+1*per_len),'g','LineWidth',3);
%   plot(t, movmean(ser5_traj(1,1*1100:1*1100+1*per_len)./bef_ser5,3),'k','LineWidth',2);
plot([0,0],[-10,150],'k--','LineWidth',1);

ylim([-10,110])
xlim([-5,per_len-10])
plot([-10,per_len-10],[0,0],'k--','LineWidth',1);
%ylim([-1,2])
xlim([-5,per_len-10])
plot([-5,0],[19.8575,19.8575],'k-','LineWidth',2)

if i == 1
    
   ylabel('Ser5ph Count','FontSize',fntsize,'FontWeight','bold') 
end

subplot(3,4,i+8)

 plot(t, ts_traj(1:n,1*1100:1*1100+1*per_len),'b','LineWidth',1)
 hold on;
 plot(tode,ts_ode,'k-','LineWidth',2)
 plot([-5,0],[15.2925,15.2925],'k-','LineWidth',2)
%plot(t, avts(1*1100:1*1100+1*per_len),'b','LineWidth',3)
% 
% plot(t, movmean(ts_traj(1,1*1100:1*1100+1*per_len)./bef_ts,3),'k','LineWidth',2)
plot([0,0],[-10,110],'k--','LineWidth',1);
plot([-10,per_len-10],[0,0],'k--','LineWidth',1);
xlim([-5,per_len-10])
ylim([-10,75])
if i == 1
    
   ylabel('mRNA Count','FontSize',fntsize,'FontWeight','bold') 
end

% ylim([-1,2])
% xlim([-10,per_len-10])

 xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
 


if i == 5
    xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
end
if i == 6
    xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
end







end
saveas(gca,'perturbs_sep1.epsc')




%% Simulated chip with SSA

kon = parameters(1);
koff = 1000;%parameters(2);
kesc = parameters(3);
kproc = parameters(4);
kin =  1000*parameters(5);
kout = parameters(6);
frac = parameters(7);
eta_rnap = parameters(8);
eta_ser5 = parameters(9);
eta_ts = parameters(10);
npts_pol2 = [];
npts_ts = [];
figure(20)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
xh = 2*3*6;
yh = 2*14.3/3;
fig1.PaperPosition = [0, 0, 1*xh, yh]; % x,y, width, height

%while length(npts) < 2000
   
pol2_on = [];
pol2_off = [];
pol2_transient = [];

for j = 1:3
    T_array = [0:1:40000];
    sol = run_single_SSA_linda(x0,S,W,T_array,time_var,signal_update_rate);

    pol2_ssa = sol(2,:)';
    ser5_ssa = sol(2,:)+ sol(3,:)';
    ts_ssa = sol(3,:)';
    
    
    
    decorr_ts = ts_ssa(400:40:40000);
    decorr_pol2 = pol2_ssa(400:40:40000);
    
    npts_pol2 = [npts_pol2, decorr_pol2'];
    npts_ts = [npts_ts, decorr_ts'];
    
end





X = [mean(npts_pol2),mean(npts_ts)];


time_elongating = 5.2/4.1;     
residence = 1/kproc;
time_processing = residence - time_elongating;

frac_proc = time_processing/residence;
processing = mean(npts_ts)*frac_proc;
elongating = mean(npts_ts)*(1-frac_proc)/6;
mids = ones(1,6)*elongating;
X = [mean(npts_pol2),mids,processing];


bar(X)
%set(gca,'xticklabel',{'1})
ylabel('Molecule counts')
hold on;

text(1-.2,6,string(round(X(1),1)),'FontSize',fntsize-4,'FontWeight','bold') 
text(8-.2,13,string(round(X(8),1)),'FontSize',fntsize-4,'FontWeight','bold') 
n = length(npts_pol2);
title(strcat('Average All nsim = ',string(n)),'FontSize',fntsize,'FontWeight','bold')
xlabel({'Location' },'FontSize',fntsize,'FontWeight','bold')
ylabel({'Simulated molecule counts'},'FontSize',fntsize,'FontWeight','bold')
ylim([0,35])

saveas(gca,'ss.epsc')


%end

%%
npts_pol2 = [];
npts_ts = [];
figure(21)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
xh = 2*3*6;
yh = 2*14.3/3;
fig1.PaperPosition = [0, 0, 1*xh, yh]; % x,y, width, height

%while length(npts) < 2000
while length(npts_ts) < 500
    length(npts_ts) 
    T_array = [0:1:40000];
    sol = run_single_SSA_linda(x0,S,W,T_array,time_var,signal_update_rate);

    pol2_ssa = sol(2,:)';
    ser5_ssa = sol(2,:)+ sol(3,:)';
    ts_ssa = sol(3,:)';
    
    decorr_ts = ts_ssa(400:40:40000);
    decorr_pol2 = pol2_ssa(400:40:40000);
    
    for i = 1:length(decorr_ts)
        if decorr_ts(i)+decorr_pol2(i) > 50
            npts_pol2 = [npts_pol2, decorr_pol2(i)];
            npts_ts = [npts_ts, decorr_ts(i)];
        end
        
    end
    
end
    

time_elongating = 5.2/4.1;     
residence = 1/kproc;
time_processing = residence - time_elongating;

frac_proc = time_processing/residence;
processing = mean(npts_ts)*frac_proc;
elongating = mean(npts_ts)*(1-frac_proc)/6;
mids = ones(1,6)*elongating;
X2 = [mean(npts_pol2),mids,processing];

X3 = [X,X2];
bar(X2)
%set(gca,'xticklabel',{'1'})
ylabel('Molecule counts')
hold on;

text(1-.2,30,string(round(X2(1),1)),'FontSize',fntsize-4,'FontWeight','bold') 
text(8-.2,30,string(round(X2(8),1)),'FontSize',fntsize-4,'FontWeight','bold') 
n = length(npts_pol2);
title(strcat('Average On (Total CTD > 50) nsim = ',string(n)),'FontSize',fntsize,'FontWeight','bold')
xlabel({'Location' },'FontSize',fntsize,'FontWeight','bold')
ylabel({'Simulated molecule counts'},'FontSize',fntsize,'FontWeight','bold')
ylim([0,35])
saveas(gca,'bursting_chip.epsc')

%%

npts_pol2 = [];
npts_ts = [];
figure(22)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
xh = 2*3*6;
yh = 2*14.3/3;
fig1.PaperPosition = [0, 0, 1*xh, yh]; % x,y, width, height

%while length(npts) < 2000
    
for j = 1:3
    T_array = [0:1:40000];
    sol = run_single_SSA_linda(x0,S,W,T_array,time_var,signal_update_rate);

    pol2_ssa = sol(2,:)';
    ser5_ssa = sol(2,:)+ sol(3,:)';
    ts_ssa = sol(3,:)';
    
    decorr_ts = ts_ssa(400:40:40000);
    decorr_pol2 = pol2_ssa(400:40:40000);
    
    for i = 1:length(decorr_ts)
        if decorr_pol2(i) < 5
            npts_pol2 = [npts_pol2, decorr_pol2(i)];
            npts_ts = [npts_ts, decorr_ts(i)];
        end
        
    end
end

time_elongating = 5.2/4.1;     
residence = 1/kproc;
time_processing = residence - time_elongating;

frac_proc = time_processing/residence;
processing = mean(npts_ts)*frac_proc;
elongating = mean(npts_ts)*(1-frac_proc)/6;
mids = ones(1,6)*elongating;
X2 = [mean(npts_pol2),mids,processing];

X3 = [X,X2];
bar(X2)
%set(gca,'xticklabel',{'1'})
ylabel('Molecule counts')
hold on;

text(1-.2,2,string(round(X2(1),1)),'FontSize',fntsize-4,'FontWeight','bold') 
text(8-.2,13,string(round(X2(8),1)),'FontSize',fntsize-4,'FontWeight','bold') 

n = length(npts_ts)

title(strcat('Average Off (Unescaped CTD < 5) nsim = ',string(n)),'FontSize',fntsize,'FontWeight','bold')
xlabel({'Location' },'FontSize',fntsize,'FontWeight','bold')
ylabel({'Simulated molecule counts'},'FontSize',fntsize,'FontWeight','bold')
ylim([0,35])
saveas(gca,'off_chip.epsc')

%%
npts_pol2 = [];
npts_ts = [];
figure(23)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
xh = 2*3*6;
yh = 2*14.3/3;
fig1.PaperPosition = [0, 0, 1*xh, yh]; % x,y, width, height

%while length(npts) < 2000
    
for j = 1:3
    T_array = [0:1:40000];
    sol = run_single_SSA_linda(x0,S,W,T_array,time_var,signal_update_rate);

    pol2_ssa = sol(2,:)';
    ser5_ssa = sol(2,:)+ sol(3,:)';
    ts_ssa = sol(3,:)';
    
    decorr_ts = ts_ssa(400:40:40000);
    decorr_pol2 = pol2_ssa(400:40:40000);
    
    for i = 1:length(decorr_ts)
        if decorr_pol2(i) >= 5 && decorr_pol2(i) < 15
            npts_pol2 = [npts_pol2, decorr_pol2(i)];
            npts_ts = [npts_ts, decorr_ts(i)];
        end
        
    end
end

time_elongating = 5.2/4.1;     
residence = 1/kproc;
time_processing = residence - time_elongating;

frac_proc = time_processing/residence;
processing = mean(npts_ts)*frac_proc;
elongating = mean(npts_ts)*(1-frac_proc)/6;
mids = ones(1,6)*elongating;
X2 = [mean(npts_pol2),mids,processing];

X3 = [X,X2];
bar(X2)
%set(gca,'xticklabel',{'1'})
ylabel('Molecule counts')
hold on;

text(1-.2,10,string(round(X2(1),1)),'FontSize',fntsize-4,'FontWeight','bold') 
text(8-.2,17,string(round(X2(8),1)),'FontSize',fntsize-4,'FontWeight','bold') 

n = length(npts_ts)

title(strcat('Average Transient (5 < Unescaped CTD < 15) nsim = ',string(n)),'FontSize',fntsize,'FontWeight','bold')
xlabel({'Location' },'FontSize',fntsize,'FontWeight','bold')
ylabel({'Simulated molecule counts'},'FontSize',fntsize,'FontWeight','bold')
ylim([0,35])
saveas(gca,'trans_chip.epsc')


%% Calculate analytical parameters


Nstates = 3;
b = zeros(Nstates,1);
b(1) = kon;
b(2) = 0;
b(3) = 0;per_len = 40;

c = zeros(3,3);
c(1,1:3)=[0,1,1];
c(2,1:3)=[0,frac,1];
c(3,3)=1;

S = zeros(3,6);
W1 = zeros(6,3);
W0 = zeros(6,1);
S(1,1) = 1;  W1(1,1) = -kon; W0(1,1) = kon;
S(1,2) = -1; W1(2,1) = koff; 

S(2,3) = 1; W1(3,1) = kin; 
S(2,4) = -1; W1(4,2) = kout; 

S(2:3,5) = [-1;1]; W1(5,2) = kesc*frac; 
S(3,6) = -1; W1(6,3) = kproc; 

A=S*W1;
tode = 0:.01:per_len-10;

W = @(x) W1*x + W0;


burst_sizes = [];
burst_freqs = [];
p2_a = [];
ab = [];
es = [];
x0 = [0,0,0]';
nbursts = [];
burst_sizes = [];
for j = 1:5
    T_array = [0:1:40000];
    [sol,pol2_arrivals,aborted,escaped,p2_arrive_times] = run_single_SSA_linda_recorder(x0,S,W,T_array,0,0);
    p2_a = [p2_a,pol2_arrivals];
    ab = [ab,aborted];
    es = [es,escaped];
    pol2_ssa = sol(2,:)';
    ser5_ssa = sol(2,:)+ sol(3,:)';
    ts_ssa = sol(3,:)';
    
    diffs = p2_arrive_times(2:end) - p2_arrive_times(1:end-1); 

    nbursts = [nbursts, sum(diffs > .01)./40000    ];
    decorr_ts = ts_ssa(400:40:40000);
    decorr_pol2 = pol2_ssa(400:40:40000);

    inds = (diffs > .001);
 
    cuminds = cumsum(inds);
    for i = 1:sum(inds)
        burst_sizes = [burst_sizes,sum((cuminds == i))];

    end    
    
    for i = 1:length(decorr_ts)
        if decorr_pol2(i) > 5 
            burst_sizes = [burst_sizes, decorr_pol2(i)];
            burst_freqs = [burst_freqs, decorr_ts(i)];
        end
    end
        
end

arrival_rate = mean(p2_a)./40000;
prob_abort = mean(ab)/ (mean(ab) + mean(es));


%%
Nstates = 3;
b = zeros(Nstates,1);
b(1) = kon;
b(2) = 0;
b(3) = 0;per_len = 40;

c = zeros(3,3);
c(1,1:3)=[0,1,1];
c(2,1:3)=[0,frac,1];
c(3,3)=1;

S = zeros(3,6);
W1 = zeros(6,3);
W0 = zeros(6,1);
S(1,1) = 1;  W1(1,1) = -kon; W0(1,1) = kon;
S(1,2) = -1; W1(2,1) = koff; 

S(2,3) = 1; W1(3,1) = kin; 
S(2,4) = -1; W1(4,2) = kout; 

S(2:3,5) = [-1;1]; W1(5,2) = kesc*frac; 
S(3,6) = -1; W1(6,3) = kproc; 

A=S*W1;
tode = 0:.01:per_len-10;
W0 = zeros(6,1);
W = @(x) W1*x + W0;

x0 = [0,1,0]'
pol2_leaves = [];
mrna_times = [];
for j = 1:10000
    T_array = [0:.01:50];
    [sol] = run_single_SSA_linda_one_mol(x0,S,W,T_array,time_var,signal_update_rate);
   
    t_left = sum(sol(2,:)~= 0)*.01;
    t_transcribed = sum(sol(3,:)~= 0)*.01;
    pol2_leaves = [pol2_leaves, t_left];
    if t_transcribed ~= 0
        mrna_times = [mrna_times,t_transcribed];
    end
        
end

mean(pol2_leaves)
mean(mrna_times)

%%

m_cluster = (15.405*parameters(1) / (kesc + kout));


disp('sim arrive, ana arrive')
[arrival_rate, 15.4*.43]
disp('sim pesc, ana pesc')
[1-prob_abort, kesc / (kesc + kout)]
disp('sim av mRNA, ana av mRNA')
[mean(ts_ssa),  m_cluster*(kesc/kproc)]
disp('sim av Pol2, ana av Pol2')
[mean(ts_ssa+pol2_ssa), m_cluster + m_cluster*(kesc/kproc)]

disp('sim av Pol2 leave, ana av Pol2 leave')
[mean(pol2_leaves),1/(kesc+ kout)]

disp('sim av  mrna completion, ana av mrna completion')
[mean(mrna_times),1/(kproc)]

disp('sim av mrna production, sim av mrna production')
[mean(es)./40000, m_cluster*kesc]

disp('sim av rnap in cluster, sim av rnap in cluster')
[mean(pol2_ssa), m_cluster]



%% CI's



addpath ../Data_files/
addpath ../

parnames = {'kon','kesc','kproc','beta','kout'};
parnames = {'beta','omega','k out','k esc','k complete'};

par_changed = [1:5];

mh_pars = [];
mh_vals = [];
ikeep = 1;
i=1;
while ikeep==1
    try
        fn = ['met_hast_pars_2x_',num2str(i),'.mat'];
        load(fn)
        mh_pars = [mh_pars;mh_smpl];
        mh_vals = [mh_vals;mh_value];
        i=i+1;
    catch
        ikeep=0;
    end    
end

Np = 5;

mh_pars = mh_pars(:,[4,1,5,2,3]);

beta = mh_pars(:,1);
omega = mh_pars(:,2);
kout = mh_pars(:,3);
kesc = mh_pars(:,4);
kcomp = mh_pars(:,5);




r = beta.*omega;

r =  sort(r);
low = r(floor(length(r)/10));
high = r(ceil(9*length(r)/10));
disp([low, high])

t_clust = 1./(kesc + kout);

t_clust =  sort(t_clust);
low = t_clust(floor(length(t_clust)/10));
high = t_clust(ceil(9*length(t_clust)/10));
disp([low, high])


u_clust = (beta.*omega)./(kesc + kout);

u_clust =  sort(u_clust);
low = u_clust(floor(length(u_clust)/10));
high = u_clust(ceil(9*length(u_clust)/10));
disp([low, high])


f = kesc./(kesc + kout);

f =  sort(f);
low = f(floor(length(f)/10));
high = f(ceil(9*length(f)/10));
disp([low, high])


bmrna = kesc./(kesc + kout).*beta;

bmrna =  sort(bmrna);
low = bmrna(floor(length(bmrna)/10));
high = bmrna(ceil(9*length(bmrna)/10));
disp([low, high])


rmrna = ((beta.*omega)./(kesc + kout)).*kesc;

rmrna =  sort(rmrna);
low = rmrna(floor(length(rmrna)/10));
high = rmrna(ceil(9*length(rmrna)/10));
disp([low, high])

rmrna = ((beta.*omega)./(kesc + kout)).*(kesc ./kcomp)  ;

rmrna =  sort(rmrna);
low = rmrna(floor(length(rmrna)/10));
high = rmrna(ceil(9*length(rmrna)/10));
disp([low, high])


totalrnap = ((beta.*omega)./(kesc + kout)).*(kesc ./kcomp) + ((beta.*omega)./(kesc + kout))  ;

totalrnap =  sort(totalrnap);
low = totalrnap(floor(length(totalrnap)/10));
high = totalrnap(ceil(9*length(totalrnap)/10));
disp([low, high])


mrnat = 1./kcomp;
mrnat =  sort(mrnat);
low = mrnat(floor(length(mrnat)/10));
high = mrnat(ceil(9*length(mrnat)/10));
disp([low, high])



%% Bursting analytical parameters

Nstates = 4;
b = zeros(Nstates,1);
b(1) = 0;
b(2) = 0;
b(3) = 0;per_len = 40;

c = zeros(3,3);
c(1,1:3)=[0,1,1];
c(2,1:3)=[0,frac,1];
c(3,3)=1;

S = zeros(4,6);
W1 = zeros(6,4);
W0 = zeros(6,1);
S(2,1) = 1;  W1(1,1) = kon; W0(1,1) = 0;
S(2,2) = -1; W1(2,2) = koff;  
S(1,1) = -1; S(1,2) = 1;

S(3,3) = 1; W1(3,2) = kin;
S(3,4) = -1; W1(4,3) = kout; 

S(3:4,5) = [-1;1]; W1(5,3) = kesc*frac; 
S(4,6) = -1; W1(6,4) = kproc; 

A=S*W1;
tode = 0:.01:per_len-10;
W0 = zeros(6,1);
W = @(x) W1*x + W0;

T_array = [0:1:40000];
[sol,bursts,p2_per_burst,burst_t,mrna_cnts] = run_single_SSA_linda_recorder_bursting([1,0,0,0]',S,W,T_array,0,0);
figure


plot(sol(2,:))
ylim([0,1.5])

figure

plot(sol(3,:))
hold on;
plot(sol(4,:))








%%

n = 0;
figure(16)
titles = {'Blocking k_{esc} and 30% k_{proc}','Blocking k_{on} and 30% k_{proc} ','Blocking \beta and 30% k_{proc}'};
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height
for i = 1:3


inhibs = ones(size(parameters));
n = n+1;
switch n

    case 1 %reduce esc and proc
        inhibs(4) = .3;
        inhibs(3) = .01;
        
    case 2%reduce on and kproc
        inhibs(1) = .01;
        inhibs(4) = .3;
        
    case 3 %reduce kproc and beta
        inhibs(4) = .3;
        inhibs(6) = 1000;
end
inhibs;
pol2_traj = [];
ser5_traj = [];
ts_traj = [];

pol2_traji = [];
ser5_traji = [];
ts_traji = [];

for j = 1:200
    sol = run_single_SSA_linda_inhib(x0,S,W,[0:1:1200],time_var,signal_update_rate,parameters,inhibs,1110);
    
        
    
    pol2_ssa = sol(2,:)';
    ser5_ssa = sol(2,:)';
    ts_ssa = sol(3,:)';
    
    
 pol2_ssa = sol(2,:)' + sol(3,:)';
ser5_ssa = pol2_ssa;
ts_ssa = sol(3,:)';

pol2_ssa = add_shot_std(pol2_ssa,eta_rnap,30/200, std(pol2_ssa(1:1100)));
ser5_ssa = add_shot_std(ser5_ssa,eta_ser5,23/200, std(ser5_ssa(1:1100)));
ts_ssa = add_shot_std(ts_ssa,eta_ts,5/200, std(ts_ssa(1:1100))  );



    %[pol2_ssa,ser5_ssa,ts_ssa,~] = get_model_intesities(sol,eta_rnap,eta_ser5,eta_ts); 
    
    top_pol2 = quantile(pol2_ssa(1:1*1100),.95); 
    top_ser5 = quantile(ser5_ssa(1:1*1100),.95); 
    top_ts = quantile(ts_ssa(1:1*1100),.95); 

    pol2_ssa = pol2_ssa./top_pol2;
    ser5_ssa = ser5_ssa./top_ser5;
    ts_ssa = ts_ssa./top_ts;
    
    pol2_ssa = min(pol2_ssa,1.5);
    ser5_ssa = min(ser5_ssa,1.5);
    ts_ssa = min(ts_ssa,1.5);

    pol2_traj = [pol2_traj; pol2_ssa'];
    ser5_traj = [ser5_traj; ser5_ssa'];
    ts_traj = [ts_traj; ts_ssa'];
    


end
avpol2 = mean(pol2_traj,1);
avser5 = mean(ser5_traj,1);
avts = mean(ts_traj,1);

bef_pol2 = mean((avpol2(1*1100:1110)));
bef_ser5 = mean((avser5(1*1100:1110)));
bef_ts = mean((avts(1*1100:1110)));
avpol2 = avpol2./mean((avpol2(1*1100:1110)));
avser5 = avser5./mean((avser5(1*1100:1110)));
avts = avts./mean((avts(1*1100:1110)));
per_len = 40;
t = [-10:1:per_len-10];

subplot(3,3,i)
plot(t, pol2_traj(1,1*1100:1*1100+1*per_len)./bef_pol2,'r','LineWidth',1); 
hold on;
plot(t, pol2_traj(2:5,1*1100:1*1100+1*per_len)./bef_pol2,'r','LineWidth',1); 
plot(t, avpol2(1*1100:1*1100+1*per_len),'r','LineWidth',3); 
plot([0,0],[-2,2],'k--','LineWidth',2);

ylim([-1,2])
title(titles{i},'FontSize',fntsize,'FontWeight','bold')


if i == 1
    
   ylabel('CTD Norm Signal','FontSize',fntsize,'FontWeight','bold') 
end

subplot(3,3,i+3)
plot(t, ser5_traj(1,1*1100:1*1100+1*per_len)./bef_ser5,'g','LineWidth',1);
hold on;
 plot(t, ser5_traj(2:5,1*1100:1*1100+1*per_len)./bef_ser5,'g','LineWidth',1); 
 plot(t, avser5(1*1100:1*1100+1*per_len),'g','LineWidth',3);
plot([0,0],[-2,2],'k--','LineWidth',2);
ylim([-1,2])

if i == 1
    
   ylabel('Ser5ph Norm Signal','FontSize',fntsize,'FontWeight','bold') 
end

subplot(3,3,i+6)
plot(t, ts_traj(1,1*1100:1*1100+1*per_len)./bef_ts,'b','LineWidth',1)
hold on;
 plot(t, ts_traj(2:5,1*1100:1*1100+1*per_len)./bef_ts,'b','LineWidth',1)
plot(t, avts(1*1100:1*1100+1*per_len),'b','LineWidth',3)
plot([0,0],[-2,2],'k--','LineWidth',2);

if i == 1
    
   ylabel('mRNA Norm Signal','FontSize',fntsize,'FontWeight','bold') 
end

ylim([-1,2])


 xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
 


if i == 5
    xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
end
if i == 6
    xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
end






end
saveas(gca,'perturbs_sep2.epsc')


%%

model_options.trypt = [0,0,0,0,0,0,0];
model_options.freepars
real_valscl3 = real_valscl;
real_valscl3(11) = real_valscl3(11)*15.8*.908;
real_valscl3(12) = real_valscl3(12)*15.8*.992;
real_valscl3(13) = real_valscl3(13)*15.8*.93;

kon = parameters(1);
koff = 1000;%parameters(2);
kesc = parameters(3);
kproc = parameters(4);
kin =  1000*parameters(5);
kout = parameters(6);
frac = parameters(7);
eta_rnap = parameters(8);
eta_ser5 = parameters(9);
eta_ts = parameters(10);

kelong = 4.1;
kdephos = 0;
kinit = 10;
kabort = 10;

rv = [kin, kout, kinit, kabort, kesc,kdephos, kelong, kon, koff, kproc, eta_rnap, eta_ser5, eta_ts];

[Ci,~,~,POL2,SER5,TS,POL2i,SER5i,TSi] = get_stoch_corrs_Nov25(rv,[],1,model_options, 40000 );

norm_signal = 0
real_valscl2 = real_valscl;
real_valscl2(11:13) = 0;
%[C_nn,~,~,POL2_nn,SER5_nn,TS_nn,POL2i_nn,SER5i_nn,TSi_nn] = get_stoch_corrs_Nov25(real_valscl2(1:end-4),[],1,model_options, 100*200  );

            i = 1;
            mod_rnap = Ci(:,(i-1)*3+i);
            mod_rnap = mod_rnap(((length(mod_rnap)+1)/2):((length(mod_rnap)+1)/2)+31);
            %mod_rnap = mod_rnap(2:16);

            
            i = 2;
            mod_ser5 = Ci(:,(i-1)*3+i);
            mod_ser5 = mod_ser5(((length(mod_ser5)+1)/2):((length(mod_ser5)+1)/2)+31);
            
            i = 3;
            mod_mrna = Ci(:,(i-1)*3+i);
            mod_mrna = mod_mrna(((length(mod_mrna)+1)/2):((length(mod_mrna)+1)/2)+31);
            
            [G1_rnap,G1_ser5,G1_mrna] = get_g0(mod_rnap,mod_ser5,mod_mrna,'G0_intp');
            mod_rnap = mod_rnap(1:30)/G1_rnap;
         
            mod_mrna = mod_mrna(1:30)/G1_mrna;

            mod_ser5 = mod_ser5(1:30)/G1_ser5;
            
            
            
            mod_mrna_rnap = Ci(:,7);
            mod_mrna_rnap = mod_mrna_rnap(((length(mod_mrna_rnap)+1)/2-10):((length(mod_mrna_rnap)+1)/2)+10);
            
            if norm_signal == 0
            mod_mrna_rnap = mod_mrna_rnap/(sqrt(G1_rnap)*sqrt(G1_mrna));
            end
            %mod_mrna_rnap = mod_mrna_rnap/mod_mrna_rnap(11);
            
            mod_ser5_rnap = Ci(:,4);
            mod_ser5_rnap = mod_ser5_rnap(((length(mod_ser5_rnap)+1)/2-10):((length(mod_ser5_rnap)+1)/2)+10);
            
            if norm_signal == 0
            mod_ser5_rnap = mod_ser5_rnap/(sqrt(G1_rnap)*sqrt(G1_ser5));
            end
            %mod_ser5_rnap = mod_ser5_rnap/mod_ser5_rnap(11);
            
            mod_mrna_ser5 = Ci(:,8);
            mod_mrna_ser5 = mod_mrna_ser5(((length(mod_mrna_ser5)+1)/2-10):((length(mod_mrna_ser5)+1)/2)+10);
            
            if norm_signal == 0
            mod_mrna_ser5 = mod_mrna_ser5/(sqrt(G1_mrna)*sqrt(G1_ser5));
            end%mod_mrna_ser5 = mod_mrna_ser5/mod_mrna_ser5(11);
            
            
figure
plot(mod_rnap); hold on; plot(mod_ser5); plot(mod_mrna);
figure
plot(mod_mrna_rnap);hold on; plot(mod_ser5_rnap); plot(mod_mrna_ser5);
            
%%          
% figure
% 
% plot([0:31],mod_rnap/G1s_rnap)
% hold on;
% b1=errorbar(lags1(1:20), rnap.mn_ac(1:20), rnap.sem_ac(1:20),'o','MarkerSize',5,'MarkerFaceColor',[138, 25, 38]./256,'Color',[138, 25, 38]./256);hold on;
%   
% 
% 
% 
% figure
% 
% plot([0:31],mod_ser5/G1s_ser5)
% hold on;
% b1=errorbar(lags1(1:20), ser5.mn_ac(1:20), ser5.sem_ac(1:20),'o','MarkerSize',5,'MarkerFaceColor',[138, 25, 38]./256,'Color',[138, 25, 38]./256);hold on;
%   
% 
% figure
% 
% plot([0:31],mod_mrna/G1s_mrna)
% hold on;
% b1=errorbar(lags1(1:20), mrna.mn_ac(1:20), mrna.sem_ac(1:20),'o','MarkerSize',5,'MarkerFaceColor',[138, 25, 38]./256,'Color',[138, 25, 38]./256);hold on;
%   


            
pol2_traj= [];
ser5_traj= [];
mrna_traj= [];

pol2i_traj= [];
ser5i_traj= [];
mrnai_traj= [];
for k = 1:100
    pol2_traj = [pol2_traj,POL2((k-1)*200+1: (k)*200)];
    ser5_traj = [ser5_traj,SER5((k-1)*200+1: (k)*200)];
    mrna_traj = [mrna_traj,TS((k-1)*200+1: (k)*200)];

    pol2i_traj = [pol2i_traj,POL2i((k-1)*200+1: (k)*200)];
    ser5i_traj = [ser5i_traj,SER5i((k-1)*200+1: (k)*200)];
    mrnai_traj = [mrnai_traj,TSi((k-1)*200+1: (k)*200)];

end




%%
[rnap_I, ser5_I, mrna_I] = Normalize_raw_intensities(.95);   %get the normalized intensities from the data 95 percentile
rnap_data = reshape(rnap_I(1:20:end,:),1,200);  %rnap_data sorted by 20 time incrementes (decorrelated) this was for histograms
ser5_data = reshape(ser5_I(1:20:end,:),1,200);
mrna_data = reshape(mrna_I(1:20:end,:),1,200);

xh = 6.3;
yh = 6.4;
fntsize = 12;


subplot(4,3,3)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, yh]; % x,y, width, height
[x,n] = hist(reshape(pol2_traj.',1,[]),30);
[11, 252, 3]./256
[218, 51, 255]./256
histogram(reshape(pol2_traj.',1,[]),n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256,'DisplayStyle','stairs')
hold on;
histogram(rnap_data ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')

a = legend(gca, {'Model','Data'},'Location','Best')
set(a, 'Box', 'off');

set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
set(a,'TextColor',global_color);

title({'CTD Normalized Intensity'},'FontSize',fntsize,'FontWeight','bold','Color',global_color)
xlabel({'Intensity (Norm)'},'FontSize',fntsize,'FontWeight','bold')
ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')

%ylim([minYlim maxYlim])
set (gca ,'TickLength',[.01,.3],'LineWidth',1);
set (gca ,'FontSize',fntsize,'FontName', 'Arial');

%set(gca,'Color','k')
%set(gcf,'Color','k')
saveas(gca,'CTD_dist.epsc') 

subplot(4,3,6)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, yh]; % x,y, width, height
[x,n] = hist(reshape(ser5_traj.',1,[]),30);
histogram(reshape(ser5_traj.',1,[]),n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256,'DisplayStyle','stairs')
hold on;
histogram(ser5_data ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')

a = legend(gca, {'Model','Data'},'Location','Best')
set(a, 'Box', 'off');

set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
set(a,'TextColor',global_color);
legend('Model','Data')
title({'SER5 Normalized Intensity'},'FontSize',fntsize,'FontWeight','bold','Color',global_color)
xlabel({'Intensity (Norm)'},'FontSize',fntsize,'FontWeight','bold')
ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
%ylim([minYlim maxYlim])
ylim([0,.13])
set (gca ,'TickLength',[.01,.3],'LineWidth',1);
set (gca ,'FontSize',fntsize,'FontName', 'Arial');
saveas(gca,'SER5_dist.epsc') 
%set(gca,'Color','k')
%set(gcf,'Color','k')

subplot(4,3,9)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, yh]; % x,y, width, height

[x,n] = hist(mrnai_traj,30)

histogram(mrnai_traj,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256,'DisplayStyle','stairs')
hold on;
histogram(mrna_hist ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')

a = legend(gca, {'Model','Data'},'Location','Best')
set(a, 'Box', 'off');

set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
set(a,'TextColor',global_color);

title({'mRNA count'},'FontSize',fntsize,'FontWeight','bold','Color',global_color)
xlabel({'N molecules'},'FontSize',fntsize,'FontWeight','bold')
ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
%ylim([minYlim maxYlim])
set (gca ,'TickLength',[.01,.3],'LineWidth',1);
set (gca ,'FontSize',fntsize,'FontName', 'Arial');
%set(gca,'Color','k')
%set(gcf,'Color','k')
saveas(gca,'TS_dist.epsc')   


subplot(4,3,[10,11,12])
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, yh]; % x,y, width, height
avpol2 = movmean(POL2,3);
avser5 = movmean(SER5,3);
avts = movmean(TS,3);
plot(avpol2(1:201),'r','linewidth',2); hold on; plot(avser5(1:201),'g','linewidth',2); plot(avts(1:201),'b','linewidth',2 )
legend('POL2','SER5ph','mRNA')
title('Representative trace','FontSize',fntsize,'FontWeight','bold')
ylabel('Normalized intensity','FontSize',fntsize,'FontWeight','bold')
xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
xlim([0,200])

saveas(gca,'ccs.epsc')


%%


chip_ctd = [];
chip_ser5 = [];
chip_ts = [];

for i = 1:80
    
   [~,~,~,~,~,~,~,~,~,tsvec,ser5vec,pol2vec] = get_stoch_corrs_Nov25_chip(real_valscl,[],1,model_options, 1000  );
    chip_ctd = [chip_ctd, pol2vec];
    chip_ser5 = [chip_ser5, ser5vec];
    chip_ts = [chip_ts, tsvec];
    
 
    
end
   figure(14) 
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, 3*xh, yh]; % x,y, width, height


subplot(3,1,1)
histogram(chip_ctd,8,'Normalization','Probability','FaceColor',[256, 0, 0]./256,'FaceAlpha',1)
axis([0 ,5.3 ,0,1])
title('CTD Signal location')
ylabel('Probability')
set (gca ,'FontSize',fntsize-4,'FontName', 'Arial');
subplot(3,1,2)
histogram(chip_ser5,8,'Normalization','Probability','FaceColor',[0, 256, 0]./256,'FaceAlpha',1)
title('ser5ph Signal location')

axis([0 ,5.3 ,0,1])
ylabel('Probability')
set (gca ,'FontSize',fntsize-4,'FontName', 'Arial');
subplot(3,1,3)
histogram(chip_ts,8,'Normalization','Probability','FaceColor',[0, 0, 256]./256,'FaceAlpha',1)
title('TS Signal location')
ylabel('Probability')
xlabel('kBP')
axis([0 ,5.3 ,0,1])
set (gca ,'FontSize',fntsize-4,'FontName', 'Arial');

saveas(gca,'simChip.epsc')      



function [thresh,std_mod] = compute_thresh(ssa,eta,tr)
std_mod = std(ssa);
shot_mod = std_mod*eta;
ssa_w_shot = ssa + randn(size(ssa))*shot_mod;
% Define fraction tr as zero
tmp = sort(ssa_w_shot);
thresh = tmp(ceil(length(tmp)*tr));

end



function ssa_w_shot = add_shot_std(ssa,eta,stdev,thresh)
std_mod = stdev;
shot_mod = std_mod*eta;
ssa_w_shot = ssa + randn(size(ssa))*shot_mod;

% Define fraction tr as zero

ssa_w_shot = ssa_w_shot-thresh;
end


