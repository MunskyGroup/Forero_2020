%%
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

addpath('./Data_files')
addpath('./Model_files')
addpath('./Parameter_files')

colors = [[239, 71, 111];
           [255, 209, 102];
           [6, 214, 160];
           [17, 138, 178];
           [86, 11, 173];
           [7, 59, 76];]./256 ;




%% DATA

%Load Cross Correlations with G0 normalization
[~,~,~,mrna_rnap,mrna_ser5,ser5_rnap] = load_normalization_variance_norm_cc(0);  %load normalization 21 pts for cross correlations
%Load Autocorrelations with G0_intp normalization and no rezeroing 
[mrna,ser5,rnap,~,~,~] = load_normalization_variance(1,'G0_intp','none',10);


%% Simple model

load('best_simple_pars.mat')  %load parameter file
real_valscl = parameters;
[sigred,TT,means,sig,TT2,mins] = get_ac_and_cc_mod_simplified(parameters,[0:.1:30]);


%%  Ctd --> Ser5 model

cd Model_fits\ctd_to_ser5_model
load('best_simple_pars_cs.mat')  %load parameter file
[sigred_cs,TT_cs,means_cs,sig_cs,TT2_cs,mins_cs] = get_ac_and_cc_mod_cs(parameters,[0:.1:30]);
cd ..\..

%% mRNA solo model

cd Model_fits\mrna_retention_model
load('best_simple_pars_mrna_solo.mat')  %load parameter file
[sigred_ms,TT_ms,means_ms,sig_ms,TT2_ms,mins_ms] = get_ac_and_cc_mod_mrnasolo(parameters,[0:.1:30]);
cd ..\..

%% fractional model
load('best_simple_pars.mat')  %load parameter file
real_valscl = parameters;
parameters(7) = 1;
[sigred_frac,TT_frac,means_frac,sig_frac,TT2_frac,mins_frac] = get_ac_and_cc_mod_simplified(parameters,[0:.1:30]);



%% 2 state

cd Model_fits\2state_model
load('best_simple_pars_2s.mat')  %load parameter file
[sigred_2s,TT_2s,means_2s,sig_2s,TT2_2s,mins_2s] = get_ac_and_cc_mod_2s(parameters,[0:.1:30]);
cd ..\..


%% 2 state unfixed
cd Model_fits\2state_model_unfixed
load('best_simple_pars_2s_freekoff.mat')  %load parameter file
[sigred_2sf,TT_2sf,means_2sf,sig_2sf,TT2_2sf,mins_2sf] = get_ac_and_cc_mod_2s(parameters,[0:.1:30]);
cd ..\..


%%


sizes = [8,6,4,2,2,2];
loc = [1,5,9];
plts = [1,3,5];
lags1 = [0:31]; 
titles = {'CTD','Ser5ph','mRNA'};
fntsize = 18;   %set up figure 3
close all
figure(10)

for i = 1:3
    
sigs = [sigred(loc(i),:); sigred_cs(loc(i),:); sigred_ms(loc(i),:); sigred_frac(loc(i),:); sigred_2s(loc(i),:); sigred_2sf(loc(i),:)  ];
subplot(3,2,plts(i))
fig1= gcf;
%set(gcf,'color','k');
fig1.PaperUnits = 'centimeters';
xh = 2*3*6;
yh = 2*14.3/3;
fig1.PaperPosition = [0, 0, 2/3*xh, 2*yh]; % x,y, width, height

for j = 1:6
    plot(TT, sigs(j,:),'-','linewidth',sizes(j),'Color',colors(j,:))  %plot model rnap_acc
    if j == 1
        hold on;
    end
end

if i == 1
b1=errorbar(lags1(1:25), rnap.mn_ac(1:25), rnap.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
end

if i == 2
b1=errorbar(lags1(1:25), ser5.mn_ac(1:25), ser5.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
    
end

if i == 3
b1=errorbar(lags1(1:25), mrna.mn_ac(1:25), mrna.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;    
end

if i == 2
legend('Chosen','Phosphorylation','mRNA retention','Fractional','two state', 'two state unfixed')
end
title(titles(i))
ylim([0,1.2])

end




sizes = [8,6,4,2,2,2];
loc1 = [4,7,8];
loc2 = [2,3,6];

plts = [2,4,6];
lags1 = [-10:10]; 
titles = {'Ser5ph-CTD','mRNA-CTD','mRNA-Ser5ph'};
fntsize = 18;   %set up figure 3

figure(10)


for i = 1:3
    
sigs = [  [sigred(loc1(i),end:-1:2),sigred(loc2(i),:)]  ; [sigred_cs(loc1(i),end:-1:2),sigred_cs(loc2(i),:)] ; [sigred_ms(loc1(i),end:-1:2),sigred_ms(loc2(i),:)] ; [sigred_frac(loc1(i),end:-1:2),sigred_frac(loc2(i),:)] ; [sigred_2s(loc1(i),end:-1:2),sigred_2s(loc2(i),:)] ; [sigred_2sf(loc1(i),end:-1:2),sigred_2sf(loc2(i),:)]   ];
subplot(3,2,plts(i))

for j = 1:6
    plot([-TT(end:-1:2);TT], sigs(j,:),'-','linewidth',sizes(j),'Color',colors(j,:))  %plot model rnap_acc
    if j == 1
        hold on;
    end
end

if i == 1
b1=errorbar(lags1, ser5_rnap.mn_cc, ser5_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;end

if i == 2
b1=errorbar(lags1, mrna_rnap.mn_cc, mrna_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;    
end

if i == 3
b1=errorbar(lags1, mrna_ser5.mn_cc, mrna_ser5.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;end


title(titles(i))
ylim([0,1.2])
xlim([-11,11])

end
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
saveas(gcf, 'comp_models','epsc')
