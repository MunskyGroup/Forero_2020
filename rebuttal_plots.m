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

cd ctd_to_ser5_model
load('best_simple_pars_cs.mat')  %load parameter file
[sigred_cs,TT_cs,means_cs,sig_cs,TT2_cs,mins_cs] = get_ac_and_cc_mod_cs(parameters,[0:.1:30]);
cd ..

%% mRNA solo model

cd mrna_solo_model
load('best_simple_pars_mrna_solo.mat')  %load parameter file
[sigred_ms,TT_ms,means_ms,sig_ms,TT2_ms,mins_ms] = get_ac_and_cc_mod_mrnasolo(parameters,[0:.1:30]);
cd ..

%% fractional model
load('best_simple_pars.mat')  %load parameter file
real_valscl = parameters;
parameters(7) = 1;
[sigred_frac,TT_frac,means_frac,sig_frac,TT2_frac,mins_frac] = get_ac_and_cc_mod_simplified(parameters,[0:.1:30]);



%% 2 state

cd 2state_model
load('best_simple_pars_2s.mat')  %load parameter file
[sigred_2s,TT_2s,means_2s,sig_2s,TT2_2s,mins_2s] = get_ac_and_cc_mod_2s(parameters,[0:.1:30]);
cd ..

%%
lags1 = [0:31]; 
fntsize = 18;   %set up figure 3
close all
figure(1)
subplot(3,2,1)
fig1= gcf;
%set(gcf,'color','k');
fig1.PaperUnits = 'centimeters';
xh = 2*3*6;
yh = 2*14.3/3;
fig1.PaperPosition = [0, 0, 2/3*xh, 2*yh]; % x,y, width, height
%sem_rnap = std(Xdata_sem(:,1:30),1)/(sqrt(20));

plot(TT, sigred(1,:),'g-','linewidth',6)  %plot model rnap_acc
hold on;
plot(TT, sigred_cs(1,:),'b-','linewidth',4)  %plot model rnap_acc
plot(TT, sigred_ms(1,:),'-','linewidth',2, 'Color', [242, 180,7]./256)  %plot model rnap_acc
plot(TT, sigred_frac(1,:),'y-','linewidth',1)  %plot model rnap_acc

plot(TT, sigred_2s(1,:),'m-','linewidth',1)  %plot model rnap_acc

%plot the data rnap_acc
b1=errorbar(lags1(1:25), rnap.mn_ac(1:25), rnap.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
ylim([0,1.2])
legend('Chosen','Phosphorylation','mRNA retention','Fractional','two state')

title('CTD')

%%

% lags1 = [0:31]; 
% fntsize = 18;   %set up figure 3
% 
% 
% subplot(3,3,1)
% fig1= gcf;
% %set(gcf,'color','k');
% fig1.PaperUnits = 'centimeters';
% xh = 2*3*6;
% yh = 2*14.3/3;
% fig1.PaperPosition = [0, 0, 2/3*xh, 2*yh]; % x,y, width, height
% %sem_rnap = std(Xdata_sem(:,1:30),1)/(sqrt(20));
% 
% plot(TT, sigred(1,:),'g-','linewidth',6)  %plot model rnap_acc
% hold on;
% plot(TT, sigred_cs(1,:),'b-','linewidth',4)  %plot model rnap_acc
% plot(TT, sigred_ms(1,:),'-','linewidth',2, 'Color', [242, 180,7]./256)  %plot model rnap_acc
% plot(TT, sigred_frac(1,:),'y-','linewidth',1)  %plot model rnap_acc
% 
% %plot(TT, sigred_2s(1,:),'m-','linewidth',1)  %plot model rnap_acc
% 
% %plot the data rnap_acc
% b1=errorbar(lags1(1:25), rnap.mn_ac(1:25), rnap.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
% ylim([.3852,.3856])
% xlim([3.7016,3.702])
% legend('simple','ctd>ser5','mrna solo','two state')
% 


%%
subplot(3,2,3)      %SER5 ACC figure
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

plot(TT, sigred(5,:),'g-','linewidth',6)  %plot model rnap_acc
hold on;

plot(TT, sigred_cs(5,:),'b-','linewidth',4)  %plot model rnap_acc

plot(TT, sigred_ms(5,:),'-','linewidth',2, 'Color', [242, 180,7]./256)  %plot model rnap_acc
plot(TT, sigred_frac(5,:),'y-','linewidth',1)  %plot model rnap_acc


plot(TT, sigred_2s(5,:),'m-','linewidth',1)  %plot model rnap_acc

b1=errorbar(lags1(1:25), ser5.mn_ac(1:25), ser5.sem_ac(1:25),'s','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
ylim([0,1.2])

title('Ser5ph')
ylabel('Correlation Function','FontSize',12,'FontWeight','bold')

%%
subplot(3,2,5)   %mRNA ACC figure
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

plot(TT, sigred(9,:),'g-','linewidth',6)  %plot model rnap_acc
hold on;
plot(TT, sigred_cs(9,:),'b-','linewidth',4)  %plot model rnap_acc
plot(TT, sigred_ms(9,:),'-','linewidth',2, 'Color', [242, 180,7]./256)  %plot model rnap_acc
plot(TT, sigred_frac(9,:),'y-','linewidth',1)  %plot model rnap_acc

plot(TT, sigred_2s(9,:),'m-','linewidth',1)  %plot model rnap_acc

b1=errorbar(lags1(1:25), mrna.mn_ac(1:25), mrna.sem_ac(1:25),'d','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
ylim([0,1.2])

title('mRNA')
xlabel('Lag Time (min) ','FontSize',12,'FontWeight','bold')

%%

lags2 = [-10:10]; %time vector for CCs

subplot(3,2,2)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

% Ser5 RNAP cross correlation
plot([-TT(end:-1:2);TT],[sigred(4,end:-1:2),sigred(2,:)],'g-','linewidth',6) %model
hold on;
plot([-TT(end:-1:2);TT],[sigred_cs(4,end:-1:2),sigred_cs(2,:)],'b-','linewidth',4) %model
plot([-TT(end:-1:2);TT],[sigred_ms(4,end:-1:2),sigred_ms(2,:)],'-','linewidth',2, 'Color', [242, 180,7]./256 ) %model
plot([-TT(end:-1:2);TT],[sigred_frac(4,end:-1:2),sigred_frac(2,:)],'y-','linewidth',1) %model


plot([-TT(end:-1:2);TT],[sigred_2s(4,end:-1:2),sigred_2s(2,:)],'m-','linewidth',1) %model

%data 
b1=errorbar(lags2, ser5_rnap.mn_cc, ser5_rnap.sem_cc,'s','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
xlim([-11,11])
ylim([0,1.2])

title('Ser5ph-CTD')


%%


lags2 = [-10:10]; %time vector for CCs
% 
% subplot(4,3,3)
% fig1= gcf;
% fig1.PaperUnits = 'centimeters';
% fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height
% 
% % Ser5 RNAP cross correlation
% plot([-TT(end:-1:2);TT],[sigred(4,end:-1:2),sigred(2,:)],'r-','linewidth',1) %model
% hold on;
% plot([-TT(end:-1:2);TT],[sigred_cs(4,end:-1:2),sigred_cs(2,:)],'g-.','linewidth',1) %model
% plot([-TT(end:-1:2);TT],[sigred_ms(4,end:-1:2),sigred_ms(2,:)],'b--','linewidth',1) %model
% plot([-TT(end:-1:2);TT],[sigred_2s(4,end:-1:2),sigred_2s(2,:)],'m-','linewidth',1) %model
% 
% %data 
% %b1=errorbar(lags2, mrna_rnap.mn_cc, mrna_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
% b1=errorbar(lags2, ser5_rnap.mn_cc, ser5_rnap.sem_cc,'s','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
% 
% xlim([-11,11])
% ylim([0,1.2])
% title('mRNA ACC')
% xlabel('tau')

%%

lags2 = [-10:10]; %time vector for CCs

subplot(3,2,4)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

plot([-TT(end:-1:2);TT],[sigred(7,end:-1:2),sigred(3,:)],'g-','linewidth',6) % mrna_rnap CC model
hold on;  %mrna_rnap CC data

plot([-TT(end:-1:2);TT],[sigred_cs(7,end:-1:2),sigred_cs(3,:)],'b-','linewidth',4) %model
plot([-TT(end:-1:2);TT],[sigred_ms(7,end:-1:2),sigred_ms(3,:)],'-','linewidth',2, 'Color', [242, 180,7]./256 ) %model

plot([-TT(end:-1:2);TT],[sigred_frac(7,end:-1:2),sigred_frac(3,:)],'y-', 'linewidth',1) %model



plot([-TT(end:-1:2);TT],[sigred_2s(7,end:-1:2),sigred_2s(3,:)],'m-','linewidth',1) %model

%data 
b1=errorbar(lags2, mrna_rnap.mn_cc, mrna_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
xlim([-11,11])
ylim([0,1.2])

title('mRNA-CTD')


%%
lags2 = [-10:10]; %time vector for CCs

subplot(3,2,6) %mrna_ser5 signal
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height


plot([-TT(end:-1:2);TT],[sigred(8,end:-1:2),sigred(6,:)],'g-','linewidth',6)
hold on;
plot([-TT(end:-1:2);TT],[sigred_cs(8,end:-1:2),sigred_cs(6,:)],'b-','linewidth',4) %model
plot([-TT(end:-1:2);TT],[sigred_ms(8,end:-1:2),sigred_ms(6,:)],'-','linewidth',2, 'Color', [242, 180,7]./256 ) %model
plot([-TT(end:-1:2);TT],[sigred_frac(8,end:-1:2),sigred_frac(6,:)],'y-','linewidth',1 ) %model


plot([-TT(end:-1:2);TT],[sigred_2s(8,end:-1:2),sigred_2s(6,:)],'m-','linewidth',1) %model

%data 
b1=errorbar(lags2, mrna_ser5.mn_cc, mrna_ser5.sem_cc,'d','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
xlim([-11,11])
ylim([0,1.2])
title('mRNA-Ser5ph')
xlabel('Lag Time (min)','FontSize',12,'FontWeight','bold')

annotation('textbox', [0.25, 0.9, 0.1, 0.1], 'String', "Auto-Correlation",'FontSize',12,'FontWeight','bold' )
annotation('textbox', [0.7, 0.9, 0.1, 0.1], 'String', "Cross-Correlation",'FontSize',12,'FontWeight','bold')

pause

%%
lags1 = [0:31]; 
fntsize = 18;   %set up figure 3
close all
figure(1)
subplot(3,2,1)
fig1= gcf;
%set(gcf,'color','k');
fig1.PaperUnits = 'centimeters';
xh = 2*3*6;
yh = 2*14.3/3;
fig1.PaperPosition = [0, 0, 2/3*xh, 2*yh]; % x,y, width, height
%sem_rnap = std(Xdata_sem(:,1:30),1)/(sqrt(20));

plot(TT, sigred(1,:),'g-','linewidth',6)  %plot model rnap_acc
hold on;
%plot(TT, sigred_cs(1,:),'b-','linewidth',4)  %plot model rnap_acc
plot(TT, sigred_ms(1,:),'-','linewidth',2, 'Color', [8, 16, 161]./256)  %plot model rnap_acc
%plot(TT, sigred_frac(1,:),'y-','linewidth',1)  %plot model rnap_acc

%plot(TT, sigred_2s(1,:),'m-','linewidth',1)  %plot model rnap_acc

%plot the data rnap_acc
b1=errorbar(lags1(1:25), rnap.mn_ac(1:25), rnap.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
ylim([0,1.2])
legend('Chosen','mRNA retention')

title('CTD')

%%

% lags1 = [0:31]; 
% fntsize = 18;   %set up figure 3
% 
% 
% subplot(3,3,1)
% fig1= gcf;
% %set(gcf,'color','k');
% fig1.PaperUnits = 'centimeters';
% xh = 2*3*6;
% yh = 2*14.3/3;
% fig1.PaperPosition = [0, 0, 2/3*xh, 2*yh]; % x,y, width, height
% %sem_rnap = std(Xdata_sem(:,1:30),1)/(sqrt(20));
% 
% plot(TT, sigred(1,:),'g-','linewidth',6)  %plot model rnap_acc
% hold on;
% plot(TT, sigred_cs(1,:),'b-','linewidth',4)  %plot model rnap_acc
% plot(TT, sigred_ms(1,:),'-','linewidth',2, 'Color', [242, 180,7]./256)  %plot model rnap_acc
% plot(TT, sigred_frac(1,:),'y-','linewidth',1)  %plot model rnap_acc
% 
% %plot(TT, sigred_2s(1,:),'m-','linewidth',1)  %plot model rnap_acc
% 
% %plot the data rnap_acc
% b1=errorbar(lags1(1:25), rnap.mn_ac(1:25), rnap.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
% ylim([.3852,.3856])
% xlim([3.7016,3.702])
% legend('simple','ctd>ser5','mrna solo','two state')
% 


%%
subplot(3,2,3)      %SER5 ACC figure
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

plot(TT, sigred(5,:),'g-','linewidth',6)  %plot model rnap_acc
hold on;

%plot(TT, sigred_cs(5,:),'b-','linewidth',4)  %plot model rnap_acc

plot(TT, sigred_ms(5,:),'-','linewidth',2, 'Color', [8, 16, 161]./256)  %plot model rnap_acc
%plot(TT, sigred_frac(5,:),'y-','linewidth',1)  %plot model rnap_acc


%plot(TT, sigred_2s(5,:),'m-','linewidth',1)  %plot model rnap_acc

b1=errorbar(lags1(1:25), ser5.mn_ac(1:25), ser5.sem_ac(1:25),'s','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
ylim([0,1.2])

title('Ser5ph')
ylabel('Correlation Function','FontSize',12,'FontWeight','bold')

%%
subplot(3,2,5)   %mRNA ACC figure
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

plot(TT, sigred(9,:),'g-','linewidth',6)  %plot model rnap_acc
hold on;
%plot(TT, sigred_cs(9,:),'b-','linewidth',4)  %plot model rnap_acc
plot(TT, sigred_ms(9,:),'-','linewidth',2, 'Color', [8, 16, 161]./256)  %plot model rnap_acc
%plot(TT, sigred_frac(9,:),'y-','linewidth',1)  %plot model rnap_acc

%plot(TT, sigred_2s(9,:),'m-','linewidth',1)  %plot model rnap_acc

b1=errorbar(lags1(1:25), mrna.mn_ac(1:25), mrna.sem_ac(1:25),'d','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
ylim([0,1.2])

title('mRNA')
xlabel('Lag Time (min) ','FontSize',12,'FontWeight','bold')

%%

lags2 = [-10:10]; %time vector for CCs

subplot(3,2,2)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

% Ser5 RNAP cross correlation
plot([-TT(end:-1:2);TT],[sigred(4,end:-1:2),sigred(2,:)],'g-','linewidth',6) %model
hold on;
%plot([-TT(end:-1:2);TT],[sigred_cs(4,end:-1:2),sigred_cs(2,:)],'b-','linewidth',4) %model
plot([-TT(end:-1:2);TT],[sigred_ms(4,end:-1:2),sigred_ms(2,:)],'-','linewidth',2, 'Color', [8, 16, 161]./256 ) %model
%plot([-TT(end:-1:2);TT],[sigred_frac(4,end:-1:2),sigred_frac(2,:)],'y-','linewidth',1) %model


%plot([-TT(end:-1:2);TT],[sigred_2s(4,end:-1:2),sigred_2s(2,:)],'m-','linewidth',1) %model

%data 
b1=errorbar(lags2, ser5_rnap.mn_cc, ser5_rnap.sem_cc,'s','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
xlim([-11,11])
ylim([0,1.2])

title('Ser5ph-CTD')


%%


lags2 = [-10:10]; %time vector for CCs
% 
% subplot(4,3,3)
% fig1= gcf;
% fig1.PaperUnits = 'centimeters';
% fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height
% 
% % Ser5 RNAP cross correlation
% plot([-TT(end:-1:2);TT],[sigred(4,end:-1:2),sigred(2,:)],'r-','linewidth',1) %model
% hold on;
% plot([-TT(end:-1:2);TT],[sigred_cs(4,end:-1:2),sigred_cs(2,:)],'g-.','linewidth',1) %model
% plot([-TT(end:-1:2);TT],[sigred_ms(4,end:-1:2),sigred_ms(2,:)],'b--','linewidth',1) %model
% plot([-TT(end:-1:2);TT],[sigred_2s(4,end:-1:2),sigred_2s(2,:)],'m-','linewidth',1) %model
% 
% %data 
% %b1=errorbar(lags2, mrna_rnap.mn_cc, mrna_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
% b1=errorbar(lags2, ser5_rnap.mn_cc, ser5_rnap.sem_cc,'s','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
% 
% xlim([-11,11])
% ylim([0,1.2])
% title('mRNA ACC')
% xlabel('tau')

%%

lags2 = [-10:10]; %time vector for CCs

subplot(3,2,4)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

plot([-TT(end:-1:2);TT],[sigred(7,end:-1:2),sigred(3,:)],'g-','linewidth',6) % mrna_rnap CC model
hold on;  %mrna_rnap CC data

%plot([-TT(end:-1:2);TT],[sigred_cs(7,end:-1:2),sigred_cs(3,:)],'b-','linewidth',4) %model
plot([-TT(end:-1:2);TT],[sigred_ms(7,end:-1:2),sigred_ms(3,:)],'-','linewidth',2, 'Color', [8, 16, 161]./256 ) %model

%plot([-TT(end:-1:2);TT],[sigred_frac(7,end:-1:2),sigred_frac(3,:)],'y-', 'linewidth',1) %model



%plot([-TT(end:-1:2);TT],[sigred_2s(7,end:-1:2),sigred_2s(3,:)],'m-','linewidth',1) %model

%data 
b1=errorbar(lags2, mrna_rnap.mn_cc, mrna_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
xlim([-11,11])
ylim([0,1.2])

title('mRNA-CTD')


%%
lags2 = [-10:10]; %time vector for CCs

subplot(3,2,6) %mrna_ser5 signal
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height


plot([-TT(end:-1:2);TT],[sigred(8,end:-1:2),sigred(6,:)],'g-','linewidth',6)
hold on;
%plot([-TT(end:-1:2);TT],[sigred_cs(8,end:-1:2),sigred_cs(6,:)],'b-','linewidth',4) %model
plot([-TT(end:-1:2);TT],[sigred_ms(8,end:-1:2),sigred_ms(6,:)],'-','linewidth',2, 'Color', [8, 16, 161]./256 ) %model
%plot([-TT(end:-1:2);TT],[sigred_frac(8,end:-1:2),sigred_frac(6,:)],'y-','linewidth',1 ) %model


%plot([-TT(end:-1:2);TT],[sigred_2s(8,end:-1:2),sigred_2s(6,:)],'m-','linewidth',1) %model

%data 
b1=errorbar(lags2, mrna_ser5.mn_cc, mrna_ser5.sem_cc,'d','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
xlim([-11,11])
ylim([0,1.2])
title('mRNA-Ser5ph')
xlabel('Lag Time (min)','FontSize',12,'FontWeight','bold')

annotation('textbox', [0.25, 0.9, 0.1, 0.1], 'String', "Auto-Correlation",'FontSize',12,'FontWeight','bold' )
annotation('textbox', [0.7, 0.9, 0.1, 0.1], 'String', "Cross-Correlation",'FontSize',12,'FontWeight','bold')
pause 


%%
% 
% f = figure(1);
% 
% fig1= gcf;
% pos = get(subplot(3,3,7),'position'); %mrna_ser5 signal
% delete(subplot(3,3,7))
% 
% %set(f,'Position',[0, 3000, xh, 3300]);
% dat =  {'        Simple', 39.95, '        --';...
%         '        CTD>Ser5', 40.86, '        --';...   
%         '        mRNA solo', 41.76, '        --';...
%         };
% columnname =   {'Bayes Information Criterion', 'Value', 'Units'};
% columnformat = {'char', 'numeric', 'char'}; 
% t = uitable('Data', dat,... 
%             'ColumnName', columnname,...
%             'ColumnFormat', columnformat,...
%             'RowName',[]);
%   
% set(t,'units','normalized')
% set(t,'position',pos)       
%         
%%


X_SIZE = 13; Y_SIZE = 15;
figure(10);clf;
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

cd ctd_to_ser5_model
load('best_simple_pars_cs.mat')  %load parameter file
[sigred_cs,TT_cs,means_cs,sig_cs,TT2_cs,mins_cs] = get_ac_and_cc_mod_cs(parameters,[0:.1:30]);
cd ..

%% mRNA solo model

cd mrna_solo_model
load('best_simple_pars_mrna_solo.mat')  %load parameter file
[sigred_ms,TT_ms,means_ms,sig_ms,TT2_ms,mins_ms] = get_ac_and_cc_mod_mrnasolo(parameters,[0:.1:30]);
cd ..

%% 2 state

cd 2state_model
load('best_simple_pars_2s.mat')  %load parameter file
[sigred_2s,TT_2s,means_2s,sig_2s,TT2_2s,mins_2s] = get_ac_and_cc_mod_2s(parameters,[0:.1:30]);
cd ..

%%
lags1 = [0:31]; 
fntsize = 18;   %set up figure 3

figure(10)
subplot(3,3,2)
fig1= gcf;
%set(gcf,'color','k');
fig1.PaperUnits = 'centimeters';
xh = 2*3*6;
yh = 2*14.3/3;
fig1.PaperPosition = [0, 0, 2/3*xh, 2*yh]; % x,y, width, height
%sem_rnap = std(Xdata_sem(:,1:30),1)/(sqrt(20));

plot(TT, sigred(1,:),'r-','linewidth',1)  %plot model rnap_acc
hold on;
%plot(TT, sigred_cs(1,:),'g-.','linewidth',1)  %plot model rnap_acc
%plot(TT, sigred_ms(1,:),'b--','linewidth',1)  %plot model rnap_acc
plot(TT, sigred_2s(1,:),'m-','linewidth',1)  %plot model rnap_acc

%plot the data rnap_acc
b1=errorbar(lags1(1:25), rnap.mn_ac(1:25), rnap.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
ylim([0,1.2])
legend('simple','ctd>ser5','mrna solo','two state')

title('RNAP ACC')

%%

lags1 = [0:31]; 
fntsize = 18;   %set up figure 3


subplot(3,3,1)
fig1= gcf;
%set(gcf,'color','k');
fig1.PaperUnits = 'centimeters';
xh = 2*3*6;
yh = 2*14.3/3;
fig1.PaperPosition = [0, 0, 2/3*xh, 2*yh]; % x,y, width, height
%sem_rnap = std(Xdata_sem(:,1:30),1)/(sqrt(20));

plot(TT, sigred(1,:),'b-','linewidth',1)  %plot model rnap_acc
hold on;
%plot(TT, sigred_cs(1,:),'g-.','linewidth',1)  %plot model rnap_acc
%plot(TT, sigred_ms(1,:),'b--','linewidth',1)  %plot model rnap_acc
plot(TT, sigred_2s(1,:),'m-','linewidth',1)  %plot model rnap_acc

%plot the data rnap_acc
b1=errorbar(lags1(1:25), rnap.mn_ac(1:25), rnap.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
ylim([.3852,.3856])
xlim([3.7016,3.702])
legend('simple','ctd>ser5','mrna solo','two state')



%%
subplot(3,3,5)      %SER5 ACC figure
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

plot(TT, sigred(5,:),'b-','linewidth',1)  %plot model rnap_acc
hold on;
%plot(TT, sigred_cs(5,:),'g-.','linewidth',1)  %plot model rnap_acc
%plot(TT, sigred_ms(5,:),'b--','linewidth',1)  %plot model rnap_acc
plot(TT, sigred_2s(5,:),'m-','linewidth',1)  %plot model rnap_acc

b1=errorbar(lags1(1:25), ser5.mn_ac(1:25), ser5.sem_ac(1:25),'s','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
ylim([0,1.2])

title('Ser5ph ACC')


%%
subplot(3,3,8)   %mRNA ACC figure
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

plot(TT, sigred(9,:),'b-','linewidth',1)  %plot model rnap_acc
hold on;
%plot(TT, sigred_cs(9,:),'g-.','linewidth',1)  %plot model rnap_acc
%plot(TT, sigred_ms(9,:),'b--','linewidth',1)  %plot model rnap_acc
plot(TT, sigred_2s(9,:),'m-','linewidth',1)  %plot model rnap_acc

b1=errorbar(lags1(1:25), mrna.mn_ac(1:25), mrna.sem_ac(1:25),'d','MarkerSize',5,'MarkerFaceColor',[0, 0, 0]./256,'Color',[0, 0, 0]./256);hold on;
ylim([0,1.2])

title('mRNA ACC')
xlabel('tau')

%%

lags2 = [-10:10]; %time vector for CCs

subplot(3,3,3)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

% Ser5 RNAP cross correlation
plot([-TT(end:-1:2);TT],[sigred(4,end:-1:2),sigred(2,:)],'b-','linewidth',1) %model
hold on;
%plot([-TT(end:-1:2);TT],[sigred_cs(4,end:-1:2),sigred_cs(2,:)],'g-.','linewidth',1) %model
%plot([-TT(end:-1:2);TT],[sigred_ms(4,end:-1:2),sigred_ms(2,:)],'b--','linewidth',1) %model
plot([-TT(end:-1:2);TT],[sigred_2s(4,end:-1:2),sigred_2s(2,:)],'m-','linewidth',1) %model

%data 
b1=errorbar(lags2, ser5_rnap.mn_cc, ser5_rnap.sem_cc,'s','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
xlim([-11,11])
ylim([0,1.2])

title('Ser5ph-Rnap CC')


%%


lags2 = [-10:10]; %time vector for CCs
% 
% subplot(4,3,3)
% fig1= gcf;
% fig1.PaperUnits = 'centimeters';
% fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height
% 
% % Ser5 RNAP cross correlation
% plot([-TT(end:-1:2);TT],[sigred(4,end:-1:2),sigred(2,:)],'r-','linewidth',1) %model
% hold on;
% plot([-TT(end:-1:2);TT],[sigred_cs(4,end:-1:2),sigred_cs(2,:)],'g-.','linewidth',1) %model
% plot([-TT(end:-1:2);TT],[sigred_ms(4,end:-1:2),sigred_ms(2,:)],'b--','linewidth',1) %model
% plot([-TT(end:-1:2);TT],[sigred_2s(4,end:-1:2),sigred_2s(2,:)],'m-','linewidth',1) %model
% 
% %data 
% %b1=errorbar(lags2, mrna_rnap.mn_cc, mrna_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
% b1=errorbar(lags2, ser5_rnap.mn_cc, ser5_rnap.sem_cc,'s','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
% 
% xlim([-11,11])
% ylim([0,1.2])
% title('mRNA ACC')
% xlabel('tau')

%%

lags2 = [-10:10]; %time vector for CCs

subplot(3,3,6)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

plot([-TT(end:-1:2);TT],[sigred(7,end:-1:2),sigred(3,:)],'b-','linewidth',1) % mrna_rnap CC model
hold on;  %mrna_rnap CC data

%plot([-TT(end:-1:2);TT],[sigred_cs(7,end:-1:2),sigred_cs(3,:)],'g-.','linewidth',1) %model
%plot([-TT(end:-1:2);TT],[sigred_ms(7,end:-1:2),sigred_ms(3,:)],'b--','linewidth',1) %model
plot([-TT(end:-1:2);TT],[sigred_2s(7,end:-1:2),sigred_2s(3,:)],'m-','linewidth',1) %model

%data 
b1=errorbar(lags2, mrna_rnap.mn_cc, mrna_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
xlim([-11,11])
ylim([0,1.2])

title('mRNA-RNAP CC')


%%
lags2 = [-10:10]; %time vector for CCs

subplot(3,3,9) %mrna_ser5 signal
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height


plot([-TT(end:-1:2);TT],[sigred(8,end:-1:2),sigred(6,:)],'b-','linewidth',1)
hold on;
%plot([-TT(end:-1:2);TT],[sigred_cs(8,end:-1:2),sigred_cs(6,:)],'g-.','linewidth',1) %model
%plot([-TT(end:-1:2);TT],[sigred_ms(8,end:-1:2),sigred_ms(6,:)],'b--','linewidth',1) %model
plot([-TT(end:-1:2);TT],[sigred_2s(8,end:-1:2),sigred_2s(6,:)],'m-','linewidth',1) %model

%data 
b1=errorbar(lags2, mrna_ser5.mn_cc, mrna_ser5.sem_cc,'d','MarkerSize',5,'MarkerFaceColor',[0, 0, 0],'Color',[0, 0, 0]);hold on;
xlim([-11,11])
ylim([0,1.2])
title('mRNA-Ser5ph CC')
xlabel('tau')

%%

f = figure(10);

fig1= gcf;
pos = get(subplot(3,3,7),'position'); %mrna_ser5 signal
delete(subplot(3,3,7))

%set(f,'Position',[0, 3000, xh, 3300]);
dat =  {'        Paper Model', 39.95, '        --';...
        '        Two State Bursting',  59.51, '        --';...
        };
columnname =   {'Bayes Information Criterion', 'Value', 'Units'};
columnformat = {'char', 'numeric', 'char'}; 
t = uitable('Data', dat,... 
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'RowName',[]);
  
set(t,'units','normalized')
set(t,'position',pos)       
        














%% delay times plot 
figure(2);

cd ctd_to_ser5_model

load('sub_min_data.mat')
load('best_simple_pars_cs.mat')
rnap_ser5_submin = flip(movmean(sub_min_data(:,2),60 ));

%sub_norm = rnap_ser5_submin / max(rnap_ser5_submin );


sub_norm = rescale(rnap_ser5_submin,.66,1);
TTs = linspace(-1,1,1198).*60;

valid_delays = TTs(find(sub_norm > .99));
highlight = [min(valid_delays),max(valid_delays)];

load('best_simple_pars_cs.mat')  %load parameter file
[sigred_cs,TT_cs,means_cs,sig_cs,TT2_cs,mins_cs] = get_ac_and_cc_mod_cs(parameters,[0:.1:30]);

par_sweep = parameters;
tpars = linspace(1,100,200);
delay_times = [];

for i = 1:length(tpars)
    par_sweep(4) = tpars(i);
    [sigred,TT,means,sig] = get_ac_and_cc_mod_cs(par_sweep,[0:1/600:15]);
    TT = 60.*TT;
    t = [-TT(end:-1:2);TT];
    delay_time = t(find([sigred(4,end:-1:2),sigred(2,:)]/max([sigred(4,end:-1:2),sigred(2,:)]) == 1));
    delay_times = [delay_times, delay_time];

end


plot(tpars,delay_times,'k','LineWidth',2);
xlabel('ctd --> ctd+ser5ph rate (1/min)')
ylabel('delay time (s)')
hold on;
h = fill([tpars(1),tpars(end),tpars(end),tpars(1)],[highlight(1),highlight(1),highlight(2),highlight(2)],'green','FaceAlpha',.3, 'LineStyle','none');


valid_pars = tpars(find(delay_times <= highlight(2) & delay_times >= highlight(1)   ));
par_highlight = [min(valid_pars),max(valid_pars)];

h = fill([par_highlight(1),par_highlight(end),par_highlight(end),par_highlight(1)],[highlight(1),highlight(1),highlight(2),highlight(2)],'blue','FaceAlpha',.3, 'LineStyle','none');

legend('Model Ser5ph-CTD Delay', 'Experimental delay range','Overlap')

cd ..


%%
figure(3)



load('best_simple_pars.mat')  %load parameter file
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

figure(30); 
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, X_SIZE, Y_SIZE]; % x,y, width, height
global_color = 'k';
cla;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 1.1*yh]; % x,y, width, height

for i = 1:4
    
    sol = run_single_SSA_linda(x0,S,W,T_array,time_var,signal_update_rate);  
    subplot(4,2,2*i)
    pol2_ssa = sol(2,:)';
    ser5_ssa = sol(2,:)';
    ts_ssa = sol(3,:)';

    [pol2_ssa,ser5_ssa,ts_ssa,~] = get_model_intesities(sol,eta_rnap,eta_ser5,eta_ts); %convert molecules to signal
    [pol2norm,ser5norm,tsnorm] = Normalize_simulated_intensities(.95,pol2_ssa,ser5_ssa,ts_ssa);

    X_SIZE = 13; Y_SIZE = 15;

    avpol2 =movmean(pol2norm(end-200:end),3);
    avser5 = movmean(ser5norm(end-200:end),3);
    avts = movmean(tsnorm(end-200:end),3);
    plot(avpol2(1:201),'r','linewidth',2); hold on; plot(avser5(1:201),'g','linewidth',2); plot(avts(1:201),'b','linewidth',2 )
    if i == 1
        legend('POL2','SER5ph','mRNA')

        title('Simulated trace','FontSize',fntsize,'FontWeight','bold')
        rng(88)
    end
    
    
    ylabel('Normalized Int','FontSize',fntsize-4,'FontWeight','bold')
    if i == 4
        xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
    end
    xlim([0,200])
    ylim([-.5,1.5])
end


[POL2_I_Norm, SER5_I_Norm, MRNA_I_Norm] = Normalize_raw_intensities(.95);

toplot = [1,14,7,19];
for i = 1:1
    subplot(4,2,2*(i-1)+1)
    avpol2 =movmean(POL2_I_Norm(:,toplot(i)) ,3);
    avser5 = movmean(SER5_I_Norm(:,toplot(i)) ,3);
    avts = movmean(MRNA_I_Norm(:,toplot(i)) ,3);
    plot(avpol2,'r','linewidth',2); hold on; plot(avser5,'g','linewidth',2); plot(avts,'b','linewidth',2 )
     if i == 1
        legend('POL2','SER5ph','mRNA')

        title('Experimental trace (Shown in Fig. 2c)','FontSize',fntsize,'FontWeight','bold')
    end
    ylabel('Normalized Int','FontSize',fntsize-4,'FontWeight','bold')
    if i == 4
        xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
    end
    xlim([0,200])
    ylim([-.5,1.5])  
%     if i == 3:
%     text(10,1.3,string(toplot(i)),'Color','black','FontSize',14,'FontWeight','bold')
%     
%     if i==1
%         text(35,1.3,'Figure 2c','Color','black','FontSize',14,'FontWeight','bold')
%         
%     end
     if i==2
        text(10,1.3,'Sup. Figure 2b','Color','black','FontSize',14,'FontWeight','bold')
        
    end
    
end

%% Minima Plots

figure(4)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, X_SIZE, Y_SIZE]; % x,y, width, height
global_color = 'k';
sm = SimpleModel;
load('best_simple_pars.mat')  %load parameter file
real_valscl = parameters;


pol2_traj = [];
ser5_traj = [];
mrna_traj = [];

for i = 1:20
    [~,~,norm_ints] = sm.ssa_traj([200:399]);
    pol2_traj = [pol2_traj; norm_ints(1,:)];
    ser5_traj = [ser5_traj; norm_ints(2,:)];
    mrna_traj = [mrna_traj; norm_ints(3,:)];
end

pol2_traj = pol2_traj';
ser5_traj = ser5_traj';
mrna_traj = mrna_traj';

%plot(mrna_traj)
PlusMinusRange = 10; % renge of values to select after finding a local minima.
rangeToDisplay = 5; % seconds to display in the plot it should be smaller than PlusMinusRange.
fileDataComplete = 'Data_with_Oscillations.xlsx'; % experimental data with cells that show oscillatory behaviour.
detectionThreshold = 0.8; % minimal threshold to detect a minimum in the signal.
withRespectToMinValue = 'mRNA'; % 'mRNA' 'Unph'  %% Sets time interval with respect to the minim of selected signal.


UnphFluc=xlsread('Data_with_Oscillations.xlsx','Unph');
Ser5phFluc=xlsread('Data_with_Oscillations.xlsx','Ser5ph');
mRNAFluc=xlsread('Data_with_Oscillations.xlsx','mRNA');

[POL2_I_Norm, SER5_I_Norm, MRNA_I_Norm] = Normalize_raw_intensities(.95);


[mRNAFluc_Rave,Ser5phFluc_Rave, UnphFluc_Rave, vector_minima_mRNA,lmin_mRNA] = minimaLocation(pol2_traj,ser5_traj,mrna_traj ,detectionThreshold);
minima_sum = 0;
for n = 1:length(vector_minima_mRNA)
    minima_sum = minima_sum + length(vector_minima_mRNA{n});
end
[mean_mRNA, mean_Ser5ph, mean_Unph, sem_mRNA, sem_Ser5ph, sem_Unph,mRNA_mRNA] = intensitySelected_FromMinimum_InARange (mRNAFluc_Rave,UnphFluc_Rave,Ser5phFluc_Rave,vector_minima_mRNA, withRespectToMinValue, PlusMinusRange);
[fit_mRNA, mu_mRNA, sem_mRNA,fit_Ser5ph, mu_Ser5ph, sem_Ser5ph, fit_Unph, mu_Unph, sem_Unph]= fit_Gaussian_MainCode_2 (mean_mRNA, mean_Ser5ph, mean_Unph, sem_mRNA, sem_Ser5ph, sem_Unph, rangeToDisplay);




[mRNAFluc_data,Ser5phFluc_data, UnphFluc_data, vector_minima_mRNA_data,lmin_mRNA_data] = minimaLocation(POL2_I_Norm,SER5_I_Norm,MRNA_I_Norm ,detectionThreshold);
[mean_mRNA_data, mean_Ser5ph_data, mean_Unph_data, sem_mRNA_data, sem_Ser5ph_data, sem_Unph_data,mRNA_mRNA_data] = intensitySelected_FromMinimum_InARange (mRNAFluc_data,UnphFluc_data,Ser5phFluc_data,vector_minima_mRNA_data, withRespectToMinValue, PlusMinusRange);


subplot(2,2,1)
x= -rangeToDisplay:rangeToDisplay;
b1=errorbar ( x, mu_Unph, sem_Unph, 'o','MarkerSize',5,'MarkerFaceColor','r','Color','r');
hold on;
axis([-5.5 5.5 -0.4 0.7])
x2 = linspace(-5.5,5,100);


b1=errorbar (x, mean_Unph_data(6:16), sem_Unph_data(6:16), 'o','MarkerSize',5,'MarkerFaceColor',[.2,.2,.2],'Color',[.2,.2,.2]);
plot(x2,fit_Unph(x2),'-','Color',[1,0,0]);


title('CTD minima')
ylabel('Normalized Intensity')
legend('Sim','Data')
ylim([0,.7])


subplot(2,2,2)
b2=errorbar( x, mu_Ser5ph, sem_Ser5ph,'s','MarkerSize',5,'MarkerFaceColor','g','Color','g');
hold on;
axis([-5.5 5.5 -0.4 0.7])

b1=errorbar ( x, mean_Ser5ph_data(6:16), sem_Ser5ph_data(6:16), 'o','MarkerSize',5,'MarkerFaceColor',[.2,.2,.2],'Color',[.2,.2,.2]);
plot( x2, fit_Ser5ph(x2),'-','Color',[0,1,0]);

title('Ser5ph minima')
ylim([0,.7])
legend('Sim','Data')


subplot(2,2,3);
b2=errorbar( x, mu_mRNA, sem_mRNA,'s','MarkerSize',5,'MarkerFaceColor','b','Color','b');
hold on;
axis([-5.5 5.5 -0.4 0.7])

b1=errorbar (x, mean_mRNA_data(6:16), sem_mRNA_data(6:16), 'o','MarkerSize',5,'MarkerFaceColor',[.2,.2,.2],'Color',[.2,.2,.2]);
plot( x2, fit_mRNA(x2),'-','Color',[0,0,1]);
hold('off')

title('mRNA minima')
xlabel('Tau (min)')
ylabel('Normalized Intensity')
legend('Sim','Data')
ylim([0,.7])

subplot(2,2,4);
b2=errorbar( x, mu_Ser5ph, sem_Ser5ph,'s','MarkerSize',5,'MarkerFaceColor','g','Color','g');
hold on;
b1=errorbar ( x, mu_Unph, sem_Unph, 'o','MarkerSize',5,'MarkerFaceColor','r','Color','r');
b2=errorbar( x, mu_mRNA, sem_mRNA,'s','MarkerSize',5,'MarkerFaceColor','b','Color','b');

plot([0,0],[0,1],'k--')
ylim([0,.7])


axis([-5.5 5.5 -0.4 0.7])
h2 = plot(x2, fit_mRNA(x2),'-','Color',[0,0,1]);
h2 = plot(x2, fit_Ser5ph(x2),'-','Color',[0,1,0]);
h1 = plot(x2, fit_Unph(x2),'-','Color',[1,0,0]);

legend('Ser5ph','CTD','mRNA')
title('Simulation minima')
xlabel('Tau (min)')
ylim([0,.7])

figure(5)

plot(cumsum(mean_Unph_data)./sum(mean_Unph_data),'r'); hold on; plot(cumsum(mean_Ser5ph_data)./sum(mean_Ser5ph_data),'g');plot(cumsum(mean_mRNA_data)./sum(mean_mRNA_data),'b');


figure(6)
mrna_data_cdf = [];
mrna_cdf = [];
for i = 0:1/25:1
    mrna_data_cdf = [mrna_data_cdf, sum(mean_mRNA_data <= i )/length(mean_mRNA_data)];
    mrna_cdf = [mrna_cdf, sum(mu_mRNA <= i )/length(mu_mRNA)];
    
    
end

plot(mrna_data_cdf); hold on; plot(mrna_cdf);

%%

sm = SimpleModel;
load('best_simple_pars.mat')  %load parameter file
real_valscl = parameters;

mrna_sim_cdf = [];
ser5_sim_cdf = [];
ctd_sim_cdf = [];

mu_sim_mrna = [];
mu_sim_ser5 = [];
mu_sim_ctd = [];

for k = .6:.05:.95

    pol2_traj = [];
    ser5_traj = [];
    mrna_traj = [];

    for i = 1:20
        [~,~,norm_ints] = sm.ssa_traj([200:399]);
        pol2_traj = [pol2_traj; norm_ints(1,:)];
        ser5_traj = [ser5_traj; norm_ints(2,:)];
        mrna_traj = [mrna_traj; norm_ints(3,:)];
    end

    pol2_traj = pol2_traj';
    ser5_traj = ser5_traj';
    mrna_traj = mrna_traj';

    %plot(mrna_traj)
    PlusMinusRange = 10; % renge of values to select after finding a local minima.
    rangeToDisplay = 5; % seconds to display in the plot it should be smaller than PlusMinusRange.
    fileDataComplete = 'Data_with_Oscillations.xlsx'; % experimental data with cells that show oscillatory behaviour.
    detectionThreshold = k; % minimal threshold to detect a minimum in the signal.
    withRespectToMinValue = 'mRNA'; % 'mRNA' 'Unph'  %% Sets time interval with respect to the minim of selected signal.


    UnphFluc=xlsread('Data_with_Oscillations.xlsx','Unph');
    Ser5phFluc=xlsread('Data_with_Oscillations.xlsx','Ser5ph');
    mRNAFluc=xlsread('Data_with_Oscillations.xlsx','mRNA');

    [POL2_I_Norm, SER5_I_Norm, MRNA_I_Norm] = Normalize_raw_intensities(.95);


    [mRNAFluc_Rave,Ser5phFluc_Rave, UnphFluc_Rave, vector_minima_mRNA,lmin_mRNA] = minimaLocation(pol2_traj,ser5_traj,mrna_traj ,detectionThreshold);
    minima_sum = 0;
    for n = 1:length(vector_minima_mRNA)
        minima_sum = minima_sum + length(vector_minima_mRNA{n});
    end
    [mean_mRNA, mean_Ser5ph, mean_Unph, sem_mRNA, sem_Ser5ph, sem_Unph,mRNA_mRNA] = intensitySelected_FromMinimum_InARange (mRNAFluc_Rave,UnphFluc_Rave,Ser5phFluc_Rave,vector_minima_mRNA, withRespectToMinValue, PlusMinusRange);
    [fit_mRNA, mu_mRNA, sem_mRNA,fit_Ser5ph, mu_Ser5ph, sem_Ser5ph, fit_Unph, mu_Unph, sem_Unph]= fit_Gaussian_MainCode_2 (mean_mRNA, mean_Ser5ph, mean_Unph, sem_mRNA, sem_Ser5ph, sem_Unph, rangeToDisplay);




    [mRNAFluc_data,Ser5phFluc_data, UnphFluc_data, vector_minima_mRNA_data,lmin_mRNA_data] = minimaLocation(POL2_I_Norm,SER5_I_Norm,MRNA_I_Norm ,detectionThreshold);
    [mean_mRNA_data, mean_Ser5ph_data, mean_Unph_data, sem_mRNA_data, sem_Ser5ph_data, sem_Unph_data,mRNA_mRNA_data] = intensitySelected_FromMinimum_InARange (mRNAFluc_data,UnphFluc_data,Ser5phFluc_data,vector_minima_mRNA_data, withRespectToMinValue, PlusMinusRange);


    subplot(2,2,1)
    x= -rangeToDisplay:rangeToDisplay;
    b1=errorbar ( x, mu_Unph, sem_Unph, 'o','MarkerSize',5,'MarkerFaceColor','r','Color','r');
    hold on;
    axis([-5.5 5.5 -0.4 0.7])
    x2 = linspace(-5.5,5,100);
    plot(x2,fit_Unph(x2),'-','Color',[1,0,0]);

    b1=errorbar (x, mean_Unph_data(6:16), sem_Unph_data(6:16), 'o','MarkerSize',5,'MarkerFaceColor',[.2,.2,.2],'Color',[.2,.2,.2]);

    
    subplot(2,2,2)
    b2=errorbar( x, mu_Ser5ph, sem_Ser5ph,'s','MarkerSize',5,'MarkerFaceColor','g','Color','g');
    hold on;
    axis([-5.5 5.5 -0.4 0.7])
    plot( x2, fit_Ser5ph(x2),'-','Color',[0,1,0]);
    b1=errorbar ( x, mean_Ser5ph_data(6:16), sem_Ser5ph_data(6:16), 'o','MarkerSize',5,'MarkerFaceColor',[.2,.2,.2],'Color',[.2,.2,.2]);

    subplot(2,2,3);
    b2=errorbar( x, mu_mRNA, sem_mRNA,'s','MarkerSize',5,'MarkerFaceColor','b','Color','b');
    hold on;
    axis([-5.5 5.5 -0.4 0.7])
    plot( x2, fit_mRNA(x2),'-','Color',[0,0,1]);
    b1=errorbar (x, mean_mRNA_data(6:16), sem_mRNA_data(6:16), 'o','MarkerSize',5,'MarkerFaceColor',[.2,.2,.2],'Color',[.2,.2,.2]);
    hold('off')


    subplot(2,2,4);
    b2=errorbar( x, mu_Ser5ph, sem_Ser5ph,'s','MarkerSize',5,'MarkerFaceColor','g','Color','g');
    hold on;
    b1=errorbar ( x, mu_Unph, sem_Unph, 'o','MarkerSize',5,'MarkerFaceColor','r','Color','r');
    b2=errorbar( x, mu_mRNA, sem_mRNA,'s','MarkerSize',5,'MarkerFaceColor','b','Color','b');

    axis([-5.5 5.5 -0.4 0.7])
    h2 = plot(x2, fit_mRNA(x2),'-','Color',[0,0,1]);
    h2 = plot(x2, fit_Ser5ph(x2),'-','Color',[0,1,0]);
    h1 = plot(x2, fit_Unph(x2),'-','Color',[1,0,0]);



    %figure(5)

    %plot(cumsum(mean_Unph_data)./sum(mean_Unph_data),'r'); hold on; plot(cumsum(mean_Ser5ph_data)./sum(mean_Ser5ph_data),'g');plot(cumsum(mean_mRNA_data)./sum(mean_mRNA_data),'b');


    %figure(6)
    mrna_data_cdf = [];
    mrna_cdf = [];
    ser5_cdf = [];
    ctd_cdf = [];
    
    for i = 0:1/25:1
        
        mrna_cdf = [mrna_cdf, sum(mu_mRNA <= i )/length(mu_mRNA)];
        ser5_cdf = [ser5_cdf, sum(mu_Ser5ph <= i )/length(mu_Ser5ph)];
        ctd_cdf = [ctd_cdf, sum(mu_Unph <= i )/length(mu_Unph)];

    end

    mrna_sim_cdf = [mrna_sim_cdf; mrna_cdf];
    ser5_sim_cdf = [ser5_sim_cdf; ser5_cdf];
    ctd_sim_cdf = [ctd_sim_cdf; ctd_cdf];
    
    mu_sim_mrna = [mu_sim_mrna, mu_mRNA];
    mu_sim_ser5 = [mu_sim_ser5, mu_Ser5ph];
    mu_sim_ctd = [mu_sim_ctd, mu_Unph];
    
    
end
%%
figure(11)
subplot(3,2,1)
colors = viridis(length([.6:.05:.95]));
for j = 1:length([.6:.05:.95])
    plot([0:1/25:1],mrna_sim_cdf(j,:), 'Color', colors(j,:))
    if j == 1
        hold on;
    end
end
title('mRNA CDF')

subplot(3,2,2)
for j = 1:length([.6:.05:.95])
    plot([-5:5],mu_sim_mrna(:,j), 'Color', colors(j,:))
    if j == 1
        hold on;
    end
end
title('mRNA minima')

subplot(3,2,3)
colors = viridis(length([.6:.05:.95]));
for j = 1:length([.6:.05:.95])
    plot([0:1/25:1],ser5_sim_cdf(j,:), 'Color', colors(j,:))
    if j == 1
        hold on;
    end
end
title('Ser5ph CDF')
subplot(3,2,4)
for j = 1:length([.6:.05:.95])
    plot([-5:5],mu_sim_ser5(:,j), 'Color', colors(j,:))
    if j == 1
        hold on;
    end
end
title('Ser5ph minima')

subplot(3,2,5)
colors = viridis(length([.6:.05:.95]));
for j = 1:length([.6:.05:.95])
    plot([0:1/25:1],ctd_sim_cdf(j,:), 'Color', colors(j,:))
    if j == 1
        hold on;
    end
end
xlabel('+- Tau')
title('CTD CDF')
subplot(3,2,6)
for j = 1:length([.6:.05:.95])
    plot([-5:5],mu_sim_ctd(:,j), 'Color', colors(j,:))
    if j == 1
        hold on;
    end
end

xlabel('+- Tau')
title('CTD minima')
legend(string([.6:.05:.95]))
