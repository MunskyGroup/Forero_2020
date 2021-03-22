
load('best_simple_pars_2s_justmRNA.mat')
kon = parameters(1);
koff = 1000;%parameters(2);
k_mrna = 1000*parameters(3);

k_decay = parameters(4);
frac = parameters(5);


eta_rnap =  parameters(6);
eta_ser5 = parameters(7);
eta_ts = parameters(8);

Nstates = 2;
Nrxns = 4;
b = zeros(Nstates,1);
b(1) = kon;

c = [[0,1];[0,1];[0,1]];
x0 = [0,0]';

S =   [[ 1    -1     0     0    ];
     [0     0     1    -1    ]];
 
 W1 =[[-kon         0       ];
    [koff         0         ];
    [k_mrna         0         ];
     [    0    k_decay    ];];

W0 = zeros(Nrxns,1);
W0(1,1) = kon;
% W1 = zeros(Nrxns,Nstates);
% W0 = zeros(Nrxns,1);
% S(1,1) = 1;  W1(1,1) = -kon; W0(1,1) = kon;
% S(1,2) = -1; W1(2,1) = koff; 
% 
% S(2,3) = 1; W1(3,1) = k_mrna; 
% S(2,4) = -1;     S(3,4) = 1; W1(4,2) = k_decay; 
% 
% S(3,5) = -1; W1(5,3) = kout; 
% 
% S(3,6) = -1; S(4,6) = 1; W1(6,3) = kesc*frac; 
% 
% 
% S(4,7) =  -1; W1(7,4) = kproc;


A=S*W1;

W = @(x) W1*x + W0;
time_var = 0;
signal_update_rate = 0;


data = xlsread('Nasc_TranscriptsCount.xlsx');
mrna_hist = data(:,1);
%pol2_hist = data(:,1);

mrna_hist = mrna_hist(1:end-5);
mrnai_traj = [];
for i = 1:1000
    sol = run_single_SSA(x0,S,W,[0:1:2000],time_var,signal_update_rate);
    
    ts_ssa = sol(2,:)';
    
    mrnai = ts_ssa;
    
    mrnai_traj = [mrnai_traj, mrnai(1:200:end)];
    
  

end


%%



fntsize=12
xh = 2*3*6;
yh = 2*14.3/3;
figure(12)
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh/3, 3*yh/3]; % x,y, width, height
global_color = 'k'
[x,n] = hist(mrnai_traj,30)

x1 = histogram(mrnai_traj,n,'Normalization','probability','FaceColor',[11, 252, 3]./256,'FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256)%'DisplayStyle','stairs')
hold on;
%histogram(mrna_hist ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')
x2 = histogram(mrnai_traj,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256,'DisplayStyle','stairs','HandleVisibility','off')

x3 = histogram(mrna_hist ,n,'Normalization','probability','FaceColor',[218, 51, 255]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[218, 51, 255]./256)%,'DisplayStyle','stairs')

x4 = histogram(mrna_hist ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs','HandleVisibility','off')

%a = legend(gca,[x1,x3], {'Model','Data'},'Location','Best')
%set(a, 'Box', 'off');

set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
%set(gca,'TextColor',global_color);
set(gca,'linewidth',2)
saveas(gca, 'tshist_labels.epsc') 
%set(gca,'YTickLabel',[],'XtickLabel',[])
saveas(gca, 'tshist_nolabels.epsc') 

% title({'mRNA count'},'FontSize',fntsize,'FontWeight','bold','Color',global_color)
% xlabel({'N molecules'},'FontSize',fntsize,'FontWeight','bold')
% ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
% %ylim([minYlim maxYlim])
% set (gca ,'TickLength',[.01,.3],'LineWidth',1);
% set (gca ,'FontSize',fntsize,'FontName', 'Arial');
%set(gca,'Color','k')
legend('Model', 'Data')
xlabel('mRNA molecules')
ylabel('Probability')
