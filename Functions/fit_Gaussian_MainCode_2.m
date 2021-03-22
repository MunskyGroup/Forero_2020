function [fit_mRNA, mu_mRNA, sem_mRNA,fit_Ser5ph, mu_Ser5ph, sem_Ser5ph, fit_Unph, mu_Unph, sem_Unph]= fit_Gaussian_MainCode_2 (mean_mRNA, mean_Ser5ph, mean_Unph, sem_mRNA, sem_Ser5ph, sem_Unph, rangeToDisplay)

[fit_mRNA, mu_mRNA, sem_mRNA] = fit_Gaussian_MainCode (mean_mRNA, sem_mRNA, rangeToDisplay);
[fit_Ser5ph , mu_Ser5ph, sem_Ser5ph] = fit_Gaussian_MainCode (mean_Ser5ph, sem_Ser5ph, rangeToDisplay);
[fit_Unph, mu_Unph, sem_Unph]= fit_Gaussian_MainCode (mean_Unph, sem_Unph, rangeToDisplay);

x= -rangeToDisplay:rangeToDisplay;

%% Plotting valleys to Unph with their Gaussian fits
%figure('visible', 'off');
return
h =  findobj('type','figure');
n_Fig = length(h);
n_Fig = n_Fig+1;

figure(n_Fig);clf;
fig12= gcf;
fig12.PaperUnits = 'centimeters';
fig12.PaperPosition = [0, 0, 9.5, 8]; % x,y, width, height
hold on;
b1=errorbar (x, mu_Unph, sem_Unph, 'o','MarkerSize',5,'MarkerFaceColor','r','Color','r');
h1 = plot( fit_Unph,'-');
h1.Color=[1 0 0];
b2=errorbar(x, mu_Ser5ph, sem_Ser5ph,'s','MarkerSize',5,'MarkerFaceColor','g','Color','g');
h2 = plot( fit_Ser5ph,'-');
h2.Color=[0 1 0];
b3=errorbar(x, mu_mRNA, sem_mRNA,'d','MarkerSize',5,'MarkerFaceColor','b','Color','b');hold on;
h3 = plot( fit_mRNA,'-');
h3.Color=[0 0 1];
legend([b1,h1,b2,h2,b3,h3],{'Unph','Fit Data','Ser5ph','Fit Data','mRNA','Fit Data'},'FontSize',10,'Location','SouthEast')
axis([-5.5 5.5 -0.4 0.7])
xlabel({'time (min)'},'FontSize',12,'FontWeight','bold')
ylabel({'Norm. Int (a.u)'},'FontSize',12,'FontWeight','bold')
set (gca ,'FontSize',12,'FontName', 'Arial');
%print('ValleystoUnph','-dpng')
%saveas(gcf,'ValleystoUnph.epsc')%Save the fig in eps format for editing later in Corel


end