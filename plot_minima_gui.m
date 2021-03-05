
x= -rangeToDisplay:rangeToDisplay;
b1=errorbar (app.minima_plot, x, mu_Unph, sem_Unph, 'o','MarkerSize',5,'MarkerFaceColor','r','Color','r');
hold(app.minima_plot,'on')

%axis(app.minima_plot,[-5.5 5.5 -0.4 0.7])
x2 = linspace(-5.5,5,100);
plot(app.minima_plot,x2,fit_Unph(x2),'-','Color',[1,0,0]);

b1=errorbar (app.minima_plot, x, mean_Unph_data(6:16), sem_Unph_data(6:16), 'o','MarkerSize',5,'MarkerFaceColor',[.2,.2,.2],'Color',[.2,.2,.2]);

hold(app.minima_plot,'off')

b2=errorbar(app.minima_plot_2, x, mu_Ser5ph, sem_Ser5ph,'s','MarkerSize',5,'MarkerFaceColor','g','Color','g');
hold(app.minima_plot_2,'on')
%axis(app.minima_plot_2,[-5.5 5.5 -0.4 0.7])
plot( app.minima_plot_2,x2, fit_Ser5ph(x2),'-','Color',[0,1,0]);
b1=errorbar (app.minima_plot_2, x, mean_Ser5ph_data(6:16), sem_Ser5ph_data(6:16), 'o','MarkerSize',5,'MarkerFaceColor',[.2,.2,.2],'Color',[.2,.2,.2]);
hold(app.minima_plot_2,'off')

b2=errorbar(app.minima_plot_3, x, mu_mRNA, sem_mRNA,'s','MarkerSize',5,'MarkerFaceColor','b','Color','b');
hold(app.minima_plot_3,'on')
%axis(app.minima_plot_3,[-5.5 5.5 -0.4 0.7])
plot(app.minima_plot_3, x2, fit_mRNA(x2),'-','Color',[0,0,1]);
b1=errorbar (app.minima_plot_3, x, mean_mRNA_data(6:16), sem_mRNA_data(6:16), 'o','MarkerSize',5,'MarkerFaceColor',[.2,.2,.2],'Color',[.2,.2,.2]);
hold(app.minima_plot_3,'off')



b2=errorbar(app.minima_plot_5, x, mu_Ser5ph, sem_Ser5ph,'s','MarkerSize',5,'MarkerFaceColor','g','Color','g');
hold(app.minima_plot_5,'on')
b1=errorbar (app.minima_plot_5, x, mu_Unph, sem_Unph, 'o','MarkerSize',5,'MarkerFaceColor','r','Color','r');
b2=errorbar(app.minima_plot_5, x, mu_mRNA, sem_mRNA,'s','MarkerSize',5,'MarkerFaceColor','b','Color','b');

%axis(app.minima_plot_5, [-5.5 5.5 -0.4 0.7])
h2 = plot(app.minima_plot_5,x2, fit_mRNA(x2),'-','Color',[0,0,1]);
h2 = plot(app.minima_plot_5,x2, fit_Ser5ph(x2),'-','Color',[0,1,0]);
h1 = plot(app.minima_plot_5,x2, fit_Unph(x2),'-','Color',[1,0,0]);
hold(app.minima_plot_5,'off')




% legend([b1,h1,b2,h2,b3,h3],{'Unph','Fit Data','Ser5ph','Fit Data','mRNA','Fit Data'},'FontSize',10,'Location','SouthEast')
% axis([-5.5 5.5 -0.4 0.7])
% xlabel({'time (min)'},'FontSize',12,'FontWeight','bold')
% ylabel({'Norm. Int (a.u)'},'FontSize',12,'FontWeight','bold')
% set (gca ,'FontSize',12,'FontName', 'Arial');


timeFluc = linspace(0,size(mRNAFluc_Rave,1),size(mRNAFluc_Rave,1));


SelectedCell =1;
for n = length(vector_minima_mRNA)
   if length(vector_minima_mRNA{n})>0
       SelectedCell =n;
   end
    
end


hold(app.minima_plot_4,'off')
plot(app.minima_plot_4, timeFluc, UnphFluc_Rave(:,SelectedCell),'ro','markersize', 3,'linestyle','-','linewidth',1.5);
hold(app.minima_plot_4,'on')

xlabel(app.minima_plot_4,{'time (min)'},'FontSize',8,'FontWeight','bold')
ylabel(app.minima_plot_4, {'Norm. Int (a.u)'},'FontSize',8,'FontWeight','bold')

plot(app.minima_plot_4, timeFluc, Ser5phFluc_Rave(:,SelectedCell),'gs','markersize', 3,'linestyle','-','linewidth',1.5);
plot(app.minima_plot_4, timeFluc, mRNAFluc_Rave(:,SelectedCell),'bd','markersize', 3,'linestyle','-','linewidth',1.5);
plot(app.minima_plot_4, [0,200],[1-app.detectionthresholdEditField.Value,1-app.detectionthresholdEditField.Value],'--')

plot(app.minima_plot_4, vector_minima_mRNA{SelectedCell},mRNAFluc_Rave(vector_minima_mRNA{SelectedCell},SelectedCell), 'ko','markersize',10,'MarkerFaceColor','k');


hold(app.minima_plot_4,'off')



