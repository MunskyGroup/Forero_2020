function make_plots(parameters,mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap,sig_dat,sig_dat_sem,err,col)
[sigred,TT,means,sig] = get_ac_and_cc_mod_simplified(parameters,[0:.1:30]);

% means
% vars = diag(sig)
% sum(err)

figure(1);
subplot(3,3,1);
errorbar([0:31],rnap.mn_ac,rnap.sem_ac,'rs','linewidth',2); hold on
plot(TT,sigred(1,:),'Color',col,'linewidth',3); hold on
set(gca,'xlim',[-1 25])
set(gca,'fontsize',14);title('RNAP');
grid on

subplot(3,3,4);
errorbar([0:31],ser5.mn_ac,ser5.sem_ac,'gs','linewidth',2); hold on
plot(TT,sigred(5,:),'Color',col,'linewidth',3); hold on
set(gca,'xlim',[-1 25])
set(gca,'fontsize',14);title('SER5');
grid on

subplot(3,3,7);
errorbar([0:31],mrna.mn_ac,mrna.sem_ac,'bs','linewidth',2); hold on
plot(TT,sigred(9,:),'Color',col,'linewidth',3); hold on
set(gca,'xlim',[-1 25])
set(gca,'fontsize',14);title('mRNA');
grid on

subplot(3,3,2);
plot_cc(ser5_rnap.mn_cc,ser5_rnap.sem_cc,'c'); hold on
plot([-TT(end:-1:2);TT],[sigred(4,end:-1:2),sigred(2,:)],'Color',col,'linewidth',3); hold on
plot([0,0],[0 0.04],'k--')
set(gca,'fontsize',14);title('RNAP-SER5');
grid on
set(gca,'xlim',[-8 8])

subplot(3,3,5);
plot_cc(mrna_rnap.mn_cc,mrna_rnap.sem_cc,'k'); hold on
plot([-TT(end:-1:2);TT],[sigred(7,end:-1:2),sigred(3,:)],'Color',col,'linewidth',3); hold on
set(gca,'xlim',[-8 8])
plot([0,0],[0 0.02],'k--')
set(gca,'fontsize',14);title('RNAP-mRNA');
grid on

subplot(3,3,8);
plot_cc(mrna_ser5.mn_cc,mrna_ser5.sem_cc,'m'); hold on
plot([-TT(end:-1:2);TT],[sigred(8,end:-1:2),sigred(6,:)],'Color',col,'linewidth',3); hold on
set(gca,'xlim',[-8 8])
plot([0,0],[0 0.03],'k--')
set(gca,'fontsize',14);title('SER5-mRNA');
grid on

% subplot(3,3,6)
% cv12 = 1/parameters(11)*sigred(2,1);
% cv13 = 1/parameters(12)*sigred(3,1);
% cv23 = 1/parameters(13)*sigred(6,1);
% 
% bar([1 2.4 3.8],[sig_dat([2,3,6])-sig_dat_sem([2,3,6]);...
%     2*sig_dat_sem([2,3,6])]',.4,'stacked');hold on;
% bar([1.6 3.0 4.4],[cv12,cv13,cv23]',0.4);
% title(['Correl Coeffs, err = ',num2str(err(13))]);
% set(gca,'xtick',[1.4 2.8 4.2],'xticklabel',{'\sigma_{CS}','\sigma_{CR}','\sigma_{SR}'},...
%     'xlim',[0.5 5])
% legend({'Data-SEM','Data+SEM','Model'})

% subplot(3,3,3)
% mean_mrna = 15.5;
% sem_mean_mrna = 0.93;
% var_dat = sem_mean_mrna^2*130;
% 
% bar([1 2.4],[mean_mrna-sem_mean_mrna,var_dat-10;...
%     2*sem_mean_mrna,20]',.4,'stacked');hold on;
% bar([1.6 3.0],[means(3),sig(3,3)]',0.4);
% title(['mRNA Stats, err = ',num2str(err(7)+err(8))]);
% set(gca,'xticklabel',{'\mu','\sigma^2'})
% legend({'Data-SEM','Data+SEM','Model'})

% subplot(3,3,9)
% bar([9 10 11 12],err([9 10 11 12]))
% set(gca,'xticklabel',{'\lambda_i<100','k_{e}','k_p','prior'})
% title(['prior constaints, ', num2str(sum(err([9 10 11 12])))]);
% 
% for j=[1,2,4,5,7,8]
%     subplot(3,3,j)
%     set(gca,'ylim',[-0.05 1.4])
% end
end

function plot_cc(v,e,c)
errorbar([-10:10],v/v(11),e/v(11),['s',c],'linewidth',2)
end