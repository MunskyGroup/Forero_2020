function err = get_log_l_simplified_gui(parameters,mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap)
%parameters = par_fixed
%parameters

%parameters(par_change) = pars;
[sigred,~,means,SIGc] = get_ac_and_cc_mod_simplified(parameters);
err(1) =sum((sigred(1,1:21)'-rnap.mn_ac(1:21)).^2./rnap.sem_ac(1:21).^2)/2;
err(2) =sum((sigred(5,1:21)'-ser5.mn_ac(1:21)).^2./ser5.sem_ac(1:21).^2)/2;
err(3) =sum((sigred(9,1:21)'-mrna.mn_ac(1:21)).^2./mrna.sem_ac(1:21).^2)/2;

err(4) = sum(([sigred(4,8:-1:2),sigred(2,1:8)]'-ser5_rnap.mn_cc(4:18)).^2./ser5_rnap.sem_cc(4:18).^2)/2;
err(5) = sum(([sigred(7,8:-1:2),sigred(3,1:8)]'-mrna_rnap.mn_cc(4:18)).^2./mrna_rnap.sem_cc(4:18).^2)/2;
err(6) = sum(([sigred(8,8:-1:2),sigred(6,1:8)]'-mrna_ser5.mn_cc(4:18)).^2./mrna_ser5.sem_cc(4:18).^2)/2;

% match average number of mRNA
mean_mrna = 15.5;
sem_mean_mrna = 0.93;
err(7) = (means(3)-mean_mrna)^2/sem_mean_mrna^2/2;

% match variance in number of mRNA
var_dat = sem_mean_mrna^2*130;
err(8) = (SIGc(3,3)-var_dat)^2/10^2/2;

% frac should be between 0 and 1
% err(9) = 100*max(0,log(parameters(7)));

% err(9) = sum(max(0,log(parameters)-3));
% covariances
% cv12 = 1/parameters(11);
% cv13 = 1/parameters(12);
% cv23 = 1/parameters(13);
cv12 = 1/parameters(11)*sigred(2,1);
cv13 = 1/parameters(12)*sigred(3,1);
cv23 = 1/parameters(13)*sigred(6,1);
sig_d = [9.7339e-01   3.4716e-01   5.0248e-01];
sig_d_sem = [1.4262e-01   8.6887e-02   1.1456e-01];
%     err(13) = sum(abs(log([cv12,cv13,cv23]./sig_d)));
err(13) = 0*sum(([cv12,cv13,cv23]-sig_d).^2./sig_d_sem.^2)/2;

end