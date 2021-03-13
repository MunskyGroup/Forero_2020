clear all
close all
clc
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

%%
load best_simple_pars

par_changed = [1,3:6];
par_fixed = parameters;
par_opt = (parameters(par_changed));
get_err = @(pars)-2*sum(get_log_l_simplified(pars,mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap,par_fixed,par_changed));

delta = 0.1*par_opt;
proprnd = @(x)(x+delta.*randn(size(x)).*(randn(size(x))>0.5));
% delta = 0.5;
% load R_sim_for_Met_Hast
% proprnd = @(x)(x+delta*randn(size(x))*R);
parnames = {'kon','na','kesc','kproc','beta','kout','na','eta_ctd','eta_ser5','eta_mrna',...
    'sc_ctd','sc_ser5','sc_mrna'};
% parnames = {'kin','kout','kinit','kabort','kesc','kdephos','ke','kon',...
%     'koff','kproc'};
nsamples = 1000;l
thin = 10;

for i=1:10
    x0 = par_opt;
    sv_file = ['met_hast_pars_2x_',num2str(i)];
    mh_smpl=[x0];
    mh_value=[get_err(par_opt)];
    
    for iq = 1:100
        
        [mh_smpli,accepti,mh_valuei] = MetHast(mh_smpl(end,:),nsamples,'logpdf',get_err,'proprnd', ...
            proprnd,'symmetric',1,'thin',thin);
        mh_smpl = [mh_smpl;mh_smpli];
        mh_value = [mh_value;mh_valuei];
        accepti
        
        save(sv_file,'mh_smpl','mh_value')
        
    end
end
