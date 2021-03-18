clear all
close all
clc
%%

addpath Data_files/ Parameter_files/

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
get_err = @(pars)-sum(get_log_l_simplified(pars,mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap,par_fixed,par_changed));

delta = 0.05*par_opt;
proprnd = @(x)(x+delta.*randn(size(x)).*(randn(size(x))>0.5));
% delta = 0.5;
% load R_sim_for_Met_Hast
% proprnd = @(x)(x+delta*randn(size(x))*R);
parnames = {'kon','na','kesc','kproc','beta','kout','na','eta_ctd','eta_ser5','eta_mrna',...
    'sc_ctd','sc_ser5','sc_mrna'};
% parnames = {'kin','kout','kinit','kabort','kesc','kdephos','ke','kon',...
%     'koff','kproc'};
nsamples = 5001;
nchains = 40;
nsegments = 20;
thin = 20;

seeds = ceil(linspace(0,100,10));
parfor i=1:nchains    
    sv_file = ['met_hast_pars_splitC_',num2str(i)];
    run_chain(i,par_opt,sv_file,get_err,proprnd,thin,nsamples,nsegments) 
end

return
%%
par_changed = [1,3:6];
parnames = {'kon','na','kesc','kproc','beta','kout','na','eta_ctd','eta_ser5','eta_mrna',...
    'sc_ctd','sc_ser5','sc_mrna'};
nchains = 40;
nsegments = 40;

close all
figure(1)
SPS =[];
for i=1:nchains 
    try
        mh_smpl{i} = [];
        mh_value{i} = [];
        sv_file = ['met_hast_pars_splitC_',num2str(i)];
        clear mh_smpl_*  mh_value_*
        load(sv_file)
        for j=1:nsegments
            try
                eval(['mh_smpl{i} = [mh_smpl{i};mh_smpl_',num2str(j),'];']);
                eval(['mh_value{i} = [mh_value{i};mh_value_',num2str(j),'];']);
            catch
            end
        end
        figure(1)
        subplot(4,5,i)
        plot(mh_value{i}); hold on
        set(gca,'ylim',[-30,-14])
        figure(3)
        subplot(4,5,i)
        histogram(mh_value{i}); hold on        
        set(gca,'xlim',[-30,-14])
        
        figure(4)
        subplot(4,5,i)
        C = xcorr(mh_value{i}-mean(mh_value{i}),mh_value{i}-mean(mh_value{i}),'coeff');
        plot(C(floor(length(C)/2):end));
        set(gca,'xlim',[0,10000],'ylim',[-0.1 1])
      
    catch
    end
    %     figure(1+i)
%     for j = 1:5
%         subplot(2,3,j); histogram(mh_smpl{i}(:,j));
%     end
    SPS = [SPS;mh_smpl{i}];
end
figure(2)
for i=1:5  
    subplot(5,5,(i-1)*5+i); 
    histogram(SPS(:,i));
    xlabel(parnames(par_changed(i)));
    ylabel(parnames(par_changed(i)));
    
  
for j = i+1:5   
        subplot(5,5,(i-1)*5+j); 
        scatter(SPS(:,j),SPS(:,i));
end
end




function run_chain(i,x0,sv_file,get_err,proprnd,thin,nsamples,nsegments)
pause(i)
rng('shuffle')
for iq = 1:nsegments
    
    [mh_smpli,accepti,mh_valuei] = MetHast(x0,nsamples,'logpdf',get_err,'proprnd', ...
        proprnd,'symmetric',1,'thin',thin);
    
    eval(['mh_smpl_',num2str(iq),'= mh_smpli(1:end-1,:);']);
    eval(['mh_value_',num2str(iq),'= mh_valuei(1:end-1,:);']);
    
    x0 = mh_smpli(end,:);
    
    accepti
    if iq==1
        save(sv_file,['mh_value_',num2str(iq)],['mh_smpl_',num2str(iq)]);
    else
        save(sv_file,['mh_value_',num2str(iq)],['mh_smpl_',num2str(iq)],'-append');
    end
end
end

