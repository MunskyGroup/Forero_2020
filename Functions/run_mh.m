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

delta = 0.03*par_opt;
proprnd = @(x)(x+delta.*randn(size(x)).*(randn(size(x))>0.5));
% delta = 0.5;
% load R_sim_for_Met_Hast
% proprnd = @(x)(x+delta*randn(size(x))*R);
parnames = {'kon','na','kesc','kproc','beta','kout','na','eta_ctd','eta_ser5','eta_mrna',...
    'sc_ctd','sc_ser5','sc_mrna'};
% parnames = {'kin','kout','kinit','kabort','kesc','kdephos','ke','kon',...
%     'koff','kproc'};
nsamples = 5001;
nchains = 20;
nsegments = 50;
thin = 20;

seeds = ceil(linspace(0,100,10));
parfor i=1:nchains    
    sv_file = ['met_hast_pars_splitE_',num2str(i)];
    run_chain(i,par_opt,sv_file,get_err,proprnd,thin,nsamples,nsegments) 
end

return
%%
load best_simple_pars
par_changed = [1,3:6];
parnames = {'omega','na','k esc','k complete','beta','k out','na','eta_ctd','eta_ser5','eta_mrna',...
    'sc_ctd','sc_ser5','sc_mrna'};
par_opt = (parameters(par_changed));
nchains = 20;
nsegments = 100;

figure(1);clf;
figure(2);clf;
figure(3);clf;
figure(4);clf;
SPS =[];
SPv =[];
for i=1:nchains 
    try
        mh_smpl{i} = [];
        mh_value{i} = [];
        sv_file = ['met_hast_pars_splitE_',num2str(i)];
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
        title(num2str(mean(mh_value{i})))
       
        figure(3)
        subplot(4,5,i)
        histogram(mh_value{i}); hold on        
        set(gca,'xlim',[-30,-14])
        title(num2str(mean(mh_value{i})))
        
        figure(4)
        subplot(4,5,i); hold off
        C = xcorr(mh_value{i}-mean(mh_value{i}),mh_value{i}-mean(mh_value{i}),'coeff');
        plot(C(floor(length(C)/2):end)); hold on
        plot([0,50000],[0,0],'k--')
        set(gca,'xlim',[0,50000],'ylim',[-0.3 1])
      
    catch
    end
    %     figure(1+i)
    %     for j = 1:5
    %         subplot(2,3,j); histogram(mh_smpl{i}(:,j));
    %     end
    if mean(mh_value{i})>-19.5
        SPS = [SPS;mh_smpl{i}];
        SPv = [SPv;mh_value{i}];
    end
end
figure(2); clf
for i=1:5  
    subplot(5,5,(i-1)*5+i); 
    histogram(SPS(:,i),20);
    ylabel(parnames(par_changed(i)));
    A = sort(SPS(:,i)); B = [A(ceil(length(A)/10)),A(floor(9*length(A)/10))];
    title({parnames{par_changed(i)},['[',num2str(B),']']});

  
for j = 1:i-1  
        subplot(5,5,(i-1)*5+j); hold off
        scatter(SPS(1:500:end,j),SPS(1:500:end,i),4,SPv(1:500:end));
        hold on
        plot(par_opt(j),par_opt(i),'ko','markersize',10,'MarkerFaceColor','k')
end
end
subplot(5,5,5);
histogram(SPv,20);        
        set(gca,'xlim',[-30,-14])




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

