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
% 
% seeds = ceil(linspace(0,100,10));
% parfor i=1:nchains    
%     sv_file = ['met_hast_pars_splitF_',num2str(i)];
%     run_chain(i,par_opt,sv_file,get_err,proprnd,thin,nsamples,nsegments) 
% end


%%
par_changed = [1,3:6];
parnames = {'kon','na','kesc','kproc','beta','kout','na','eta_ctd','eta_ser5','eta_mrna',...
    'sc_ctd','sc_ser5','sc_mrna'};
nchains = 40;
nsegments = 40;

close all
figure(1)
SPS =[];
mh_val_all = [];



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
    mh_val_all = [mh_val_all; mh_value{i}]; 
end

average_LLs = [];
for i = 1:20
    average_LLs = [average_LLs; mean(mh_value{i})];
    
end



[~,inds] = sort(-average_LLs)

% for i = 1:10
%     
%     SPS = [SPS;mh_smpl{inds(i)}];
%     mh_val_all = [mh_val_all; mh_value{inds(i)}];  
%     
% end




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


%%
% plot for figure + 80% CI


Np = 5;
mh_vals = mh_val_all;
mh_pars = SPS(:,[4,1,5,2,3]);
par_changed = [5,1,6,3,4];
%mh_pars = SPS;
sz = 2*(1+max(-mh_vals)-(-mh_vals));
figure(8)
xh = 2*3*6;
yh = 2*14.3/3;
fntsize = 18;
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

for i=1:Np
    subplot(Np,Np,(i-1)*Np+i)

    H = histogram(mh_pars(:,i),20,'Normalization','Probability'); hold on;
    %H = histogram(mh_pars(:,i),20,'Normalization','Probability'); hold on;

    X =  sort(mh_pars(:,i));
    low = X(floor(length(X)/10));
    high = X(ceil(9*length(X)/10));
    
    plot(mh_pars(1,i)*[1,1],get(gca,'ylim'),'r--','linewidth',3)

    title([parnames{par_changed(i)}])

    
    if i==1
        ylabel({parnames{par_changed(i)}},'FontSize',fntsize,'FontWeight','bold');
    elseif i==Np
        xlabel({parnames{par_changed(i)}},'FontSize',fntsize,'FontWeight','bold');
    end
    %     hold on
    %     pr = max(0,log(H.BinEdges/100));
    %     switch parnames{par_changed(i)}
    %         case 'kproc'
    %             pr = pr + log(H.BinEdges/(1/(103/60)));
    %         case 'ke'
    %             pr = pr + log(H.BinEdges/4.1);
    %     end
    
    %     [qw,qe] =max(H.Values);
    
    %     pr = max(H.Values)*exp(-pr)/exp(-pr(qe));
    %     plot(H.BinEdges,pr,'r','linewidth',2)

        set (gca ,'TickLength',[.01,.3],'LineWidth',2);
    set (gca ,'FontSize',fntsize,'FontName', 'Arial');
    for j=i+1:Np
        subplot(Np,Np,(j-1)*Np+i)
        %scatter(mh_pars(1:100:end,i),mh_pars(1:100:end,j),sz(1:100:end),mh_vals(1:100:end),'filled'); hold on
        
         scatter(mh_pars(1:100:end,i),mh_pars(1:100:end,j),sz(1:100:end),mh_vals(1:100:end),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerFaceColor',[50, 139, 191]./256,'MarkerEdgeColor',[50, 139, 191]./256); hold on
         scatter(mh_pars(1:100:end,i),mh_pars(1:100:end,j),sz(1:100:end),mh_vals(1:100:end),'filled','MarkerFaceAlpha',.04,'MarkerEdgeAlpha',.04,'MarkerFaceColor',[0,1,0]); hold on
         scatter(mh_pars(1:100:end,i),mh_pars(1:100:end,j),sz(1:100:end),mh_vals(1:100:end),'filled','MarkerFaceAlpha',.007,'MarkerEdgeAlpha',.007,'MarkerFaceColor',[1,1,0]); hold on

        plot(mh_pars(1,i),mh_pars(1,j),'ko','markersize',8,'markerfacecolor','k')
        if i==1
            ylabel({parnames{par_changed(j)}},'FontSize',fntsize,'FontWeight','bold');
        end
        if j==Np
            xlabel({parnames{par_changed(i)}} ,'FontSize',fntsize,'FontWeight','bold');
        end
        set (gca ,'TickLength',[.01,.3],'LineWidth',2);
        set (gca ,'FontSize',fntsize,'FontName', 'Arial');
    end
    
    
end






Np = 5;
mh_vals = mh_val_all;
mh_pars = SPS(:,[4,1,5,2,3]);
par_changed = [5,1,6,3,4];
%mh_pars = SPS;
sz = 2*(1+max(-mh_vals)-(-mh_vals));
figure(9)
xh = 2*3*6;
yh = 2*14.3/3;
fntsize = 18;
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height


[sorted_mh_vals, mh_val_inds] = sort(mh_vals); 
mh_par_sort = mh_pars(mh_val_inds,:);
colormap(flipud(parula))

mh_val_new = -sort(mh_vals);
mh_val_new(mh_val_new > mh_val_new(.05*(length(mh_val_new)))) = -min(mh_vals);


for i=1:Np
    subplot(Np,Np,(i-1)*Np+i)

    H = histogram(mh_pars(:,i),20,'Normalization','Probability'); hold on;
    %H = histogram(mh_pars(:,i),20,'Normalization','Probability'); hold on;

    X =  sort(mh_pars(:,i));
    low = X(floor(length(X)/10));
    high = X(ceil(9*length(X)/10));
    
    plot(mh_pars(1,i)*[1,1],get(gca,'ylim'),'r--','linewidth',3)

    title([parnames{par_changed(i)}])

    
    if i==1
        ylabel({parnames{par_changed(i)}},'FontSize',fntsize,'FontWeight','bold');
    elseif i==Np
        xlabel({parnames{par_changed(i)}},'FontSize',fntsize,'FontWeight','bold');
    end
    %     hold on
    %     pr = max(0,log(H.BinEdges/100));
    %     switch parnames{par_changed(i)}
    %         case 'kproc'
    %             pr = pr + log(H.BinEdges/(1/(103/60)));
    %         case 'ke'
    %             pr = pr + log(H.BinEdges/4.1);
    %     end
    
    %     [qw,qe] =max(H.Values);
    
    %     pr = max(H.Values)*exp(-pr)/exp(-pr(qe));
    %     plot(H.BinEdges,pr,'r','linewidth',2)

        set (gca ,'TickLength',[.01,.3],'LineWidth',2);
    set (gca ,'FontSize',fntsize,'FontName', 'Arial');
    for j=i+1:Np
        subplot(Np,Np,(j-1)*Np+i)
         scatter(mh_par_sort(1:100:end,i),mh_par_sort(1:100:end,j),sz(1:100:end),mh_val_new(1:100:end),'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on
        
%          scatter(mh_pars(1:100:end,i),mh_pars(1:100:end,j),sz(1:100:end),mh_vals(1:100:end),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerFaceColor',[50, 139, 191]./256,'MarkerEdgeColor',[50, 139, 191]./256); hold on
%          scatter(mh_pars(1:100:end,i),mh_pars(1:100:end,j),sz(1:100:end),mh_vals(1:100:end),'filled','MarkerFaceAlpha',.04,'MarkerEdgeAlpha',.04,'MarkerFaceColor',[0,1,0]); hold on
%          scatter(mh_pars(1:100:end,i),mh_pars(1:100:end,j),sz(1:100:end),mh_vals(1:100:end),'filled','MarkerFaceAlpha',.007,'MarkerEdgeAlpha',.007,'MarkerFaceColor',[1,1,0]); hold on

        plot(mh_pars(1,i),mh_pars(1,j),'ko','markersize',8,'markerfacecolor','k')
        if i==1
            ylabel({parnames{par_changed(j)}},'FontSize',fntsize,'FontWeight','bold');
            colorbar
        end
        if j==Np
            xlabel({parnames{par_changed(i)}} ,'FontSize',fntsize,'FontWeight','bold');
        end
        set (gca ,'TickLength',[.01,.3],'LineWidth',2);
        set (gca ,'FontSize',fntsize,'FontName', 'Arial');
    end
    
    
end


[4,1,5,2,3];
kon = mh_pars(:,2);
koff = 1000;
kesc = mh_pars(:,3);
kproc = mh_pars(:,5);
kin = mh_pars(:,1)*koff;
kout = mh_pars(:,4);

mu_ctd = kon./(kon+koff).*kin./(kout+kesc);
mu_ctd(1)
% mu_ctd_on = kin./(kout+kesc);

mu_rna = kon./(kon+koff).*kin./(kout+kesc).*kesc./kproc;
mu_rna(1)
% mu_rna_on = kin./(kout+kesc).*kesc./kproc;

rate_prod = mu_ctd.*kesc;
% rate_prod_on = mu_ctd_on.*kesc;

burst_refractory_period = 1./kon;
burst_duration = 1./koff;
pol2_burst_size = kin./koff;

burst_efficiency = 100*kesc./(kesc+kout);
mrna_burst_size = pol2_burst_size.*burst_efficiency/100;

figure(7);


fntsize = 18;
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height
sbplots = {'mu_ctd','mu_rna','rate_prod',...
    'burst_refractory_period',...
    'pol2_burst_size','burst_efficiency','mrna_burst_size'};
sbplots_leg = {'\mu ctd (mol)','\mu rna (mol)',...
    'rate prod (mol/min)',...
    'burst refractory period (min)',...
    'pol2 burst size (mol)','burst efficiency (%)','mrna burst size (mol)'}
for i=1:7
    subplot(3,3,i); hold off
    histogram(eval(sbplots{i}),'Normalization','probability','FaceColor',[11, 252, 3]./256,'FaceAlpha',.5,'linewidth',.5,'edgecolor',[11, 252, 3]./256 ); hold on
    plot(eval([sbplots{i},'(1)'])*[1,1],get(gca,'ylim'),'k--','linewidth',3)
    
    X =  sort(eval(sbplots{i}));
    low = X(floor(length(X)/10));
    high = X(ceil(9*length(X)/10));
    
    title([sbplots_leg{i},'; ',num2str(eval([sbplots{i},'(1)']),3),...
        ' (',num2str(low,3),', ',num2str(high,3),')'],'FontSize',fntsize-5,'FontWeight','bold')
    if mod(i,3) == 1 
        ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
    end
    if ismember(i,[7,6,5])
        xlabel({'Value'},'FontSize',fntsize,'FontWeight','bold')
    end
    
end


%%

highv = 90/100
lowv = 10/100
beta = mh_pars(:,1);
omega = mh_pars(:,2);
kout = mh_pars(:,3);
kesc = mh_pars(:,4);
kcomp = mh_pars(:,5);

tbl = [];
beta_sort = sort(beta);
low = beta_sort(floor(length(beta_sort)*lowv));
high = beta_sort(ceil(highv* length(beta_sort)));
disp([low, high])
tbl = [tbl; [low,high]];

omega_sort = sort(omega(1:100:end));
low = omega_sort(floor(length(omega_sort)*lowv));
high = omega_sort(ceil(highv* length(omega_sort)));
disp([low, high])
tbl = [tbl; [low,high]];

kout_sort = sort(kout(1:100:end));
low = kout_sort(floor(length(kout_sort)*lowv));
high = kout_sort(ceil(highv* length(kout_sort)));
disp([low, high])
tbl = [tbl; [low,high]];


kproc_sort = sort(kcomp);
low = kproc_sort(floor(length(kproc_sort)*lowv));
high = kproc_sort(ceil(highv* length(kproc_sort)));
disp([low, high])
tbl = [tbl; [low,high]];

kesc_sort = sort(kesc);
low = kesc_sort(floor(length(kesc_sort)*lowv));
high = kesc_sort(ceil(highv* length(kesc_sort)));
disp([low, high])
tbl = [tbl; [low,high]];



%%

tbl = [];
r = beta.*omega;

r =  sort(r);
low = r(floor(length(r)*lowv));
high = r(ceil(highv* length(r)));
disp([low, high])
tbl = [tbl; [low,high]];

t_clust = 1./(kesc + kout);

t_clust =  sort(t_clust);
low = t_clust(floor(length(t_clust)*lowv));
high = t_clust(ceil(highv* length(t_clust)));
disp([low, high])
tbl = [tbl; [low,high]];

u_clust = (beta.*omega)./(kesc + kout);

u_clust =  sort(u_clust);
low = u_clust(floor(length(u_clust)*lowv));
high = u_clust(ceil(highv* length(u_clust)));
disp([low, high])
tbl = [tbl; [low,high]];

f = kesc./(kesc + kout);

f =  sort(f);
low = f(floor(length(f)*lowv));
high = f(ceil(highv* length(f)));
disp([low, high])
tbl = [tbl; [low,high]];

bmrna = kesc./(kesc + kout).*beta;

bmrna =  sort(bmrna);
low = bmrna(floor(length(bmrna)*lowv));
high = bmrna(ceil(highv* length(bmrna)));
disp([low, high])
tbl = [tbl; [low,high]];

rmrna = ((beta.*omega)./(kesc + kout)).*kesc;

rmrna =  sort(rmrna);
low = rmrna(floor(length(rmrna)*lowv));
high = rmrna(ceil(highv* length(rmrna)));
disp([low, high])
tbl = [tbl; [low,high]];

rmrna = ((beta.*omega)./(kesc + kout)).*(kesc ./kcomp)  ;

rmrna =  sort(rmrna);
low = rmrna(floor(length(rmrna)*lowv));
high = rmrna(ceil(highv* length(rmrna)));
disp([low, high])
tbl = [tbl; [low,high]];


totalrnap = ((beta.*omega)./(kesc + kout)).*(kesc ./kcomp) + ((beta.*omega)./(kesc + kout))  ;

totalrnap =  sort(totalrnap);
low = totalrnap(floor(length(totalrnap)*lowv));
high = totalrnap(ceil(highv* length(totalrnap)));
disp([low, high])
tbl = [tbl; [low,high]];

mrnat = 1./kcomp;
mrnat =  sort(mrnat);
low = mrnat(floor(length(mrnat)*lowv));
high = mrnat(ceil(highv* length(mrnat)));
disp([low, high])
tbl = [tbl; [low,high]];


%%
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

