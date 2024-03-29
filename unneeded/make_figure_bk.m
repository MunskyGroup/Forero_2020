
%% Make all figures for Forero 2020



function make_figure(paper_fig_num, save_eps)
figure_names = {'3a','3d','3f','sup7','sup10','sup8a','sup8b'};

if ~any(strcmp(figure_names,paper_fig_num))
    disp('not a valid figure number, use one of the following')
    disp(figure_names)
    return
end


%Final figure code that generates all computational model figures used for
%Forero2020 

%Set up figure sizes
%clc; close all; clear all;
X_SIZE = 13; Y_SIZE = 15;
figure(1);clf;
fig1= gcf;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0, 0, X_SIZE, Y_SIZE]; % x,y, width, height
global_color = 'k';
    fntsize = 18;
    xh = 2*3*6;
    yh = 2*14.3/3;
    
subplot(3,1,1)
norm_intensity = 0;
minYlim = -0.1;
maxYlim = 5;
yyaxis left
addpath('./Data_files')
addpath('./Model_files')
addpath('./Parameter_files')
addpath('./Functions')
addpath('./Model_fits/simple_MH_results/3.17.2021')
rng(0)

%% Load in data from excel files and normalize them

norm_signal = 1; %normalized signal is true
seednum = floor(rand*1000);  %get random seed

%Load Cross Correlations with G0 normalization
[~,~,~,mrna_rnap,mrna_ser5,ser5_rnap] = load_normalization_variance_norm_cc(0);  %load normalization 21 pts for cross correlations
%Load Autocorrelations with G0_intp normalization and no rezeroing 
[mrna,ser5,rnap,~,~,~] = load_normalization_variance(1,'G0_intp','none',10);


load('best_simple_pars.mat','parameters')  %load parameter file
real_valscl = parameters;




%% Generate analytical correlation signals

[sigred,TT,means,sig,TT2,mins] = get_ac_and_cc_mod_simplified(parameters,[0:.1:30]);

parnames = {'kon','kesc','kproc','beta','kout'};
par_changed = [1:5];

mh_pars = [];   %load in the methaste sampling data for plots
mh_vals = [];
ikeep = 1;
i=1; 


% 
% cd Model_fits\simple_MH_results\3.24.2020\
% while ikeep==1
%     try
%         fn = ['met_hast_pars_2x_',num2str(i),'.mat'];
%         load(fn)
%         mh_pars = [mh_pars;mh_smpl];
%         mh_vals = [mh_vals;mh_value];
%         i=i+1;
%     catch
%         ikeep=0;
%     end    
% end




% cd ..\..\..

rng(0)
par_changed = [1,3:6];
parnames = {'kon','na','kesc','kproc','beta','kout','na','eta_ctd','eta_ser5','eta_mrna',...
    'sc_ctd','sc_ser5','sc_mrna'};
nchains = 20;
nsegments = 50;

close all
figure(1)
SPS =[];
mh_val_all = [];



if strcmp(paper_fig_num,'3a')
    for i=1:nchains 
        %try
            mh_smpl{i} = [];
            mh_value{i} = [];
            sv_file = ['./Model_fits/simple_MH_results/3.17.2021/met_hast_pars_splitE_',num2str(i)];
            clear mh_smpl_*  mh_value_*
            tmp_dat = load(sv_file);



            for j=1:nsegments
                %try
                    eval(['mh_smpl{i} = [mh_smpl{i};tmp_dat.mh_smpl_',num2str(j),'];']);
                    eval(['mh_value{i} = [mh_value{i};tmp_dat.mh_value_',num2str(j),'];']);
                %catch
                %end
            end
    %         figure(1)
    %         subplot(4,5,i)
    %         plot(mh_value{i}); hold on
    %         set(gca,'ylim',[-30,-14])
    %         figure(3)
    %         subplot(4,5,i)
    %         histogram(mh_value{i}); hold on        
    %         set(gca,'xlim',[-30,-14])
    %         
    %         figure(4)
    %         subplot(4,5,i)
    %         C = xcorr(mh_value{i}-mean(mh_value{i}),mh_value{i}-mean(mh_value{i}),'coeff');
    %         plot(C(floor(length(C)/2):end));
    %         set(gca,'xlim',[0,10000],'ylim',[-0.1 1])
    %       
%         catch
%         end
        %     figure(1+i)
    %     for j = 1:5
    %         subplot(2,3,j); histogram(mh_smpl{i}(:,j));
    %     end
        SPS = [SPS;mh_smpl{i}];
        mh_val_all = [mh_val_all; mh_value{i}]; 
    end

    mh_pars = SPS;

    par_fixed = parameters;
    par_changed = [1,3:6];
    parameters_samp = parameters;

    lags1 = [0:31];   
    %%

    fntsize = 18;   %set up figure 3
    close all
    figure(1)
    subplot(4,3,2)
    fig1= gcf;
    %set(gcf,'color','k');
    fig1.PaperUnits = 'centimeters';
    xh = 2*3*6;
    yh = 2*14.3/3;
    fig1.PaperPosition = [0, 0, 2/3*xh, 2*yh]; % x,y, width, height
    %sem_rnap = std(Xdata_sem(:,1:30),1)/(sqrt(20));




    plot(TT, sigred(1,:),'k-','linewidth',2)  %plot model rnap_acc
    hold on;
    %plot the data rnap_acc
    b1=errorbar(lags1(1:25), rnap.mn_ac(1:25), rnap.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[256, 0, 0]./256,'Color',[256, 0, 0]./256);hold on;

    par_fixed = parameters;
    par_changed = [1,3:6];
    parameters_samp = parameters;

    minvec = sigred(1,:);
    maxvec = sigred(1,:);   %randomly sample 50 mh par lines to generate gray lines
    rng(0)
    for i=1:1000
        n = ceil(rand*size(mh_pars,1));


        par_samp = mh_pars(n,:);
         parameters_samp(par_changed) = par_samp;
        [sigred_samp,TT_samp,means_samp,sig_samp,TT2_samp,mins_samp] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
        %plot(TT_samp, sigred_samp(1,:),'-','Color',[.2,.2,.2],'linewidth',2)

        for j = 1:length(minvec)               %keep the minimum and maximum points
            if sigred_samp(1,j) > maxvec(1,j)
                maxvec(1,j) = sigred_samp(1,j);

            end    
             if sigred_samp(1,j) < minvec(1,j)
                minvec(1,j) = sigred_samp(1,j);

            end           
        end

    end


    %fill in the points between min and max MH sampled signals with gray
    fill([TT',fliplr(TT')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);
    b1=errorbar(lags1(1:25), rnap.mn_ac(1:25), rnap.sem_ac(1:25),'o','MarkerSize',5,'MarkerFaceColor',[256, 0, 0]./256,'Color',[256, 0, 0]./256);hold on;

    plot(TT, sigred(1,:),'k-','linewidth',2) %replot so model acc signal is on top

    % Setup labels and limits 
    xlim([-.5,25])
    plot([0,25],[0,0],'k--')
    ylim([-.1,1.3])
    set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color','Ycolor',global_color);
    title({'\color{black}Autocorrelations'},'FontSize',fntsize,'FontWeight','bold')
    set (gca ,'TickLength',[.01,.3],'LineWidth',1);
    title({'CTD'},'FontSize',fntsize,'FontWeight','bold')


    %%


    subplot(4,3,5)      %SER5 ACC figure
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

    rnap_color = [138, 25, 38]./256;
    ser5_color = [20, 122, 47]./256;

    plot(TT,sigred(5,:),'-k','linewidth',2)  %plot model ser5 acc
    hold on;
    % plot data ser5 acc
    b1=errorbar(lags1(1:25), ser5.mn_ac(1:25), ser5.sem_ac(1:25),'s','MarkerSize',5,'MarkerFaceColor',[0, 256, 0]./256,'Color',[0, 256, 0]./256);hold on;


    par_fixed = parameters;
    par_changed = [1,3:6];
    parameters_samp = parameters;

    pt = 5;
    minvec = sigred(pt,:);
    maxvec = sigred(pt,:);
    rng(0)
    for i=1:1000   %sample and get min / max from 50 random mh signals
        n = ceil(rand*size(mh_pars,1));

        parameters_samp(par_changed) = par_samp;
        par_samp = mh_pars(n,:);
         parameters_samp(par_changed) = par_samp;
        [sigred_samp,TT_samp,means_samp,sig_samp,TT2_samp,mins_samp] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
        %plot(TT_samp, sigred_samp(5,:),'-','Color',[.8,.8,.8],'linewidth',2)

        for j = 1:length(sigred_samp)
            if sigred_samp(pt,j) > maxvec(1,j)
                maxvec(1,j) = sigred_samp(pt,j);

            end    
             if sigred_samp(pt,j) < minvec(1,j)
                minvec(1,j) = sigred_samp(pt,j);

            end           
        end

    end
    %fill in the gray area and replot other signals on top
    fill([TT',fliplr(TT')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);
    b1=errorbar(lags1(1:25), ser5.mn_ac(1:25), ser5.sem_ac(1:25),'s','MarkerSize',5,'MarkerFaceColor',[0, 256, 0]./256,'Color',[0, 256, 0]./256);hold on;
    plot(TT,sigred(5,:),'-k','linewidth',2)

    xlim([-.5,25])   %figure labels and limits
    plot([0,25],[0,0],'k--')
    ylim([-.1,1.3])
    title({'Ser5ph'},'FontSize',fntsize,'FontWeight','bold')
    set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
    set (gca ,'TickLength',[.01,.3],'LineWidth',1);

    %%
    subplot(4,3,8)   %mRNA ACC figure
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

    rnap_color = [138, 25, 38]./256;
    ser5_color = [20, 122, 47]./256;
    mrna_color = [22, 68, 148]./256;

    plot(TT,sigred(9,:),'-k','linewidth',2)  %model mrna ACC
    hold on;
    %data mrna ACC
    b1=errorbar(lags1(1:25), mrna.mn_ac(1:25), mrna.sem_ac(1:25),'d','MarkerSize',5,'MarkerFaceColor',[0, 0, 256]./256,'Color',[0, 0, 256]./256);hold on;

    par_fixed = parameters;
    par_changed = [1,3:6];
    parameters_samp = parameters;

    pt = 9;
    minvec = sigred(pt,:);
    maxvec = sigred(pt,:);
    rng(0)
    for i=1:1000   %get mh min and max for gray fill
        n = ceil(rand*size(mh_pars,1));

        parameters_samp(par_changed) = par_samp;
        par_samp = mh_pars(n,:);
         parameters_samp(par_changed) = par_samp;
        [sigred_samp,TT_samp,means_samp,sig_samp,TT2_samp,mins_samp] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
        %plot(TT_samp, sigred_samp(9,:),'-','Color',[.8,.8,.8],'linewidth',2)
        for j = 1:length(sigred_samp)
            if sigred_samp(pt,j) > maxvec(1,j)
                maxvec(1,j) = sigred_samp(pt,j);

            end    
             if sigred_samp(pt,j) < minvec(1,j)
                minvec(1,j) = sigred_samp(pt,j);

            end           
        end

    end

    %fill the min max mh pars, and replot signals on top
    fill([TT',fliplr(TT')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);
    b1=errorbar(lags1(1:25), mrna.mn_ac(1:25), mrna.sem_ac(1:25),'d','MarkerSize',5,'MarkerFaceColor',[0, 0, 256]./256,'Color',[0, 0, 256]./256);hold on;
    plot(TT,sigred(9,:),'-k','linewidth',2)

    xlim([-.5,25])  %Limits and labels
    plot([0,25],[0,0],'k--')
    ylim([-.1,1.3])
    title({'mRNA'},'FontSize',fntsize,'FontWeight','bold')
    set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
    set (gca ,'TickLength',[.01,.3],'LineWidth',1);

    saveas(gca,'./Figures/Accs.epsc')  %save the autocorrelations

    %%
    lags2 = [-10:10]; %time vector for CCs

    subplot(4,3,3)
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

    % Ser5 RNAP cross correlation
    plot([-TT(end:-1:2);TT],[sigred(4,end:-1:2),sigred(2,:)],'k-','linewidth',2) %model
    hold on;
    %data 
    b1=errorbar(lags2, ser5_rnap.mn_cc, ser5_rnap.sem_cc,'s','MarkerSize',5,'MarkerFaceColor',[0, 0.67, 1],'Color',[0, 0.67, 1]);hold on;

    par_fixed = parameters;
    par_changed = [1,3:6];
    parameters_samp = parameters;


    minvec = [sigred(4,end:-1:2),sigred(2,:)];  %pull the mh min max for fill of gray 
    maxvec = [sigred(4,end:-1:2),sigred(2,:)];
    rng(0)

    for i=1:1000
        n = ceil(rand*size(mh_pars,1));

        parameters_samp(par_changed) = par_samp;
        par_samp = mh_pars(n,:);
         parameters_samp(par_changed) = par_samp;
        [sigred_samp,TT_samp,means_samp,sig_samp,TT2_samp,mins_samp] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
        %plot([-TT_samp(end:-1:2);TT_samp],[sigred_samp(4,end:-1:2),sigred_samp(2,:)],'-','Color',[.8,.8,.8],'linewidth',2)

        svec = [sigred_samp(4,end:-1:2),sigred_samp(2,:)];
        for j = 1:length(svec)
            if svec(1,j) > maxvec(1,j)
                maxvec(1,j) = svec(1,j);

            end    
             if svec(1,j) < minvec(1,j)
                minvec(1,j) = svec(1,j);

            end           
        end
    end

    fill([[-TT(end:-1:2);TT]',fliplr([-TT(end:-1:2);TT]')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);

    b1=errorbar(lags2, ser5_rnap.mn_cc, ser5_rnap.sem_cc,'s','MarkerSize',5,'MarkerFaceColor',[0, 0.67, 1],'Color',[0, 0.67, 1]);hold on;

    plot([-TT(end:-1:2);TT],[sigred(4,end:-1:2),sigred(2,:)],'k-','linewidth',2) %replot signals
    plot([0,0],[0,1.3],'--','Color',global_color)

    %plotting labels and legends and limits
    set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
    title({'\color{black}Cross-correlations'},'FontSize',20,'FontWeight','bold')
    title({'Ser5ph-CTD'},'FontSize',fntsize,'FontWeight','bold')
    xlim([-10,10])
    ylim([-.1,1.3]) 
    set (gca ,'TickLength',[.01,.3],'LineWidth',1);

    %%
    subplot(4,3,6)
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

    plot([-TT(end:-1:2);TT],[sigred(7,end:-1:2),sigred(3,:)],'k-','linewidth',2) % mrna_rnap CC model
    hold on;  %mrna_rnap CC data
    b1=errorbar(lags2, mrna_rnap.mn_cc, mrna_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[1, .5, 0],'Color',[1, .5, 0]);hold on;

    par_fixed = parameters;
    par_changed = [1,3:6];
    parameters_samp = parameters;

    minvec = [sigred(7,end:-1:2),sigred(3,:)];
    maxvec = [sigred(7,end:-1:2),sigred(3,:)];

    rng(0)

    for i=1:1000  %load min and max mh_pars for gray bars
        n = ceil(rand*size(mh_pars,1));

        parameters_samp(par_changed) = par_samp;
        par_samp = mh_pars(n,:);
         parameters_samp(par_changed) = par_samp;
        [sigred_samp,TT_samp,means_samp,sig_samp,TT2,mins] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
       % plot([-TT_samp(end:-1:2);TT_samp],[sigred_samp(7,end:-1:2),sigred_samp(3,:)],'-','Color',[.8,.8,.8],'linewidth',2)

        svec = [sigred_samp(7,end:-1:2),sigred_samp(3,:)];
        for j = 1:length(svec)
            if svec(1,j) > maxvec(1,j)
                maxvec(1,j) = svec(1,j);

            end    
             if svec(1,j) < minvec(1,j)
                minvec(1,j) = svec(1,j);

            end           
        end
    end

    fill([[-TT(end:-1:2);TT]',fliplr([-TT(end:-1:2);TT]')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);
    b1=errorbar(lags2, mrna_rnap.mn_cc, mrna_rnap.sem_cc,'o','MarkerSize',5,'MarkerFaceColor',[1, .5, 0],'Color',[1, .5, 0]);hold on;
    plot([-TT(end:-1:2);TT],[sigred(7,end:-1:2),sigred(3,:)],'k-','linewidth',2)
    plot([0,0],[0,1.3],'--','Color',global_color)

    set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
    set (gca ,'TickLength',[.01,.3],'LineWidth',1);
    title({'mRNA-CTD'},'FontSize',fntsize,'FontWeight','bold')
    xlim([-10,10])
    ylim([-.1,1.3])

    %%

    subplot(4,3,9) %mrna_ser5 signal
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

    plot([-TT(end:-1:2);TT],[sigred(8,end:-1:2),sigred(6,:)],'k-','linewidth',2)
    hold on;
    b1=errorbar(lags2, mrna_ser5.mn_cc, mrna_ser5.sem_cc,'d','MarkerSize',5,'MarkerFaceColor',[.5, 0, 1],'Color',[.5, 0, 1]);hold on;



    par_fixed = parameters;
    par_changed = [1,3:6];
    parameters_samp = parameters;

    minvec = [sigred(8,end:-1:2),sigred(6,:)];
    maxvec = [sigred(8,end:-1:2),sigred(6,:)];
    rng(0)

    for i=1:1000
        n = ceil(rand*size(mh_pars,1));

        parameters_samp(par_changed) = par_samp;
        par_samp = mh_pars(n,:);
         parameters_samp(par_changed) = par_samp;
        [sigred_samp,TT_samp,means_samp,sig_samp,TT2_samp,mins_samp] = get_ac_and_cc_mod_simplified(parameters_samp,[0:.1:30]);   
        %plot([-TT_samp(end:-1:2);TT_samp],[sigred_samp(8,end:-1:2),sigred_samp(6,:)],'-','Color',[.8,.8,.8],'linewidth',2)
        svec = [sigred_samp(8,end:-1:2),sigred_samp(6,:)];
        for j = 1:length(svec)
            if svec(1,j) > maxvec(1,j)
                maxvec(1,j) = svec(1,j);

            end    
             if svec(1,j) < minvec(1,j)
                minvec(1,j) = svec(1,j);

            end           
        end
    end

    fill([[-TT(end:-1:2);TT]',fliplr([-TT(end:-1:2);TT]')], [minvec, fliplr(maxvec)],[.8,.8,.8], 'EdgeColor', [.8,.8,.8]);
    b1=errorbar(lags2, mrna_ser5.mn_cc, mrna_ser5.sem_cc,'d','MarkerSize',5,'MarkerFaceColor',[.5, 0, 1],'Color',[.5, 0, 1]);hold on;
    plot([-TT(end:-1:2);TT],[sigred(8,end:-1:2),sigred(6,:)],'k-','linewidth',2)
    plot([0,0],[0,1.3],'--','Color',global_color)

    set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
    set (gca ,'TickLength',[.01,.3],'LineWidth',1);
    title({'mRNA-Ser5ph'},'FontSize',fntsize,'FontWeight','bold')

    xlim([-10,10])
    ylim([-.1,1.3])
    
    if save_eps
        saveas(gca,'./Figures/ccs_nodist.epsc')
    end




    %% Acov taus, load in the bootstrapped dwell time calculations for speed 
    % and  reproducability 
    T_array = [0:1:1000];

    mrna_bootstrapped_tau = [];
    ser5_bootstrapped_tau = [];
    rnap_bootstrapped_tau = [];

    mrna_sim_acov = zeros(400,32);
    ser5_sim_acov = zeros(400,32);
    rnap_sim_acov = zeros(400,32);

    k = 0;
    
    %% Code to boot strap taus if needed
    % while k < 400
    %     t = [0:1:199];
    % 
    % 
    %     simulated_pol2 = zeros(20,200);
    %     simulated_ser5 = zeros(20,200);
    %     simulated_mrna = zeros(20,200);
    % 
    %     for i = 1:20
    %         sol = run_single_SSA(x0,S,W,T_array,time_var,signal_update_rate);  
    %         [pol2_ssa,ser5_ssa,ts_ssa,~] = get_model_intesities(sol,0,0,0); %convert molecules to signal
    %         [pol2norm,ser5norm,tsnorm] = Normalize_simulated_intensities(.95,pol2_ssa,ser5_ssa,ts_ssa);
    %         simulated_pol2(i,:) = pol2norm(end-199:end);
    %         simulated_ser5(i,:) = ser5norm(end-199:end);
    %         simulated_mrna(i,:) = tsnorm(end-199:end);
    %     end
    % 
    %     [mrna_sim,ser5_sim,rnap_sim,~,~,~,G1_rnap_sim,G1_ser5_sim,G1_mrna_sim] = get_simulated_cc([0:1:199],simulated_pol2', simulated_ser5', simulated_mrna',0);
    % 
    % 
    %     decorr_mrna = t(mrna_sim.mn_ac < .01);
    %     if length(decorr_mrna) == 0
    %         continue
    %     end
    %     
    %     
    %     decorr_ser5 = t(ser5_sim.mn_ac < .01);
    % 
    %     if length(decorr_ser5) == 0
    %         continue
    %     end
    %     decorr_rnap = t(rnap_sim.mn_ac < .01);
    % 
    %     if length(decorr_rnap) == 0
    %         continue
    %     end
    %     
    %     k = k + 1
    %     mrna_bootstrapped_tau = [mrna_bootstrapped_tau, decorr_mrna(1)];
    %     ser5_bootstrapped_tau = [ser5_bootstrapped_tau, decorr_ser5(1)];
    %     rnap_bootstrapped_tau = [rnap_bootstrapped_tau, decorr_rnap(1)];
    %     
    %     mrna_sim_acov(k,:) = mrna_sim.mn_ac;
    %     ser5_sim_acov(k,:) = ser5_sim.mn_ac;
    %     rnap_sim_acov(k,:) = rnap_sim.mn_ac;
    % 
    % end
    
    %%
    cd taus
    load('rnap_sim_acov_nonoise.mat')
    load('ser5_sim_acov_nonoise.mat')
    load('mrna_sim_acov_nonoise.mat')
    cd ..

    %%

    thresh = .2;
    t = [0:1:199];
    mrna_bootstrapped_tau = [];
    ser5_bootstrapped_tau = [];
    rnap_bootstrapped_tau = [];
    for i = 1:400
       tmpacov = t(mrna_sim_acov(i,:) < thresh);

       mrna_bootstrapped_tau = [mrna_bootstrapped_tau, tmpacov(1)];

       tmpacov = t(ser5_sim_acov(i,:) < thresh);
       ser5_bootstrapped_tau = [ser5_bootstrapped_tau, tmpacov(1)];

       tmpacov = t(rnap_sim_acov(i,:) < thresh);
       rnap_bootstrapped_tau = [rnap_bootstrapped_tau, tmpacov(1)];

    end

    % figure(32)
    % histogram(mrna_bootstrapped_tau,'FaceColor','b')
    % hold on;
    % histogram(ser5_bootstrapped_tau,'FaceColor','g')
    % histogram(rnap_bootstrapped_tau,'FaceColor','r')



    %%
    figure(1)
    % subplot(4,3,2)
    % plot( [mean(rnap_bootstrapped_tau),mean(rnap_bootstrapped_tau)], [-1,6],'r', 'LineWidth',2)
    % x1 = mean(rnap_bootstrapped_tau)-std(rnap_bootstrapped_tau);
    % x2 = mean(rnap_bootstrapped_tau)+std(rnap_bootstrapped_tau);
    % fill( [x1,x2, x2,x1], [-1,-1,6,6],'r', 'FaceAlpha',.2,'LineStyle','none')
    % 
    % subplot(4,3,5)
    % plot( [mean(ser5_bootstrapped_tau),mean(ser5_bootstrapped_tau)], [-1,6],'g', 'LineWidth',2)
    % x1 = mean(ser5_bootstrapped_tau)-std(ser5_bootstrapped_tau);
    % x2 = mean(ser5_bootstrapped_tau)+std(ser5_bootstrapped_tau);
    % fill( [x1,x2, x2,x1], [-1,-1,6,6],'g', 'FaceAlpha',.2,'LineStyle','none')
    % 
    % subplot(4,3,8)
    % plot( [mean(mrna_bootstrapped_tau),mean(mrna_bootstrapped_tau)], [-1,6],'b', 'LineWidth',2)
    % 
    % 
    % x1 = mean(mrna_bootstrapped_tau)-std(mrna_bootstrapped_tau);
    % x2 = mean(mrna_bootstrapped_tau)+std(mrna_bootstrapped_tau);
    % fill( [x1,x2, x2,x1], [-1,-1,6,6],'b', 'FaceAlpha',.2,'LineStyle','none')
    % 
    % 




    cd taus
    load('rnap_sim_acov.mat')
    load('ser5_sim_acov.mat')
    load('mrna_sim_acov.mat')
    cd ..

    t = [0:1:199];
    mrna_bootstrapped_tau = [];
    ser5_bootstrapped_tau = [];
    rnap_bootstrapped_tau = [];
    for i = 1:400
       tmpacov = t(mrna_sim_acov(i,:) < thresh);

       mrna_bootstrapped_tau = [mrna_bootstrapped_tau, tmpacov(1)];

       tmpacov = t(ser5_sim_acov(i,:) < thresh);
       ser5_bootstrapped_tau = [ser5_bootstrapped_tau, tmpacov(1)];

       tmpacov = t(rnap_sim_acov(i,:) < thresh);
       rnap_bootstrapped_tau = [rnap_bootstrapped_tau, tmpacov(1)];

    end


    figure(1)
    subplot(4,3,2)
    plot( [mean(rnap_bootstrapped_tau),mean(rnap_bootstrapped_tau)], [-1,6],'r', 'LineWidth',2); hold on;
    x1 = mean(rnap_bootstrapped_tau)-std(rnap_bootstrapped_tau);
    x2 = mean(rnap_bootstrapped_tau)+std(rnap_bootstrapped_tau);
    fill( [x1,x2, x2,x1], [-1,-1,6,6],'r', 'FaceAlpha',.2,'LineStyle','--','edgecolor','r')
    xlim([-.5,25])
    ylim([-.1,1.3])

    subplot(4,3,5)
    plot( [mean(ser5_bootstrapped_tau),mean(ser5_bootstrapped_tau)], [-1,6],'g', 'LineWidth',2); hold on;
    x1 = mean(ser5_bootstrapped_tau)-std(ser5_bootstrapped_tau);
    x2 = mean(ser5_bootstrapped_tau)+std(ser5_bootstrapped_tau);
    fill( [x1,x2, x2,x1], [-1,-1,6,6],'g', 'FaceAlpha',.2,'LineStyle','--','edgecolor','g')
    xlim([-.5,25])
    ylim([-.1,1.3])

    subplot(4,3,8)
    plot( [mean(mrna_bootstrapped_tau),mean(mrna_bootstrapped_tau)], [-1,6],'b', 'LineWidth',2); hold on;


    x1 = mean(mrna_bootstrapped_tau)-std(mrna_bootstrapped_tau);
    x2 = mean(mrna_bootstrapped_tau)+std(mrna_bootstrapped_tau);
    fill( [x1,x2, x2,x1], [-1,-1,6,6],'b', 'FaceAlpha',.2,'LineStyle','--','edgecolor','b')
    xlim([-.5,25])
    ylim([-.1,1.3])
    
    if save_eps
        saveas(gca, './Figures/corrs_with_dwell.epsc')  %save the whole figure
    end
    return
end





if strcmp(paper_fig_num,'3f')
    
    %% SSA trajectory for bottom of figure 3
    figure(1)
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
    sol = run_single_SSA(x0,S,W,T_array,time_var,signal_update_rate);  

    pol2_ssa = sol(2,:)';
    ser5_ssa = sol(2,:)';
    ts_ssa = sol(3,:)';

    [pol2_ssa,ser5_ssa,ts_ssa,~] = get_model_intesities(sol,eta_rnap,eta_ser5,eta_ts); %convert molecules to signal
    [pol2norm,ser5norm,tsnorm] = Normalize_simulated_intensities(.95,pol2_ssa,ser5_ssa,ts_ssa);

    X_SIZE = 13; Y_SIZE = 15;
    figure(1);clf;
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, X_SIZE, Y_SIZE]; % x,y, width, height
    global_color = 'k';

    cla;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh, 1.1*yh]; % x,y, width, height


    avpol2 =movmean(pol2norm(end-200:end),3);
    avser5 = movmean(ser5norm(end-200:end),3);
    avts = movmean(tsnorm(end-200:end),3);
    plot(avpol2(1:201),'r','linewidth',2); hold on; plot(avser5(1:201),'g','linewidth',2); plot(avts(1:201),'b','linewidth',2 )
    legend('POL2','SER5ph','mRNA')
    title('Representative trace','FontSize',fntsize,'FontWeight','bold')
    ylabel('Normalized intensity','FontSize',fntsize,'FontWeight','bold')
    xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
    xlim([0,200])
    ylim([-.5,1.5])




    if save_eps
        saveas(gca, './Figures/trace.epsc')  %save the whole figure
    end
    
    
    return
end

%% Supplemental Figure
% Plot the molecule signals with bursting highlighted and filled in 
if strcmp(paper_fig_num,'sup8a')
    
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
    figure(1)
    rng(1095)
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh, yh]; % x,y, width, height
    T_array = [0:.1:200];
    sol = run_single_SSA(x0,S,W,T_array,time_var,signal_update_rate);

    pol2_ssa = (sol(2,:) + sol(3,:))';
    ser5_ssa = (sol(2,:)+ sol(3,:))';
    ts_ssa = sol(3,:)';


    bubble_top = pol2_ssa(end-500:end);
    bubble_bottom = ts_ssa(end-500:end);

    p2 = pol2_ssa(end-500:end)-ts_ssa(end-500:end);
    on = (p2 > 5);
    fill([T_array(1:501),fliplr(T_array(1:501))], [(on*500)'-500, fliplr((on*500)')],[220,220,220]./256,'linestyle','None','HandleVisibility','off' )
    hold on;
    fill([T_array(1:501),fliplr(T_array(1:501))], [bubble_top', fliplr(bubble_bottom')], [255, 212, 212]./256, 'linestyle','None')


    plot(T_array(1:501),pol2_ssa(end-500:end),'r','linewidth',2); hold on;  plot(T_array(1:501),ts_ssa(end-500:end),'b','linewidth',2 )

    ylim([0,150])
    xlabel({'Time (min)'},'FontSize',fntsize,'FontWeight','bold')
    ylabel({'Simulated molecule counts'},'FontSize',fntsize,'FontWeight','bold')
    a = legend({'Burst','CTD','mRNA'},'FontSize',fntsize,'FontWeight','bold');
    %set (gca ,'TickLength',[.01,.3],'LineWidth',2);
    set(a, 'Box', 'off');
    if save_eps
        saveas(gca,'./Figures/Burst_traj.epsc')
    end
    return
end




%% histograms

if strcmp(paper_fig_num, '3d')
    rng(0)
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

    time_var = 0;
    signal_update_rate = 0;

    W = @(x) W1*x + W0;
    x0 = [0,0,0]';
    data = xlsread('Nasc_TranscriptsCount.xlsx');
    mrna_hist = data(:,1);
    %pol2_hist = data(:,1);

    mrna_hist = mrna_hist(1:end-5);
    %pol2_hist = pol2_hist(1:end-5);



    pol2_ssa_traj = [];
    ser5_ssa_traj = [];
    ts_ssa_traj = [];
    mrnai_traj = [];
    rnap_c = [];
    ser5_c = [];
    mrna_c = [];
    for i = 1:1000
        sol = run_single_SSA(x0,S,W,[0:1:2000],time_var,signal_update_rate);
        pol2_ssa = sol(2,:)';
        ser5_ssa = sol(2,:)';
        ts_ssa = sol(3,:)';
        mrnai = ts_ssa;

     [pol2_ssa,ser5_ssa,ts_ssa,~] = get_model_intesities(sol,eta_rnap,eta_ser5,eta_ts);   

        %mrnai = ts_ssa;



        [pol2norm,ser5norm,tsnorm] = Normalize_simulated_intensities(.95,pol2_ssa,ser5_ssa,ts_ssa);


        pol2_ssa_traj = [pol2_ssa_traj, pol2norm(1:200:end)];
        ser5_ssa_traj = [ser5_ssa_traj, ser5norm(1:200:end)];
        ts_ssa_traj = [ts_ssa_traj, ts_ssa(1:200:end)];
        mrnai_traj = [mrnai_traj, mrnai(1:200:end)];



    end


    pol2_traj = pol2_ssa_traj';
    ser5_traj = ser5_ssa_traj';
    ts_traj = ts_ssa_traj';

    [rnap_I, ser5_I, mrna_I] = Normalize_raw_intensities(.95);   %get the normalized intensities from the data 95 percentile
    rnap_data = reshape(rnap_I(1:20:end,:),1,200);  %rnap_data sorted by 20 time incrementes (decorrelated) this was for histograms
    ser5_data = reshape(ser5_I(1:20:end,:),1,200);
    mrna_data = reshape(mrna_I(1:20:end,:),1,200);

    %xh = 6.3;
    %yh = 6.4;

    fntsize = 18;


    figure(11)
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0,xh/3, 3*yh/3]; % x,y, width, height



    [x,n] = hist(reshape(pol2_traj.',1,[]),30);
    [11, 252, 3]./256
    [218, 51, 255]./256
    x1 = histogram(reshape(pol2_traj.',1,[]),n,'Normalization','probability','FaceColor',[11, 252, 3]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[11, 252, 3]./256)%,'DisplayStyle','stairs')
    hold on;
    x2 = histogram(reshape(pol2_traj.',1,[]),n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256,'DisplayStyle','stairs')

    x3 = histogram(rnap_data ,n,'Normalization','probability','FaceColor',[218, 51, 255]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[218, 51, 255]./256)%,'DisplayStyle','stairs')

    x4 = histogram(rnap_data ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')

    a = legend(gca, [x1,x3],{'Model','Data'},'Location','Best')
    set(a, 'Box', 'off');

    set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
    set(a,'TextColor',global_color);

    ylim([0,.15])
    xlim([-1.3,1.5])
    set(gca,'linewidth',2)
    if save_eps
        saveas(gca,'./Figures/pol2_hist.epsc')
    end
    set(gca,'YTickLabel',[],'XtickLabel',[])
    if save_eps
        saveas(gca,'./Figures/pol2_hist_nolabels.epsc')    
    end

%%
% 
    % title({'CTD Normalized Intensity'},'FontSize',fntsize,'FontWeight','bold','Color',global_color)
    % xlabel({'Intensity (Norm)'},'FontSize',fntsize,'FontWeight','bold')
    % ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
    % 
    % ylim([minYlim maxYlim])
    % set (gca ,'TickLength',[.01,.3],'LineWidth',1);
    % set (gca ,'FontSize',fntsize,'FontName', 'Arial');
    % 
    % %set(gca,'Color','k')
    % %set(gcf,'Color','k')
    % set(gca,'YTickLabel',[],'XtickLabel',[])

    % saveas(gca, 'pol2hist.epsc')
    rng(0)
    figure(10)
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh/3, 3*yh/3]; % x,y, width, height
    [x,n] = hist(reshape(ser5_traj.',1,[]),30);
    x1 = histogram(reshape(ser5_traj.',1,[]),n,'Normalization','probability','FaceColor',[11, 252, 3]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[11, 252, 3]./256)%,'DisplayStyle','stairs')
    hold on;
    x2 = histogram(reshape(ser5_traj.',1,[]),n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256,'DisplayStyle','stairs')

    x3 = histogram(ser5_data ,n,'Normalization','probability','FaceColor',[218, 51, 255]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[218, 51, 255]./256)%,'DisplayStyle','stairs')

    x4 = histogram(ser5_data ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')
    %a = legend(gca,[x1,x3], {'Model','Data'},'Location','Best')
    %set(a, 'Box', 'off');
    ylim([0,.15])
    xlim([-1.3,1.5])
    set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
    set(a,'TextColor',global_color);
    set(gca,'linewidth',2)
    if save_eps
        saveas(gca, './Figures/ser5hist.epsc')
    end
    set(gca,'YTickLabel',[],'XtickLabel',[])
    if save_eps
        saveas(gca,'./Figures/SER5_dist_nolabels.epsc') 
    end

    %%

    % %legend('Model','Data')
    % title({'SER5 Normalized Intensity'},'FontSize',fntsize,'FontWeight','bold','Color',global_color)
    % xlabel({'Intensity (Norm)'},'FontSize',fntsize,'FontWeight','bold')
    % ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
    % %ylim([minYlim maxYlim])
    % %ylim([0,.2])
    % set (gca ,'TickLength',[.01,.3],'LineWidth',1);
    % set (gca ,'FontSize',fntsize,'FontName', 'Arial');


    %set(gca,'Color','k')
    %set(gcf,'Color','k')
    rng(0)
    figure(12)
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh/3, 3*yh/3]; % x,y, width, height

    [x,n] = hist(mrnai_traj,30)

    x1 = histogram(mrnai_traj,n,'Normalization','probability','FaceColor',[11, 252, 3]./256,'FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256)%'DisplayStyle','stairs')
    hold on;
    %histogram(mrna_hist ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')
    x2 = histogram(mrnai_traj,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[11, 252, 3]./256,'DisplayStyle','stairs')

    x3 = histogram(mrna_hist ,n,'Normalization','probability','FaceColor',[218, 51, 255]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[218, 51, 255]./256)%,'DisplayStyle','stairs')

    x4 = histogram(mrna_hist ,n,'Normalization','probability','FaceColor','none','FaceAlpha',.2,'linewidth',3,'edgecolor',[218, 51, 255]./256,'DisplayStyle','stairs')

    %a = legend(gca,[x1,x3], {'Model','Data'},'Location','Best')
    %set(a, 'Box', 'off');

    set (gca ,'FontSize',fntsize,'FontName', 'Arial','Xcolor',global_color,'Ycolor',global_color);
    set(a,'TextColor',global_color);
    set(gca,'linewidth',2)
    if save_eps
        saveas(gca, './Figures/tshist_labels.epsc') 
    end
    set(gca,'YTickLabel',[],'XtickLabel',[])
    if save_eps
        saveas(gca, './Figures/tshist_nolabels.epsc') 
    end

    % title({'mRNA count'},'FontSize',fntsize,'FontWeight','bold','Color',global_color)
    % xlabel({'N molecules'},'FontSize',fntsize,'FontWeight','bold')
    % ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
    % %ylim([minYlim maxYlim])
    % set (gca ,'TickLength',[.01,.3],'LineWidth',1);
    % set (gca ,'FontSize',fntsize,'FontName', 'Arial');
    %set(gca,'Color','k')
    return
end

%% MH plots
if strcmp(paper_fig_num,'sup7')
    % addpath ../Data_files/
    % addpath ../

    %cd ./Model_fits/simple_MH_results/3.24.2020/
    parnames = {'kon','kesc','kproc','beta','kout'};
    parnames = {'beta','omega','k out','k esc','k complete'};

    par_changed = [1:5];
    % 
    % mh_pars = [];
    % mh_vals = [];
    % ikeep = 1;
    % i=1;
    % while ikeep==1
    %     try
    %         fn = ['met_hast_pars_2x_',num2str(i),'.mat'];
    %         load(fn)
    %         mh_pars = [mh_pars;mh_smpl];
    %         mh_vals = [mh_vals;mh_value];
    %         i=i+1;
    %         if i == 3
    %             i = i +1;
    %         end
    %     catch
    %         ikeep=0;
    %     end    
    % end
    % cd ../../..
    % Np = 5;
    % 
    % mh_pars = mh_pars(:,[4,1,5,2,3]);



    par_changed = [1,3:6];
    parnames = {'ω','na','k_{esc}','k_c','β','k_{ab}','na','eta_ctd','eta_ser5','eta_mrna',...
        'sc_ctd','sc_ser5','sc_mrna'};
    nchains = 20;
    nsegments = 50;

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
            tmp_dat = load(sv_file);



            for j=1:nsegments
                try
                    eval(['mh_smpl{i} = [mh_smpl{i};tmp_dat.mh_smpl_',num2str(j),'];']);
                    eval(['mh_value{i} = [mh_value{i};tmp_dat.mh_value_',num2str(j),'];']);
                catch
                end
            end
            figure(10)
            subplot(4,5,i)
            plot(mh_value{i}); hold on
            set(gca,'ylim',[-30,-14])
            figure(3)
            subplot(4,5,i)
            histogram(mh_value{i}); hold on        
            set(gca,'xlim',[-30,-14])
            
            figure(11)
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


    mh_vals = mh_val_all;
    sz = 2*(1+max(mh_vals)-mh_vals);
    figure(1)
    xh = 2*3*6;
    yh = 2*14.3/3;
    fntsize = 18;
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height


    mh_pars = SPS(:,[4,1,5,2,3]);
    par_changed = [5,1,6,3,4];

    [sorted_mh_vals, mh_val_inds] = sort(mh_vals); 
    mh_par_sort = mh_pars(mh_val_inds,:);
    colormap(flipud(parula))

    mh_val_new = -sort(mh_vals);
    mh_val_new(mh_val_new > mh_val_new(.005*(length(mh_val_new)))) = mh_val_new(.005*(length(mh_val_new)));

    Np = 5;
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
       subplot(5,5,2)
       scatter(mh_par_sort(1:100:end,1),mh_par_sort(1:100:end,2),sz(1:100:end),mh_val_new(1:100:end),'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on



       colorbar
        %delete(subplot(5,5,2))
    saveas(gca,'./Figures/sensitivity.epsc')

    %% table of derived parameters with 95% confidence interval
 
    highv = 97.5/100
    lowv = 2.5/100
    beta = mh_pars(:,1);
    omega = mh_pars(:,2);
    kout = mh_pars(:,3);
    kesc = mh_pars(:,4);
    kcomp = mh_pars(:,5);

    table_names = {'beta','omega','kab','kproc','kesc','r','time_cluster','mean_cluster','frac','burst_mrna','rate_mrna','total_rnap','mrna_dwell'};
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

    
    

%     beta = mh_pars(:,1);
%     omega = mh_pars(:,2);
%     kout = mh_pars(:,3);
%     kesc = mh_pars(:,4);
%     kcomp = mh_pars(:,5);
% 
% 
%     mu_ctd = beta.*omega./(kout+kesc);
%     mu_ctd(1)
%     % mu_ctd_on = kin./(kout+kesc);
% 
%     mu_rna = beta.*omega./(kout+kesc).*kesc./kcomp;
%     mu_rna(1)
%     % mu_rna_on = kin./(kout+kesc).*kesc./kproc;
% 
%     rate_prod = mu_ctd.*kesc;
%     % rate_prod_on = mu_ctd_on.*kesc;
% 
%     burst_refractory_period = 1./omega;
%     burst_duration = 1./koff;
%     pol2_burst_size = kin./koff;
% 
%     burst_efficiency = 100*kesc./(kesc+kout);
%     mrna_burst_size = pol2_burst_size.*burst_efficiency/100;
% 
%     figure(7);
% 
% 
%     fntsize = 18;
%     fig1= gcf;
%     fig1.PaperUnits = 'centimeters';
%     fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height
%     sbplots = {'mu_ctd','mu_rna','rate_prod',...
%         'burst_refractory_period',...
%         'pol2_burst_size','burst_efficiency','mrna_burst_size'};
%     sbplots_leg = {'\mu ctd (mol)','\mu rna (mol)',...
%         'rate prod (mol/min)',...
%         'burst refractory period (min)',...
%         'pol2 burst size (mol)','burst efficiency (%)','mrna burst size (mol)'}
%     for i=1:7
%         subplot(3,3,i); hold off
%         histogram(eval(sbplots{i}),'Normalization','probability','FaceColor',[11, 252, 3]./256,'FaceAlpha',.5,'linewidth',.5,'edgecolor',[11, 252, 3]./256 ); hold on
%         plot(eval([sbplots{i},'(1)'])*[1,1],get(gca,'ylim'),'k--','linewidth',3)
% 
%         X =  sort(eval(sbplots{i}));
%         low = X(floor(length(X)/10));
%         high = X(ceil(9*length(X)/10));
% 
%         title([sbplots_leg{i},'; ',num2str(eval([sbplots{i},'(1)']),3),...
%             ' (',num2str(low,3),', ',num2str(high,3),')'],'FontSize',fntsize-5,'FontWeight','bold')
%         if mod(i,3) == 1 
%             ylabel({'Probability'},'FontSize',fntsize,'FontWeight','bold')
%         end
%         if ismember(i,[7,6,5])
%             xlabel({'Value'},'FontSize',fntsize,'FontWeight','bold')
%         end
% 
%     end
%     if save_eps
%         saveas(gca,'./Figures/derived_par_dists.epsc')
%     end
    return
end


%%
if strcmp(paper_fig_num, 'sup10')
    
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
    n = 0;
    rng(20)
    figure(1)
    titles = {'Blocking k_{on}','Blocking k_{escape}','Blocking k_{esc} and 30% k_{proc}','Blocking k_{on} and 30% k_{proc} '};
        fig1= gcf;
        fig1.PaperUnits = 'centimeters';
        fig1.PaperPosition = [0, 0, xh, 3*yh]; % x,y, width, height

    k = 0;
    per_len = 40;
    for i = 1:4


    inhibs = ones(size(parameters));
    k = k+1;
    switch k
        case 1  %inhibit kon
            inhibs(1) = .000001;

        case 2 %inhibit kescape
            inhibs(3) = .01;
        case 3 %reduce esc and proc
            inhibs(4) = .3;
            inhibs(3) = .01;

        case 4%reduce on and kproc
            inhibs(1) = .01;
            inhibs(4) = .3;


    end
    inhibs;
    pol2_traj = [];
    ser5_traj = [];
    ts_traj = [];

    pol2_traji = [];
    ser5_traji = [];
    ts_traji = [];
    n = 0;
    while n < 7

        sol = run_single_SSA_inhib(x0,S,W,[0:1:1200],time_var,signal_update_rate,parameters,inhibs,1110);
        sol(2,1110);
        if sol(3,1110) > means(3)
            n = n+1;
            pol2_ssa = sol(2,:)' + sol(3,:)';
            ser5_ssa = pol2_ssa;
            ts_ssa = sol(3,:)';        

            pol2_ssa = sol(2,:)' + sol(3,:)';
            ser5_ssa = pol2_ssa;
            ts_ssa = sol(3,:)';       


            pol2_traj = [pol2_traj; pol2_ssa'];
            ser5_traj = [ser5_traj; ser5_ssa'];
            ts_traj = [ts_traj; ts_ssa'];

        end




    end


    newpars = parameters;
    newpars = newpars.*inhibs;

    inhibs 
    kon = newpars(1);
    koff = 1000;%parameters(2);
    kesc = newpars(3);
    kproc = newpars(4);
    kin =  1000*newpars(5);
    kout = newpars(6);
    frac = newpars(7);
    eta_rnap = newpars(8);
    eta_ser5 = newpars(9);
    eta_ts = newpars(10);

    Nstates = 3;
    b = zeros(Nstates,1);
    b(1) = kon;
    b(2) = 0;
    b(3) = 0;

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

    A=S*W1;
    tode = 0:.1:per_len-10;
    on = [];
    ctd_ode = [];
    ts_ode = [];
    for j = 1:length(tode)

        %P = expm(A*tode(j)+ W0)*b;
    %      if j == 1
    %          b(2) = means(2);
    %          b(3) = means(3);
    %      else
    %          b(2) = 0;
    %          b(3) = 0;         
    %      end

         EX = -A\b;
         means_new = c*EX;
         %P =  A^-1*(-eye(3) + expm(A*tode(j)))*b - x0;
         P= A^-1*(-eye(3) + expm(A*tode(j)))*b  + expm(A*tode(j))*[0;4.5649;15.2925];
        on = [on, P(1)];
        ctd_ode = [ctd_ode,P(2)];
        ts_ode = [ts_ode,P(3)];
    end

    t = [-10:1:per_len-10];

    subplot(3,4,i)


    plot(t, pol2_traj(1:n,1*1100:1*1100+1*per_len),'r','LineWidth',1); 
    hold on;
    plot(tode,ctd_ode+ts_ode,'k-','LineWidth',2)
    % 
    % plot(t, avpol2(1*1100:1*1100+1*per_len),'r','LineWidth',3); 
    % 
    % fill([t, fliplr(t)], [avpol2(1*1100:1*1100+1*per_len)+stdpol2(1*1100:1*1100+1*per_len),fliplr(avpol2(1*1100:1*1100+1*per_len)-stdpol2(1*1100:1*1100+1*per_len))],'r', 'FaceAlpha',.5)
    % 
    % plot(t, avpol2(1*1100:1*1100+1*per_len)-stdpol2(1*1100:1*1100+1*per_len),'r','LineWidth',1); 
    % plot(t, avpol2(1*1100:1*1100+1*per_len)+stdpol2(1*1100:1*1100+1*per_len),'r','LineWidth',1); 
    % 
    % plot(t, movmean(pol2_traj(1,1*1100:1*1100+1*per_len)./bef_pol2,3),'k','LineWidth',2); 

    plot([0,0],[-10,150],'k--','LineWidth',1);
    plot([-10,per_len-10],[0,0],'k--','LineWidth',1);
    plot([-5,0],[19.8575,19.8575],'k-','LineWidth',2)
    ylim([-10,110])
    xlim([-5,per_len-10])


    %ylim([-1,2])
    title(titles{i},'FontSize',fntsize,'FontWeight','bold')


    if i == 1

       ylabel('CTD Count','FontSize',fntsize,'FontWeight','bold') 
    end

    subplot(3,4,i+4)


     plot(t,ser5_traj(1:n,1*1100:1*1100+1*per_len),'g','LineWidth',1); 
     hold on;
     plot(tode,ctd_ode+ts_ode,'k-','LineWidth',2)

    %  plot(t, avser5(1*1100:1*1100+1*per_len),'g','LineWidth',3);
    %   plot(t, movmean(ser5_traj(1,1*1100:1*1100+1*per_len)./bef_ser5,3),'k','LineWidth',2);
    plot([0,0],[-10,150],'k--','LineWidth',1);

    ylim([-10,110])
    xlim([-5,per_len-10])
    plot([-10,per_len-10],[0,0],'k--','LineWidth',1);
    %ylim([-1,2])
    xlim([-5,per_len-10])
    plot([-5,0],[19.8575,19.8575],'k-','LineWidth',2)

    if i == 1

       ylabel('Ser5ph Count','FontSize',fntsize,'FontWeight','bold') 
    end

    subplot(3,4,i+8)

     plot(t, ts_traj(1:n,1*1100:1*1100+1*per_len),'b','LineWidth',1)
     hold on;
     plot(tode,ts_ode,'k-','LineWidth',2)
     plot([-5,0],[15.2925,15.2925],'k-','LineWidth',2)
    %plot(t, avts(1*1100:1*1100+1*per_len),'b','LineWidth',3)
    % 
    % plot(t, movmean(ts_traj(1,1*1100:1*1100+1*per_len)./bef_ts,3),'k','LineWidth',2)
    plot([0,0],[-10,110],'k--','LineWidth',1);
    plot([-10,per_len-10],[0,0],'k--','LineWidth',1);
    xlim([-5,per_len-10])
    ylim([-10,75])
    if i == 1

       ylabel('mRNA Count','FontSize',fntsize,'FontWeight','bold') 
    end

    % ylim([-1,2])
    % xlim([-10,per_len-10])

     xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')



    if i == 5
        xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
    end
    if i == 6
        xlabel('Time (min)','FontSize',fntsize,'FontWeight','bold')
    end







    end
    if save_eps
        saveas(gca,'./Figures/perturbs_sep1.epsc')
    end
    return
end


%% Simulated chip with SSA

if strcmp(paper_fig_num, 'sup8b')
    rng(0)
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
    npts_pol2 = [];
    npts_ts = [];
    figure(20)
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    xh = 2*3*6;
    yh = 2*14.3/3;
    fig1.PaperPosition = [0, 0, 1*xh, yh]; % x,y, width, height


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
    %while length(npts) < 2000

    pol2_on = [];
    pol2_off = [];
    pol2_transient = [];

    for j = 1:3
        T_array = [0:1:40000];
        sol = run_single_SSA(x0,S,W,T_array,time_var,signal_update_rate);

        pol2_ssa = sol(2,:)';
        ser5_ssa = sol(2,:)+ sol(3,:)';
        ts_ssa = sol(3,:)';



        decorr_ts = ts_ssa(400:40:40000);
        decorr_pol2 = pol2_ssa(400:40:40000);

        npts_pol2 = [npts_pol2, decorr_pol2'];
        npts_ts = [npts_ts, decorr_ts'];

    end





    X = [mean(npts_pol2),mean(npts_ts)];


    time_elongating = 5.2/4.1;     
    residence = 1/kproc;
    time_processing = residence - time_elongating;

    frac_proc = time_processing/residence;
    processing = mean(npts_ts)*frac_proc;
    elongating = mean(npts_ts)*(1-frac_proc)/6;
    mids = ones(1,6)*elongating;
    X = [mean(npts_pol2),mids,processing];


    bar(X)
    %set(gca,'xticklabel',{'1})
    ylabel('Molecule counts')
    hold on;

    text(1-.2,6,string(round(X(1),1)),'FontSize',fntsize-4,'FontWeight','bold') 
    text(8-.2,13,string(round(X(8),1)),'FontSize',fntsize-4,'FontWeight','bold') 
    n = length(npts_pol2);
    title(strcat('Average All nsim = ',string(n)),'FontSize',fntsize,'FontWeight','bold')
    xlabel({'Location' },'FontSize',fntsize,'FontWeight','bold')
    ylabel({'Simulated molecule counts'},'FontSize',fntsize,'FontWeight','bold')
    ylim([0,35])
    
    if save_eps
    saveas(gca,'./Figures/ss.epsc')
    end


    %end

    %%
   
    rng(0)
    npts_pol2 = [];
    npts_ts = [];
    figure(21)
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    xh = 2*3*6;
    yh = 2*14.3/3;
    fig1.PaperPosition = [0, 0, 1*xh, yh]; % x,y, width, height

    %while length(npts) < 2000
    while length(npts_ts) < 500
        length(npts_ts) 
        T_array = [0:1:40000];
        sol = run_single_SSA(x0,S,W,T_array,time_var,signal_update_rate);

        pol2_ssa = sol(2,:)';
        ser5_ssa = sol(2,:)+ sol(3,:)';
        ts_ssa = sol(3,:)';

        decorr_ts = ts_ssa(400:40:40000);
        decorr_pol2 = pol2_ssa(400:40:40000);

        for i = 1:length(decorr_ts)
            if decorr_ts(i)+decorr_pol2(i) > 50
                npts_pol2 = [npts_pol2, decorr_pol2(i)];
                npts_ts = [npts_ts, decorr_ts(i)];
            end

        end

    end


    time_elongating = 5.2/4.1;     
    residence = 1/kproc;
    time_processing = residence - time_elongating;

    frac_proc = time_processing/residence;
    processing = mean(npts_ts)*frac_proc;
    elongating = mean(npts_ts)*(1-frac_proc)/6;
    mids = ones(1,6)*elongating;
    X2 = [mean(npts_pol2),mids,processing];

 
    bar(X2)
    %set(gca,'xticklabel',{'1'})
    ylabel('Molecule counts')
    hold on;

    text(1-.2,30,string(round(X2(1),1)),'FontSize',fntsize-4,'FontWeight','bold') 
    text(8-.2,30,string(round(X2(8),1)),'FontSize',fntsize-4,'FontWeight','bold') 
    n = length(npts_pol2);
    title(strcat('Average On (Total CTD > 50) nsim = ',string(n)),'FontSize',fntsize,'FontWeight','bold')
    xlabel({'Location' },'FontSize',fntsize,'FontWeight','bold')
    ylabel({'Simulated molecule counts'},'FontSize',fntsize,'FontWeight','bold')
    ylim([0,35])
    if save_eps
    saveas(gca,'./Figures/bursting_chip.epsc')
    end

    %%
    
    rng(0)
    npts_pol2 = [];
    npts_ts = [];
    figure(22)
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    xh = 2*3*6;
    yh = 2*14.3/3;
    fig1.PaperPosition = [0, 0, 1*xh, yh]; % x,y, width, height

    %while length(npts) < 2000

    for j = 1:3
        T_array = [0:1:40000];
        sol = run_single_SSA(x0,S,W,T_array,time_var,signal_update_rate);

        pol2_ssa = sol(2,:)';
        ser5_ssa = sol(2,:)+ sol(3,:)';
        ts_ssa = sol(3,:)';

        decorr_ts = ts_ssa(400:40:40000);
        decorr_pol2 = pol2_ssa(400:40:40000);

        for i = 1:length(decorr_ts)
            if decorr_pol2(i) < 5
                npts_pol2 = [npts_pol2, decorr_pol2(i)];
                npts_ts = [npts_ts, decorr_ts(i)];
            end

        end
    end

    time_elongating = 5.2/4.1;     
    residence = 1/kproc;
    time_processing = residence - time_elongating;

    frac_proc = time_processing/residence;
    processing = mean(npts_ts)*frac_proc;
    elongating = mean(npts_ts)*(1-frac_proc)/6;
    mids = ones(1,6)*elongating;
    X2 = [mean(npts_pol2),mids,processing];

 
    bar(X2)
    %set(gca,'xticklabel',{'1'})
    ylabel('Molecule counts')
    hold on;

    text(1-.2,2,string(round(X2(1),1)),'FontSize',fntsize-4,'FontWeight','bold') 
    text(8-.2,13,string(round(X2(8),1)),'FontSize',fntsize-4,'FontWeight','bold') 

    n = length(npts_ts)

    title(strcat('Average Off (Unescaped CTD < 5) nsim = ',string(n)),'FontSize',fntsize,'FontWeight','bold')
    xlabel({'Location' },'FontSize',fntsize,'FontWeight','bold')
    ylabel({'Simulated molecule counts'},'FontSize',fntsize,'FontWeight','bold')
    ylim([0,35])
    if save_eps
    saveas(gca,'./Figures/off_chip.epsc')
    end

    %%
   
    rng(0)
    npts_pol2 = [];
    npts_ts = [];
    figure(23)
    fig1= gcf;
    fig1.PaperUnits = 'centimeters';
    xh = 2*3*6;
    yh = 2*14.3/3;
    fig1.PaperPosition = [0, 0, 1*xh, yh]; % x,y, width, height

    %while length(npts) < 2000

    for j = 1:3
        T_array = [0:1:40000];
        sol = run_single_SSA(x0,S,W,T_array,time_var,signal_update_rate);

        pol2_ssa = sol(2,:)';
        ser5_ssa = sol(2,:)+ sol(3,:)';
        ts_ssa = sol(3,:)';

        decorr_ts = ts_ssa(400:40:40000);
        decorr_pol2 = pol2_ssa(400:40:40000);

        for i = 1:length(decorr_ts)
            if decorr_pol2(i) >= 5 && decorr_pol2(i) < 15
                npts_pol2 = [npts_pol2, decorr_pol2(i)];
                npts_ts = [npts_ts, decorr_ts(i)];
            end

        end
    end

    time_elongating = 5.2/4.1;     
    residence = 1/kproc;
    time_processing = residence - time_elongating;

    frac_proc = time_processing/residence;
    processing = mean(npts_ts)*frac_proc;
    elongating = mean(npts_ts)*(1-frac_proc)/6;
    mids = ones(1,6)*elongating;
    X2 = [mean(npts_pol2),mids,processing];


    bar(X2)
    %set(gca,'xticklabel',{'1'})
    ylabel('Molecule counts')
    hold on;

    text(1-.2,10,string(round(X2(1),1)),'FontSize',fntsize-4,'FontWeight','bold') 
    text(8-.2,17,string(round(X2(8),1)),'FontSize',fntsize-4,'FontWeight','bold') 

    n = length(npts_ts)

    title(strcat('Average Transient (5 < Unescaped CTD < 15) nsim = ',string(n)),'FontSize',fntsize,'FontWeight','bold')
    xlabel({'Location' },'FontSize',fntsize,'FontWeight','bold')
    ylabel({'Simulated molecule counts'},'FontSize',fntsize,'FontWeight','bold')
    ylim([0,35])
    
    if save_eps
        saveas(gca,'./Figures/trans_chip.epsc')
    end
    return

end




function [thresh,std_mod] = compute_thresh(ssa,eta,tr)
std_mod = std(ssa);
shot_mod = std_mod*eta;
ssa_w_shot = ssa + randn(size(ssa))*shot_mod;
% Define fraction tr as zero
tmp = sort(ssa_w_shot);
thresh = tmp(ceil(length(tmp)*tr));

end



function ssa_w_shot = add_shot_std(ssa,eta,stdev,thresh)
std_mod = stdev;
shot_mod = std_mod*eta;
ssa_w_shot = ssa + randn(size(ssa))*shot_mod;

% Define fraction tr as zero

ssa_w_shot = ssa_w_shot-thresh;
end



end