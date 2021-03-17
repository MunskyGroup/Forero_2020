%% MH plots

addpath ../Data_files/
addpath ../

cd ./Model_fits/simple_MH_results/3.24.2020/
parnames = {'kon','kesc','kproc','beta','kout'};
parnames = {'beta','omega','k out','k esc','k complete'};

par_changed = [1:5];

mh_pars = [];
mh_vals = [];
ikeep = 1;
i=1;
while ikeep==1
    try
        fn = ['met_hast_pars_2x_',num2str(i),'.mat'];
        load(fn)
        mh_pars = [mh_pars;mh_smpl];
        mh_vals = [mh_vals;mh_value];
        i=i+1;
        if i == 3
            i = i +1;
        end
    catch
        ikeep=0;
    end    
end
cd ../../..
Np = 5;

mh_pars = mh_pars(:,[4,1,5,2,3]);

sz = 2*(1+max(mh_vals)-mh_vals);
figure(6)
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

        scatter(mh_pars(1:100:end,i),mh_pars(1:100:end,j),sz(1:100:end),'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerFaceColor',[50, 139, 191]./256,'MarkerEdgeColor',[50, 139, 191]./256); hold on
        scatter(mh_pars(1:100:end,i),mh_pars(1:100:end,j),sz(1:100:end),'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'MarkerFaceColor',[1,.5,0]); hold on

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



kon = mh_pars(:,1);
koff = 1000;
kesc = mh_pars(:,2);
kproc = mh_pars(:,3);
kin = mh_pars(:,4)*koff;
kout = mh_pars(:,5);

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
cd ./Model_fits/simple_MH_results/3.24.2020/
figure(8) 
i = 1
ikeep = 1;
while ikeep==1
    try
        fn = ['met_hast_pars_2x_',num2str(i),'.mat'];
        load(fn)
        mh_pars = [mh_pars;mh_smpl];
        mh_vals = [mh_vals;mh_value];
        i=i+1
        if i == 3
            i = i +1;
        end
        
    catch
        ikeep=0;
    end    

    for j = 1:Np
        subplot(Np,Np,(j-1)*Np+j)
        
        plot(mh_smpl(1:100:end,j)); hold on;
        order = [4,1,5,2,3];
        title({parnames{par_changed(order(j))}} ,'FontSize',fntsize,'FontWeight','bold');
    end
end
cd ../../..