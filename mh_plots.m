parnames = {'kon','kesc','kproc','beta','kout','k_ctd_to_ser5'};
parnames = {'beta','omega','k out','k esc','k_cs','k comp'};

par_changed = [1:6];

mh_pars = [];
mh_vals = [];
ikeep = 1;
i=1;
while ikeep==1
    try
        fn = ['met_hast_pars_2x_cs_',num2str(i),'.mat'];
        load(fn)
        mh_pars = [mh_pars;mh_smpl];
        mh_vals = [mh_vals;mh_value];
        i=i+1;
    catch
        ikeep=0;
    end    
end

Np = 6;

mh_pars = mh_pars(:,[4,1,5,2,3,6]);

sz = 4*real(1+max(mh_vals)-mh_vals);
%sz = ones(size(mh_vals));

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
