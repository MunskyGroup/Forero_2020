function plot_histograms(ModelIntensities, DataIntensities, ChannelNames, nbins)

    % Function that visualizes the model intensity histograms vs the data intensity histograms

    % ModelIntensities - Intensity structure from a model
    % DataIntensities - Intensity structure from the data
    % ChannelNames - char array of channel names
    % nbins - number of bins to bin the histograms into

    figure();
    k = 1;
    
    for i = 1:length(ChannelNames)
        subplot(1,length(ChannelNames),k)
    
        model_data = ModelIntensities.(ChannelNames{i}).normalized_data;
        data_data = DataIntensities.(ChannelNames{i}).normalized_data;
        [~,n] = hist(reshape(model_data.',1,[]),nbins);
    
        x1 = histogram(reshape(model_data.',1,[]), n,'Normalization','probability','FaceColor',[11, 252, 3]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[11, 252, 3]./256);%,'DisplayStyle','stairs')
        hold on;
        x3 = histogram(reshape(data_data.',1,[]) ,n,'Normalization','probability','FaceColor',[218, 51, 255]./256,'FaceAlpha',.2,'linewidth',.01,'edgecolor',[218, 51, 255]./256);%,'DisplayStyle','stairs')
    
        a = legend(gca, [x1,x3],{'Model','Data'},'Location','Best');
        set(a, 'Box', 'off');
        xlabel([ChannelNames{i}, ' ', 'Normalized Intensity'], 'Interpreter','none')
        k = k + 1;
    end


end

