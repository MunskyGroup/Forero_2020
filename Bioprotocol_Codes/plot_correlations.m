function plot_correlations(ModelCorrelations, DataCorrelations, ChannelNames)

    % Function that visualizes the model correlations vs the data correlations

    % ModelCorrelations - Correlation structure from a model
    % DataCorrelations - Correlation structure from the data
    % ChannelNames - char array of channel names
    
    figure(1);
    k = 1;
    
    % get correlation related indices
    Nchannels = length(ChannelNames);
    diagonal_inds = zeros(Nchannels); %indices of autocorrelations
    for i = 1:Nchannels
        diagonal_inds(i) = 1 + (i-1)*Nchannels + (i-1)*1;
    end
    
    all_inds = zeros(Nchannels^2,2);
    m = 1;
    for i = 1:Nchannels
        for j = 1:Nchannels
            all_inds(m,:) = [i,j];
            m = m+1;
        end
    end
  
    % loop to plot the correlations
    for i = 1:length(ChannelNames)
        for j = 1:length(ChannelNames)
            if i >= j
                subplot(Nchannels,Nchannels,k);
                
                tmp_dat = DataCorrelations.(strcat(ChannelNames{i}, '_',ChannelNames{j})).mean_corr; 
                tmp_sem_dat = DataCorrelations.(strcat(ChannelNames{i}, '_',ChannelNames{j})).sem_corr; 
                
                if i == j
                    errorbar([0:length(tmp_dat)-1],tmp_dat',tmp_sem_dat','rs','linewidth',2); hold on
                    plot([0:length(ModelCorrelations(k,:))-1],ModelCorrelations(k,:),'Color','k','linewidth',3); hold on
                    
                else            
                    [~, idx_1] = ismember(all_inds, [i,j], 'rows');
                    [~, idx_2] = ismember(all_inds, [j,i], 'rows');
                    idx_1 = find(idx_1);
                    idx_2 = find(idx_2);
                    errorbar([-ceil(length(tmp_dat)/2)+1:ceil(length(tmp_dat)/2)-1],tmp_dat',tmp_sem_dat','rs','linewidth',2); hold on
                    plot([-ceil(length(tmp_dat)/2)+1:ceil(length(tmp_dat)/2)-1], [ModelCorrelations(idx_1,ceil(length(tmp_dat)/2):-1:2),ModelCorrelations(idx_2,1:ceil(length(tmp_dat)/2))],'Color','k','linewidth',3); hold on;
                end
                set(gca,'fontsize',14);title(strcat(ChannelNames{i}, '_',ChannelNames{j}), 'Interpreter','none');
                grid on
                if i == Nchannels
                    xlabel('Delays')
                end
    
                if j == 1
                    ylabel('Norm Corr')
                end
                if i == j
                   legend('Data','Model')
                end
            end
             k = k + 1;
        end
    end

