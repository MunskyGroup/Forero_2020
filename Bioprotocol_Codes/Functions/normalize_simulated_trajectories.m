function [IntensityData] = normalize_simulated_trajectories(IntensityTrajectories, NChannels, ...
                                                        ChannelNames, NormalizationQuantile, ...
                                                        MaxValue)

    % Function that normalizes simulated intensity trajectories from SSAs

    % IntensityTrajectories - Intensity data structure
    % NChannels - number of channels
    % ChannelNames - character array of channel names
    % NormalizationQuantile - percentile to set as 1, default is .95
    % MaxValue - maximum value after normalization, default is 1.5

    IntensityData = struct();
    for i = 1:NChannels
        IntensityData.(ChannelNames{i}) = struct();
        tmp_data = squeeze(IntensityTrajectories(i,:,:)); % load data in 
        IntensityData.(ChannelNames{i}).raw_data = tmp_data; % save raw_data 
    
        for j = 1:size(tmp_data,2) %Normalize by % quantile
            top_percentile_value = quantile(tmp_data(:,j),NormalizationQuantile);
            tmp_data(:,j) = (tmp_data(:,j))./(top_percentile_value);
        end
        tmp_data = min(tmp_data,MaxValue); %set to max value
        IntensityData.(ChannelNames{i}).normalized_data = tmp_data; %store normalized data   
    end

