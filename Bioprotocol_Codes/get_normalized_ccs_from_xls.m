function [IntensityData, CorrelationData] = get_normalized_ccs_from_xls(FileName, DataCells, NChannels, ...
                                                        ChannelNames, NCells, LTrajectory, NormalizationQuantile, ...
                                                        MaxValue, RezeroAutoCorrelation, RezeroPoints, G0_type, CC_G0, CC_Norm, CC_delays, ACC_delays)

%Make Intensity data struct

% Intensity.Name.raw_data = raw data loaded in
% Intensity.Name.normalized_data = normalized data by quantile and max
% value
IntensityData = struct();
for i = 1:NChannels
    IntensityData.(ChannelNames{i}) = struct();
    tmp_data = xlsread(FileName,'Sheet1',DataCells{i}); % load data in 
    tmp_data = reshape(tmp_data, LTrajectory, NCells);  % reshape to time x cells
    IntensityData.(ChannelNames{i}).raw_data = tmp_data; % save raw_data  

    for j = 1:NCells %Normalize by % quantile
        top_percentile_value = quantile(tmp_data(:,j),NormalizationQuantile);
        tmp_data(:,j) = (tmp_data(:,j))./(top_percentile_value);
    end

    tmp_data = min(tmp_data,MaxValue); %set to max value
    IntensityData.(ChannelNames{i}).normalized_data = tmp_data; %store normalized data
end


% Calculate all autocorrelation pairs first so we can use their g0 info for
% ccs

CorrelationData = struct(); % Set up structure

for i = 1:NChannels
    for j = 1:NChannels
        signal1 = IntensityData.(ChannelNames{i}).normalized_data; %get channeli intensities
        signal2 = IntensityData.(ChannelNames{j}).normalized_data; %get channeli intensities
        
        if i == j
            ccs = zeros(NCells, ACC_delays+1);  %preallocate matrices to fill fo
            G0s = zeros(1,NCells);
            for k = 1:NCells
                V = [signal1(:,k), signal2(:,k)];
                Nt = LTrajectory;
                V = (V - repmat(mean(V),LTrajectory,1));
                [C,~] = xcorr(V);
              
                C(round(length(C)/2):end,:) = C(round(length(C)/2):end,:)./(Nt:-1:1)'; %Xcorr scale fix
                C(1:round(length(C)/2)-1,:) = C(1:round(length(C)/2)-1,:)./(1:1:Nt-1)'; 
    
                corr = C(:,1); %get correlation
                corr = corr(((length(corr)+1)/2):((length(corr)+1)/2)+ACC_delays); %adjust and snip N delays
                G0s(k) = get_g0_single(corr,G0_type); %get g0
                ccs(k,:) = corr; %store correlation
            end
            G0s_mean = mean(G0s); %take the mean of all trajectory g0s
            ccs = ccs./G0s_mean;

            CorrelationData.(strcat(ChannelNames{i}, '_', ChannelNames{j})).mean_corr = mean(ccs,1);
            CorrelationData.(strcat(ChannelNames{i}, '_', ChannelNames{j})).sem_corr = std(ccs,0,1)./sqrt(NCells);            
            CorrelationData.(strcat(ChannelNames{i}, '_', ChannelNames{j})).normalized_corrs = ccs; %store normalized corr
            CorrelationData.(strcat(ChannelNames{i}, '_', ChannelNames{j})).G0_mean = G0s_mean; %store g0 mean

        end
    end
end

% calculate all correlation pairs
for i = 1:NChannels
    for j = 1:NChannels
        signal1 = IntensityData.(ChannelNames{i}).normalized_data;
        signal2 = IntensityData.(ChannelNames{j}).normalized_data;
        
        if i ~= j
            ccs = zeros(NCells, CC_delays*2 + 1);

            for k = 1:NCells
                V = [signal1(:,k), signal2(:,k)];
                Nt = LTrajectory;
                V = (V - repmat(mean(V),LTrajectory,1));
                [C,~] = xcorr(V);
                
                C(round(length(C)/2):end,:) = C(round(length(C)/2):end,:)./(Nt:-1:1)'; %Xcorr scale fix
                C(1:round(length(C)/2)-1,:) = C(1:round(length(C)/2)-1,:)./(1:1:Nt-1)'; 
                corr = C(:,2);
                corr = corr(((length(corr)+1)/2-CC_delays):((length(corr)+1)/2)+CC_delays);
                ccs(k,:) = corr;

            end
            G0_sig1 = CorrelationData.(strcat(ChannelNames{i},'_',ChannelNames{i})).G0_mean;
            G0_sig2 = CorrelationData.(strcat(ChannelNames{j},'_',ChannelNames{j})).G0_mean;
            
            if strcmp(CC_Norm, 'ACC_SQRT')
                ccs = ccs./(sqrt(G0_sig1)*sqrt(G0_sig2));
            end

            G0_cc  = get_g0_cc_single(mean(ccs,1),CC_G0); %get g0 for cc
            

            CorrelationData.(strcat(ChannelNames{i}, '_', ChannelNames{j})).mean_corr = mean(ccs,1)./G0_cc;
            CorrelationData.(strcat(ChannelNames{i}, '_', ChannelNames{j})).sem_corr = std(ccs,0,1)./sqrt(NCells)./G0_cc;
            CorrelationData.(strcat(ChannelNames{i}, '_', ChannelNames{j})).normalized_corrs = ccs;
            CorrelationData.(strcat(ChannelNames{i}, '_', ChannelNames{j})).G0_mean = G0_cc;

        end
    end
end

if RezeroAutoCorrelation == true
    for i = 1:NChannels
        cc = CorrelationData.(strcat(ChannelNames{i}, '_', ChannelNames{i})).mean_corr;
        rezero_value = mean(cc(end-RezeroPoints:end));
        scale_factor = cc(1)/(cc(1)-rezero_value);
        CorrelationData.(strcat(ChannelNames{i}, '_', ChannelNames{i})).mean_corr = (cc-rezero_value)*scale_factor;
    end
end