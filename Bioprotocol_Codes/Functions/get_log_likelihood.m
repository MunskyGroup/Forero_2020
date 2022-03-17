function err = get_log_likelihood(pars, NormalizedCorrelationData, DataMeans, DataSEMs, S, W1, W0, c, noise_parameters, time_vec, cc_range, acc_range, channels_quantified, Nspots, ChannelNames, par_fixed, par_changed)

    % Function that calculates the log likliehood from data correlation and a set of parameters
    % pars - full list of parameters
    % NormalizedCorrelationData - Correlation data structure [Nchannels^2]
    % DataMeans - data spot means [1xNchannels]
    % DataSEMs - data spot standard error of the means [1xNchannels]
    % S, W1, W0, c, b - model matrices
    % noise parameters (etas, scales)
    % time_vec - time vector to solve the model for
    
    % cc_range - the cc points to consider in the SSE calculation, 8 would consider -+8 from center delay
    % acc_range - the acc points to consider in their SSE calculation, 31 would consider delays 0 - 31 (32 total)

    % channels quantified - which channels true or false have quantifications to calculate [1xNchannels]
    % Nspots - n total spots quantified for each channel [1xNchannels]
    % ChannelNames - char array of names for each channel [1xNchannels]

    % par_fixed - a backup copy of parameters so fixed parameters stays constant
    % par_changed - indices of parameters that are allowed to change for MLE estimation.
    
    parameters = par_fixed;
    parameters(par_changed) = pars;
    
    [ModelCorrelations, ModelMeans, ModelVariances] = solve_model_from_mats(S, W1, W0, c, noise_parameters, time_vec, parameters);
    
    Nchannels = length(ChannelNames);
    acc_points = length(NormalizedCorrelationData.(strcat(ChannelNames{1},'_',ChannelNames{1})).mean_corr); 
    cc_points = length(NormalizedCorrelationData.(strcat(ChannelNames{1},'_',ChannelNames{2})).mean_corr); 
    sigred_data = zeros(length(ChannelNames)^2, max(acc_points,cc_points) );
    sem_sigred_data = zeros(length(ChannelNames)^2, max(acc_points,cc_points));
    
    acc_range = acc_range + 1; %adjust the range + 1
    cc_range = ceil(cc_range/2);
    k = 1;
    
    %Put correlation data in a single matrix
    for i = 1:Nchannels
        for j = 1:Nchannels
            if i ==j
                tmp_dat = NormalizedCorrelationData.(strcat(ChannelNames{i},'_',ChannelNames{j})).mean_corr;
                sigred_data(k,1:acc_points) = tmp_dat;
                sem_sigred_data(k,1:acc_points) =  NormalizedCorrelationData.(strcat(ChannelNames{i},'_',ChannelNames{j})).sem_corr;
            else
                tmp_dat = NormalizedCorrelationData.(strcat(ChannelNames{i},'_',ChannelNames{j})).mean_corr;
                sigred_data(k,1:cc_points) = tmp_dat;
                sem_sigred_data(k,1:cc_points) =  NormalizedCorrelationData.(strcat(ChannelNames{i},'_',ChannelNames{j})).sem_corr;
            end
            k = k + 1;
        end
    end
    
    % Calculate the indices for correlations
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
    sigred = ModelCorrelations;
    
    % Calculate the errors for auto correlations
    k = 1;
    for i = 1:Nchannels
        err(k) = sum((sigred(diagonal_inds(i),1:acc_range)-sigred_data(diagonal_inds(i), 1:acc_range)).^2 ./ sem_sigred_data(diagonal_inds(i) , 1:acc_range).^2)/2;
        k = k+1;
    end
    
    
    % Calculate the errors for cross correlations
    n = 1;
    for i = 1:Nchannels
        for j = 1:Nchannels
            if j < i  %only do cross correlations on one triangle of the Corr mat
    
                [~, idx_1] = ismember(all_inds, [i,j], 'rows');
                [~, idx_2] = ismember(all_inds, [j,i], 'rows');
                idx_1 = find(idx_1);
                idx_2 = find(idx_2);
                tmp_dat = sigred_data(n, 1:cc_points);
                tmp_dat_sem = sem_sigred_data(n, 1:cc_points);
                center = ceil(length(tmp_dat)/2);
                err(k) = sum(([sigred(idx_2,cc_range:-1:2),sigred(idx_1,1:cc_range)]-flip(tmp_dat(center - cc_range+1 : center + cc_range-1))).^2./flip(tmp_dat_sem(center - cc_range+1 : center + cc_range-1)).^2)/2;
                k = k+1;
            end
            n = n + 1;
        end
        
    end
    
    % match average number of mRNA
    mean_dat = DataMeans(channels_quantified);
    var_dat = DataSEMs(channels_quantified).^2.*Nspots(channels_quantified);
    
    %Error for model / data means for channels quantified
    err(k) = sum((ModelMeans(channels_quantified) - mean_dat)^.2 ./ var_dat.^2 ./ 2);
    
    % match variance in number of mRNA
    err(k+1) = sum((diag(ModelVariances(channels_quantified,channels_quantified)) - var_dat).^2 ./ 10^2 ./2);
    
    
    %% Optional priors here if needed
    % Prior - parameter values should be positive.
    %err(k+3) = 1000*sum(max(-parameters+1e-8,0));
    
    % All other remaining priors have been ignored.
    % Uncomment if needed
    
    % frac should be between 0 and 1
    % err(9) = 100*max(0,log(parameters(7)));
    
    % err(9) = sum(max(0,log(parameters)-3));
    % covariances
    % cv12 = 1/parameters(11);
    % cv13 = 1/parameters(12);
    % cv23 = 1/parameters(13);
    % cv12 = 1/parameters(11)*sigred(2,1);
    % cv13 = 1/parameters(12)*sigred(3,1);
    % cv23 = 1/parameters(13)*sigred(6,1);
    % sig_d = [9.7339e-01   3.4716e-01   5.0248e-01];
    % sig_d_sem = [1.4262e-01   8.6887e-02   1.1456e-01];
    % %err(13) = sum(abs(log([cv12,cv13,cv23]./sig_d)));
    % err(13) = 0*sum(([cv12,cv13,cv23]-sig_d).^2./sig_d_sem.^2)/2;

end