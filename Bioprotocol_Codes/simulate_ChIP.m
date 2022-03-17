function [SimulatedChIP, SimulatedChIP_std] = simulate_ChIP(S, W1, W0, c, ...
                      time_vec, parameters,  gene_size, elongation_rate, total_residence_time, ...
                      number_trajectories, number_samples, n_sites)

    % Function that simulates chromatin immunoprecipitation assay of a gene

    % S, W1, W0, c - model matrices
    % time_vec - time vector to run the assay over
    % gene_size - kilobases of the gene 
    % elognation_rate - kb/min rate of RNAP elongation
    % total_residence_time - min, how long does a RNAP reside after initation? calculated from 1/kprocessing
    % number_trajectories - number of trajectories to sample a chip from
    % n_samples - number of samples to take of n_trajectories
    % n_sites - number of internal ChIP sites to split elongation over
    
    pol2_samples = zeros(number_samples, number_trajectories);
    ts_samples = zeros(number_samples, number_trajectories);
    for i = 1:number_samples
        [MoleculeTrajectories] = sample_model_from_mats(S, W1, W0, c, time_vec, parameters, ...
                              number_trajectories);
        pol2_samples(i,:) = squeeze(MoleculeTrajectories(1,end,:));
        ts_samples(i,:) = squeeze(MoleculeTrajectories(3,end,:));
    
    end
    
    ts_samples = reshape(squeeze(ts_samples).',1,[]);
    pol2_samples = reshape(squeeze(pol2_samples).',1,[]);
    
    time_elongating = gene_size /elongation_rate;     
    residence = total_residence_time;
    time_processing = residence - time_elongating;
    
    frac_proc = time_processing/residence;
    processing = mean(ts_samples)*frac_proc;
    elongating = mean(ts_samples)*(1-frac_proc)/n_sites;
    mids = ones(1,n_sites)*elongating;
    SimulatedChIP = [mean(pol2_samples),mids,processing];
    SimulatedChIP_std = [std(pol2_samples), ones(1,n_sites)*(std(ts_samples)*(1-frac_proc)/n_sites),std(ts_samples)*frac_proc];

end