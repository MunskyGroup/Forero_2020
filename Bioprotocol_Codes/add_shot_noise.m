function [IntensityTrajectories] = add_shot_noise(MoleculeTrajectories, parameters, noise_parameters, zero_values)
    
    % Function that takes model molecule counts from SSA trajectories and converts them to 
    % intensities with shot noise based on the noise parameters

    % Molecule Trajectories  - Nchannels x Ntime molecule counts
    % parameters - Parameters used to make the molecule trajectories -  1xNpar
    % noise_parameters - the specific noise parameters calculated from ACC G0 for each channel - 1xNchannel*2
    % zeros_value - approximation of zero bg value (counts of times from real data trajectories that dip below zero)

    np_c = noise_parameters(parameters);
    IntensityTrajectories = zeros(size(MoleculeTrajectories));
    for i = 1:size(MoleculeTrajectories,3)
        for j = 1:size(MoleculeTrajectories,1)
            MoleculeTrajectories(j,:,i);
            IntensityTrajectories(j,:,i) = add_individual_noise(MoleculeTrajectories(j,:,i),np_c(j),zero_values(j));
        end
    end
end


function ssa_w_shot = add_individual_noise(ssa,eta,tr)
    std_mod = std(ssa);
    shot_mod = std_mod*eta;
    ssa_w_shot = ssa + randn(size(ssa))*shot_mod;

    % Define fraction tr as zero
    tmp = sort(ssa_w_shot);
    thresh = tmp(ceil(length(tmp)*tr));
    ssa_w_shot = ssa_w_shot-thresh;
 end