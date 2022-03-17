
function [PerturbedMoleculeTrajectories] = simulated_perturbation(S, W1, W0, c, ...
                      time_vec, parameters, perturbation_vector, ...
                      perturbation_time, number_trajectories);


    % Sample the model for molecular SSA counts with a perturbation

    % x' = S*W1(x) + W0

    % S - stoichiometry matrix 
    % W1 - Linear Propensity matrix(x)
    % W0 - Independent Propensity matrix
    % c - intensity transformation matrix, ie molecule to intensity matrix
    % time_vec - time vector of each sample
    % parameters - parameters to solve for
    % number of trajectories - number of trajectories to return
    % perturbation_time - time to apply the perturbation
    % perturbation_vector - change value of each paramter at perturbation time [1xNparameters]

    x0 = zeros(size(S,1),1);
    
    time_var = 0;
    signal_update_rate = 0;
    tmp_W1 = W1(parameters);
    tmp_W0 = W0(parameters);

    W = @(x) tmp_W1*x + tmp_W0;

    tmp_W0_after = W0(parameters.*perturbation_vector);
    tmp_W1_after = W1(parameters.*perturbation_vector);

    Wafter = @(x) tmp_W1_after*x + tmp_W0_after;

    c_t = c(parameters);
    PerturbedMoleculeTrajectories = zeros( size(S,1), length(time_vec), number_trajectories);
    for i = 1:number_trajectories

        sol = run_single_SSA_generic_inhib(x0,S,W, Wafter,time_vec,time_var,signal_update_rate, perturbation_time);  
        sol = c_t*sol; %convert using fraction matrix 
        PerturbedMoleculeTrajectories(:,:,i) = sol;
    end




