function [MoleculeTrajectories] = sample_model_from_mats(S, W1, W0, c, ...
                       time_vec, parameters, ...
                      number_trajectories);

    % Sample the model for molecular SSA counts

    % x' = S*W1(x) + W0

    % S - stoichiometry matrix 
    % W1 - Linear Propensity matrix(x)
    % W0 - Independent Propensity matrix
    % c - intensity transformation matrix, ie molecule to intensity matrix
    % time_vec - time vector of each sample
    % parameters - parameters to solve for
    % number of trajectories - number of trajectories to return

    x0 = zeros(size(S,1),1);
    
    time_var = 0;
    signal_update_rate = 0;
    tmp_W1 = W1(parameters);
    tmp_W0 = W0(parameters);
    W = @(x) tmp_W1*x + tmp_W0;
    c_t = c(parameters);
    MoleculeTrajectories = zeros( size(S,1), length(time_vec), number_trajectories);
    for i = 1:number_trajectories
        sol = run_single_SSA(x0,S,W,time_vec,time_var,signal_update_rate);  
    
        sol = c_t*sol; %convert using fraction matrix 
    
        MoleculeTrajectories(:,:,i) = sol;
    end


     

     