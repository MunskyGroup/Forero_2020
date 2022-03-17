function [PerturbationAnalyticalSolution] = solve_perturbed_analytical_solution(S, W1, W0, c, b, parameters, perturbation_vector)
        
        % Function that solves analytical means solution 

        % S, W1, W0, c, b - Model Matrices
        % parameters - parameters to solve 
        % perturbation_vector - perturbuation multiplication for each parameter

        per_len = 40;

        Nchannels = size(c(parameters),1);

        tmp_W1 = W1(parameters);
        tmp_W0 = W0(parameters);

        W = @(x) tmp_W1*x + tmp_W0;

        tmp_W1_after = W1(parameters.*perturbation_vector);


        c_t = c(parameters);

        A=S*W1(parameters);  %get A matrix for dp/dt = Ap
        EX_original =  -A\b(parameters);

        A_new = S*tmp_W1_after;
        tode = 0:.1:per_len-10;
        ode_sol = zeros(Nchannels,length(tode) );
        ode_sol_before = zeros(Nchannels,30 );
        for j = 1:30
            ode_sol_before(:,j) = EX_original;
        end

        for j = 1:length(tode)   
             EX = -A_new\b(parameters.*perturbation_vector);
             %P =  A^-1*(-eye(3) + expm(A*tode(j)))*b - x0;
             P= A_new^-1*(-eye(Nchannels) + expm(A_new*tode(j)))*b(parameters.*perturbation_vector)  + expm(A_new*tode(j))*EX_original; %solve the ODE
             ode_sol(:,j) = P;

        end
        PerturbationAnalyticalSolution = horzcat(c_t*ode_sol_before, c_t*ode_sol);

     