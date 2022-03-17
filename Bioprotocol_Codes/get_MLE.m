function [best_par, min_err] = get_MLE(get_LL_fun, parameter_guess, par_changed, search_chains, Constraints, GA_pop, save_file_name)
    
    % Function that finds the MLE from subsequent rounds of GA and fminsearch

    % get_LL_fun - autonomous function to get the LL that only varies on logspace parameters that change
    % parameter_guess - initial full parameter guess 
    % par_changed - indicies of parameters that are allowed to change
    % search_chains - number of search iterations to run 
    % Constraints - constraints upper and lower bound for the GA sweep
    % GA_pop - Genetic algorithm population
    % save_file_name - file name to save the best found parameters as


    par_opt = log10(parameter_guess(par_changed));
    min_err = 1e7;
    
    
    for i=1:search_chains

        try
            parpool
        catch
        end
        x = par_opt;
        size(x)
        
        % Setup and run the GA
        pctRunOnAll warning('off', 'all')
        Mut_Fun = @(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation)Mutation(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation,Constraints);
        Npop = GA_pop;
        gaopt = gaoptimset('Display','final','useparallel',1,...
            'MutationFcn',Mut_Fun,...
            'PopulationSize',Npop,...
            'EliteCount',1,'Generations',200,'CrossoverFraction',0);  %add output function
        gaopt.InitialPopulation = [repmat(x,Npop,1)+randn(Npop,length(x))];
    
        [x,~,~,~,~,~] = ga(get_LL_fun,length(x),gaopt); % Run the G.A.
        par_opt = x;
        
        % Setup and run the fminsearch
        options = optimset('display','iter','MaxIter',5000);
        par_opt = fminsearch(get_LL_fun,par_opt,options);

        % save the parameters if the min error was less than the best one found so far
        err = get_LL_fun(par_opt);
        if sum(err)<min_err 
            min_err=sum(err)
            parameter_guess(par_changed) = 10.^par_opt;
            save(save_file_name,'parameter_guess')
            best_par = parameter_guess;
           
        end
    end
    
    end
    
    function mut_Chil = Mutation(parents,~,~,~,~,~,this_Pop,Constrs)
    % Custom Mutation function for G.A. search
    % Choose which parameters to permute from the parents.
    PTB = rand(size(this_Pop(parents,:)))>0.3;
    while min(max(PTB,[],2))==0
        J = find(max(PTB,[],2)==0);
        PTB(J,:) = rand(size(PTB(J,:)))>0.3;
    end
    
    % Make the permutation by adding random numbers to the perturbed
    % parameters.
    mut_Chil = this_Pop(parents,:)+...
        PTB.*randn(size(this_Pop(parents,:)))/10^(5*rand)+...
        10^(-3*rand)*randn*(rand(size(this_Pop(parents,:)))>0.95);
    % mut_Chil = this_Pop(parents,:)+...
    %     0.08*PTB.*randn(size(this_Pop(parents,:)));
    
    % This mutation function takes the original parents, then chooses 80% of
    % the number to mutate, then mutates these by a normally distributed random
    % variable (multiplicatively). Finally, we add another small normally distributed
    % random variable (again to 50%) in order to push the values away from
    % zero.
    
    % Flip the signs for some of the random parameter values.
    % FLP = ones(size(mut_Chil))-2*(rand(size(mut_Chil))>0.99);
    % mut_Chil=mut_Chil.*FLP;
    for i=1:size(mut_Chil,1)
        mut_Chil(i,:) = max([mut_Chil(i,:);Constrs.LB]);
        mut_Chil(i,:) = min([mut_Chil(i,:);Constrs.UB]);
    end

end