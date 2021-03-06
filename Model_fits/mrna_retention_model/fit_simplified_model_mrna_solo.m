clear all
close all
clc
addpath ../Data_files/
addpath ../
%%
[mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap] = load_normalization_variance_gui(1,'G0_intp');

sig_dat = [rnap.mn_ac(1),ser5_rnap.mn_cc(11),mrna_rnap.mn_cc(11);...
    ser5_rnap.mn_cc(11),ser5.mn_ac(1),mrna_ser5.mn_cc(11);...
    mrna_rnap.mn_cc(11),mrna_ser5.mn_cc(11),mrna.mn_ac(1)];
sig_dat_sem = [rnap.sem_ac(1),ser5_rnap.sem_cc(11),mrna_rnap.sem_cc(11);...
    ser5_rnap.sem_cc(11),ser5.sem_ac(1),mrna_ser5.sem_cc(11);...
    mrna_rnap.sem_cc(11),mrna_ser5.sem_cc(11),mrna.sem_ac(1)];

mrna_rnap.sem_cc = mrna_rnap.sem_cc/mrna_rnap.mn_cc(11);
mrna_rnap.mn_cc = mrna_rnap.mn_cc/mrna_rnap.mn_cc(11);

mrna_ser5.sem_cc = mrna_ser5.sem_cc/mrna_ser5.mn_cc(11);
mrna_ser5.mn_cc = mrna_ser5.mn_cc/mrna_ser5.mn_cc(11);

ser5_rnap.sem_cc = ser5_rnap.sem_cc/ser5_rnap.mn_cc(11);
ser5_rnap.mn_cc = ser5_rnap.mn_cc/ser5_rnap.mn_cc(11);
%%
load best_simple_pars
parameters = ones(1,13)*2;
parameters(2) = 1;
parameters(7) = 1;
parameters(14) = .3;
parameters(15) = .3;
par_fixed = parameters;
par_changed = [1,3:6,8:15];
par_opt = log10(par_fixed(par_changed));
get_err = @(pars)sum(get_log_l_simplified_mrna(10.^pars,mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap,par_fixed,par_changed));
min_err=inf;
%%
for i=1:20
%         try
%             parpool
%         catch
%         end
%         x = par_opt;
%         Constrs.LB=-5*ones(size(x));
%         Constrs.UB=5*ones(size(x));
%         % Constrs.UB(2) = 0;
%     
%         pctRunOnAll warning('off', 'all')
%     
%         Mut_Fun = @(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation)Mutation(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation,Constrs);
%         Npop = 200;
%         gaopt = gaoptimset('Display','final','useparallel',1,...
%             'MutationFcn',Mut_Fun,...
%             'PopulationSize',Npop,...
%             'EliteCount',1,'Generations',200,'CrossoverFraction',0);  %add output function
%         gaopt.InitialPopulation = [repmat(x,Npop,1)+randn(Npop,length(x))];
%     
%         [x,~,~,~,~,~] = ga(get_err,length(x),gaopt); % Run the G.A.
%         par_opt = x;
    
    options = optimset('display','iter','MaxIter',5000);
    par_opt = fminsearch(get_err,par_opt,options);
    err = get_log_l_simplified_mrna(10.^par_opt,mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap,par_fixed,par_changed);
    if sum(err)<min_err
        min_err=sum(err)
        parameters(par_changed) = 10.^par_opt;
        save best_simple_pars_mrna_solo parameters
        
        make_plots(parameters,mrna,ser5,rnap,mrna_rnap,mrna_ser5,ser5_rnap,sig_dat,sig_dat_sem,err,[0 0 0])
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
