function [X_array,burst_times,pol2_per_burst,burst_t,mrna_cnts] = run_single_SSA_recorder_bursting(x0,S,W,T_array,time_var,signal_update_rate)
% Start the simulation.
t = 0;   % initial time of simulation.
x = x0;
iprint = 1;  % The next time at which to record results.
Nsp = size(S,1);  %Number of species.
Nt = length(T_array);
X_array = zeros(Nsp,Nt);
recorded_pol2_arrivals = 0;
aborted = 0;
escaped = 0;
escaped_times = [];
if time_var
        S = [zeros(Nsp,1),S];
        props(1) = signal_update_rate;
        jt = 1;
else
    jt=0;
end
   
burst_times = 0;
pol2_per_burst = [];
turned_on = 0;
pol2_cnt = 0;
burst_t = [];
mrna_cnts = [];
mrna_cnt = 0;
while t<max(T_array)

    %% Choose time of reaction
   % for i=length(W):-1:1
    %    props(i+jt) = W{i}(x(1),x(2),x(3),t);     % evaluate the propensity functions at the current state.
   % end
    
    props = W(x);
    w0 = sum(props);  % sum the propensity functions (inverse of ave. waiting time).
    tau = -1/w0*log(rand); % The time until the next reaction.
   
    %% update time
    t = t+tau;  % The time of the next reaction.
    
    while iprint<=Nt&&t>T_array(iprint)
        X_array(:,iprint) = x;
        iprint=iprint+1;
    end  
        
    if t>max(T_array)
        break;
    end
    
    %% Choose which reaction.
    r2 = w0*rand;    % this is a uniform r.v. betweeen zero and w0;
    j = 1;
    while sum(props(1:j))<r2
        j=j+1;
    end
    % At this point j is the chosen reaction.
    
    %% Update state
    
    
    if j == 1
        burst_times = burst_times+1;
        turned_on = 1;
        burst_t = [burst_t, t];
        mrna_cnts = [mrna_cnts, mrna_cnt];
        mrna_cnt = 0;
    end
    
    if j == 2
        turned_on = 0;
        pol2_per_burst = [pol2_per_burst,pol2_cnt];
        pol2_cnt = 0;
    end
    
    if turned_on == 1
        if j == 3
           pol2_cnt = pol2_cnt+1; 
        end
        
    end
    
    if j == 5
       mrna_cnt = mrna_cnt+1; 
    end

    
    x = x + S(:,j);
    
    

end

end
