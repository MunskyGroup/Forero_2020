function [X_array,recorded_pol2_arrivals,aborted,escaped,escaped_times] = run_single_SSA_linda_recorder(x0,S,W,T_array,time_var,signal_update_rate)
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
    
    
    if j == 3 
        recorded_pol2_arrivals = recorded_pol2_arrivals + 1;
        escaped_times = [escaped_times, t];
    end
    
    if j == 4
        aborted = aborted + 1;
    end
    if j == 5
       escaped = escaped + 1; 
    end

    
    x = x + S(:,j);
    
    

end

end
