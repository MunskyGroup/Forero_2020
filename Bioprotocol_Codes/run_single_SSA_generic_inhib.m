function [X_array] = run_single_SSA_generic_inhib(x0,S,W,Wafter, T_array,time_var,signal_update_rate,inhib_t)
% Start the simulation.
t = 0;   % initial time of simulation.
x = x0;
iprint = 1;  % The next time at which to record results.
Nsp = size(S,1);  %Number of species.
Nt = length(T_array);
X_array = zeros(Nsp,Nt);

if time_var
        S = [zeros(Nsp,1),S];
        props(1) = signal_update_rate;
        jt = 1;
else
    jt=0;
end
            
perturbed = 0;
while t<max(T_array)

    %% Choose time of reaction
   % for i=length(W):-1:1
    %    props(i+jt) = W{i}(x(1),x(2),x(3),t);     % evaluate the propensity functions at the current state.
   % end
    
   if t >= inhib_t  && perturbed == 0
       
        W = Wafter;
        
%         kon = parameters(1)*inhib_per(1);
%         koff = 1000;%parameters(2);
%         kesc = parameters(3)*inhib_per(3);
%         kproc = parameters(4)*inhib_per(4);
%         kin =  1000*parameters(5)*inhib_per(5);
%         kout = parameters(6)*inhib_per(6);
%         frac = parameters(7)*inhib_per(7);
%         eta_rnap = parameters(8)*inhib_per(8);
%         eta_ser5 = parameters(9)*inhib_per(9);
%         eta_ts = parameters(10)*inhib_per(10);
% 
%        W1 = zeros(6,3);
%        W0 = zeros(6,1);
%        W1(1,1) = -kon; W0(1,1) = kon;
%        W1(2,1) = koff;    
%        W1(3,1) = kin; 
%        W1(4,2) = kout;
%        W1(6,3) = kproc; 
%        W1(5,2) = kesc*frac; 
%        W = @(x) W1*x + W0;
%        perturbed = 1;
   end
   
   
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
    x = x + S(:,j);
end

