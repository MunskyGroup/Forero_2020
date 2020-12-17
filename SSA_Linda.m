%% Run a SSA for Lindas model

load('best_simple_pars.mat')  %load parameter file
real_valscl = parameters;

kon = parameters(1);
koff = 1000;%parameters(2);
kesc = parameters(3);
kproc = parameters(4);
kin =  1000*parameters(5);
kout = parameters(6);
frac = parameters(7);
eta_rnap = parameters(8);
eta_ser5 = parameters(9);
eta_ts = parameters(10);

%%

% State 1: On/ Off
% State 2: Pol2/Ser5 these are merged for the purpose of this mode, the
%          Resolution wasnt enough to tell these two signals apart in model
%          fitting
% State 3: TS+Pol2+Ser5 


% [ off ]  <-koff -kon-> [on] --kin--> [Pol2+ser5] --kesc--> [Pol2+Ser5+Ts] -kproc-> out
%                                         \/                   
%                                        kout

Nstates = 3;  %set up propensity and stoich for lindas model
b = zeros(Nstates,1);
b(1) = kon;

c = zeros(3,3);
c(1,1:3)=[0,1,1];
c(2,1:3)=[0,frac,1];
c(3,3)=1;

S = zeros(3,6);
W1 = zeros(6,3);
W0 = zeros(6,1);
S(1,1) = 1;  W1(1,1) = -kon; W0(1,1) = kon;
S(1,2) = -1; W1(2,1) = koff; 

S(2,3) = 1; W1(3,1) = kin; 
S(2,4) = -1; W1(4,2) = kout; 

S(2:3,5) = [-1;1]; W1(5,2) = kesc*frac; 
S(3,6) = -1; W1(6,3) = kproc; 


x0 = [0,0,0]';
T_array = [0:1:1000];
time_var = 0;
signal_update_rate = 0;

W = @(x) W1*x + W0;


%%
%Solve a singal trajectory 
sol = run_single_SSA_linda(x0,S,W,T_array,time_var,signal_update_rate);  

pol2_ssa = sol(2,:)';
ser5_ssa = sol(2,:)';
ts_ssa = sol(3,:)';

[pol2_ssa_intensity,ser5_ssa_intensity,ts_ssa_intensity,~] = get_model_intesities(sol,eta_rnap,eta_ser5,eta_ts); %convert molecules to signal
[pol2norm,ser5norm,tsnorm] = Normalize_simulated_intensities(.95,pol2_ssa_intensity,ser5_ssa_intensity,ts_ssa_intensity);

close all

figure
plot(pol2_ssa+ts_ssa,'r'); hold on; plot(ts_ssa,'b'); 
xlabel('time')
ylabel('molecules')
legend('pol2','ts')

figure

plot(pol2norm,'r'); hold on; plot(ser5norm,'g'); plot(tsnorm,'b'); 
xlabel('time')
ylabel('norm int')
legend('pol2','ser5','ts')
figure

plot(pol2_ssa_intensity,'r'); hold on; plot(ser5_ssa_intensity,'g'); plot(ts_ssa_intensity,'b'); 
xlabel('time')
ylabel('int')
legend('pol2','ser5','ts')



