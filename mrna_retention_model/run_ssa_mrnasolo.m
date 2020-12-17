

load('best_simple_pars.mat');

parameters(14) = .3;
parameters(15) = .5;

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

k_release_mrna = parameters(14);
k_proc_mRNAsolo = parameters(15);

%%

Nstates = 4;
Nrxns = 8;
b = zeros(Nstates,1);
b(1) = kon;

c = zeros(Nstates,Nstates);
c(1,1:Nstates)=[0,1,1,1];
c(2,1:Nstates)=[0,frac,1,1];
c(Nstates,Nstates)=1;

S = zeros(Nstates,Nrxns);
W1 = zeros(Nrxns,Nstates);
W0 = zeros(Nrxns,1);
S(1,1) = 1;  W1(1,1) = -kon; W0(1,1) = kon;
S(1,2) = -1; W1(2,1) = koff; 

S(2,3) = 1; W1(3,1) = kin; 
S(2,4) = -1; W1(4,2) = kout; 

S(2:3,5) = [-1;1]; W1(5,2) = kesc*frac; 
S(3,6) = -1; W1(6,3) = kproc; 


S(3,7) =  -1; W1(7,3) = k_release_mrna;

S(4,7) = 1;

S(4,8) = -1; W1(8,4) = k_proc_mRNAsolo;


%%



x0 = [0,0,0,0]';
T_array = [0:1:1000];
time_var = 0;
signal_update_rate = 0;

W = @(x) W1*x + W0;

%Solve a singal trajectory 
sol = run_single_SSA_linda(x0,S,W,T_array,time_var,signal_update_rate);  

pol2_ssa = sol(2,:)';
ser5_ssa = sol(2,:)';
ts_ssa = sol(3,:)' + sol(4,:)';

figure
plot(sol(2,:)')
hold on;
plot(sol(3,:)')
plot(sol(4,:)')
legend('ctdserp5ph','tricolor','mrnasolo')
[pol2_ssa,ser5_ssa,ts_ssa,~] = get_model_intesities_mrnasolo(sol,eta_rnap,eta_ser5,eta_ts); %convert molecules to signal
[pol2norm,ser5norm,tsnorm] = Normalize_simulated_intensities(.95,pol2_ssa,ser5_ssa,ts_ssa);

figure

plot(pol2norm(500:700), 'r')
hold on
plot(ser5norm(500:700), 'g')
plot(tsnorm(500:700), 'b')


