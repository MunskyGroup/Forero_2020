function [sigred,TT,means,SIGc,TT2,mins] = get_ac_and_cc_mod_simplified(parameters,TT)
if nargin<2
    TT = [0:30];
end

parameters = max(1e-8,parameters);  % enforce positive parameter constraint
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

if frac > 1
    frac = 1;
end

%%

Nstates = 3;
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

A=S*W1;

%% Find steady state covariance:
EX = -A\b;
means = c*EX;
Q = S*diag(W1*EX+W0)*S';
phi = S*W1;
try
SIG = lyap(phi,Q);
catch 
    SIG=eye(length(phi));
end

%%
x0 = zeros(length(A)^2,1);
x0(:) = SIG;
SIGc = c*SIG*c';
%%
n = length(phi);
PHI = spalloc(n^2,n^2,n^3);
for i=1:length(phi)
    PHI((i-1)*n+1:i*n,(i-1)*n+1:i*n)=phi;
end
%%
fun = @(t,x)PHI*x;
options = odeset('jacobian',PHI);
[TT,YY] = ode23s(fun,TT,x0,options);

%%
SIGt = 0*SIG;
sigred = zeros(size(c,1)^2,length(TT));
for i=1:length(TT)
    SIGt(:) = YY(i,:);
    tmp = c*SIGt*c';
    sigred(:,i) = tmp(:);
end

%% Scale
sgr = sigred;
sc1 = parameters(11);
sc2 = parameters(12);
sc3 = parameters(13);
sigred(1,:) = sgr(1,:)/sgr(1,1);
sigred(2,:) = sc1*sgr(2,:)/sqrt(sgr(1,1)*sgr(5,1));
sigred(3,:) = sc2*sgr(3,:)/sqrt(sgr(1,1)*sgr(9,1));
sigred(4,:) = sc1*sgr(4,:)/sqrt(sgr(5,1)*sgr(1,1));
sigred(5,:) = sgr(5,:)/sgr(5,1);
sigred(6,:) = sc3*sgr(6,:)/sqrt(sgr(5,1)*sgr(9,1));
sigred(7,:) = sc2*sgr(7,:)/sqrt(sgr(9,1)*sgr(1,1));
sigred(8,:) = sc3*sgr(8,:)/sqrt(sgr(9,1)*sgr(5,1));
sigred(9,:) = sgr(9,:)/sgr(9,1);

%% Add shot noise
sigred([1 5 9],1) = sigred([1 5 9],1)+([eta_rnap,eta_ser5,eta_ts].^2)';

%%
if nargout>4
    x = EX;
    x(1) = 0;
    fun = @(t,x)A*x+b;
    options = odeset('jacobian',A);
    [TT2,YY] = ode23s(fun,[0:0.1:20],x,options);
    mins = c*YY';
    TT2 = [-5:0.1:0,TT2'];
    mins = [repmat(c*EX,1,51),mins];
end
