function [pol2_ssa_w_shot,ser5_ssa_w_shot,ts_ssa_w_shot,ts_ssa] = ...
    get_model_intesities(sol,eta_rnap,eta_ser5,eta_ts)

% Get integer numbers from simulations
pol2_ssa = sol(2,:)' + sol(3,:)';
ser5_ssa = pol2_ssa;
ts_ssa = sol(3,:)';

pol2_ssa_w_shot = add_shot(pol2_ssa,eta_rnap,30/200);
ser5_ssa_w_shot = add_shot(ser5_ssa,eta_ser5,23/200);
ts_ssa_w_shot = add_shot(ts_ssa,eta_ts,5/200);
% The third inputs are the fraction of data points that were assigned as
% negative values.

end

function ssa_w_shot = add_shot(ssa,eta,tr)
std_mod = std(ssa);
shot_mod = std_mod*eta;
ssa_w_shot = ssa + randn(size(ssa))*shot_mod;

% Define fraction tr as zero
tmp = sort(ssa_w_shot);
thresh = tmp(ceil(length(tmp)*tr));
ssa_w_shot = ssa_w_shot-thresh;
end
