function [r_model, local_auc] = model_rotvel_from_dna02(mean_pulse_srRL, time_pulse, k, f, w)
% model_rotvel_from_dna02
% CREATED: 04/17/2025 - MC
%
% Computes causal moving window area under the curve (AUC) of R-L firing rate,
% scaled by gain factor k, thresholded by friction f, and delayed by 200 ms.
%
% Inputs:
%   mean_pulse_srRL : time x 1 x sweep
%   time_pulse      : time vector in ms
%   k               : gain multiplier
%   f               : friction threshold
%   w               : window size in timepoints (e.g. 44 for ~50ms)
%
% Outputs:
%   r_model         : thresholded, delayed scaled AUC (time x sweep)
%   local_auc       : raw causal AUC (time x sweep)

% Convert ms to seconds
time_sec = time_pulse / 1000;
dt = mean(diff(time_sec));

% Delay in samples (200 ms)
n_delay = round(0.200 / dt);

% Extract firing rate
rl_FR = squeeze(mean_pulse_srRL(:,1,:));  % time x sweep

% Instantaneous area
inst_area = rl_FR * dt;

% Causal moving sum (past window)
local_auc = movsum(inst_area, [w-1, 0], 1);

% Apply gain and friction
r_model = k * local_auc - f;
r_model(r_model < 0) = 0;

% Apply delay (pad with zeros at the start, trim from end)
r_model = [zeros(n_delay, size(r_model,2)); r_model(1:end-n_delay,:)];

end
