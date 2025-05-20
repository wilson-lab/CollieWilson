% auto_xcorr
% 
% Calculates the autocorrelation for each directional velocity (forward, angular, and 
% sideways) over a specified time window.
%
% INPUTS:
% forward - Array of forward velocities (units: m/s).
% angular - Array of angular velocities (units: deg/s).
% sideway - Array of sideways velocities (units: m/s).
% ttime   - Time array (units: seconds).
%
% OUTPUTS:
% r_val   - Struct containing r values across the autocorrelation window for 
%            forward, angular, and sideways velocities.
% lag_t   - Array of times (in milliseconds) corresponding to the autocorrelation 
%            window.
%
% PROCESS:
% The function initializes by determining the number of trials and setting the size 
% of the autocorrelation window. It prepares the data by adding a buffer window to 
% the velocity arrays and reshaping them for autocorrelation calculations. The 
% autocorrelation is then computed for each directional velocity using a specified 
% function (xcorrWGaps) that accommodates gaps in the data. The lag indices are 
% converted to time in milliseconds for output.
%
% CREATED: 11/17/23 by MC
%
function [r_val,lag_t] = auto_xcorr(forward,angular,sideway,ttime)
%% initialize

% pull number of trials
nTrials = size(forward,2);

% set size of xcorr window
xc_t = 2; %s, total (e.g., 1/2 on each side of 0)
[xc_window] = fetchTimeIdx(ttime,xc_t);

%% prepare data

% add buffer window
forward_buff = [forward ; nan(xc_window*2,nTrials)];
angular_buff = [angular ; nan(xc_window*2,nTrials)];
sideway_buff = [sideway ; nan(xc_window*2,nTrials)];

% concatenate
forward_cc = reshape(forward_buff,[],1);
angular_cc = reshape(angular_buff,[],1);
sideway_cc = reshape(sideway_buff,[],1);

%% run xcorr

% for each directional velocity
disp('Running auto corr...')
[r_val.fwd, lag_fwd] = xcorrWGaps(forward_cc,forward_cc,xc_window);
[r_val.ang, lag_ang] = xcorrWGaps(angular_cc,angular_cc,xc_window);
[r_val.sid, lag_sid] = xcorrWGaps(sideway_cc,sideway_cc,xc_window);

% convert lag index to time (msec)
lag_t = (lag_fwd.*ttime(2))*1000;


end