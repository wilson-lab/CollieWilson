% spikert_xcorr
% determine optimal lag between spike rate and visual target
%
% INPUT
% spikert - array of spike rates
% panelps - array of panel positions
%
% OUTPUT
% r_val - r values across xcorr window
% lag_t - times (msec) of xcorr window
% lag_opt - optimal lag time (msec) based on r peak
%
% CREATED 04/23/24 MC from spikert_xcorr
%

function [r_val,lag_t] = spikert_xcorr2(spikert,panelps,ttime)
%% initialize

% pull number of trials
nTrials = size(spikert,2);

% set size of xcorr window
xc_window = 4000; %index

%% prepare data

% add buffer window
panelps_buff = [panelps ; nan(xc_window*2,nTrials)];
spikert_buff = [spikert ; nan(xc_window*2,nTrials)];

% concatenate
panelps_cc = reshape(panelps_buff,[],1);
spikert_cc = reshape(spikert_buff,[],1);

%% run xcorr

% for each directional velocity
disp('Running xcorr...')
[r_val, lag_obj] = xcorrWGaps(panelps_cc,spikert_cc,xc_window);
disp('Complete.')

% convert lag index to time (msec)
lag_t = (lag_obj.*ttime(2))*1000;


end