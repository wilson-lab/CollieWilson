% computePanelVelocity
% This function calculates the velocity of the panel position data by taking the 
% gradient of the smoothed panel positions. It unwraps the angles to avoid discontinuities 
% and applies a Gaussian smoothing filter to the position data before calculating the 
% velocity. The function outputs a 3D array of panel velocities for each trial and condition.
%
% INPUTS:
%   panelps - Matrix of panel positions (time x trials x conditions) in degrees
%   ttime   - Time vector (in seconds) corresponding to the panel positions
%
% OUTPUTS:
%   panelvel - Matrix of calculated panel velocities (time x trials x conditions) in degrees/second
%
% CREATED: [Date] MC
%
function panelvel = computePanelVelocity(panelps,ttime)
%% initialize
nCond = size(panelps,3);
nTrial = size(panelps,2);

% fetch sample rate
sampRate = fetchTimeIdx(ttime,1)-1;

%% compute velocity and smooth
% initialize
panelvel = nan(size(panelps));
% for each condition
for c = 1:nCond
    % for each trial
    for t = 1:nTrial
        % unwrap and smooth
        thispos = unwrap(deg2rad(panelps(:,t,c)));
        smoothpos = gaussSmooth(thispos, 500, 100);
        smoothpos_deg = rad2deg(smoothpos);
        % calculate velocity
        thisvel = gradient(smoothpos_deg) .* sampRate;
        % store
        panelvel(:,t,c) = thisvel;
    end
end
% test
% close all;plot(smoothpos_deg); hold on; plot(thisvel)
end