% SETPOINT_VELVTURN - This function analyzes the relationship between the fly's turn behavior (angular
% velocity) and panel velocity during setpoint fixation. It outputs the combined (RL) binned angular 
% velocity for all velocities and for the front field of view (±30°), across different conditions.

% INPUTS:
%   panelps     - 3D array of panel position data (yaw gain-modified), where each slice represents 
%                 a different condition and each column represents a trial.
%   panelvel    - 3D array of panel velocity data corresponding to panel position.
%   angular     - 3D array of angular velocity data representing the fly's turn behavior.
%   ttime       - Time vector corresponding to the temporal resolution of the data.
%   settings    - Structure containing parameters like visuomotor lag, pursuit gain, and plot settings.
%   optLag      - Optional flag (1/0) to shift angular data based on visuomotor lag.
%   optPlot     - Optional flag (1/0) to plot the results (1 for yes, 0 for no).

% OUTPUTS:
%   velvangRL          - Combined (RL) binned angular velocity for all velocity bins.
%   velvangFrontFOVRL  - Combined (RL) binned angular velocity for velocities within the front field of view (±30°).
%   velBins            - Velocity bins (center values) used for binning panel velocity data.

% Created: 09/03/24 by MC
% Updated: N/A

% The function:
% - Bins panel velocity data and computes the mean angular velocity for each bin.
% - Outputs binned angular velocity for all velocities and for the front field of view (±30°).
% - Optionally applies a visuomotor lag shift and plots the results across conditions.
%
function [velvangRL, velvangFrontFOVRL, velBins] = setpoint_velvturn(panelps, panelvel, angular, ttime, settings, optLag, optPlot)
%% initialize
% fetch number of conditions
nCond = size(panelvel,3);

% set velocity bin parameters
velMax = 500; %+/- deg/s
velBin = 20; %deg/s

% create velocity bins
velEdge = -velMax-velBin/2:velBin:velMax+velBin/2; % bin edges
velBins = -velMax:velBin:velMax; % bin labels (center)
nVelBins = length(velBins);

% Define front field of view limits
frontFOV = 20;

%% (optional) shift according to lag estimates
% if lag estimates were provided, shift
if optLag
    % fetch shift indices for each lag
    [idx_vm] = fetchTimeIdx(ttime,settings.visuomotorLag);
    idx_vm = idx_vm-1;
    
    % shift and exclude data at start/stop of trial
    angular = circshift(angular,-idx_vm,1); %shift
    angular(end-idx_vm+1:end,:) = nan; %exclude shifts
end

%% bin turn magnitude according to panel velocity

% initialize
velBinnedAngular = nan(nVelBins,nCond);
velBinnedAngularFrontFOV = nan(nVelBins,nCond); % For front field of view
minBin = 500;

% for each condition
for c = 1:nCond
    % fetch and reshape data
    thisPanelps = reshape(panelps(:,:,c),[],1);
    thisPanelvel = reshape(panelvel(:,:,c),[],1);
    thisAngular = reshape(angular(:,:,c),[],1);

    % discretize panel velocity data
    discPanelvel = discretize(thisPanelvel,velEdge,velBins);

    % calculate mean for each velocity bin
    for v = 1:nVelBins
        thisBin = velBins(v);
        thisBinIdx = find(discPanelvel==thisBin);

        % only include the bin if it has at least 3 points
        if length(thisBinIdx) >= minBin
            velBinnedAngular(v,c) = mean(thisAngular(thisBinIdx),'omitnan');
        else
            velBinnedAngular(v,c) = NaN; % exclude bin with fewer than 3 points
        end

        % Check if the panel position is within the front field of view (±30°)
        frontFOVIdx = thisBinIdx(thisPanelps(thisBinIdx) >= -frontFOV & thisPanelps(thisBinIdx) <= frontFOV);

        % Only include front FOV points if there are enough points
        if length(frontFOVIdx) >= minBin
            velBinnedAngularFrontFOV(v,c) = mean(thisAngular(frontFOVIdx),'omitnan');
        else
            velBinnedAngularFrontFOV(v,c) = NaN; % exclude bin with fewer than 3 points
        end
    end
end

% combine R and L to account for bias
RL_velvang = (velBinnedAngular + flip(-velBinnedAngular,1))./2;
RL_velvangFrontFOV = (velBinnedAngularFrontFOV + flip(-velBinnedAngularFrontFOV,1))./2;

% clean isolated points
minSequenceLength = 5;
RLclean_velvang = removeIsolatedPoints(RL_velvang, minSequenceLength);
RLclean_velvangFrontFOV = removeIsolatedPoints(RL_velvangFrontFOV, minSequenceLength);

% check to make sure no errors were made
for c = 1:nCond
    if isnan(RLclean_velvang(find(velBins == 0),c))
        RLclean_velvang(:,c) = nan;
    end
    if isnan(RLclean_velvangFrontFOV(find(velBins == 0),c))
        RLclean_velvangFrontFOV(:,c) = nan;
    end
end

% store
velvangRL = RLclean_velvang;
velvangFrontFOVRL = RLclean_velvangFrontFOV;

% optional plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1500 600]) % Adjust figure size for 3 rows
    tiledlayout(2,nCond,'TileSpacing','compact')  % 3 rows now

    % Plot all velocities for each condition
    for c = 1:nCond
        nexttile
        plot(velBins,velvangRL(:,c),'Color',settings.HDColor)
        xline(0); yline(0); ylim([-100 100])
        title([num2str(settings.pursuitGain(c)) 'X - All Velocities'])
    end

    % Plot front FOV velocities for each condition
    for c = 1:nCond
        nexttile
        plot(velBins,velvangFrontFOVRL(:,c),'Color',settings.HDColor)
        xline(0); yline(0); ylim([-100 100])
        title([num2str(settings.pursuitGain(c)) 'X - Front FOV'])
    end
end


end