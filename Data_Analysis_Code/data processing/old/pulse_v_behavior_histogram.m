% pulse_v_spikert_histogram
%
% data processing function, for generating velocity distribution for each
% pulse separately
%
% INPUT
% panelps
% angular
% forward
% sideways
%
% OUTPUT
% pulseHist - binned histogram for directional velocities
%
% CREATED   05/23/2024 - MC
%
function pulseHist = pulse_v_behavior_histogram(panelps,forward,angular,sideways)
%% initialize
nTrial = size(panelps,2);
thisSpeed = 1;
nSpeeds = 1;

%% pull motion pulses from pseudorandomized dataset

% initialize
bufferOptions = [1650 1000];
bufferWindow = bufferOptions(thisSpeed);

% for each trial
for t = 1:nTrial
    % pull this panel data
    thisTrial = panelps(:,t);
    try
        % run ordering function to find motion pulses in order
        [thisOrderR,thisOrderL] = order_motion_pulse(thisTrial,nSpeeds,thisSpeed,bufferWindow);
        % if any indices fall outside of trial duration, fix
        thisOrderR(thisOrderR>size(panelps,1)) = size(panelps,1);
        thisOrderL(thisOrderL>size(panelps,1)) = size(panelps,1);

        % pull and store data based on ordered motion pulse indices
        % rows = data, columns = trials, z = sweeps
        for p = 1:size(thisOrderR,2)
            pulsePanelpsR(:,t,p) = panelps(thisOrderR(:,p),t);
            pulsePanelpsL(:,t,p) = panelps(thisOrderL(:,p),t);

            pulseForwardR(:,t,p) = forward(thisOrderR(:,p),t);
            pulseForwardL(:,t,p) = forward(thisOrderL(:,p),t);

            pulseAngularR(:,t,p) = angular(thisOrderR(:,p),t);
            pulseAngularL(:,t,p) = angular(thisOrderL(:,p),t);

            pulseSidewayR(:,t,p) = sideways(thisOrderR(:,p),t);
            pulseSidewayL(:,t,p) = sideways(thisOrderL(:,p),t);
        end
    catch
        % pull and store data based on ordered motion pulse indices
        % rows = data, columns = trials, z = sweeps
        for p = 1:size(thisOrderR,2)
            pulsePanelpsR(:,t,p) = nan;
            pulsePanelpsL(:,t,p) = nan;

            pulseForwardR(:,t,p) = nan;
            pulseForwardL(:,t,p) = nan;

            pulseAngularR(:,t,p) = nan;
            pulseAngularL(:,t,p) = nan;

            pulseSidewayR(:,t,p) = nan;
            pulseSidewayL(:,t,p) = nan;
        end
    end
end
nPulse = size(pulseAngularL,3);

%% calculate mean

% omit trials where fly was largely stationary
for p = 1:nPulse
    omitTrialsR = find(sum(isnan(pulseAngularR(:,:,p)))>200);
    omitTrialsL = find(sum(isnan(pulseAngularL(:,:,p)))>200);
    pulseForwardR(:,omitTrialsR,p) = nan;
    pulseForwardL(:,omitTrialsL,p) = nan;
    pulseAngularR(:,omitTrialsR,p) = nan;
    pulseAngularL(:,omitTrialsL,p) = nan;
    pulseSidewayR(:,omitTrialsR,p) = nan;
    pulseSidewayL(:,omitTrialsL,p) = nan;
end

% calculate means for each
meanPanelps(:,1,:) = mean(pulsePanelpsR,2,"omitnan");
meanPanelps(:,2,:) = mean(pulsePanelpsL,2,"omitnan");


%% analyze distribution for each sweep position separately

for p = 1
    % fetch velocities for this panel position
    this_forward = 
% analyze velocity distributions
[fwdHist,angHist,sidHist] = velocity_histogram(int_forward,int_angular,int_sideway,hNorm);
end

end