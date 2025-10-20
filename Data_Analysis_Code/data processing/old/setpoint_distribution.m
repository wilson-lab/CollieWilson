function [hBins,hData,setpointProb] = setpoint_distribution(panelps)
%% initialize
% pull trial info
nTrials = size(panelps,2); %number of trials
nTypes = size(panelps,3); %number of trial types

% generate histogram bins
psBin = 5;
psEdges = -180:psBin:180;

%% analyze distribution around 0

for t = 1:nTypes
    % bin paneldata
    [h,e] = histcounts(panelps(:,:,t),psEdges,'Normalization','probability');
    % pull probability of being within +/- 5 degrees of setpoint
    ec = e(1:end-1)+psBin/2; %center bins
    setpointProb(t,1) = sum(h(abs(ec)<=5));
    % store histogram data
    hData(:,t) = h';
end
% store bins
hBins = ec';

end