% function that converts panel data into panel direction data
%
% INPUT
% panelps - panel position data
%
% OUTPUT
% paneldir - panel direction data
%
% Created 4/19/24 MC
%

function [paneldir] = panelDirection(panelps)
%% initialize
trialInfo = size(panelps);

% round panelps to reduce noise
panelps = round((panelps*2))/2;

%% pull sweep directions
paneldir = ones(trialInfo);

% for each trial
for t = 1:trialInfo(2)
    try
        % fetch data
        this_panelps = panelps(:,t);
        % calculate max/min
        maxPos = max(this_panelps,[],'all');
        minPos = min(this_panelps,[],'all');

        % pull max/min indices
        maxIdx = find(this_panelps==maxPos);
        minIdx = find(this_panelps==minPos);

        % indices repeat due to framerate, so find changes instead
        maxIdx = [maxIdx(diff(maxIdx)>5); maxIdx(end)];
        minIdx = [minIdx(diff(minIdx)>5); minIdx(end)];
        minIdx = [minIdx; trialInfo(1)]; %add until end for last sweep

        % for each max (right) to min (left), set direction to -1 (left)
        for m = 1:length(maxIdx)
            leftSweep = maxIdx(m):minIdx(m);
            paneldir(leftSweep,t) = -1;
        end
    catch
        paneldir(:,t) = nan;
    end
end

end