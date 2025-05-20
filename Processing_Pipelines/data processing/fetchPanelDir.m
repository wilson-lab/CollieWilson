% fetchPanelDir
% This function determines the direction of the panel movement (left or right) for each trial
% based on the panel position data. It assigns a value of 1 for rightward movement and -1 
% for leftward movement. The function identifies motion pulses by detecting the start and stop 
% of sweeps and categorizes the direction of each sweep accordingly.
%
% INPUTS:
%   panelps   - Matrix of panel positions (time x trials)
%
% OUTPUTS:
%   paneldir  - Matrix of panel directions (1 for right, -1 for left; same size as panelps)
%
% CREATED: [Date] MC
%
function paneldir = fetchPanelDir(panelps)
%% initialize
% fetch trial number
nTrial = size(panelps,2);

% determine what left and right should be assigned to (e.g., 1 and -1)
setL = -1;
setR = 1;

%% determine panel directions for motion pulse experiment

% initialize
paneldir = zeros(size(panelps));

% for each trial
for t = 1:nTrial
    thispos = panelps(:,t);

    % find when motion pulse based on where sweeps start/stop
    findSweeps = ~isnan(thispos); %find sweeps (1) vs blanks (0)
    findStartStop = diff(findSweeps); %find both sweep starts (1) and stops (-1)
    findStarts = find(findStartStop>0)+1;
    findStops = find(findStartStop<0);
    % ensure there are no incomplete sweeps
    if length(findStarts)>length(findStops)
        findStarts = findStarts(1:length(findStops));
    end
    % store
    idxStartStop=[];
    idxStartStop(:,1) = findStarts;
    idxStartStop(:,2) = findStops;

    % for each pulse
    for p = 1:size(idxStartStop,1)
        thisSweep = thispos(idxStartStop(p,2)) - thispos(idxStartStop(p,1));
        if thisSweep>0 %right
            paneldir(idxStartStop(p,1):idxStartStop(p,2),t) = setR;
        else %left
            paneldir(idxStartStop(p,1):idxStartStop(p,2),t) = setL;
        end
    end

end
end