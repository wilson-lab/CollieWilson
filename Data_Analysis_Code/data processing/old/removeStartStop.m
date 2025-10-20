% removeStartStop
% Analysis function that uses the animal's total speed to find and blank
% all start/stops using a flexible window provided by the user
%
% INPUT:
% forward - array containing forward speeds (mm/s)
% angular - array containing angular speeds (deg/s)
% sideway - array containing sidewyas speeds (mm/s)
% spikerate
% removeWindow - size of start/stop window to blank out
%
% OUTPUT:
% totalSpeed
%
% CREATED: 12/2/2022 MC
%

function [forwardSS,angularSS,sidewaySS,spikertSS] = removeStartStop(forward,angular,sideway,spikert,removeWindow)
%% set parameters
% set threshold for denoting when the animal is moving vs not moving
%minSpeedThresh = 2; %deg/s
minSpeedThresh = 0.1; %mm/s

%% calculate total speed
%totalSpeed = calculateTotalSpeed(forward,angular,sideway);
totalSpeed=abs(forward);

%% find start/stops
% set all rest points to 0 and move points to 1
[m,n] = size(totalSpeed);
ssMatrix = zeros(m,n); %initialize
ssMatrix(totalSpeed>=minSpeedThresh)=1;
% find changes between rest and move (+1 = start, -1 = stop)
ssDiff = diff(ssMatrix,1,1);
ssDiff = [ssDiff; zeros(1,n)]; %match dimensions by adding buffer

% create tracker for finding start/stops
ssTracker = zeros(m,n); %initialize

% for each column (aka trial)
for nt = 1:n
    % find starts in this trial (+1)
    startIdx = find(ssDiff(:,nt)==1);
    for s=1:length(startIdx)
        % select this start index
        thisIdx = startIdx(s);
        % define this start window
        winStart = thisIdx-removeWindow;
        winEnd = thisIdx+removeWindow;
        % ensure window does not pass trial boundaries
        if winStart<1
            winStart=1;
        end
        if winEnd>m
            winEnd=m;
        end
        % blank out start through window
        ssTracker(winStart:winEnd,nt)=1;
    end

    % find stops in this trial (-1)
    stopIdx = find(ssDiff(:,nt)==-1);
    for s=1:length(stopIdx)
        % select this stop index
        thisIdx = stopIdx(s);
        % define this window
        winStart = thisIdx-removeWindow;
        winEnd = thisIdx+removeWindow;
        % ensure window does not pass trial boundaries
        if winStart<1
            winStart=1;
        end
        if winEnd>m
            winEnd=m;
        end
        % blank out start through window
        ssTracker(winStart:winEnd,nt)=1;
    end
end

ssIdx = find(ssTracker==1);

%% extract directional velocities w/o start/stops

% initialize 
forwardSS = forward;
angularSS = angular;
sidewaySS = sideway;
spikertSS = spikert;

% remove start/stops
forwardSS(ssIdx)=NaN;
angularSS(ssIdx)=NaN;
sidewaySS(ssIdx)=NaN;
spikertSS(ssIdx)=NaN;


end