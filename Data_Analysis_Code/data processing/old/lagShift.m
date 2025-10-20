% lagShift
% simple function for shifting variable according to calculated lag
%
% INPUT
% thisVar - any experimental variable (e.g., velocity, object pos, etc)
% t - time
% lag - +/- amount of lag to be shifted by
% OUTPUT
% shiftVar - experimental variable shifted by lag of interest
%
% CREATED 4/23/24 MC
%

function [shiftVar] = lagShift(thisVar,t,lag)
%% initialize

% fetch number of trials
nTrial = size(thisVar,2);

% fetch lag shift index
[shiftIdx] = fetchTimeIdx(t,lag);
shiftIdx = shiftIdx-1;

%% shift variable

% initialize output matrix
shiftVar = zeros(size(thisVar));
% for each trial
for n = 1:nTrial
    thisShift = circshift(thisVar(:,n),-shiftIdx,1);
    % nan out shifted data
    thisShift(end-shiftIdx+1:end) = nan;
    % store in output matrix
    shiftVar(:,n) = thisShift;
end

end