% calculateTotalSpeed
% analysis function that is fed the animal's directional speeds (forward,
% angular, and sideways), converts to common units, and then takes the
% absolute sum to generate an estimate of "total speed". Note that this
% does not account for the directional vector, but is more useful for
% estimating start/stop epochs.
%
% INPUT:
% forward - array containing forward speeds (mm/s)
% angular - array containing angular speeds (deg/s)
% sideway - array containing sidewyas speeds (mm/s)
%
% OUTPUT:
% totalSpeed
%
% CREATED: 12/2/2022 MC
%

function totalSpeed = calculateTotalSpeed(forward,angular,sideway)

% convert all directional speeds to degrees/sec
ballCirc = 2*pi*(9/2);

forward_d = (forward ./ ballCirc) .*360;
sideway_d = (sideway ./ ballCirc) .*360;
angular_d = angular;


% sum absolute speed
totalSpeed = abs(forward_d)+abs(sideway_d)+abs(angular_d);

end