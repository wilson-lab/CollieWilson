% velocity_histogram
% This function bins forward, angular, and sideways velocity data 
% and returns histograms of each based on user-specified normalization.
%
% INPUT
% fwdVel - forward velocity (vector)
% angVel - angular velocity (vector)
% sidVel - sideways velocity (vector)
% norm   - 0 for counts, 1 for normalized probability
%
% OUTPUT
% fwdHist - binned forward velocity [bin centers, counts/probabilities]
% angHist - binned angular velocity [bin centers, counts/probabilities]
% sidHist - binned sideways velocity [bin centers, counts/probabilities]
%
% 05/09/24 - MC created
% 10/23/24 - updated (added description and output details)
%
function [fwdHist,angHist,sidHist] = velocity_histogram(fwdVel,angVel,sidVel,norm)
%% initialize bin variables
% set bin variables
fwdMax = 30;
angMax = 250;
sidMax = 10;
fwdSize = 0.5;
angSize = 5;
sidSize = 0.25;

% set bin edges
fwdEdges = -5+fwdSize/2:fwdSize:fwdMax;
angEdges = -angMax+angSize/2:angSize:angMax;
sidEdges = -sidMax+sidSize/2:sidSize:sidMax;

% set bin centers
fwdBin = fwdEdges(1:end-1)+fwdSize/2;
angBin = angEdges(1:end-1)+angSize/2;
sidBin = sidEdges(1:end-1)+sidSize/2;

%% reshape and remove nans
fwdVel_all = reshape(fwdVel,[],1);
angVel_all = reshape(angVel,[],1);
sidVel_all = reshape(sidVel,[],1);

fwdVel_all(isnan(fwdVel_all)) = [];
angVel_all(isnan(angVel_all)) = [];
sidVel_all(isnan(sidVel_all)) = [];

%% invert turn velocities to account for R-L biases

angVel_all = [angVel_all -angVel_all];
sidVel_all = [sidVel_all -sidVel_all];

%% bin velocity data
if norm
    [nFwd,~] = histcounts(fwdVel_all,fwdEdges,'Normalization','probability');
    [nAng,~] = histcounts(angVel_all,angEdges,'Normalization','probability');
    [nSid,~] = histcounts(sidVel_all,sidEdges,'Normalization','probability');
else
    [nFwd,~] = histcounts(fwdVel_all,fwdEdges);
    [nAng,~] = histcounts(angVel_all,angEdges);
    [nSid,~] = histcounts(sidVel_all,sidEdges);
end

% store output arrays
fwdHist = [fwdBin' , nFwd'];
angHist = [angBin' , nAng'];
sidHist = [sidBin' , nSid'];

end