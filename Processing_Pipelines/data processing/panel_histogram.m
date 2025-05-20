% panel_histogram
% This function bins panel position and velocity data, generating histograms 
% for both variables. It can return normalized histograms based on user input.
%
% INPUTS:
%   panelpos - Array of panel position data (degrees)
%   panelvel - Array of panel velocity data (degrees/second)
%   norm     - Flag indicating whether to normalize the histograms (1 for normalized, 0 for raw counts)
%
% OUTPUTS:
%   posHist  - Binned panel position data (matrix with bin centers and counts/probabilities)
%   velHist  - Binned panel velocity data (matrix with bin centers and counts/probabilities)
%
% CREATED: 09/15/24 - MC
%
function [posHist,velHist] = panel_histogram(panelpos,panelvel,norm)
%% initialize bin variables
% set bin variables
posMax = 180;
velMax = 600;
posBinSize = 9;
velBinSize = 20;

% set bin edges
posEdges = -posMax+posBinSize/2:posBinSize:posMax;
velEdges = -velMax+velBinSize/2:velBinSize:velMax;

% set bin centers
posBins = posEdges(1:end-1)+posBinSize/2;
velBins = velEdges(1:end-1)+velBinSize/2;

%% reshape and remove nans
panelpos_all = reshape(panelpos,[],1);
panelvel_all = reshape(panelvel,[],1);

panelpos_all(isnan(panelpos_all)) = [];
panelvel_all(isnan(panelvel_all)) = [];

%% adjust for bias in fixation
biasHD = mean(panelpos_all,'all','omitnan');
panelpos_all = panelpos_all - biasHD;

%% invert to account for R-L biases

panelpos_all = [panelpos_all -panelpos_all];
panelvel_all = [panelvel_all -panelvel_all];

%% bin velocity data
if norm
    [nPos,~] = histcounts(panelpos_all,posEdges,'Normalization','probability');
    [nVel,~] = histcounts(panelvel_all,velEdges,'Normalization','probability');
else
    [nPos,~] = histcounts(panelpos_all,posEdges);
    [nVel,~] = histcounts(panelvel_all,velEdges);
end

% store output arrays
posHist = [posBins' , nPos'];
velHist = [velBins' , nVel'];

end