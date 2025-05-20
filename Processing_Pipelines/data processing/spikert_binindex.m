% spikert_binindex
% This function generates a summary plot of spike rate versus the pursuit index,
% which is calculated as fidelity x vigor. The spike rate is averaged for each 
% pursuit index bin.
%
% INPUTS:
% spikert - spike rate data (vector or matrix)
% index   - pursuit index (vector or matrix)
%
% OUTPUTS:
% srvidx  - binned spike rate vs. pursuit index (matrix with two columns: index bins and average spike rate)
%
% ORIGINAL: 07/16/2024 - MC (created from binvelocity)
%
function srvidx = spikert_binindex(spikert,index)
%% set analysis parameters

% set bin size
i_bin = 0.05;

% reshape all datasets
r_spikert = reshape(spikert,[],1);
r_index = reshape(index,[],1);

%% discretize datasets

% discretize pursuit index
i_edge = -1-i_bin/2:i_bin:1+i_bin/2; %bin edges
i_bins = -1:i_bin:1; %bin labels (center)
i_disc=discretize(r_index,i_edge,i_bins);

%% for ALL behavior, calculate firing rate averages

% initialize
sr_mean=[];

% calculate mean and sem for each angular bin
for i = 1:length(i_bins)
    thisBin = i_bins(i);
    sr_mean(i,1) = mean(r_spikert(i_disc==thisBin),'omitnan');
end

srvidx = [i_bins' sr_mean];

end

