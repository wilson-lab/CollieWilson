% spikert_heatmapvelocity
% This function generates a summary plot of spike rate versus directional velocities
% (forward, angular, and sideways) as heatmaps. The function also bins the data for analysis,
% with optional plotting of the heatmaps.
%
% INPUTS:
% forward      - downsampled forward velocities
% angular      - downsampled angular velocities
% sideway      - downsampled sideways velocities
% spikert      - downsampled spike rates
% ttime        - time vector (in seconds)
% lagOpt       - 0 or 1 to apply lag shifts to velocity data
% optPlot      - 0 to skip plotting, 1 to generate heatmap plots
%
% OUTPUT:
% binnedData   - binned spike rate and directional velocity data (table)
%
% ORIGINAL: 06/24/2022 - MC
% UPDATED:  03/14/2024 - MC added optional plotting, output variables
%
function binnedData = spikert_heatmapvelocity(forward,angular,sideway,spikert,ttime,lagOpt,optPlot)
%% settings
% set bin size
fs = 0.5;
as = 20;
ss = 0.25;

% fetch settings
settings = processSettings();
% fetch trial duration
[trialDur,~] = size(forward);

%% (optional) exclude start/stop transitions
ex_startstop = 1;
if ex_startstop
    % set transition window to exclude
    [tidx] = fetchTimeIdx(ttime,settings.ssExclude);
    tidxh = floor(tidx/2); %half before/after each start/stop transition

    % find all transitions
    forward_startstop = zeros(size(forward));
    forward_startstop(forward>settings.ssThresh) = 1; %0 not moving, 1 moving
    forward_trans = find(abs(diff(forward_startstop))); %find all start/stop transitions 0 -> 1

    % for each start/stop transition
    for ft = 1:length(forward_trans)
        t1 = forward_trans(ft)-tidxh;
        t2 = forward_trans(ft)+tidxh;
        % ensure start/stop does not exceed trial start/stop
        if t1<1
            t1 = 1;
        elseif t2>trialDur
            t2 = trialDur;
        end

        % remove start/stop transition
        spikert(t1:t2) = nan;
    end

    % test plot
    %figure(1); clf(1);plot(cellact_ss(:,1)); hold on;plot(forward_startstop(:,1)*30);plot(forward(:,1))
end

%% (optional) shift according to lag estimates
% if lag estimates were provided, shift
if lagOpt
    % fetch shift indices for each lag
    [idxf] = fetchTimeIdx(ttime,settings.fwdLag);
    [idxa] = fetchTimeIdx(ttime,settings.angLag);
    [idxs] = fetchTimeIdx(ttime,settings.sidLag);
    idxf = idxf-1; idxa = idxa-1; idxs = idxs-1;
    
    % shift and exclude data at start/stop of trial
    forward = circshift(forward,-idxf,1); %shift
    forward(end-idxf+1:end,:) = nan; %exclude shifts
    angular = circshift(angular,-idxa,1); %shift
    angular(end-idxa+1:end,:) = nan; %exclude shifts
    sideway = circshift(sideway,-idxs,1); %shift
    sideway(end-idxs+1:end,:) = nan; %exclude shifts
end

%% reshape and discretize the data set
% reshape each to single column
rsp_forward = reshape(forward,[],1);
rsp_angular = reshape(angular,[],1);
rsp_sideway = reshape(sideway,[],1);
rsp_spikerate = reshape(spikert,[],1);

% discretize each variable
f_max = 10;
f_edge = -fs/2:fs:f_max+fs/2; %bin edges
f_bins = 0:fs:f_max; %bin labels (center)
disc_forward=discretize(rsp_forward,f_edge,f_bins);

a_max = 200;
a_edge = -a_max-as/2:as:a_max+as/2; %bin edges
a_bins = -a_max:as:a_max; %bin labels (center)
disc_angular=discretize(rsp_angular,a_edge,a_bins);

s_max = 2.5;
s_edge = -s_max-ss/2:ss:s_max+ss/2; %bin edges
s_bins = -s_max:ss:s_max; %bin labels (center)
disc_sideway=discretize(rsp_sideway,s_edge,s_bins);

rs = 1; %bin size
r_max = round(max(rsp_spikerate),-1); %max
r_edge = -rs/2:rs:r_max+rs/2; %bin edges
r_bins = 0:rs:r_max; %bin labels (center)
disc_spikerate=discretize(rsp_spikerate,r_edge,r_bins);

% compile data and remove NaN values
trialData_pre = [disc_spikerate disc_forward disc_angular disc_sideway];
[rNan, ~] = find(isnan(trialData_pre));
trialData_pre(rNan,:)=[];

% convert to table
colNames = {'SpikeRate','Forward','Angular','Sideway'};
binnedData = array2table(trialData_pre,'VariableNames',colNames);

%% plot

if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1000 500])
    tiledlayout(1,2,'TileSpacing','compact')
    
    % plot forward and angular
    nexttile
    h(1) = heatmap(binnedData,'Angular','Forward','ColorVariable','SpikeRate');
    h(1).YDisplayData = flipud(h(1).YDisplayData);
    h(1).CellLabelColor = 'None';
    h(1).GridVisible = 'off';
    
    % plot forward and sideways
    nexttile
    h(2) = heatmap(binnedData,'Sideway','Forward','ColorVariable','SpikeRate');
    h(2).YDisplayData = flipud(h(2).YDisplayData);
    h(2).CellLabelColor = 'None';
    h(2).GridVisible = 'off';
end

end

