% spikert_binacceleration
% This function generates a summary plot of spike rate binned according to 
% directional acceleration (forward, angular, and sideways). It calculates the 
% acceleration from the velocity data and bins the spike rate for each acceleration bin.
% The function also provides optional exclusion of start/stop transitions and lag shifts.
%
% INPUTS:
% forward   - forward velocities
% angular   - angular velocities
% sideway   - sideways velocities
% spikert   - spike rate or vm data
% ttime     - trial time (in seconds)
% lagOpt    - 1 to apply lag shifts, 0 to omit
%
% OUTPUTS:
% summaryData - structure containing binned means for each directional acceleration
%               (fwdBin, angBin, sidBin, fwdMean, angMean, sidMean)
%
% ORIGINAL: 08/07/2024 - MC created from binvelocities
%
function summaryData = spikert_binacceleration(forward,angular,sideway,spikert,ttime,lagOpt)
%% set analysis parameters
% fetch settings
settings = processSettings();
% fetch trial duration
[trialDur,~] = size(forward);

% set bin max
fwdMax = 200; %mm/s2
angMax = 1000; %deg/s2
sidMax = 200; %mm/s2
% set bin size
fs = 20;
as = 50;
ss = 20;

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

%% calculate acceleration
[accel] = velocity2acceleration(forward,angular,sideway,ttime);

%% reshape all dataset
spikert_r = reshape(spikert,[],1);
forward_r = reshape(accel.forward,[],1);
angular_r = reshape(accel.angular,[],1);
sideway_r = reshape(accel.sideway,[],1);

%% discretize behavioral datasets
% forward
f_edge = -fwdMax-fs/2:fs:fwdMax+fs/2; %bin edges
f_bins = -fwdMax:fs:fwdMax; %bin labels (center)
forward_disc=discretize(forward_r,f_edge,f_bins);
% angular
a_edge = -angMax-as/2:as:angMax+as/2; %bin edges
a_bins = -angMax:as:angMax; %bin labels (center)
angular_disc=discretize(angular_r,a_edge,a_bins);
% sideway
s_edge = -sidMax-ss/2:ss:sidMax+ss/2; %bin edges
s_bins = -sidMax:ss:sidMax; %bin labels (center)
sideway_disc=discretize(sideway_r,s_edge,s_bins);

%% for ALL behavior, calculate cell activity averages

% initialize
fwdAllMean=[];
angAllMean=[];
sidAllMean=[];

% calculate mean and sem for each forward bin
for b = 1:length(f_bins)
    thisBin = f_bins(b);
    thisFwdIdx = find(forward_disc==thisBin);
    fwdAllMean(b) = mean(spikert_r(thisFwdIdx),'omitnan');
end
% calculate mean and sem for each angular bin
for b = 1:length(a_bins)
    thisBin = a_bins(b);
    thisAngIdx = find(angular_disc==thisBin);
    angAllMean(b) = mean(spikert_r(thisAngIdx),'omitnan');
end
% calculate mean and sem for each sideways bin
for b = 1:length(s_bins)
    thisBin = s_bins(b);
    thisSidIdx = find(sideway_disc==thisBin);
    sidAllMean(b) = mean(spikert_r(thisSidIdx),'omitnan');
end

% store distribution data
summaryData.fwdBin = f_bins;
summaryData.angBin = a_bins;
summaryData.sidBin = s_bins;
summaryData.fwdMean = fwdAllMean;
summaryData.angMean = angAllMean;
summaryData.sidMean = sidAllMean;

end

