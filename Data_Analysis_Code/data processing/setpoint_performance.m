% SETPOINT_PERFORMANCE - This function assesses the steering performance of flies during pursuit
% behavior by analyzing various metrics in relation to a fixed setpoint (assumed to be 0 or fixation).
% The analysis excludes timepoints immediately following a bar jump (saccade) to focus on continuous 
% pursuit behavior. Several metrics, such as integral absolute error, integral squared error, 
% circular standard deviation, and setpoint proximity, are computed to evaluate the fly's performance.
%
% INPUTS:
%   panelps  - 3D array containing the panel position (yaw gain-modified), where the 3rd dimension 
%              represents different experimental conditions, and the 2nd dimension represents trials.
%   barjump  - 3D array indicating bar jump triggers, used to identify saccadic events.
%   forward  - 3D array representing the forward velocity of the fly during trials.
%   ttime    - Time vector representing the temporal resolution of the data.
%   optPlot  - Optional flag (1/0) indicating whether to plot the results (1 for yes, 0 for no).
%
% OUTPUTS:
%   sp_out   - Structure containing various performance metrics:
%              - sp_out.prob  : Probability of being near the setpoint (fixation).
%              - sp_out.runt  : Median run time for each condition.
%              - sp_out.iae   : Integral of Absolute Error (IAE), representing total deviation.
%              - sp_out.ise   : Integral of Squared Error (ISE), penalizing large deviations.
%              - sp_out.cstd  : Circular standard deviation of the panel position.
%              - sp_out.cvar  : Circular variance of the panel position.
%
% Created: 08/21/24 by MC
%
function [sp_out] = setpoint_performance(panelps,barjump,forward,ttime,optPlot)
%% initialize
% dataset info
nCond = size(panelps,3);
nTrials = size(panelps,2);

% load processing settings
settings = processSettings();

% omit timepoints right after a bar jump
tomit = 2; % s
iomit = fetchTimeIdx(ttime,tomit); %idx

%% analyze setpoint
% initialize
probSP = [];
runT = [];
iae = [];
ise = [];
cSTD = [];
cVar = [];

% for each condition
for c = 1:nCond
    % estimate HD bias
    biasHD = mean(panelps(:,:,c),'all','omitnan');
    % for each trial
    for t = 1:nTrials
        % fetch data
        thisPanelps = panelps(:,t,c);
        thisBarJump = barjump(:,t,c);
        thisForward = forward(:,t,c);
        thisTime = ttime;
        
        % fictrac occassionally breaks, ensure this trial was run properly else omit
        if sum(~isnan(thisPanelps))
            % post-process HD data
            % remove bar jumps
            idxJumps = find(diff(thisBarJump)>0);
            for j = 1:length(idxJumps)
                thisPanelps(idxJumps(j):idxJumps(j)+iomit,1) = nan;
            end
            % omit nans from calculations
            thisTime(isnan(thisPanelps)) = [];
            thisPanelps(isnan(thisPanelps)) = [];
            
            % adjust for bias
            thisPanelps = thisPanelps-biasHD;

            % calculate time spent running
            runT(c,t) = ttime(sum(~isnan(thisForward))+1);

            % calculate probability of being near setpoint (+/-)
            [h,e] = histcounts(thisPanelps,'BinWidth',settings.spBin,'Normalization','probability');
            e = e(1:end-1)+settings.spBin/2; %center bins
            probSP(c,t) = sum(h(abs(e)<=settings.spHD));

            % calculate IAE (integral of the absolute error)
            % less aggressive, treats small and large errors equally
            abs_error = abs(thisPanelps-0);
            iae(c,t) = trapz(thisTime,abs_error);
            % calculate ISE (integral of the squared error)
            % more aggressive, heavily penalizes large errors
            sqd_error = (thisPanelps-0).^2;
            ise(c,t) = trapz(thisTime,sqd_error);

            % calculate circular standard deviation
            cSTD(c,t) = circ_std(deg2rad(thisPanelps),[],[],1);
            % calculate circular variance
            cVar(c,t) = 1-circ_var(deg2rad(thisPanelps),[],[],1);
        else %omit
            probSP(c,t) = nan;
            runT(c,t) = nan;
            iae(c,t) = nan;
            ise(c,t) = nan;
            cSTD(c,t) = nan;
            cVar(c,t) = nan;
        end
    end
end

sp_out.prob = median(probSP,2,'omitnan');
sp_out.runt = median(runT,2,'omitnan');
sp_out.iae = median(iae,2,'omitnan');
sp_out.ise = median(ise,2,'omitnan');
sp_out.cstd = median(cSTD,2,'omitnan');
sp_out.cvar = median(cVar,2,'omitnan');

%% (optional) plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 600 800])
    tiledlayout(2,3,'TileSpacing','compact')
    gainOpt = settings.pursuitGain;

    % plot trials and median
    nexttile; hold on
    plot(gainOpt,runT,'.','Color',settings.trialColor,'LineWidth',settings.lwTri)
    plot(gainOpt,sp_out.runt,'o-','Color',settings.HDColor,'LineWidth',settings.lwAvg)
    axis padded; xlabel('k'); ylabel('run time (s)'); xticks(gainOpt); ylim([0 115])

    nexttile; hold on
    plot(gainOpt,iae,'.','Color',settings.trialColor,'LineWidth',settings.lwTri)
    plot(gainOpt,sp_out.iae,'o-','Color',settings.HDColor,'LineWidth',settings.lwAvg)
    axis padded; xlabel('k'); ylabel('Integral of Absolute Error'); xticks(gainOpt);

    nexttile; hold on
    plot(gainOpt,ise,'.','Color',settings.trialColor,'LineWidth',settings.lwTri)
    plot(gainOpt,sp_out.ise,'o-','Color',settings.HDColor,'LineWidth',settings.lwAvg)
    axis padded; xlabel('k'); ylabel('Integral of Squared Error'); xticks(gainOpt);

    nexttile; hold on
    plot(gainOpt,probSP,'.','Color',settings.trialColor,'LineWidth',settings.lwTri)
    plot(gainOpt,sp_out.prob,'o-','Color',settings.HDColor,'LineWidth',settings.lwAvg)
    axis padded; xlabel('k'); ylabel('probability(setpoint +/-10)'); xticks(gainOpt); ylim([0 1])

    nexttile; hold on
    plot(gainOpt,cSTD,'.','Color',settings.trialColor,'LineWidth',settings.lwTri)
    plot(gainOpt,sp_out.cstd,'o-','Color',settings.HDColor,'LineWidth',settings.lwAvg)
    axis padded; xlabel('k'); ylabel('circularSTD'); xticks(gainOpt); ylim([0 1])

    nexttile; hold on
    plot(gainOpt,cVar,'.','Color',settings.trialColor,'LineWidth',settings.lwTri)
    plot(gainOpt,sp_out.cvar,'o-','Color',settings.HDColor,'LineWidth',settings.lwAvg)
    axis padded; xlabel('k'); ylabel('circularVariance'); xticks(gainOpt); ylim([0 1])

end
end