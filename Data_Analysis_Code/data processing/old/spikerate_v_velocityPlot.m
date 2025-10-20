% spikerate_v_velocityPlot
% analysis function for generating a summary plot of spike rate versus
% directional velocity as line plot
%
% OUTPUTS:
% velocityOutput - bins, mean, and SEM for each directional velocity
%
% INPUTS:
% allForward - downsampled forward velocities
% allAngular - downsampled angular velocities
% allSideway - downsampled sideways velocities
% ttime - downsampled trial time
% allSpikeRt - downsampled spikerate
% fwdThresh - for pursuit
% nTrial - number of trials or flies included
% optPlot - 0 no plot, 1 for plot
%
% ORIGINAL: 12/12/2022 - MC
% UPDATED:  02/07/2023 - MC added output variables
% UPDATED:  06/07/2023 - MC updated to focus on pursuit distributions
%

function summaryData = spikerate_v_velocityPlot(allForward,allAngular,allSideway,ttime,allSpikeRt,fwdThresh,nTrial,optPlot)
%% set analysis parameters
% set bin max
fwdMax = 10; %mm/s
angMax = 260; %deg/s
sidMax = 3.5; %mm/s

% set bin size
fs = 0.5;
as = 20;
ss = 0.25;

% reshape all dataset
rall_spikert = reshape(allSpikeRt,[],1);
rall_forward = reshape(allForward,[],1);
rall_angular = reshape(allAngular,[],1);
rall_sideway = reshape(allSideway,[],1);

%% isolate pursuit behavior

% pull idx for only vs not pursuit behavior
[~,~,~,purSpikeRt] = pursuitFinder(allForward,0,0,allSpikeRt,ttime,fwdThresh);
purIdx = ~isnan(purSpikeRt);
notIdx = isnan(purSpikeRt);

% pull data for only vs not pursuit behavior
rpur_spikert = rall_spikert;
rpur_spikert(notIdx) = nan; %blank anything NOT pursuit
rnot_spikert = rall_spikert;
rnot_spikert(purIdx) = nan; %blank anything pursuit

% determine pursuit threshold met
minThresh = 30; %sec
timespentchasing = (sum(allForward>fwdThresh,'all')/length(ttime))*60;
chaseLog = (timespentchasing>minThresh);


%% discretize behavioral datasets
% forward
f_edge = -fs/2:fs:fwdMax+fs/2; %bin edges
f_bins = 0:fs:fwdMax; %bin labels (center)
disc_all_forward=discretize(rall_forward,f_edge,f_bins);
% angular
a_edge = -angMax-as/2:as:angMax+as/2; %bin edges
a_bins = -angMax:as:angMax; %bin labels (center)
disc_all_angular=discretize(rall_angular,a_edge,a_bins);
% sideway
s_edge = -sidMax-ss/2:ss:sidMax+ss/2; %bin edges
s_bins = -sidMax:ss:sidMax; %bin labels (center)
disc_all_sideway=discretize(rall_sideway,s_edge,s_bins);

% if sufficient pursuit, bin data
if chaseLog
    % pull discretized data for pursuit ONLY
    disc_pur_forward = disc_all_forward(purIdx);
    disc_pur_angular = disc_all_angular(purIdx);
    disc_pur_sideway = disc_all_sideway(purIdx);

    % pull discretized data for NOT pursuit only
    disc_not_forward = disc_all_forward(notIdx);
    disc_not_angular = disc_all_angular(notIdx);
    disc_not_sideway = disc_all_sideway(notIdx);
end

% store distribution data
summaryData.fwdBin = f_bins;
summaryData.angBin = a_bins;
summaryData.sidBin = s_bins;


%% for ALL behavior, calculate firing rate averages

% initialize
fwdAllMean=[];
angAllMean=[];
sidAllMean=[];

% calculate mean and sem for each forward bin
for b = 1:length(f_bins)
    thisBin = f_bins(b);
    thisFwdIdx = find(disc_all_forward==thisBin);
    fwdAllMean(b) = mean(rall_spikert(thisFwdIdx),'omitnan');
    fwdAllSEM(b) = std(rall_spikert(thisFwdIdx),'omitnan')/sqrt(nTrial);
end
% calculate mean and sem for each angular bin
for b = 1:length(a_bins)
    thisBin = a_bins(b);
    thisAngIdx = find(disc_all_angular==thisBin);
    angAllMean(b) = mean(rall_spikert(thisAngIdx),'omitnan');
    angAllSEM(b) = std(rall_spikert(thisAngIdx),'omitnan')/sqrt(nTrial);
end
% calculate mean and sem for each sideways bin
for b = 1:length(s_bins)
    thisBin = s_bins(b);
    thisSidIdx = find(disc_all_sideway==thisBin);
    sidAllMean(b) = mean(rall_spikert(thisSidIdx),'omitnan');
    sidAllSEM(b) = std(rall_spikert(thisSidIdx),'omitnan')/sqrt(nTrial);
end

summaryData.fwdAllMean2 = fwdAllMean;
summaryData.angAllMean2 = angAllMean;
summaryData.sidAllMean2 = sidAllMean;

% check for nans
f_binsA = f_bins;
nanidx = find(isnan(fwdAllMean));
if nanidx
    f_binsA(nanidx) = [];
    fwdAllMean(nanidx) = [];
    fwdAllSEM(nanidx) = [];
end
a_binsA = a_bins;
nanidx = find(isnan(angAllMean));
if nanidx
    a_binsA(nanidx) = [];
    angAllMean(nanidx) = [];
    angAllSEM(nanidx) = [];
end
s_binsA = s_bins;
nanidx = find(isnan(sidAllMean));
if nanidx
    s_binsA(nanidx) = [];
    sidAllMean(nanidx) = [];
    sidAllSEM(nanidx) = [];
end

% store data of interest
summaryData.fwdAllBins = f_binsA;
summaryData.angAllBins = a_binsA;
summaryData.sidAllBins = s_binsA;
summaryData.fwdAllMean = fwdAllMean;
summaryData.angAllMean = angAllMean;
summaryData.sidAllMean = sidAllMean;
summaryData.fwdAllSEM = fwdAllSEM;
summaryData.angAllSEM = angAllSEM;
summaryData.sidAllSEM = sidAllSEM;

%% if sufficient pursuit, calculate respective firing rate averages for BINNED behavior
if chaseLog
    % initialize
    fwdPurMean=[];
    angPurMean=[];
    sidPurMean=[];
    fwdNotMean=[];
    angNotMean=[];
    sidNotMean=[];

    % calculate mean and sem for each forward bin
    for b = 1:length(f_bins)
        thisBin = f_bins(b);
        thisFwdIdx = find(disc_all_forward==thisBin);
        fwdPurMean(b) = mean(rpur_spikert(thisFwdIdx),'omitnan');
        fwdPurSEM(b) = std(rpur_spikert(thisFwdIdx),'omitnan')/sqrt(nTrial);
        fwdNotMean(b) = mean(rnot_spikert(thisFwdIdx),'omitnan');
        fwdNotSEM(b) = std(rnot_spikert(thisFwdIdx),'omitnan')/sqrt(nTrial);
    end
    % calculate mean and sem for each angular bin
    for b = 1:length(a_bins)
        thisBin = a_bins(b);
        thisAngIdx = find(disc_all_angular==thisBin);
        angPurMean(b) = mean(rpur_spikert(thisAngIdx),'omitnan');
        angPurSEM(b) = std(rpur_spikert(thisAngIdx),'omitnan')/sqrt(nTrial);
        angNotMean(b) = mean(rnot_spikert(thisAngIdx),'omitnan');
        angNotSEM(b) = std(rnot_spikert(thisAngIdx),'omitnan')/sqrt(nTrial);
    end
    % calculate mean and sem for each sideways bin
    for b = 1:length(s_bins)
        thisBin = s_bins(b);
        thisSidIdx = find(disc_all_sideway==thisBin);
        sidPurMean(b) = mean(rpur_spikert(thisSidIdx),'omitnan');
        sidPurSEM(b) = std(rpur_spikert(thisSidIdx),'omitnan')/sqrt(nTrial);
        sidNotMean(b) = mean(rnot_spikert(thisSidIdx),'omitnan');
        sidNotSEM(b) = std(rnot_spikert(thisSidIdx),'omitnan')/sqrt(nTrial);
    end

    % store data of interest
    summaryData.fwdMean = fwdPurMean;
    summaryData.angMean = angPurMean;
    summaryData.sidMean = sidPurMean;

    % check for nans
    f_binsP = f_bins;
    f_binsN = f_bins;
    nanidx = find(isnan(fwdPurMean));
    if nanidx
        f_binsP(nanidx) = [];
        fwdPurMean(nanidx) = [];
        fwdPurSEM(nanidx) = [];
    end
    nanidx = find(isnan(fwdNotMean));
    if nanidx
        f_binsN(nanidx) = [];
        fwdNotMean(nanidx) = [];
        fwdNotSEM(nanidx) = [];
    end
    a_binsP = a_bins;
    a_binsN = a_bins;
    nanidx = find(isnan(angPurMean));
    if nanidx
        a_binsP(nanidx) = [];
        angPurMean(nanidx) = [];
        angPurSEM(nanidx) = [];
    end
    nanidx = find(isnan(angNotMean));
    if nanidx
        a_binsN(nanidx) = [];
        angNotMean(nanidx) = [];
        angNotSEM(nanidx) = [];
    end
    s_binsP = s_bins;
    s_binsN = s_bins;
    nanidx = find(isnan(sidPurMean));
    if nanidx
        s_binsP(nanidx) = [];
        sidPurMean(nanidx) = [];
        sidPurSEM(nanidx) = [];
    end
    nanidx = find(isnan(sidNotMean));
    if nanidx
        s_binsN(nanidx) = [];
        sidNotMean(nanidx) = [];
        sidNotSEM(nanidx) = [];
    end
end


%% plot

% check if spikerate data or voltage data
if min(allSpikeRt(:,1,1))<0
    srmin = floor(min(fwdAllMean,[],'all'))-5;
    srmax = ceil(max(fwdAllMean,[],'all'))+5;
    ylbl = 'Voltage (mV)';
else
    srmin = 0;
    srmax = 30;
    ylbl = 'SR (spikes/sec)';
end

if chaseLog
    plotHeight = 800; %figure height
    nR = 3; %number of subplot rows
else
    plotHeight = 300; %figure height
    nR = 1; %number of subplot rows
end
if optPlot
    figure; set(gcf,'Position',[100 100 1000 plotHeight])
    plotlabels = {'Forward Velocity (mm/s)'; 'Angular Velocity (deg/s)'; 'Sideway Velocity (mm/s)'}; %velocity names
    colorlabels = {'#D95319';'#0072BD';'#7E2F8E';"#A2142F";[.7 .7 .7]}; %velocity colors
    a = 0.3; %face alpha
    lw = 2; %linewidth
end


% first, plot firing rate relationship for ALL data
if optPlot
    % plot forward results
    subplot(nR,3,1)
    % plot SEM band
    r(1) = patch([f_binsA'; flipud(f_binsA')],[(fwdAllMean-fwdAllSEM)'; flipud((fwdAllMean+fwdAllSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
    r(1).FaceColor = colorlabels{1};
    % plot average
    hold on
    plot(f_binsA,fwdAllMean,'LineWidth',lw,'Color',colorlabels{1})
    xlabel(plotlabels{1})
    ylabel([ylbl ' for ALL'])
    ylim([srmin srmax])
    xlim([0 fwdMax])
    hold off

    % plot angular
    subplot(nR,3,2)
    % plot SEM band
    r(2) = patch([a_binsA'; flipud(a_binsA')],[(angAllMean-angAllSEM)'; flipud((angAllMean+angAllSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
    r(2).FaceColor = colorlabels{2};
    % plot average
    hold on
    plot(a_binsA,angAllMean,'LineWidth',lw,'Color',colorlabels{2})
    xlabel(plotlabels{2})
    xlim([-angMax angMax])
    ylim([srmin srmax])

    % plot sideway
    subplot(nR,3,3)
    % plot SEM band
    r(3) = patch([s_binsA'; flipud(s_binsA')],[(sidAllMean-sidAllSEM)'; flipud((sidAllMean+sidAllSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
    r(3).FaceColor = colorlabels{3};
    % plot average
    hold on
    plot(s_binsA,sidAllMean,'LineWidth',lw,'Color',colorlabels{3})
    xlabel(plotlabels{3})
    xlim([-sidMax sidMax])
    ylim([srmin srmax])


% second, if pursuit, plot firing rate relationship for BINNED data
    if chaseLog
        % plot forward results
        subplot(nR,3,4)
        % plot not SEM band
        r(1) = patch([f_binsN'; flipud(f_binsN')],[(fwdNotMean-fwdNotSEM)'; flipud((fwdNotMean+fwdNotSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
        r(1).FaceColor = colorlabels{5};
        % plot not average
        hold on
        plot(f_binsN,fwdNotMean,'LineWidth',lw,'Color',colorlabels{5})
        % plot pursuit SEM band
        r(1) = patch([f_binsP'; flipud(f_binsP')],[(fwdPurMean-fwdPurSEM)'; flipud((fwdPurMean+fwdPurSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
        r(1).FaceColor = colorlabels{4};
        % plot pursuit average
        plot(f_binsP,fwdPurMean,'LineWidth',lw,'Color',colorlabels{4})
        xlabel(plotlabels{1})
        ylabel([ylbl ' BINNED'])
        ylim([srmin srmax])
        xlim([0 fwdMax])
        hold off

        % plot angular
        subplot(nR,3,5)
        % plot not SEM band
        r(2) = patch([a_binsN'; flipud(a_binsN')],[(angNotMean-angNotSEM)'; flipud((angNotMean+angNotSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
        r(2).FaceColor = colorlabels{5};
        % plot not average
        hold on
        plot(a_binsN,angNotMean,'LineWidth',lw,'Color',colorlabels{5})
        % plot pursuit SEM band
        r(2) = patch([a_binsP'; flipud(a_binsP')],[(angPurMean-angPurSEM)'; flipud((angPurMean+angPurSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
        r(2).FaceColor = colorlabels{4};
        % plot pursuit average
        plot(a_binsP,angPurMean,'LineWidth',lw,'Color',colorlabels{4})
        xlabel(plotlabels{2})
        xlim([-angMax angMax])
        ylim([srmin srmax])

        % plot sideway
        subplot(nR,3,6)
        % plot not SEM band
        r(3) = patch([s_binsN'; flipud(s_binsN')],[(sidNotMean-sidNotSEM)'; flipud((sidNotMean+sidNotSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
        r(3).FaceColor = colorlabels{5};
        % plot not average
        hold on
        plot(s_binsN,sidNotMean,'LineWidth',lw,'Color',colorlabels{5})
        % plot pursuit SEM band
        r(3) = patch([s_binsP'; flipud(s_binsP')],[(sidPurMean-sidPurSEM)'; flipud((sidPurMean+sidPurSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
        r(3).FaceColor = colorlabels{4};
        % plot pursuit average
        plot(s_binsP,sidPurMean,'LineWidth',lw,'Color',colorlabels{4})
        xlabel(plotlabels{3})
        xlim([-sidMax sidMax])
        ylim([srmin srmax])


% last, if pursuit, plot distribution of pursuit behavior
        % forward velocity histogram
        subplot(nR,3,7)
        hf = histogram(disc_pur_forward,f_edge,'FaceColor',colorlabels{4},'FaceAlpha',a);
        xlim([0 fwdMax])
        xlabel(plotlabels{1})
        ylabel('frequency')

        % angular velocity histogram
        subplot(nR,3,8)
        a_edgeH = -angMax-as/2:as:angMax+as/2; %bin edges
        ha = histogram(disc_pur_angular,a_edgeH,'FaceColor',colorlabels{4},'FaceAlpha',a);
        xlim([-angMax angMax])
        xlabel(plotlabels{2})
        ylabel('frequency')
        xline(0)

        % sideway velocity histogram
        subplot(nR,3,9)
        s_edgeH = -sidMax-ss/2:ss:sidMax+ss/2; %bin edges
        hs = histogram(disc_pur_sideway,s_edgeH,'FaceColor',colorlabels{4},'FaceAlpha',a);
        xlim([-sidMax sidMax])
        xlabel(plotlabels{3})
        ylabel('frequency')
        xline(0)

        % store distribution data
        summaryData.fwdCounts = hf.Values;
        summaryData.angCounts = ha.Values;
        summaryData.sidCounts = hs.Values;
    end
end

end

