% spikerate_v_velocityPlot
% analysis function for generating a summary plot of spike rate versus
% directional velocity as line plot
%
% OUTPUTS:
% velocityOutput - bins, mean, and SEM for each directional velocity
%
% INPUTS:
% int_forward - downsampled forward velocities
% int_angular - downsampled angular velocities
% int_sideway - downsampled sideways velocities
% int_time - downsampled trial time
% int_spikerate - downsampled spikerate
% fwdThresh - for pursuit
% nTrial - number of trials or flies included
% optPlot - 0 no plot, 1 for plot
%
% ORIGINAL: 12/12/2022 - MC
% UPDATED:  02/07/2023 - MC added output variables
%

function velocityOutput = spikerate_v_velocityPlot(int_forward,int_angular,int_sideway,int_time,int_spikerate,fwdThresh,nTrial,optPlot)
%% set plot parameters
% avoid plotting outliers or datapoints with insufficient samples
fwdMax = 8; %mm/s
angMax = 200; %deg/s
sidMax = 3.5; %mm/s

srmax = 50;


%% optional, threshold behavior using forward velocity
%analyze pursuit and all behavior separately
% for all behavior, simply reshape
rsp_spikerate = reshape(int_spikerate,[],1);
% for pursuit behavior, threshold then reshape
[~,~,~,int_spikerate_high] = pursuitFinder(int_forward,0,0,int_spikerate,int_time,fwdThresh);
rsp_spikerate_high = reshape(int_spikerate_high,[],1);

% determine if experiment met minimum requirements to be considered "pursuit"
% min threshold set as time spent running above 5mm/s (e.g., likely pursuit)
minThresh = 30; %sec
timespentchasing = (sum(int_forward>fwdThresh,'all')/length(int_time))*60;
chaseLog = (timespentchasing>minThresh);


%% reshape and discretize the data set
% reshape each velocity to single column
rsp_forward = reshape(int_forward,[],1);
rsp_angular = reshape(int_angular,[],1);
rsp_sideway = reshape(int_sideway,[],1);

% discretize each variable
fs = 1; %bin size
%f_max = round(max(rsp_forward),0); %max
f_max = fwdMax;
f_edge = -fs/2:fs:f_max+fs/2; %bin edges
f_bins = 0:fs:f_max; %bin labels (center)
disc_forward=discretize(rsp_forward,f_edge,f_bins);

as = 20; %bin size
%a_max = round(max([rsp_angular -rsp_angular],[],'all'),-1); %max
a_max = angMax;
a_edge = -a_max-as/2:as:a_max+as/2; %bin edges
a_bins = -a_max:as:a_max; %bin labels (center)
disc_angular=discretize(rsp_angular,a_edge,a_bins);

ss = 0.5; %bin size
%s_max = round(max([rsp_sideway -rsp_sideway],[],'all'),0); %max
s_max = sidMax;
s_edge = -s_max-ss/2:ss:s_max+ss/2; %bin edges
s_bins = -s_max:ss:s_max; %bin labels (center)
disc_sideway=discretize(rsp_sideway,s_edge,s_bins);


%% initialize plot
% initialize
if chaseLog
    plotHeight = 800; %figure height
    nR = 3; %number of subplot rows
else
    plotHeight = 500; %figure height
    nR = 2; %number of subplot rows
end
if optPlot
    figure; set(gcf,'Position',[100 100 1000 plotHeight])
    plotlabels = {'fwd (mm/s)'; 'ang (deg/s)'; 'side (mm/s)'}; %velocity names
    colorlabels = {'#D95319';'#0072BD';'#7E2F8E'}; %velocity colors
    a = 0.5; %face alpha
    lw = 3; %linewidth
end

%% plot distribution for each behavior

if optPlot
    % forward velocity histogram
    subplot(nR,3,1)
    hf = histogram(disc_forward,f_edge,'FaceColor',colorlabels{1},'FaceAlpha',a);
    xlim([0 fwdMax])
    ylim([0 hf.Values(2)])
    xlabel(plotlabels{1})
    ylabel('frequency')

    % angular velocity histogram
    subplot(nR,3,2)
    a_edgeH = -a_max-as/2:as:a_max+as/2; %bin edges
    ha = histogram(disc_angular,a_edgeH,'FaceColor',colorlabels{2},'FaceAlpha',a);
    xlim([-angMax angMax])
    xlabel(plotlabels{2})
    ylabel('frequency')

    % sideway velocity histogram
    subplot(nR,3,3)
    s_edgeH = -s_max-ss/2:ss:s_max+ss/2; %bin edges
    hs = histogram(disc_sideway,s_edgeH,'FaceColor',colorlabels{3},'FaceAlpha',a);
    xlim([-sidMax sidMax])
    xlabel(plotlabels{3})
    ylabel('frequency')
end

%% plot for all behavior
fwdMean=[];
angMean=[];
sidMean=[];

% calculate mean and stdev for each forward bin
for b = 1:length(f_bins)
    thisBin = f_bins(b);
    thisFwdIdx = find(disc_forward==thisBin);
    fwdMean(b) = mean(rsp_spikerate(thisFwdIdx),'omitnan');
    fwdSEM(b) = std(rsp_spikerate(thisFwdIdx),'omitnan')/sqrt(nTrial);
end
% calculate mean and stdev for each angular bin
for b = 1:length(a_bins)
    thisBin = a_bins(b);
    thisAngIdx = find(disc_angular==thisBin);
    angMean(b) = mean(rsp_spikerate(thisAngIdx),'omitnan');
    angSEM(b) = std(rsp_spikerate(thisAngIdx),'omitnan')/sqrt(nTrial);
end
% calculate mean and stdev for each sideways bin
for b = 1:length(s_bins)
    thisBin = s_bins(b);
    thisSidIdx = find(disc_sideway==thisBin);
    sidMean(b) = mean(rsp_spikerate(thisSidIdx),'omitnan');
    sidSEM(b) = std(rsp_spikerate(thisSidIdx),'omitnan')/sqrt(nTrial);
end

% store
velocityOutput.forward = [f_bins; fwdMean; fwdSEM; hf.Values];
velocityOutput.sideway = [s_bins; sidMean; sidSEM; hs.Values];
velocityOutput.angular = [a_bins; angMean; angSEM; ha.Values];

% check for nans
nanidx = find(isnan(fwdMean));
if nanidx
    f_bins(nanidx) = [];
    fwdMean(nanidx) = [];
    fwdSEM(nanidx) = [];
end

if optPlot
    if ~chaseLog
        % plot forward results
        subplot(nR,3,4)
        % plot SEM band
        r(1) = patch([f_bins'; flipud(f_bins')],[(fwdMean-fwdSEM)'; flipud((fwdMean+fwdSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
        r(1).FaceColor = colorlabels{1};
        % plot average
        hold on
        plot(f_bins,fwdMean,'LineWidth',lw,'Color',colorlabels{1})
        xlabel(plotlabels{1})
        ylabel('spike rate')
        ylim([0 srmax])
        xlim([0 fwdMax])
        hold off
    end

    % plot angular results
    subplot(nR,3,5)
    % plot SEM band
    angMean(isnan(angMean))=0; %remove nans
    angSEM(isnan(angSEM))=0; %remove nans
    r(2) = patch([a_bins'; flipud(a_bins')],[(angMean-angSEM)'; flipud((angMean+angSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
    r(2).FaceColor = colorlabels{2};
    % plot average
    hold on
    plot(a_bins,angMean,'LineWidth',lw,'Color',colorlabels{2})
    xlabel(plotlabels{2})
    ylabel('spike rate')
    xlim([-angMax angMax])
    ylim([0 srmax])

    % plot sideway results
    subplot(nR,3,6)
    % plot SEM band
    sidMean(isnan(sidMean))=0; %remove nans
    sidSEM(isnan(sidSEM))=0; %remove nans
    r(3) = patch([s_bins'; flipud(s_bins')],[(sidMean-sidSEM)'; flipud((sidMean+sidSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
    r(3).FaceColor = colorlabels{3};
    % plot average
    hold on
    plot(s_bins,sidMean,'LineWidth',lw,'Color',colorlabels{3})
    xlabel(plotlabels{3})
    ylabel('spike rate')
    xlim([-sidMax sidMax])
    ylim([0 srmax])
end


%% plot for pursuit behavior, if sufficient behavior

if chaseLog
    angMean=[];
    sidMean=[];

    % calculate mean and stdev for each angular bin
    for b = 1:length(a_bins)
        thisBin = a_bins(b);
        thisAngIdx = find(disc_angular==thisBin);
        angMean(b) = mean(rsp_spikerate_high(thisAngIdx),'omitnan');
        angSEM(b) = std(rsp_spikerate_high(thisAngIdx),'omitnan')/sqrt(nTrial);
    end
    % calculate mean and stdev for each sideways bin
    for b = 1:length(s_bins)
        thisBin = s_bins(b);
        thisSidIdx = find(disc_sideway==thisBin);
        sidMean(b) = mean(rsp_spikerate_high(thisSidIdx),'omitnan');
        sidSEM(b) = std(rsp_spikerate_high(thisSidIdx),'omitnan')/sqrt(nTrial);
    end

    if optPlot
        % plot forward results
        subplot(nR,3,7)
        % plot SEM band
        r(1) = patch([f_bins'; flipud(f_bins')],[(fwdMean-fwdSEM)'; flipud((fwdMean+fwdSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
        r(1).FaceColor = colorlabels{1};
        % plot average
        hold on
        plot(f_bins,fwdMean,'LineWidth',lw,'Color',colorlabels{1})
        xlabel(plotlabels{1})
        ylabel('spike rate')
        ylim([0 srmax])
        xlim([0 fwdMax])
        hold off

        % plot angular
        subplot(nR,3,8)
        % plot SEM band
        angMean(isnan(angMean))=0; %remove nans
        angSEM(isnan(angSEM))=0; %remove nans
        r(2) = patch([a_bins'; flipud(a_bins')],[(angMean-angSEM)'; flipud((angMean+angSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
        r(2).FaceColor = colorlabels{2};
        % plot average
        hold on
        plot(a_bins,angMean,'LineWidth',lw,'Color',colorlabels{2})
        ylabel(['spike rate (fwd>' num2str(fwdThresh) 'mm/s)'])
        xlabel(plotlabels{2})
        xlim([-angMax angMax])
        ylim([0 srmax])

        % plot sideway
        subplot(nR,3,9)
        % plot SEM band
        sidMean(isnan(sidMean))=0; %remove nans
        sidSEM(isnan(sidSEM))=0; %remove nans
        r(3) = patch([s_bins'; flipud(s_bins')],[(sidMean-sidSEM)'; flipud((sidMean+sidSEM)')], 'r', 'FaceAlpha',a, 'EdgeColor','none');
        r(3).FaceColor = colorlabels{3};
        % plot average
        hold on
        plot(s_bins,sidMean,'LineWidth',lw,'Color',colorlabels{3})
        ylabel(['spike rate (fwd>' num2str(fwdThresh) 'mm/s)'])
        xlabel(plotlabels{3})
        xlim([-sidMax sidMax])
        ylim([0 srmax])
    end
end


end

