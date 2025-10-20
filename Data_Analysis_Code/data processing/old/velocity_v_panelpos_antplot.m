% velocity_v_panelpos_antplot
% analysis function for plotting when left and right turns cross the
% midline, in order to determine whether anticipation changes over time
%
% OUTPUTS:
%
% INPUTS:
% int_panelps - downsampled panel positions
% int_forward - downsampled forward velocities
% int_angular - downsampled angular velocities
% int_sideway - downsampled sideways velocities
% int_time - downsampled trial time
% chaseThreshold - forward velocity threshold, + for behavior >, - for behavior <
%
% ORIGINAL: 07/28/2022 - MC
%

function velocity_v_panelpos_antplot(int_panelps, int_forward,int_angular,int_sideway,int_time,chaseThreshold,optPlot)
%% find sweeps

% find right/left sweep start/stops based on min/max overlap
minPx = min(int_panelps(:,1));
maxPx = max(int_panelps(:,1));

rSweepIdx = find(int_panelps(:,1)==minPx);
lSweepIdx = find(int_panelps(:,1)==maxPx);
% plot for optional troubleshooting
%figure; plot(int_panelps(:,1),'k'); hold on; xline(rSweepIdx,'r'); xline(lSweepIdx,'b');hold off


%% optional, threshold behavior using forward velocity
if chaseThreshold>0 %pull ONLY behavior above threshold
    chase_Idx = schmittTrigger(int_forward,chaseThreshold,1);
    int_forward(~chase_Idx) = NaN;
    int_angular(~chase_Idx) = NaN;
    int_sideway(~chase_Idx) = NaN;
elseif chaseThreshold<0 %pull ONLY behavior above 0 and below threshold
    nonChase_Idx = int_forward>(-chaseThreshold) | int_forward<0;
    int_forward(nonChase_Idx) = NaN;
    int_angular(nonChase_Idx) = NaN;
    int_sideway(nonChase_Idx) = NaN;
end

%% pull directional velocity for each sweep
% initialize
nTrial = size(int_panelps,2);
nSweep = size(rSweepIdx,1);
int_angular_sweepR=[];
int_sideway_sweepR=[];
int_angular_sweepL=[];
int_sideway_sweepL=[];

% pull sweeps
for nt = 1:nTrial
    for ns = 1:nSweep
        % pull right sweeps
        int_angular_sweepR(:,end+1) = int_angular(rSweepIdx(ns):lSweepIdx(ns+1),nt);
        int_sideway_sweepR(:,end+1) = int_sideway(rSweepIdx(ns):lSweepIdx(ns+1),nt);
        % pull left sweeps
        int_angular_sweepL(:,end+1) = int_angular(lSweepIdx(ns):rSweepIdx(ns),nt);
        int_sideway_sweepL(:,end+1) = int_sideway(lSweepIdx(ns):rSweepIdx(ns),nt);
    end
end
% pull sweep px positions
pos_sweep(:,1) = int_panelps(rSweepIdx(ns):lSweepIdx(ns+1),nt);
pos_sweep(:,2) = int_panelps(lSweepIdx(ns):rSweepIdx(ns),nt);

%% calculate mean and adjust for sweep offset
% right sweep means
mean_vR(:,2) = mean(int_angular_sweepR,2,'omitnan');
mean_vR(:,3) = mean(int_sideway_sweepR,2,'omitnan');
% left sweep means
mean_vL(:,2) = mean(int_angular_sweepL,2,'omitnan');
mean_vL(:,3) = mean(int_sideway_sweepL,2,'omitnan');

% calculate individual fly bias
biasAng = mean([mean_vR(:,2) ;mean_vL(:,2)]);
biasSide = mean([mean_vR(:,3) ;mean_vL(:,3)]);

% adjust for individual fly bias
int_angular_sweepR = int_angular_sweepR - biasAng;
int_angular_sweepL = int_angular_sweepL - biasAng;
int_sideway_sweepR = int_sideway_sweepR - biasSide;
int_sideway_sweepL = int_sideway_sweepL - biasSide;

%% threshold
% remove any sweeps during which the fly was not actively chasing for the
% majority of the time

% find total NaN values for each sweep
nanTotR = sum(isnan(int_angular_sweepR),1);
nanTotL = sum(isnan(int_angular_sweepL),1);
% find each sweep with fewer NaNs than threshold
thresh = 20; %must have fewer than
chasesweep_idxR = find(nanTotR<thresh);
chasesweep_idxL = find(nanTotL<thresh);
% pull just those sweeps
int_angular_sweepRc = int_angular_sweepR(:,chasesweep_idxR);
int_angular_sweepLc = int_angular_sweepL(:,chasesweep_idxL);


%% find crossover points

% initialize
nsweepR = size(int_angular_sweepRc,2);
nsweepL = size(int_angular_sweepLc,2);
crossRight = NaN(1,nsweepR);
crossLeft = NaN(1,nsweepL);

for iR = 1:nsweepR
    crossIdx = find(int_angular_sweepRc(1:end-1,iR).*int_angular_sweepRc(2:end,iR)<0);
    if ~isempty(crossIdx)
        crossRight(iR) = pos_sweep(crossIdx(1),1); %pull first cross in panel pos
    end
end

for iL = 1:nsweepL
    crossIdx = find(int_angular_sweepLc(1:end-1,iL).*int_angular_sweepLc(2:end,iL)<0);
    if ~isempty(crossIdx)
        crossLeft(iL) = pos_sweep(crossIdx(1),2); %pull first cross in panel pos
    end
end


%% plot
% set plot variables
p1 = round(min(pos_sweep(:,1)));
p2 = round(max(pos_sweep(:,1)));
nSweeps = length(int_angular_sweepR);

if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1500 500])

    plot(chasesweep_idxR,crossRight,'o','Color', '#ff0080');
    hold on
    plot(chasesweep_idxL,crossLeft,'o','Color', '#0032A0');

    % adjust axes
    axis tight
    ylim([p1 p2])
    ylabel('crossover point (deg)')
    xlim([0 nSweeps])
    xlabel('sweep number')

    % quick linear fit
    for s = 1:2
        switch s
            case 1 %right fit
                x = chasesweep_idxR;
                xp = x;
                y = crossRight;
                linecolor = '#ff0080';
            case 2 %left fit
                x = chasesweep_idxL;
                xp = x;
                y = crossLeft;
                linecolor = '#0032A0';
        end

        % optional specify x range
        if 1
            xmin = 500;
            xmax = max(x);

            xidx = find(x>xmin & x<xmax);
            x=x(xidx);
            y=y(xidx);
        end

        % remove NaNs
        yidx = ~isnan(y);
        x=x(yidx);
        y=y(yidx);

        % generate linear fit
        coeff = polyfit(x,y,1);
        yfit = polyval(coeff,xp);
        plot(xp,yfit,'Color', linecolor,'LineWidth',2)
        hold on
    end
    hold off


end

end

