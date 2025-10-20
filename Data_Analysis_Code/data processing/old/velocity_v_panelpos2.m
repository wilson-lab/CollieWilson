% velocity_v_panelpos2
% analysis function for generating a summary plot of panel position versus
% directional velocity. plots left and right sweeps across same position
% separately.
%
% OUTPUTS:
% mean_vR - mean velocities for rightward sweep
% mean_vL - mean velocities for leftward sweep
%
% INPUTS:
% int_panelps - downsampled panel positions
% int_forward - downsampled forward velocities
% int_angular - downsampled angular velocities
% int_sideway - downsampled sideways velocities
% int_time - downsampled trial time
% chaseThreshold - forward velocity threshold, + for behavior >, - for behavior <
%
% ORIGINAL: 06/21/2022 - MC
% UPDATED:  12/13/2022 - MC improved pursuit tracker
%

function [mean_vR,mean_vL,pos_sweep] = velocity_v_panelpos2(int_panelps, int_forward,int_angular,int_sideway,int_time,chaseThreshold,optPlot)
%% optional, threshold behavior using forward velocity
if chaseThreshold>0 %pull ONLY behavior above threshold
    [int_forward,int_angular,int_sideway,~] = pursuitFinder(int_forward,int_angular,int_sideway,0,int_time,chaseThreshold);
elseif chaseThreshold<0 %pull ONLY behavior above 0 and below threshold
    nonChase_Idx = int_forward>(-chaseThreshold) | int_forward<0;
    int_forward(nonChase_Idx) = NaN;
    int_angular(nonChase_Idx) = NaN;
    int_sideway(nonChase_Idx) = NaN;
end

%% pull directional velocity for each sweep

% remove panel data noise
int_panelps_r = round((int_panelps*2))/2; %round to nearest 0.5

% for each trial
nTrial = size(int_panelps,2);
for nt = 1:nTrial
    % find start/stop of sweeps by finding max left/right pos
    [~,peakIdx_left] = findpeaks(-int_panelps_r(:,nt),'MinPeakHeight',35); %left peaks
    [~,peakIdx_right] = findpeaks(int_panelps_r(:,nt),'MinPeakHeight',35); %right peaks
    if peakIdx_left(1)<10 %occassional error at start, remove it if needed
        peakIdx_left(1)=[];
    end
    nSweeps = size(peakIdx_left,1);

    % if first trial, initialize
    if nt==1
        sweepDur = round(mean((peakIdx_left(1:end)-peakIdx_right(1:end-1))))+1;
        int_forward_rightward=NaN(sweepDur,nTrial*nSweeps);
        int_angular_rightward=NaN(sweepDur,nTrial*nSweeps);
        int_sideway_rightward=NaN(sweepDur,nTrial*nSweeps);
        int_forward_leftward=NaN(sweepDur,nTrial*nSweeps);
        int_angular_leftward=NaN(sweepDur,nTrial*nSweeps);
        int_sideway_leftward=NaN(sweepDur,nTrial*nSweeps);
        is=0; %sweep index
    end

    % pull sweeps
    for ns = 1:nSweeps
        dIdx = ns + is; %data index
        % pull rightward sweeps by taking sweep following left-most peak
        int_forward_rightward(1:sweepDur,dIdx) = int_forward(peakIdx_left(ns):peakIdx_left(ns)+sweepDur-1,nt);
        int_angular_rightward(1:sweepDur,dIdx) = int_angular(peakIdx_left(ns):peakIdx_left(ns)+sweepDur-1,nt);
        int_sideway_rightward(1:sweepDur,dIdx) = int_sideway(peakIdx_left(ns):peakIdx_left(ns)+sweepDur-1,nt);
        % pull leftward sweeps by taking sweep following right-most peak
        int_forward_leftward(1:sweepDur,dIdx) = int_forward(peakIdx_right(ns):peakIdx_right(ns)+sweepDur-1,nt);
        int_angular_leftward(1:sweepDur,dIdx) = int_angular(peakIdx_right(ns):peakIdx_right(ns)+sweepDur-1,nt);
        int_sideway_leftward(1:sweepDur,dIdx) = int_sideway(peakIdx_right(ns):peakIdx_right(ns)+sweepDur-1,nt);
    end
    is = is+nSweeps; %update sweep index
end

% pull sweep px positions
pos_sweep(:,1) = int_panelps_r(peakIdx_left(ns):peakIdx_left(ns)+sweepDur-1,nt);
pos_sweep(:,2) = int_panelps_r(peakIdx_right(ns):peakIdx_right(ns)+sweepDur-1,nt);


%% calculate mean and error (sem)
% right sweep means
mean_vR(:,1) = mean(int_forward_rightward,2,'omitnan');
mean_vR(:,2) = mean(int_angular_rightward,2,'omitnan');
mean_vR(:,3) = mean(int_sideway_rightward,2,'omitnan');
% left sweep means
mean_vL(:,1) = mean(int_forward_leftward,2,'omitnan');
mean_vL(:,2) = mean(int_angular_leftward,2,'omitnan');
mean_vL(:,3) = mean(int_sideway_leftward,2,'omitnan');

ntot = size(int_forward_rightward,2);
% right sweep SEM
sem_vR(:,1) = std(int_forward_rightward,0,2,'omitnan')/sqrt(ntot);
sem_vR(:,2) = std(int_angular_rightward,0,2,'omitnan')/sqrt(ntot);
sem_vR(:,3) = std(int_sideway_rightward,0,2,'omitnan')/sqrt(ntot);
% left sweep SEM
sem_vL(:,1) = std(int_forward_leftward,0,2,'omitnan')/sqrt(ntot);
sem_vL(:,2) = std(int_angular_leftward,0,2,'omitnan')/sqrt(ntot);
sem_vL(:,3) = std(int_sideway_leftward,0,2,'omitnan')/sqrt(ntot);


%% plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1500 800])

    % set plot variables
    p1 = round(min(pos_sweep(:,1)));
    p2 = round(max(pos_sweep(:,1)));
    pm = 0;
    plotlabels = {'fwd (mm/s)'; 'ang (deg/s)'; 'side (mm/s)'};

    for dv = 1:3
        subplot(1,3,dv)
        % plot sem as patch
        r(dv) = patch([pos_sweep(:,1); flipud(pos_sweep(:,1))],[mean_vR(:,dv)-sem_vR(:,dv); flipud(mean_vR(:,dv)+sem_vR(:,dv))], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
        r(dv).FaceColor = '#ff0080';
        hold on
        l(dv) = patch([pos_sweep(:,2); flipud(pos_sweep(:,2))],[mean_vL(:,dv)-sem_vL(:,dv); flipud(mean_vL(:,dv)+sem_vL(:,dv))], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
        l(dv).FaceColor = '#0032A0';
        hold on

        % plot mean as line
        plot(pos_sweep(:,1),mean_vR(:,dv),'Color', '#ff0080','LineWidth',2)
        hold on
        plot(pos_sweep(:,2),mean_vL(:,dv),'Color', '#0032A0','LineWidth',2)

        % adjust axes
        axis tight
        xlim([p1 p2])
        xticks([p1 0 p2])
        xline(pm,'Color','k')

        if dv~=1
            minV = min([mean_vR(:,dv) mean_vR(:,dv)],[],'all');
            maxV = max([mean_vR(:,dv) mean_vR(:,dv)],[],'all');
            cV = (maxV-minV)/2 + minV;
            if ~isnan(cV)
                yline(cV,'Color','k')
            end
        end
        ylabel(plotlabels{dv})
        hold off
    end
    xlabel('panel pos (deg)')
end

end

