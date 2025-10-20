% spikerate_v_panelpos2
% analysis function for generating a summary plot of panel position versus
% spike rate. plots left and right sweeps across same position separately.
%
% OUTPUTS:
% mean_vR - mean velocities for rightward sweep
% mean_vL - mean velocities for leftward sweep
%
% INPUTS:
% int_panelps - downsampled panel positions
% int_forward - downsampled forward velocities
% int_spikerate - downsampled spikerate
% int_time - downsampled trial time
%
% ORIGINAL: 06/21/2022 - MC
%

function [mean_srR,mean_srL,pos_sweep] = spikerate_v_panelpos2(int_panelps, int_forward,int_time,int_spikerate,chaseThreshold,lowThresh,minRun,optPlot)
%% optional, threshold behavior using forward velocity
chaseThresh = chaseThreshold; %mm/s
nonchaseMax = lowThresh; %mm/s
%pull ONLY behavior above threshold
int_spikerateH = int_spikerate;
chase_IdxH = schmittTrigger(int_forward,chaseThresh,1);
int_spikerateH(~chase_IdxH) = NaN;
%pull ONLY behavior above 0 and below threshold
int_spikerateL = int_spikerate;
nonChase_Idx = int_forward>nonchaseMax | int_forward<0;
int_spikerateL(nonChase_Idx) = NaN;

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
        int_spiker_sweepR=NaN(sweepDur,nTrial*nSweeps);
        int_spikerH_sweepR=NaN(sweepDur,nTrial*nSweeps);
        int_spikerL_sweepR=NaN(sweepDur,nTrial*nSweeps);
        int_spiker_sweepL=NaN(sweepDur,nTrial*nSweeps);
        int_spikerH_sweepL=NaN(sweepDur,nTrial*nSweeps);
        int_spikerL_sweepL=NaN(sweepDur,nTrial*nSweeps);
        is=0; %sweep index
    end

    % pull sweeps
    for ns = 1:nSweeps
        dIdx = ns + is; %data index
        % pull rightward sweeps by taking sweep following left-most peak
        int_spiker_sweepR(1:sweepDur,dIdx) = int_spikerate(peakIdx_left(ns):peakIdx_left(ns)+sweepDur-1,nt);
        int_spikerH_sweepR(1:sweepDur,dIdx) = int_spikerateH(peakIdx_left(ns):peakIdx_left(ns)+sweepDur-1,nt);
        int_spikerL_sweepR(1:sweepDur,dIdx) = int_spikerateL(peakIdx_left(ns):peakIdx_left(ns)+sweepDur-1,nt);
        % pull leftward sweeps by taking sweep following right-most peak
        int_spiker_sweepL(1:sweepDur,dIdx) = int_spikerate(peakIdx_right(ns):peakIdx_right(ns)+sweepDur-1,nt);
        int_spikerH_sweepL(1:sweepDur,dIdx) = int_spikerateH(peakIdx_right(ns):peakIdx_right(ns)+sweepDur-1,nt);
        int_spikerL_sweepL(1:sweepDur,dIdx) = int_spikerateL(peakIdx_right(ns):peakIdx_right(ns)+sweepDur-1,nt);
    end
    is = is+nSweeps; %update sweep index
end

% pull sweep px positions
pos_sweep(:,1) = int_panelps_r(peakIdx_left(ns):peakIdx_left(ns)+sweepDur-1,nt);
pos_sweep(:,2) = int_panelps_r(peakIdx_right(ns):peakIdx_right(ns)+sweepDur-1,nt);


%% calculate mean and error (sem)
% right sweep means
mean_srR(:,1) = mean(int_spiker_sweepR,2,'omitnan');
mean_srR(:,2) = mean(int_spikerH_sweepR,2,'omitnan');
mean_srR(:,3) = mean(int_spikerL_sweepR,2,'omitnan');
% left sweep means
mean_srL(:,1) = mean(int_spiker_sweepL,2,'omitnan');
mean_srL(:,2) = mean(int_spikerH_sweepL,2,'omitnan');
mean_srL(:,3) = mean(int_spikerL_sweepL,2,'omitnan');

ntot = size(int_spiker_sweepR,2);
% right sweep SEM
sem_srR(:,1) = std(int_spiker_sweepR,0,2,'omitnan')/sqrt(ntot);
sem_srR(:,2) = std(int_spikerH_sweepR,0,2,'omitnan')/sqrt(ntot);
sem_srR(:,3) = std(int_spikerL_sweepR,0,2,'omitnan')/sqrt(ntot);
% left sweep SEM
sem_srL(:,1) = std(int_spiker_sweepL,0,2,'omitnan')/sqrt(ntot);
sem_srL(:,2) = std(int_spikerH_sweepL,0,2,'omitnan')/sqrt(ntot);
sem_srL(:,3) = std(int_spikerL_sweepL,0,2,'omitnan')/sqrt(ntot);


%% determine if experiment met minimum requirements to be considered "pursuit"

% min threshold set as time spent running above 5mm/s (e.g., likely pursuit)
minThresh = minRun; %sec
timespentchasing = (sum(int_forward>chaseThresh,'all')/length(int_time))*60;
chaseLog = (timespentchasing>minThresh);


%% plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1500 800])

    % set plot variables
    srmax = ceil(max([mean_srR mean_srL],[],'all'));
    srmin = floor(min([mean_srR mean_srL],[],'all'));
    p1 = round(min(pos_sweep(:,1)));
    p2 = round(max(pos_sweep(:,1)));
    pm = 0;
    plotlabels = {'all (spikes/sec)'; 'chasing (spikes/sec)'; 'walking (spikes/sec)'};

    for dv = 1:3
        if ~chaseLog && dv==2
            %skip pursuit plot if threshold failed
        else
            subplot(1,3,dv)
            % plot sem as patch
            r(dv) = patch([pos_sweep(:,1); flipud(pos_sweep(:,1))],[mean_srR(:,dv)-sem_srR(:,dv); flipud(mean_srR(:,dv)+sem_srR(:,dv))], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
            r(dv).FaceColor = '#ff0080';
            hold on
            l(dv) = patch([pos_sweep(:,2); flipud(pos_sweep(:,2))],[mean_srL(:,dv)-sem_srL(:,dv); flipud(mean_srL(:,dv)+sem_srL(:,dv))], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
            l(dv).FaceColor = '#0032A0';
            hold on

            % plot mean as line
            plot(pos_sweep(:,1),mean_srR(:,dv),'Color', '#ff0080','LineWidth',2)
            hold on
            plot(pos_sweep(:,2),mean_srL(:,dv),'Color', '#0032A0','LineWidth',2)

            % adjust axes
            axis tight
            ylim([srmin srmax])
            xlim([p1 p2])
            xticks([p1 0 p2])
            xline(pm,'Color','k')
            ylabel(plotlabels{dv})
            hold off
        end
    end
    xlabel('panel pos (deg)')
end

end

