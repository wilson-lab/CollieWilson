% spikerate_v_fullsweep2
% analysis function for generating a summary plot of panel position versus
% spike rate. plots left and right sweeps across same position together.
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
% highThresh - pursuit threshold
% lowThresh - walk/rest threshold
% trackLR - for ephys, 0 for none, 1 for L, 2 for R
%
% ORIGINAL: 06/21/2022 - MC
% MODIFIED: 09/07/2022 - MC added laterality
% MODIFIED: 09/13/2022 - MC seperated walking and rest
% MODIFIED: 03/18/2024 - MC added optional restriction that fly must be
% running for full sweep duration
%

function [mean_sr,pos_sweep] = spikerate_v_fullsweep2(int_panelps, int_forward,int_time,int_spikerate,highThresh,lowThresh,minRun,trackLR, optPlot)
%% optional, threshold behavior using forward velocity
%pull ONLY pursuit behavior above HIGH threshold
[~,~,~,int_spikerate_high] = pursuitFinder(int_forward,0,0,int_spikerate,int_time,highThresh);

%pull ONLY walking behavior below LOW threshold
int_spikerate_low = int_spikerate;
walking_idx = int_forward<lowThresh & int_forward>0.1;
int_spikerate_low(~walking_idx) = NaN;

%pull ONLY rest behavior
int_spikerate_rest = int_spikerate;
rest_idx = int_forward<0.1 & int_forward>-0.1;
int_spikerate_rest(~rest_idx) = NaN;

%% pull spikerate for each sweep

% remove panel data noise
int_panelps_r = round((int_panelps*2))/2; %round to nearest 0.5

% for each trial
nTrial = size(int_panelps,2);
for nt = 1:nTrial
    %start sweep when object crosses midline (right-left-right)
    % find where sweep crosses midline
    midpointIdx = find(int_panelps_r(:,nt)==0);
    % ensure no oversampling has occurred
    midpointIdx = midpointIdx([1; find(diff(midpointIdx)>1)+1]);


    % depending on which side the cell is on
    if trackLR=='L'
        % for a left hemisphere cell, select right-left-right sweep
        sweepIdx = midpointIdx(1:2:end);
    else
        % for a right hemisphere cell, select left-right-left sweep
        sweepIdx = midpointIdx(2:2:end);
    end
    nSweeps = length(sweepIdx)-1; %ignore last (incomplete) sweep

    % if first trial, initialize
    if nt==1
        sweepDur = round(mean(sweepIdx(2:end)-sweepIdx(1:end-1)))+1;
        sweep_spikerate_all = NaN(sweepDur,nTrial*nSweeps);
        sweep_spikerate_high = NaN(sweepDur,nTrial*nSweeps);
        sweep_spikerate_low = NaN(sweepDur,nTrial*nSweeps);
        sweep_spikerate_rest = NaN(sweepDur,nTrial*nSweeps);
        is=0; %sweep index
    end

    % pull data by selected sweep start/stop
    for ns = 1:nSweeps
        dIdx = ns + is; %data index
        sweep_spikerate_all(1:sweepDur,dIdx) = int_spikerate(sweepIdx(ns):sweepIdx(ns)+sweepDur-1,nt);
        sweep_spikerate_high(1:sweepDur,dIdx) = int_spikerate_high(sweepIdx(ns):sweepIdx(ns)+sweepDur-1,nt);
        sweep_spikerate_low(1:sweepDur,dIdx) = int_spikerate_low(sweepIdx(ns):sweepIdx(ns)+sweepDur-1,nt);
        sweep_spikerate_rest(1:sweepDur,dIdx) = int_spikerate_rest(sweepIdx(ns):sweepIdx(ns)+sweepDur-1,nt);

        pos_singlesweep(:,ns,nt) = int_panelps_r(sweepIdx(ns):sweepIdx(ns)+sweepDur-1,nt);
    end
    is = is+nSweeps; %update sweep index

end
pos_sweep = mean(mean(pos_singlesweep,3),2);
int_time_sweep = int_time(1:sweepDur)';


%% optional restriction: omit sweeps where fly started/stopped mid sweep
% nTot = size(sweep_spikerate_all,2);
% trialLength = size(sweep_spikerate_all,1);
% minEx = 0.3; %amount of time the fly could have spent NOT in the behavior state, resulting in exclusion
% minLength = round(trialLength*minEx);
% 
% for t = 1:nTot
%     % for pursuit trials
%     if sum(isnan(sweep_spikerate_high(:,t)))>minLength
%         sweep_spikerate_high(:,t) = nan;
%     end
% 
%     % for walking trials
%     if sum(isnan(sweep_spikerate_low(:,t)))>minLength
%         sweep_spikerate_low(:,t) = nan;
%     end
% 
%     % for rest trials
%     if sum(isnan(sweep_spikerate_rest(:,t)))>minLength
%         sweep_spikerate_rest(:,t) = nan;
%     end
% end


%% calculate mean and error (sem)
% right sweep means
mean_sr(:,1) = mean(sweep_spikerate_all,2,'omitnan');
mean_sr(:,2) = mean(sweep_spikerate_high,2,'omitnan');
mean_sr(:,3) = mean(sweep_spikerate_low,2,'omitnan');
mean_sr(:,4) = mean(sweep_spikerate_rest,2,'omitnan');

ntot = size(sweep_spikerate_all,2);
% right sweep SEM
sem_sr(:,1) = std(sweep_spikerate_all,0,2,'omitnan')/sqrt(ntot);
sem_sr(:,2) = std(sweep_spikerate_high,0,2,'omitnan')/sqrt(ntot);
sem_sr(:,3) = std(sweep_spikerate_low,0,2,'omitnan')/sqrt(ntot);
sem_sr(:,4) = std(sweep_spikerate_rest,0,2,'omitnan')/sqrt(ntot);

n = size(sem_sr,2); %store number of variables


%% determine if experiment met minimum requirements to be considered "pursuit"

% min threshold set as time spent running above 5mm/s (e.g., likely pursuit)
minThresh = minRun; %sec
timespentchasing = (sum(int_forward>highThresh,'all')/length(int_time))*60;
chaseLog = (timespentchasing>minThresh);


%% plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1500 800])

    % set plot variables
    srmax = ceil(max(mean_sr,[],'all'));
    srmin = 0;
    t1 = 0;
    t2 = round(max(int_time_sweep),2);
    tm = mean(int_time_sweep);
    plotlabels = {'all (spikes/sec)'; 'pursuit (spikes/sec)'; 'walking (spikes/sec)'; 'rest (spikes/sec)'};

    for dv = 1:n
        if ~chaseLog && dv==2
            %skip pursuit plot if threshold failed
        else
            subplot(1,n,dv)
            fx = gca;

            yyaxis right
            % plot bar position
            plot(int_time_sweep,pos_sweep,'k','LineWidth',1.5)
            fx.YAxis(1).Color = 'k';
            if dv==n
                fx.YAxis(2).Color = 'k';
                ylabel('Object Position (deg)')
            else
                fx.YAxis(2).Visible = 'off';
            end
            yline(0,'Color','k')

            yyaxis left
            % plot sem as patch
            vx(dv) = patch([int_time_sweep; flipud(int_time_sweep)],[mean_sr(:,dv)-sem_sr(:,dv); flipud(mean_sr(:,dv)+sem_sr(:,dv))], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
            vx(dv).FaceColor = '#77AC30';
            hold on
            % plot mean as line
            plot(int_time_sweep,mean_sr(:,dv),'Color', '#77AC30','LineWidth',2)
            ylabel(plotlabels{dv})

            % adjust axes
            axis tight
            ylim([srmin srmax])
            xlim([t1 t2])
            xticks(t1:1:t2)
            xline(tm,'Color','k') %midline
            hold off
        end
    end
    xlabel('time (sec)')
end

end

