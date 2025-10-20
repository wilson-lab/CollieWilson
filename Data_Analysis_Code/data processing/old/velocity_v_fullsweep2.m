% velocity_v_fullsweep2
% analysis function for generating a summary plot of panel position versus
% directional velocity. plots full sweep.
%
% OUTPUTS:
% mean_v - mean velocities for full sweep
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

function [mean_v,pos_sweep] = velocity_v_fullsweep2(int_panelps, int_forward,int_angular,int_sideway,int_time,chaseThreshold,optPlot)
%% optional, threshold behavior using forward velocity
if chaseThreshold>0 %pull ONLY behavior above threshold
    [int_forward,int_angular,int_sideway,~] = pursuitFinder(int_forward,int_angular,int_sideway,0,int_time,chaseThreshold);
elseif chaseThreshold<0 %pull ONLY behavior below threshold
    nonChase_Idx = int_forward>(-chaseThreshold) | int_forward<0;
    int_forward(nonChase_Idx) = NaN;
    int_angular(nonChase_Idx) = NaN;
    int_sideway(nonChase_Idx) = NaN;
end


%% pull directional velocity for each full sweep

% remove panel data noise
int_panelps = round((int_panelps*2))/2; %round to nearest 0.5

% for each trial
nTrial = size(int_panelps,2);
for nt = 1:nTrial
    %start sweep when object crosses midline (right-left-right)
    % find where sweep crosses midline
    midpointIdx = find(int_panelps(:,nt)==0);
    % ensure no oversampling has occurred
    midpointIdx = midpointIdx([1; find(diff(midpointIdx)>1)+1]);
    % select a full right-left-right sweep by selecting every-other midpoint
    sweepIdx = midpointIdx(1:2:end);
    nSweeps = length(sweepIdx)-1; %ignore last (incomplete) sweep

    % if first trial, initialize
    if nt==1
        sweepDur = round(mean(sweepIdx(2:end)-sweepIdx(1:end-1)))+1;
        int_forward_sweep = NaN(sweepDur,nTrial*nSweeps);
        int_angular_sweep = NaN(sweepDur,nTrial*nSweeps);
        int_sideway_sweep = NaN(sweepDur,nTrial*nSweeps);
        is=0; %sweep index
    end

    % pull data by selected sweep start/stop
    for ns = 1:nSweeps
        dIdx = ns + is; %data index
        int_forward_sweep(1:sweepDur,dIdx) = int_forward(sweepIdx(ns):sweepIdx(ns)+sweepDur-1,nt);
        int_angular_sweep(1:sweepDur,dIdx) = int_angular(sweepIdx(ns):sweepIdx(ns)+sweepDur-1,nt);
        int_sideway_sweep(1:sweepDur,dIdx) = int_sideway(sweepIdx(ns):sweepIdx(ns)+sweepDur-1,nt);
    end
    is = is+nSweeps; %update sweep index

end
pos_sweep = int_panelps(sweepIdx(ns):sweepIdx(ns)+sweepDur-1,nt);
int_time_sweep = int_time(1:sweepDur)';


%% optional restriction: omit sweeps where fly started/stopped mid sweep
% nTot = size(int_forward_sweep,2);
% trialLength = size(int_forward_sweep,1);
% minEx = 0.3; %amount of time the fly could have spent NOT in the behavior state, resulting in exclusion
% minLength = round(trialLength*minEx);
% 
% for t = 1:nTot
%     % for pursuit trials
%     if sum(isnan(int_forward_sweep(:,t)))>minLength
%         int_forward_sweep(:,t) = nan;
%         int_angular_sweep(:,t) = nan;
%         int_sideway_sweep(:,t) = nan;
%     end
% end

%% calculate mean and error (sem)
mean_v(:,1) = mean(int_forward_sweep,2,'omitnan');
mean_v(:,2) = mean(int_angular_sweep,2,'omitnan');
mean_v(:,3) = mean(int_sideway_sweep,2,'omitnan');

ntot = size(int_forward_sweep,2);
sem_v(:,1) = std(int_forward_sweep,0,2,'omitnan')/sqrt(ntot);
sem_v(:,2) = std(int_angular_sweep,0,2,'omitnan')/sqrt(ntot);
sem_v(:,3) = std(int_sideway_sweep,0,2,'omitnan')/sqrt(ntot);

%% plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1500 800])

    % plot variables
    t1 = 0;
    t2 = max(int_time_sweep);
    tm = t1 + (t2-t1)/2;
    plotlabels = {'Forward Velocity (mm/sec)'; 'Angular Velocity (deg/sec)'; 'Sideways Velocity (mm/sec)'};
    colorlabels = {'#D95319';'#0072BD';'#7E2F8E'};

    for dv = 1:3
        subplot(1,3,dv)
        fx = gca;

        yyaxis right
        % plot bar position
        plot(int_time_sweep,pos_sweep,'k','LineWidth',1.5)
        fx.YAxis(1).Color = 'k';
        if dv==3
            fx.YAxis(2).Color = 'k';
            ylabel('obj pos (deg)')
        else
            fx.YAxis(2).Visible = 'off';
        end

        yyaxis left
        % plot sem as patch
        vx(dv) = patch([int_time_sweep; flipud(int_time_sweep)],[mean_v(:,dv)-sem_v(:,dv); flipud(mean_v(:,dv)+sem_v(:,dv))], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
        vx(dv).FaceColor = colorlabels{dv};
        hold on
        % plot mean as line
        plot(int_time_sweep,mean_v(:,dv),'Color', colorlabels{dv},'LineWidth',2)
        ylabel(plotlabels{dv})
        if dv~=1
            minV = min([mean_v(:,dv) mean_v(:,dv)],[],'all');
            maxV = max([mean_v(:,dv) mean_v(:,dv)],[],'all');
            cV = (maxV-minV)/2 + minV;
            if ~isnan(cV)
                yline(cV,':','Color','k')
            end
        end

        % adjust axes
        axis tight
        xlim([t1 t2])
        switch dv
            case 1
                max1 = ceil(max(mean_v(:,dv)));
                if isnan(max1)
                    ylim([0 5])
                else
                    ylim([0 max1])
                end
            case 2
                max2 = round(max(abs(mean_v(:,dv))),-1)+5;
                if isnan(max2)
                    ylim([-100 100])
                else
                    ylim([-max2 max2])
                end
            case 3
                max3 = round(max(abs(mean_v(:,dv))),1)+0.025;
                if isnan(max3)
                    ylim([-2 2])
                else
                    ylim([-max3 max3])
                end
        end
        xticks(t1:1:t2)
        xline(tm,':','Color','k')
        hold off
    end
    xlabel('Time (sec)')
end


end