% pursuit_performance
% analysis function for assessing pursuit performance across repeated
% trials
%
% OUTPUTS:
% mean_v - mean velocities for full sweep
%
% INPUTS:
% panelps - downsampled panel positions
% forward - downsampled forward velocities
% angular - downsampled angular velocities
% sideway - downsampled sideways velocities
% time - downsampled trial time
% chaseThreshold - forward velocity threshold, + for behavior >, - for behavior <
%
% ORIGINAL: 04/04/2023 - MC created from lag analysis
%

function mean_delay_across = pursuit_performance(panelps, forward,angular,sideway,time,highThreshold,lowThreshold,optPlot)
%% initialize and prepare data for processing

% remove panel position noise for easier processing
panelps_adj = round((panelps*2))/2; %round to nearest 0.5

% find pursuit indices using schmitt trigger
pursuit_idx = schmittTrigger(forward,highThreshold,lowThreshold);

% pull trial number
nTrial = size(panelps,2);

%% pull behavior for each sweep

for nt = 1:nTrial
    % find where target crosses the midline in either direction
    crosspoint = find(panelps_adj(:,nt)==0);
    % find center of each cross (as each px point sampled across time)
    crosspoint = crosspoint([1; find(diff(crosspoint)>1)+1]);
    nCross = length(crosspoint);

    % if first, initialize
    if nt ==1
        % initialize
        peak_angular = nan(nCross,nTrial);
        peak_forward = nan(nCross,nTrial);

        turn_delay = nan(nCross,nTrial);

        crossWindow = round((crosspoint(2)-crosspoint(1))/2,-1);
    end

    % for each cross point
    for c = 1:length(crosspoint)-2
        % if the fly was engaged in pursuit for both this sweep and the next
        if pursuit_idx(crosspoint(c),nt)&&pursuit_idx(crosspoint(c+1),nt)
            % pull indices for this and next sweep
            priorSweep = crosspoint(c):crosspoint(c+1);
            cross_idx_target = crosspoint(c+1);
            crossSweep = crosspoint(c+1)-crossWindow:crosspoint(c+1)+crossWindow;

            % pull behavior variables of interest
            peak_angular(c,nt) = max(abs(angular(priorSweep,nt)));
            peak_forward(c,nt) = max(forward(priorSweep,nt));

            % pull delay
            [~,cross_min_fly] = min(abs(angular(crossSweep,nt)));
            cross_idx_fly = crossSweep(cross_min_fly);
            turn_delay(c,nt) = (time(cross_idx_fly) - time(cross_idx_target))*100; %msec
        end
    end
end


%% summary results

% performance WITHIN trials
mean_delay_within = nanmean(turn_delay,2);
std_delay_within = nanstd(turn_delay,0,2);
mean_angular_within = nanmean(peak_angular,2);
std_angular_within = nanstd(peak_angular,0,2);
mean_forward_within = nanmean(peak_forward,2);
std_forward_within = nanstd(peak_forward,0,2);

% performance ACROSS trials
mean_delay_across = nanmean(turn_delay,1);
std_delay_across = nanstd(turn_delay,0,1);
mean_angular_across = nanmean(peak_angular,1);
std_angular_across = nanstd(peak_angular,0,1);
mean_forward_across = nanmean(peak_forward,1);
std_forward_across = nanstd(peak_forward,0,1);

%% plot performance

if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1500 600])
    lw = 1.5;
    gr = [0.7 0.7 0.7];
    x = repmat(1:nTrial,size(turn_delay,1),1);

    % plot performance WITHIN trials
    % plot turn delay within trials
    subplot(3,2,1)
    plot(1:nCross,turn_delay,'-o','Color',gr)
    hold on; errorbar(1:nCross,mean_delay_within,std_delay_within,'-o','Color','#77AC30','LineWidth',lw)
    yline(0)
    xlim([0 nCross])
    title('Pursuit WITHIN trials')
    ylabel('Turn delay (msec.)')

    % plot forward velocity across trials
    subplot(3,2,3)
    plot(1:nCross,peak_forward,'o','Color',gr)
    hold on; errorbar(1:nCross,mean_forward_within,std_forward_within,'-o','Color','#D95319','LineWidth',lw)
    xlim([0 nCross])
    ylabel('Peak Forward (mm/sec)')

    % plot angular velocity across trials
    subplot(3,2,5)
    plot(1:nCross,peak_angular,'o','Color',gr)
    hold on; errorbar(1:nCross,mean_angular_within,std_angular_within,'-o','Color','#0072BD','LineWidth',lw)
    xlim([0 nCross])
    ylabel('Peak Angular (deg/sec)')


    % plot performance ACROSS trials
    % plot turn delay across trials
    subplot(3,2,2)
    plot(x,turn_delay,'o','Color',gr)
    hold on; errorbar(1:nTrial,mean_delay_across,std_delay_across,'-o','Color','#77AC30','LineWidth',lw)
    yline(0)
    xlim([0 nTrial+1])
    title('Pursuit ACROSS trials')
    ylabel('Turn delay (msec.)')

    % plot forward velocity across trials
    subplot(3,2,4)
    plot(x,peak_forward,'o','Color',gr)
    hold on; errorbar(1:nTrial,mean_forward_across,std_forward_across,'-o','Color','#D95319','LineWidth',lw)
    xlim([0 nTrial+1])
    ylabel('Peak Forward (mm/sec)')

    % plot angular velocity across trials
    subplot(3,2,6)
    plot(x,peak_angular,'o','Color',gr)
    hold on; errorbar(1:nTrial,mean_angular_across,std_angular_across,'-o','Color','#0072BD','LineWidth',lw)
    xlim([0 nTrial+1])
    ylabel('Peak Angular (deg/sec)')

    xlabel('Trial Number')
end

end