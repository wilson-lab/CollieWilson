% trialtype_performance
% analysis function for assessing pursuit performance across trial types
%
% INPUTS:
% panelps - downsampled panel positions
% forward - downsampled forward velocities
% angular - downsampled angular velocities
% sideway - downsampled sideways velocities
% time - downsampled trial time
% chaseThreshold - forward velocity threshold, + for behavior >, - for behavior <
%
% ORIGINAL: 04/05/2023 - MC created velocity_performance
%

function trialtype_performance(panelps, forward,angular,sideway,time,highThreshold,lowThreshold)
%% initialize and prepare data for processing
% generate plot
nTrial = size(panelps,2);
if nTrial < 15
    figwidth = 100*size(panelps,2);
else
    figwidth=1500;
end
figure; set(gcf,'Position',[100 100 figwidth 800])

% pull number of trial types
nTypes = size(panelps,3);

%% plot average behavior over experiment duration

%for each trial type
for nt = 1:nTypes
    % calculate time spent in pursuit for each trial over time
    pursuit_idx = schmittTrigger(forward(:,:,nt),highThreshold,lowThreshold); %find pursuit indices
    trialPursuit = sum(pursuit_idx,1)*time(2); %total time spent pursuing
    % calculate average velocity for each trial over time
    trialfwd = mean(forward(:,:,nt),1);
    trialang = mean(abs(angular(:,:,nt)),1);

    % plot everything together
    subplot(3,1,1)
    hold on; plot(trialPursuit,':o')

    subplot(3,1,2)
    hold on; plot(trialfwd,':o')

    subplot(3,1,3)
    hold on; plot(trialang,':o')
end
% add labels
subplot(3,1,1); ylabel('Time Pursuing (sec.)')
subplot(3,1,2); ylabel('Average Forward (mm/s)')
subplot(3,1,3); ylabel('Average Angular (deg/s)')
xlabel('Trial Number')


end