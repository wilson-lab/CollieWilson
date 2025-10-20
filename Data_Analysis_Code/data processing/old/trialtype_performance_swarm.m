% trialtype_performance_swarm
% analysis function for assessing pursuit performance across trial types
%
% INPUTS:
% forward - downsampled forward velocities
% angular - downsampled angular velocities
% trialTypes - cell containing trial type names
%
% ORIGINAL: 04/05/2023 - MC created velocity_performance
%

function trialtype_performance_swarm(forward,angular,trialTypes)
%% initialize and prepare data for processing
% pull trial types
nTypes = size(forward,3);

% generate plot
figwidth = 150*nTypes;
figure; set(gcf,'Position',[100 100 figwidth 600])

%% run analysis

% downsample
forward_ds = interp1((1:length(forward)),forward,(1:100:length(forward)),'linear');
angular_ds = interp1((1:length(angular)),angular,(1:100:length(angular)),'linear');
xi = ones(1,size(forward_ds,1)*size(forward_ds,2));

% for each trial type
for nt = 1:nTypes
    subplot(2,1,1)
    swarmchart(xi*nt,reshape(forward_ds(:,:,nt),[],1),5)
    hold on

    subplot(2,1,2)
    swarmchart(xi*nt,reshape(angular_ds(:,:,nt),[],1),5)
    hold on
end

subplot(2,1,1)
ylabel('Forward Velocity (mm/s)')
xticks(1:nTypes)
xticklabels(trialTypes)

subplot(2,1,2)
ylabel('Angular Velocity (mm/s)')
yline(0)
xticks(1:nTypes)
xticklabels(trialTypes)


end