% spikerate_v_behavior_xcorr
% Analysis function that computes an estimate of the "ideal lag" between
% spike rate and each directional velocity by taking the mean
% cross-correlation across all trials in a data set.
%
% INPUTS:
% allForward - forward velocity array
% allAngular - angular velocity array
% allSideways - sidways velocity array
% allSpikeRate - spike rate array
% expttime - matching time array for single trial
%
% OUTPUTS:
% model - lags, and r value for forward, angular, and sideways velocities
%
% CREATED: 12/12/2022 MC
%

function model = spikerate_v_behavior_xcorr(allForward,allAngular,allSideway,allSpikeRt,expttime)
%% set parameters

% pull sample rate
sr = 20e3;
% whether to plot model output (1=yes, 0=no)
plt=1;


%% downsample

% set downsampling parameters
ds_t = 10; %msec, time to downsample to
ds_tSR = (ds_t/1000)*sr; %convert to sample rate
sz = size(allForward,1); %sample size

% downsample velocity data
int_forward = interp1((1:sz),allForward,(1:ds_tSR:sz),'linear');
int_angular = interp1((1:sz),allAngular,(1:ds_tSR:sz),'linear');
int_sideway = interp1((1:sz),allSideway,(1:ds_tSR:sz),'linear');

% downsample spikerate data
int_spikert = interp1((1:sz),allSpikeRt,(1:ds_tSR:sz),'linear');
% downsample time data
int_time = interp1((1:sz),expttime,(1:ds_tSR:sz),'linear')';


%% perform cross correlation

% initialize
f_r = [];
a_r = [];
s_r = [];
f_pk = [];
a_pk = [];
s_pk = [];

% run xcorr for each trial
c_range = 50; %xcorr window
nTrial = size(int_spikert,2);
for t = 1:nTrial
    [f_r(:,t),lags]=xcorr(int_spikert(:,t),int_forward(:,t),c_range);
    [a_r(:,t),~]=xcorr(int_spikert(:,t),int_angular(:,t),c_range);
    [s_r(:,t),~]=xcorr(int_spikert(:,t),int_sideway(:,t),c_range);
end

% store output variables
model.lags = lags'.*ds_t; %convert lag to time (msec.)
model.f_r = mean(f_r,2);
model.a_r = mean(a_r,2);
model.s_r = mean(s_r,2);


%% plot

% find peaks
[~,f_pk] = findpeaks(model.f_r);
[~,a_pk] = findpeaks(model.a_r);
[~,s_pk] = findpeaks(model.s_r);


% initialize
if plt
    figure; set(gcf,'Position',[100 100 1000 500])
    colorselect = {'#D95319';'#0072BD';'#7E2F8E'}; %velocity colors
    velocities = {'forward'; 'angular'; 'sideway'}; %velocity names
    lw = 1; %linewidth
    lw2 = 3; %linewidth

    subplot(1,3,1)
    hold on
    plot(model.lags,f_r,'Color','k','LineWidth',lw,'LineStyle',':')
    plot(model.lags,model.f_r,'Color',colorselect{1},'LineWidth',lw2)
    xline(model.lags(f_pk),'LineWidth',lw2)
    ylabel('r-value')
    xlabel([velocities{1} ' shift (msec)'])
    legend([num2str(model.lags(f_pk)) 'msec'])
    hold off

    subplot(1,3,2)
    hold on
    plot(model.lags,a_r,'Color','k','LineWidth',lw,'LineStyle',':')
    plot(model.lags,model.a_r,'Color',colorselect{2},'LineWidth',lw2)
    xline(model.lags(a_pk),'LineWidth',lw2)
    xlabel([velocities{2} ' shift (msec)'])
    legend([num2str(model.lags(a_pk)) 'msec'])
    hold off

    subplot(1,3,3)
    hold on
    plot(model.lags,s_r,'Color','k','LineWidth',lw,'LineStyle',':')
    plot(model.lags,model.s_r,'Color',colorselect{3},'LineWidth',lw2)
    xline(model.lags(s_pk),'LineWidth',lw2)
    xlabel([velocities{3} ' shift (msec)'])
    legend([num2str(model.lags(s_pk)) 'msec'])
    hold off
end


end