% spikerate_v_behaviorlag
% Analysis function that computes an estimate of the "ideal lag" between
% spike rate and each directional velocity by incrementally shifting the
% lag time between the two and calculating the r2 value for each.
%
% INPUTS:
% allForward - forward velocity array
% allAngular - angular velocity array
% allSideways - sidways velocity array
% allSpikeRate - spike rate array
% expttime - matching time array for single trial
%
% OUTPUTS:
% fwd_model - forward slope, intercept, and r-squared for each shift
% ang_model - angular slope, intercept, and r-squared for each shift
% sid_model - sideway slope, intercept, and r-squared for each shift
%
% CREATED: 12/07/2022 MC
% UPDATED: 12/08/2022 MC implemented shift to behavior instead of activity
%

function [fwd_model,ang_model,sid_model] = spikerate_v_behaviorlag(allForward,allAngular,allSideways,allSpikeRate,expttime)
%% set parameters
shift_range = 500; %range to try, in msec.
shift_incrt = 10; %range increment, in msec.

shift = -shift_range:shift_incrt:shift_range;

% pull sample rate
sr = 20e3; %expt sample rate from settings

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
int_sideway = interp1((1:sz),allSideways,(1:ds_tSR:sz),'linear');

% downsample spikerate data
int_spikert = interp1((1:sz),allSpikeRate,(1:ds_tSR:sz),'linear');
% reshape spikerate data
rsh_spikert = reshape(int_spikert,[],1);
% downsample time data
int_time = interp1((1:sz),expttime,(1:ds_tSR:sz),'linear')';

% calculat new sample rate based on downsample
sr_ds = sr/ds_tSR; %new sample rate


%% implement shifts to dataset

% convert incremental shifts to new sample rate
shift_setDS = (shift./1000).*sr_ds;
nShift = length(shift_setDS); %n
% must blank out beginning/end of trial to avoid artifacts caused by
% different amount of NaNs across different shifts
shiftBlank = max(shift_setDS);

% for each possible shift...
for sh = 1:nShift
    % select shift for this itteration
    thisShift = round(shift_setDS(sh));

    % shift data accordingly
    shift_forward = circshift(int_forward,thisShift,1);
    shift_angular = circshift(int_angular,thisShift,1);
    shift_sideway = circshift(int_sideway,thisShift,1);

    % regardless of shift direction...
    % remove start data
    shift_forward((end-shiftBlank+1):end,:) = NaN;
    shift_angular((end-shiftBlank+1):end,:) = NaN;
    shift_sideway((end-shiftBlank+1):end,:) = NaN;
    % remove stop data
    shift_forward(1:shiftBlank,:) = NaN;
    shift_angular(1:shiftBlank,:) = NaN;
    shift_sideway(1:shiftBlank,:) = NaN;
    
    % reshape into single column
    track_forward(:,sh) = reshape(shift_forward,[],1);
    track_angular(:,sh) = reshape(shift_angular,[],1);
    track_sideway(:,sh) = reshape(shift_sideway,[],1);
end


%% assess relationship b/n spike rate and behavior data across different shifts

% for each possible shift...
for sh = 1:nShift
    % forward analysis
    % fit linear regression model
    f_mdl = fitlm(rsh_spikert,track_forward(:,sh),'linear');
    % store values
    intercept(sh,1) = f_mdl.Coefficients.Estimate(1);
    coefficient(sh,1) = f_mdl.Coefficients.Estimate(2);
    rsquared(sh,1) = f_mdl.Rsquared.Adjusted;
    
    % angular analysis
    % fit linear regression model
    a_mdl = fitlm(rsh_spikert,track_angular(:,sh),'linear');
    % store values
    intercept(sh,2) = a_mdl.Coefficients.Estimate(1);
    coefficient(sh,2) = a_mdl.Coefficients.Estimate(2);
    rsquared(sh,2) = a_mdl.Rsquared.Adjusted;
    
    % sideway analysis
    % fit linear regression model
    s_mdl = fitlm(rsh_spikert,track_sideway(:,sh),'linear');
    % store values
    intercept(sh,3) = s_mdl.Coefficients.Estimate(1);
    coefficient(sh,3) = s_mdl.Coefficients.Estimate(2);
    rsquared(sh,3) = s_mdl.Rsquared.Adjusted;
    
end
% store all model outputs
varNames = {'shift','intercept','coefficient','rsquared'};
fwd_model = table(shift',intercept(:,1),coefficient(:,1),rsquared(:,1),'VariableNames',varNames);
ang_model = table(shift',intercept(:,2),coefficient(:,2),rsquared(:,2),'VariableNames',varNames);
sid_model = table(shift',intercept(:,3),coefficient(:,3),rsquared(:,3),'VariableNames',varNames);


%% plot

if plt
    % initialize
    figure; set(gcf,'Position',[100 100 800 500])
    colorselect = {'#D95319';'#0072BD';'#7E2F8E'}; %velocity colors
    ylabels = {'forward'; 'angular'; 'sideway'}; %velocity names
    lw = 3; %linewidth
    
    hold on
    plot(fwd_model.shift,fwd_model.rsquared,'Color',colorselect{1},'LineWidth',lw)
    plot(ang_model.shift,ang_model.rsquared,'Color',colorselect{2},'LineWidth',lw)
    plot(sid_model.shift,sid_model.rsquared,'Color',colorselect{3},'LineWidth',lw)
    
    legend(ylabels)
    ylabel('rsquared value')
    xlabel('shift (msec.)')
    hold off
end


end