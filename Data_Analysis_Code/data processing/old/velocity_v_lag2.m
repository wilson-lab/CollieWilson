% velocity_v_lag2
% analysis function for assessing how behavior during previous turn may
% influence lag of subsequent turn
%
% OUTPUTS:
%
% INPUTS:
% int_panelps - downsampled panel positions
% int_forward - downsampled forward velocities
% int_angular - downsampled angular velocities
% int_sideway - downsampled sideways velocities
% int_time    - downsampled trial time
% highThresh  - pursuit forward velocity threshold
%
% ORIGINAL: 03/30/2023 MC original, converted from fullsweep plot analysis
%

function velocity_v_lag2(int_panelps, int_forward,int_angular,int_sideway,int_time,highThresh,lowThresh,optPlot)
%% initialize and prepare data for processing

% remove panel position noise for easier processing
int_panelps_adj = round((int_panelps*2))/2; %round to nearest 0.5

% find pursuit indices using schmitt trigger
pursuit_idx = schmittTrigger(int_forward,highThresh,lowThresh);

% pull trial number
nTrial = size(int_panelps,2);

% plot variables
maxAng = 200;
%maxAng = ceil(max(peak_angular,[],'all'));
maxFwd = 12;
%maxFwd = ceil(max([peak_forward cross_forward],[],'all'));
ymax = 60;

%% pull behavior for each sweep

for nt = 1:nTrial
    % find where target crosses the midline in either direction
    crosspoint = find(int_panelps_adj(:,nt)==0);
    % find center of each cross (as each px point sampled across time)
    crosspoint = crosspoint([1; find(diff(crosspoint)>1)+1]);
    nCross = length(crosspoint);

    % if first, initialize
    if nt ==1
        % initialize
        peak_angular = nan(nCross,nTrial);
        peak_forward = nan(nCross,nTrial);
        cross_forward = nan(nCross,nTrial);

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
            peak_angular(c,nt) = max(abs(int_angular(priorSweep,nt)));
            peak_forward(c,nt) = max(int_forward(priorSweep,nt));
            cross_forward(c,nt) = int_forward(crosspoint(c+1),nt);

            % pull delay
            [~,cross_min_fly] = min(abs(int_angular(crossSweep,nt)));
            cross_idx_fly = crossSweep(cross_min_fly);
            turn_delay(c,nt) = (int_time(cross_idx_fly) - int_time(cross_idx_target))*100; %msec
        end
    end
end

%% analyze relationships

% omit nans from fit
omitIdx = isnan(peak_angular);

% generate linear fit models
mdl_peak_angular = fitlm(peak_angular(~omitIdx),turn_delay(~omitIdx));
mdl_peak_forward = fitlm(peak_forward(~omitIdx),turn_delay(~omitIdx));
mdl_cross_forward = fitlm(cross_forward(~omitIdx),turn_delay(~omitIdx));

%% plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1500 400])

    if mdl_peak_angular.Coefficients.Estimate(1)~=0

        % plot peak angular velocity vs lag
        subplot(1,3,1)
        x1 = plot(mdl_peak_angular);
        x1(1).Marker = 'o';
        x1(1).Color = '#0072BD';
        x1(2).Color = 'k';
        x1(3).Color = 'k';
        x1(4).Color = 'k';
        xlim([0 maxAng])
        ylim([-ymax ymax])
        title(['R2 = ' num2str(mdl_peak_angular.Rsquared.Adjusted)])
        ylabel('Turn Delay (msec.)')
        xlabel('Previous Peak Angular (deg/sec)')

        % plot peak forward velocity vs lag
        subplot(1,3,2)
        x2 = plot(mdl_peak_forward);
        x2(1).Marker = 'o';
        x2(1).Color = '#D95319';
        x2(2).Color = 'k';
        x2(3).Color = 'k';
        x2(4).Color = 'k';
        xlim([0 maxFwd])
        ylim([-ymax ymax])
        title(['R2 = ' num2str(mdl_peak_forward.Rsquared.Adjusted)])
        ylabel('Turn Delay (msec.)')
        xlabel('Previous Peak Forward (mm/sec)')

        % plot forward velocity at cross vs lag
        subplot(1,3,3)
        x3 = plot(mdl_cross_forward);
        x3(1).Marker = 'o';
        x3(1).Color = '#D95319';
        x3(2).Color = 'k';
        x3(3).Color = 'k';
        x3(4).Color = 'k';
        xlim([0 maxFwd])
        ylim([-ymax ymax])
        title(['R2 = ' num2str(mdl_cross_forward.Rsquared.Adjusted)])
        ylabel('Turn Delay (msec.)')
        xlabel('Forward at Cross (mm/sec)')
    end
end


end