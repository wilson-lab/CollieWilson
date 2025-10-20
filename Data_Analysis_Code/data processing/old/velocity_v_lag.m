% velocity_v_lag
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
% ORIGINAL: 03/29/2023 MC original, converted from fullsweep plot analysis
%

function velocity_v_lag(int_panelps, int_forward,int_angular,int_sideway,int_time,highThresh,lowThresh,optPlot)
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
ymax = 150;

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
        cross_angular = nan(nCross,nTrial);
        cross_forward = nan(nCross,nTrial);
    end

    % for each cross point
    for c = 1:length(crosspoint)-2
        % if the fly was engaged in pursuit for both this sweep and the next
        if pursuit_idx(crosspoint(c),nt)&&pursuit_idx(crosspoint(c+1),nt)
            % pull indices for this and next sweep
            thisSweep = crosspoint(c):crosspoint(c+1);
            % normalize cross angular velocity so + is ipsi and - is contra
            if median(int_panelps_adj(thisSweep,nt))>0
                a = +1;
            else
                a = -1;
            end

            % pull behavior variables of interest
            peak_angular(c,nt) = max(abs(int_angular(thisSweep,nt)));
            peak_forward(c,nt) = max(int_forward(thisSweep,nt));

            cross_angular(c,nt) = int_angular(crosspoint(c+1),nt)*a;
            cross_forward(c,nt) = int_forward(crosspoint(c+1),nt);
        end
    end
end

%% analyze relationships

% omit nans from fit
omitIdx = isnan(peak_angular);

% generate linear fit models
mdl_peak_angular = fitlm(peak_angular(~omitIdx),cross_angular(~omitIdx));
mdl_peak_forward = fitlm(peak_forward(~omitIdx),cross_angular(~omitIdx));
mdl_cross_forward = fitlm(cross_forward(~omitIdx),cross_angular(~omitIdx));

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
        ylabel('Angular Velocity at Cross (deg/sec)')
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
        ylabel('Angular Velocity at Cross (deg/sec)')
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
        ylabel('Angular Velocity at Cross (deg/sec)')
        xlabel('Forward at Cross (mm/sec)')
    end
end


end