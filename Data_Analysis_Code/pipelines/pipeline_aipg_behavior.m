% pipeline_battery_behavior
%
% This pipeline function processes and analyzes data from all flies
% in a given closed-loop experiment. It pulls all relevant processed files,
% performs the necessary analyses, and generates plots accordingly.
%
% INPUT:
%   exptFolder - String representing the path to the overarching experiment
%                folder containing processed data for all flies.
%
% Created: 07/24/2025 by MC
% Adapted from the kir closed loop pipeline
%
function pipeline_aipg_behavior(exptFolder,trialTypes)
%% Initialize
disp('STARTING ANALYSES FOR POOLED AIPG PURSUIT...')
close all  % Close all open figures

% Load processing settings
settings = processSettings();

% Replace spaces in the folder name with underscores
filebase = strrep(exptFolder, ' ', '_');

% Generate folder structure for saving outputs
folder = generateFolders(exptFolder);

% Change to the intermediate folder
cd(folder.int)

% Find all processed data files for the experiment
dataFiles = dir('*int.mat');

% Number of flies based on data files
nFlies = length(dataFiles);

% Number of visual objects
nObj = length(trialTypes);

% Minimum fixation time to include a fly (in seconds)
minFixationTime = settings.minFixationTime;

%% Load in and pool pursuit data
disp('Loading in and analyzing pursuit datasets...')
walkThresh = 1; %mm/s
arcRange = 0.5;

for e = 1:nFlies
    disp(['Processing fly ' num2str(e) '/' num2str(nFlies) '...'])
    % Load stimulation data
    cd(folder.int)
    thisTrial = dataFiles(e).name;
    load(thisTrial)

    % Find times when fly was walking
    walkIdx = schmittTrigger(int_forward(:,:,1), walkThresh, 0.1);
    walk_panelps = int_panelps(:,:,1);
    walk_panelps(~walkIdx) = nan;
    walk_forward = int_forward(:,:,1);
    walk_forward(~walkIdx) = nan;
    walk_probability(e,2) = (sum(~isnan(walk_panelps(:))) / numel(walk_panelps))*100;

    % Analyze velocity distribution
    [thisFwd, thisAng, ~] = velocity_histogram(int_forward(:,:,1), int_angular(:,:,1), int_sideway(:,:,1), 1);
    fwdHist(:,e) = thisFwd(:,2);
    angHist(:,e) = thisAng(:,2);

    % Analyze HD distributions (all and walking only)
    [posHist, ~] = panel_histogram(int_panelps(:,:,1), int_panelps(:,:,1), 1);
    HDHist(:,e) = posHist(:,2);
    [posHist, ~] = panel_histogram(walk_panelps, walk_panelps, 1);
    runHDHist(:,e) = posHist(:,2);
    cv(e) = circ_var(deg2rad(mod(walk_panelps(~isnan(walk_panelps)),360)));
    [~, width_deg(e), ~, ~] = centeredArcWidth(walk_panelps);

    % Analyze directional velocity vs object position
    [posvang_rl, ~, posBins] = setpoint_errorvturn(walk_panelps, int_angular(:,:,1), int_time, settings, 1, 0);
    [posvfwd, ~, posBins] = setpoint_errorvturn(walk_panelps, int_forward(:,:,1), int_time, settings, 1, 0);
    objVSang(:,e) = posvang_rl;
    objVSfwd(:,e) = posvfwd;

    % Load acclimitization data
    cd(folder.accl)
    load(strrep(thisTrial, '_int', '_acc'))
    % Restrict to last 5 min of trial
    restrict_accl = (length(int_accl_time)-length(int_accl_time)/4)+1:length(int_accl_time);
    int_accl_panelPs = int_accl_panelPs(restrict_accl);
    int_accl_forwardVelocity = int_accl_forwardVelocity(restrict_accl);

    % Find times when fly was walking
    walkIdx = schmittTrigger(int_accl_forwardVelocity, walkThresh, 0.1);
    walk_panelps_acc = int_accl_panelPs;
    walk_panelps_acc(~walkIdx) = nan;
    walk_forward_acc = int_accl_forwardVelocity;
    walk_forward_acc(~walkIdx) = nan;
    walk_probability(e,1) = (sum(~isnan(walk_panelps_acc(:))) / numel(walk_panelps_acc))*100;

    % Analyze HD distributions (all and walking only)
    [posHist, ~] = panel_histogram(int_accl_panelPs, int_accl_panelPs, 1);
    HDHist_acc(:,e) = posHist(:,2);
    [posHist, ~] = panel_histogram(walk_panelps_acc, walk_panelps_acc, 1);
    walkHDHist_acc(:,e) = posHist(:,2);
    cv_pre(e) = circ_var(deg2rad(mod(walk_panelps_acc(~isnan(walk_panelps_acc))',360)));
    [~, width_deg_pre(e), ~, ~] = centeredArcWidth(walk_panelps_acc');

    % Generate binarized heading
    heading_range = 35; %deg
    [this_acc_binary,this_exp_binary] = HD_binary(walk_panelps_acc,walk_panelps,heading_range);
    % Reshape and store
    acc_binary(:,e) = reshape(this_acc_binary,[],1);
    exp_binary(:,e) = reshape(this_exp_binary,[],1);
    % Fetch forward velocities when fly was oriented towards target
    fwdMean(e,1) = mean(walk_forward_acc(this_acc_binary),'omitnan');
    fwdMean(e,2) = mean(walk_forward(this_exp_binary),'omitnan');

    if contains(exptFolder,'P1')
        % Find fixation periods
        thisFixation = fixationFinder(int_panelps,int_forward,int_time,0);
        fix_panelps = thisFixation.panelps_run;
        fix_forward = thisFixation.forward_run;
        [cv_slow(e), cv_fast(e)] = cv_byfwd(fix_panelps, fix_forward);
        [arc_slow(e), arc_fast(e)] = arc_byfwd(fix_panelps, fix_forward);
    end
end

%% Plot basic parameters
% run_probability: rows = animals (e), col 1 = pre-stim, col 2 = during-stim
x = [1 2];

figure('Color','w'); set(gcf, 'Position', [100 100 1000 600])
tiledlayout(1,2,"TileSpacing","compact")

% Percentage of time spent walking
% Individual animals (gray dots connected by lines)
nexttile; hold on
for e = 1:nFlies
    plot(x, walk_probability(e,:), '.-', 'Color', [0.6 0.6 0.6], 'MarkerSize', 16, 'LineWidth', 1);
end

% Median for each condition as horizontal dash
meds = median(walk_probability,'omitnan');   % 1×2
dashHalf = 0.20;                                % half-length of the dash in x units
plot([1-dashHalf 1+dashHalf], [meds(1) meds(1)], 'k-', 'LineWidth', 3);
plot([2-dashHalf 2+dashHalf], [meds(2) meds(2)], 'r-', 'LineWidth', 3);

% Cosmetics
xlim([0.5 2.5]); ylim([0 100])
xticks(x); xticklabels({'Pre-stim','During-stim'});
ylabel(['% time spent walking >' num2str(walkThresh) 'mm/s']);
grid on; box off

% Plot average walking speeds when oriented
% Individual animals (gray dots connected by lines)
nexttile; hold on
for e = 1:nFlies
    plot(x, fwdMean(e,:), '.-', 'Color', [0.6 0.6 0.6], 'MarkerSize', 16, 'LineWidth', 1);
end

% Median for each condition as horizontal dash
meds = median(fwdMean,'omitnan');   % 1×2
dashHalf = 0.20;                                % half-length of the dash in x units
plot([1-dashHalf 1+dashHalf], [meds(1) meds(1)], 'k-', 'LineWidth', 3);
plot([2-dashHalf 2+dashHalf], [meds(2) meds(2)], 'r-', 'LineWidth', 3);

% Cosmetics
xlim([0.5 2.5]); ylim([0 18])
xticks(x); xticklabels({'Pre-stim','During-stim'});
ylabel('Avg forward speed while orienting (mm/s)');
grid on; box off

% Stats
[~,p,~,~] = ttest(fwdMean(:,1), fwdMean(:,2));
% Add p-value to top right corner
text(0.95,0.95,sprintf('p = %.4g',p), 'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');

% Save plot
cd(folder.summary)
plotname = 'summary_runprob';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.')

%% Plot HD distribution

% set whole column to NaN if it contains a 1
walkHDHist_acc(:, any(walkHDHist_acc==1,1)) = NaN;
runHDHist(:,    any(runHDHist==1,1))    = NaN;

binCenters_deg = posHist(:,1);             % Vector of bin centers in degrees
binWidth_deg = median(diff(binCenters_deg));
binEdges_deg = [binCenters_deg' - binWidth_deg/2, binCenters_deg(end) + binWidth_deg/2];
binEdges_rad = deg2rad(binEdges_deg);     % Convert to radians
binEdges_rad = binEdges_rad(:)';          % Ensure row vector

figure('Color','w'); set(gcf, 'Position', [100 100 1000 600])
tiledlayout(1,2)
rmax = 0.085;

% For pre-stimulation trials
% Across walking timepoints
avgHist = mean(walkHDHist_acc, 2, 'omitnan');
nexttile;
polarhistogram('BinEdges', binEdges_rad, ...
    'BinCounts', avgHist(:)', ...  % row vector
    'FaceColor', 'k', 'FaceAlpha', 0.6, ...
    'Normalization', 'probability');
rlim([0 rmax]);
title('Walking only')

% For stimulation trials
% Across walking timepoints
avgHist = mean(runHDHist, 2, 'omitnan');
nexttile;
polarhistogram('BinEdges', binEdges_rad, ...
    'BinCounts', avgHist(:)', ...  % row vector
    'FaceColor', 'r', 'FaceAlpha', 0.6, ...
    'Normalization', 'probability');
rlim([0 rmax]);
title('Walking only +aIPg')

% Save plot
cd(folder.summary)
plotname = 'summary_HDpolar';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.')

%% Plot circular variance
% Data: cv_pre(e) and cv(e), one value per fly
x = [1 2];
cv_data = [cv_pre(:) cv(:)];   % nFlies × 2


figure('Color','w'); set(gcf, 'Position', [100 100 800 500])
hold on

% Plot individual animals (gray lines with dots)
for e = 1:nFlies
    plot(x, cv_data(e,:), '.-', 'Color', [0.6 0.6 0.6], ...
         'MarkerSize', 16, 'LineWidth', 1);
end

% Medians for each condition as horizontal dashes
meds = median(cv_data,'omitnan');   % 1×2
dashHalf = 0.20;                    
plot([1-dashHalf 1+dashHalf], [meds(1) meds(1)], 'k-', 'LineWidth', 3);
plot([2-dashHalf 2+dashHalf], [meds(2) meds(2)], 'r-', 'LineWidth', 3);

% Cosmetics
xlim([0.5 2.5]); ylim([0 1]);   % circ variance ranges 0–1
xticks(x); xticklabels({'Pre-stim','During-stim'});
ylabel('Circular variance');
grid on; box off

% Paired t-test on circular variance
[~,p,~,~] = ttest(cv_pre, cv);

% Add p-value to top right corner
text(0.95, 0.95, sprintf('p = %.4g', p), ...
     'Units','normalized', ...
     'HorizontalAlignment','right', ...
     'VerticalAlignment','top');

% Save plot
cd(folder.summary)
plotname = 'summary_cv';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.')


%% Plot median absolute deviation from 0
x = [1 2];
arc_data = [width_deg_pre(:) width_deg(:)];   % nFlies × 2


figure('Color','w'); set(gcf, 'Position', [100 100 800 500])
hold on

% Plot individual animals (gray lines with dots)
for e = 1:nFlies
    plot(x, arc_data(e,:), '.-', 'Color', [0.6 0.6 0.6], ...
         'MarkerSize', 16, 'LineWidth', 1);
end

% Medians for each condition as horizontal dashes
meds = median(arc_data,'omitnan');   % 1×2
dashHalf = 0.20;                    
plot([1-dashHalf 1+dashHalf], [meds(1) meds(1)], 'k-', 'LineWidth', 3);
plot([2-dashHalf 2+dashHalf], [meds(2) meds(2)], 'r-', 'LineWidth', 3);

% Cosmetics
xlim([0.5 2.5]); ylim([0 360]);   % circ variance ranges 0–1
xticks(x); xticklabels({'Pre-stim','During-stim'});
ylabel('Half-containment width (°)');
grid on; box off

% Paired t-test on circular variance
[~,p,~,~] = ttest(width_deg_pre, width_deg);

% Add p-value to top right corner
text(0.95, 0.95, sprintf('p = %.4g', p), ...
     'Units','normalized', ...
     'HorizontalAlignment','right', ...
     'VerticalAlignment','top');

% Save plot
cd(folder.summary)
plotname = 'summary_arc';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.')

%% Plot HD raster

% Uniform sampling (minutes)
dt = median(diff(int_accl_time), 'omitnan')/60;
if ~isfinite(dt) || dt<=0
    dt = median(diff(int_time), 'omitnan')/60;
end
if ~isfinite(dt) || dt<=0
    error('Cannot infer sampling interval.');
end

% Time bases (minutes): pre ≤ 0, stim ≥ 0
Tpre  = size(acc_binary,1);
Tstim = size(exp_binary,1);
t_pre  = ( (0:Tpre-1)'  - (Tpre-1) ) * dt;   % [-Dpre, 0]
t_stim = (0:Tstim-1)' * dt;                  % [0, Dstim]

% Styling
colRaster = [0.1 0.1 0.1];
colLine   = 'r';
rowHalf   = 0.5; rowH = 2*rowHalf;

% Initialize
figure('Color','w'); set(gcf,'Position',[100 100 1200 520]); axes('NextPlot','add');

for e = 1:nFlies
    y0 = e - rowHalf;

    % Pre
    x = acc_binary(:,e);
    d = diff([0; x; 0]);
    st = find(d==1); en = find(d==-1)-1;
    for k = 1:numel(st)
        i1 = st(k); i2 = en(k);
        ts = t_pre(i1);
        te = t_pre(i2) + dt;
        if te > ts
            rectangle('Position',[ts, y0, te-ts, rowH], ...
                      'FaceColor', colRaster, 'EdgeColor','none');
        end
    end

    % Stim
    x = exp_binary(:,e);
    d = diff([0; x; 0]);
    st = find(d==1); en = find(d==-1)-1;
    for k = 1:numel(st)
        i1 = st(k); i2 = en(k);
        ts = t_stim(i1);
        te = t_stim(i2) + dt;
        if te > ts
            rectangle('Position',[ts, y0, te-ts, rowH], ...
                      'FaceColor', colRaster, 'EdgeColor','none');
        end
    end
end

xline(0, '-', 'Color', colLine, 'LineWidth', 1.5);

ylim([0.5, nFlies+0.5]); yticks(1:nFlies); ylabel('Fly');
xlim([min(t_pre), max(t_stim)+dt]);
xlabel('Time (min)');
title(['Orientation raster - walking >' num2str(walkThresh) 'mm/s and within +/-' num2str(heading_range) 'deg']);
box on

% Save plot
cd(folder.summary)
plotname = 'summary_headingRaster';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.')


%% Plot circular variance by forward (for P1 pursuit only)
if contains(exptFolder,'P1')
    % Data: cv_pre(e) and cv(e), one value per fly
    x = [1 2];
    arc_data = [arc_slow(:) arc_fast(:)];   % nFlies × 2

    figure('Color','w'); set(gcf, 'Position', [100 100 800 500])
    hold on

    % Plot individual animals (gray lines with dots)
    for e = 1:nFlies
        plot(x, arc_data(e,:), '.-', 'Color', [0.6 0.6 0.6], ...
            'MarkerSize', 16, 'LineWidth', 1);
    end

    % Medians for each condition as horizontal dashes
    meds = median(arc_data,'omitnan');   % 1×2
    dashHalf = 0.20;
    plot([1-dashHalf 1+dashHalf], [meds(1) meds(1)], 'k-', 'LineWidth', 3);
    plot([2-dashHalf 2+dashHalf], [meds(2) meds(2)], 'r-', 'LineWidth', 3);

    % Cosmetics
    xlim([0.5 2.5]); ylim([0 180]);   % circ variance ranges 0–1
    xticks(x); xticklabels({'Low fwd','High fwd'});
    ylabel('50% arc width');
    grid on; box off

    % Paired t-test on circular variance
    [~,p,~,~] = ttest(arc_slow, arc_fast);

    % Add p-value to top right corner
    text(0.95, 0.95, sprintf('p = %.4g', p), ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top');

    % Save plot
    cd(folder.summary)
    plotname = 'summary_arcfwd';
    saveas(gcf, [plotname '.png']);
    copyfile([plotname '.png'], folder.dropbox, 'f');

    % Save vectorized plot
    cd(folder.vector)
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, [plotname '.svg'])
    copyfile([plotname '.svg'], folder.dropbox, 'f');
    disp('Complete.')
end
%% end
disp('ALL ANALYSES COMPLETE.')
end

