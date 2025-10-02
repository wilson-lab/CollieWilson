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
function pipeline_battery_behavior(exptFolder,trialTypes)
%% Initialize
disp('STARTING ANALYSES FOR POOLED OBJECT BATTERY PURSUIT...')
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

for e = 1:nFlies
    disp(['Processing fly ' num2str(e) '/' num2str(nFlies) '...'])
    % Load current trial
    cd(folder.int)
    thisTrial = dataFiles(e).name;
    load(thisTrial)
    etime = int_time; etime(end+1) = 60;  % Append trial end time

    % Find fixation periods and bin data
    thisFixation = fixationFinder(int_panelps,int_forward,int_time,0);
    fix_panelps = thisFixation.panelps_run;
    fix_panelvel = int_panelvel;
    fix_forward = thisFixation.forward_run;
    fix_angular = int_angular;
    fix_angular(~thisFixation.idx_run) = nan;
    fix_sideway = int_sideway;
    fix_sideway(~thisFixation.idx_run) = nan;

    fixationFwdVels(e,:) = reshape(mean(fix_forward,[1 2],'omitnan'),1,nObj);  % Avg forward speed (mm/s)
    fixationAngVels(e,:) = reshape(mean(abs(fix_angular),[1 2],'omitnan'),1,nObj);  % Avg angular speed (deg/s)

    % Sample period (assumes roughly uniform sampling)
    dt = median(diff(int_time));   % seconds per sample

    fixationTimes = nan(1,nObj);
    totalTimes    = nan(1,nObj);

    for o = 1:nObj
        % Masks over time × trials for this object
        mask_total = ~isnan(int_angular(:,:,o));
        mask_fix   = ~isnan(fix_angular(:,:,o));   % you set non-running to NaN above

        % Sum valid samples and convert to seconds
        totalTimes(1,o)    = dt * sum(mask_total(:));
        fixationTimes(1,o) = dt * sum(mask_fix(:));
    end

    fixationRatios(e,:) = (fixationTimes ./ totalTimes) * 100;

    % Analyze HD distributions
    for o = 1:nObj
        % Find times when fly was running
        runIdx = schmittTrigger(int_forward(:,:,o), settings.runThreshB, 0.1);
        run_panelps = int_panelps(:,:,o);
        run_panelps(~runIdx) = nan;

        [posHist, ~] = panel_histogram(run_panelps, fix_panelvel(:,:,o), 1);
        fixationHDHist(:,o,e) = posHist(:,2);
    end

    [posvang, posvang_rl, posBins] = setpoint_errorvturn(fix_panelps, fix_angular, int_time, settings, 1, 0);
    fixationEvO(:,:,e) = posvang;
    fixationEvO_rl(:,:,e) = posvang_rl;
end

%% Plot basic pursuit parameters
figure; set(gcf, 'Position', [100 100 800 500])
tiledlayout(1,3,"TileSpacing","compact");
xvals = 1:numel(trialTypes);

for p = 1:3
    switch p
        case 1
            thisData = fixationRatios;
            nameData = '% time spent pursuing';
            yrange = [0 100];
        case 2
            thisData = fixationFwdVels;
            nameData = 'Forward Speed (mm/s)';
            yrange = [0 20];
        case 3
            thisData = fixationAngVels;
            nameData = 'Angular Speed (deg/s)';
            yrange = [0 80];
    end

    nexttile;
    hold on;

    % Scatter each fly's data (columns = object conditions)
    plot(xvals, thisData', '.-', 'Color', settings.trialColor);

    % Plot median across flies
    plot(xvals, median(thisData,1,'omitnan'), '_', 'Color', [0 0 0], 'MarkerSize', 10);

    ylabel(nameData);
    ylim(yrange); xlim([0 nObj+1])
    xticks(xvals); xticklabels(trialTypes);
    if p == 3
        xlabel('Object');
    end

    % ===== Linear Mixed-Effects Model (Condition fixed, Fly random) =====
    valid = all(isfinite(thisData),2);
    pval = NaN;
    if sum(valid) >= 2 && size(thisData,2) >= 2
        nFlies = sum(valid);
        nCond  = size(thisData,2);

        % Build long-format table: one row per (Fly × Condition)
        [condID, flyID] = ndgrid(1:nCond, 1:nFlies);
        tbl = table;
        tbl.Fly       = categorical(flyID(:));
        tbl.Condition = categorical(condID(:), 1:nCond, trialTypes(1:nCond));
        vals          = thisData(valid,:)';   % transpose so conds in rows
        tbl.Value     = vals(:);

        % Fit LME (random intercept for Fly)
        lme = fitlme(tbl, 'Value ~ Condition + (1|Fly)');

        % ANOVA table for fixed effect of Condition
        stats = anova(lme);
        pval  = stats.pValue(2);   % row 2 = Condition effect

        % Add p-value to plot
        text(0.98, 0.02, sprintf('p = %.5g', pval), ...
            'Units','normalized','HorizontalAlignment','right', ...
            'VerticalAlignment','bottom','FontSize',8,'Color',[0 0 0]);

% --------- Post-hoc comparisons (manual Tukey–Kramer from LME) ----------
if pval < 0.05 && nCond >= 2
    disp(['--- Tukey–Kramer post-hoc for ' nameData ' (manual) ---'])

    % Pull fixed-effect estimates and covariance from the LME
    beta      = lme.Coefficients.Estimate;           % column vector
    betaNames = string(lme.Coefficients.Name);        % e.g. "(Intercept)","Condition_<level>"
    covB      = lme.CoefficientCovariance;           % covariance matrix of fixed effects
    df        = lme.DFE;                              % error DF

    % Condition level names in model order
    lvl = string(categories(tbl.Condition));          % matches how you built tbl
    k   = numel(lvl);

    % All pairwise contrasts μ_a - μ_b
    pairs = nchoosek(1:k, 2);
    out   = table('Size',[size(pairs,1) 7], ...
                  'VariableTypes', {'string','string','double','double','double','double','double'}, ...
                  'VariableNames', {'Cond1','Cond2','Diff','SE','t','q','p_Tukey'});

    % Helper to find coefficient index for a given level name
    coefIdx = @(lev) find(betaNames == "Condition_" + lev, 1);

    for r = 1:size(pairs,1)
        a = pairs(r,1); b = pairs(r,2);
        c = zeros(numel(beta),1);

        if a==1 && b>=2
            % μ1 - μb = -β_b
            idxb = coefIdx(lvl(b));
            c(idxb) = -1;
        elseif b==1 && a>=2
            % μa - μ1 = +β_a
            idxa = coefIdx(lvl(a));
            c(idxa) = +1;
        else
            % μa - μb = β_a - β_b
            idxa = coefIdx(lvl(a));
            idxb = coefIdx(lvl(b));
            c(idxa) = +1;  c(idxb) = -1;
        end

        diff_ab = c' * beta;
        se_ab   = sqrt(max(0, c' * covB * c));
        t_ab    = diff_ab / se_ab;

        % Tukey–Kramer p via studentized range (q = sqrt(2)*|t|)
        q_ab = sqrt(2) * abs(t_ab);
        if exist('studrangecdf','file') == 2
            p_tukey = 1 - studrangecdf(q_ab, k, df);
        else
            % Fallback: Sidák-adjusted p from t (conservative without studrangecdf)
            p_raw   = 2 * (1 - tcdf(abs(t_ab), df));
            m       = size(pairs,1);
            p_tukey = 1 - (1 - p_raw).^m;   % Sidák familywise
        end

        out.Cond1(r)  = lvl(a);
        out.Cond2(r)  = lvl(b);
        out.Diff(r)   = diff_ab;
        out.SE(r)     = se_ab;
        out.t(r)      = t_ab;
        out.q(r)      = q_ab;
        out.p_Tukey(r)= p_tukey;
    end
    % Save posthoc stats to Excel
    fname = fullfile(folder.summary, ...
        sprintf('posthoc_%s.xlsx', regexprep(nameData,'\W','_')));
    writetable(out, fname);
    disp(['Posthoc results saved to ' fname])

    disp(out)
end

    end

end

% Save plot
cd(folder.summary)
plotname = 'summary_basics';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector)
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox, 'f');
disp('Complete.')


%% Plot HD distribution

binCenters_deg = posHist(:,1);             % Vector of bin centers in degrees
binWidth_deg = median(diff(binCenters_deg));
binEdges_deg = [binCenters_deg' - binWidth_deg/2, binCenters_deg(end) + binWidth_deg/2];
binEdges_rad = deg2rad(binEdges_deg);     % Convert to radians
binEdges_rad = binEdges_rad(:)';          % Ensure row vector

figure; set(gcf, 'Position', [100 100 1500 500])
tiledlayout(1,nObj,"TileSpacing","compact");

for o = 1:nObj
    % Average across flies
    avgHist = mean(fixationHDHist(:,o,:), 3, 'omitnan');

    nexttile;
    polarhistogram('BinEdges', binEdges_rad, ...
                   'BinCounts', avgHist(:)', ...  % row vector
                   'FaceColor', 'k', 'FaceAlpha', 0.6, ...
                   'Normalization', 'probability');

    rlim([0 0.3]);
    title(trialTypes{o});
end

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


%% Plot error vs turn response

condColors = lines(nObj);  % default colormap
figure; set(gcf, 'Position', [100 100 300 600]);
hold on;

for o = 1:nObj
    % Extract data for object condition o across all flies
    thisData = fixationEvO(:, o, :);  % (posBins × 1 × nFlies)
    thisData = squeeze(thisData);     % (posBins × nFlies)

    % Compute mean and SEM across flies
    meanEv = mean(thisData, 2, 'omitnan');
    semEv = std(thisData, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(thisData), 2));

    % Only plot positions with sufficient data
    validPts = any(~isnan(thisData), 2);

    % Plot mean
    hLine = plot(posBins(validPts), meanEv(validPts), ...
                          'Color', condColors(o,:), 'LineWidth', settings.lwAvg);
    lineHandles(o) = hLine(1);  % ensure only one handle is stored

    % Plot SEM patch (excluded from legend)
    semPatch = patch([posBins(validPts)'; flipud(posBins(validPts)')], ...
                     [meanEv(validPts) - semEv(validPts); flipud(meanEv(validPts) + semEv(validPts))], ...
                     condColors(o,:), 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    uistack(semPatch, 'bottom');  % push patch behind line
end

xline(0); yline(0);
xlabel('Object Error (deg)');
ylabel('Angular Velocity (deg/s)');
xlim([-60 60]);
ylim([-200 200]);

% Use only line handles for legend
legend(lineHandles, trialTypes, 'Location', 'southeast');

% Save combined plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_combined';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox,'f');

% Save zoomed combined plot
ylim([-130 130]); xlim([-40 40])
% Save zoomed plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_combined_20x';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');
disp('Complete.')

%% Fit slope
% Set rangeValue for plotting (e.g., -20 to 20)
rangeValue = 20;

% Fit slope
[allSlopes, ~, ~] = fitTurnVelocity_single(posBins, fixationEvO, rangeValue);

% Convert slopes cell array -> matrix (flies × conditions), padding with NaN if needed
maxFlies = max(cellfun(@numel, allSlopes));
nCond    = numel(allSlopes);
thisData = nan(maxFlies, nCond);

for cond = 1:nCond
    nFlies = numel(allSlopes{cond});
    thisData(1:nFlies, cond) = allSlopes{cond};
end

% Plot
figure; set(gcf, 'Position', [100 100 300 600]); hold on;
xvals = 1:nCond;
nameData = 'Slope of Turn vs. Object Position';
yrange = [0 5];   % adjust based on your data

% Scatter each fly’s slope across conditions
plot(xvals, thisData', '.-', 'Color', settings.trialColor);

% Plot median slope across flies
plot(xvals, median(thisData,1,'omitnan'), '_', 'Color', [0 0 0], 'MarkerSize', 10);

ylabel(nameData);
ylim(yrange); xlim([0 nCond+1])
xticks(xvals); xticklabels(trialTypes);
xlabel('Condition')

% ===== Linear Mixed-Effects Model (Condition fixed, Fly random) =====
valid = all(isfinite(thisData),2);
pval = NaN;
if sum(valid) >= 2 && nCond >= 2
    nFlies = sum(valid);

    % Long-format table
    [flyID, condID] = ndgrid(1:nFlies, 1:nCond);
    tbl = table;
    tbl.Fly       = categorical(flyID(:));
    tbl.Condition = categorical(condID(:), 1:nCond, trialTypes(1:nCond));
    tbl.Value     = reshape(thisData(valid,:)', [], 1);

    % LME: slope ~ condition + (1|fly)
    lme = fitlme(tbl, 'Value ~ Condition + (1|Fly)');

    % ANOVA table for Condition
    stats = anova(lme, 'DFMethod','Satterthwaite');
    condRow = find(strcmp(stats.Term, 'Condition'));
    if ~isempty(condRow), pval = stats.pValue(condRow); end
end

% Annotate p-value
if ~isnan(pval)
    text(0.98, 0.02, sprintf('p = %.4g', pval), ...
        'Units','normalized','HorizontalAlignment','right', ...
        'VerticalAlignment','bottom','FontSize',8,'Color',[0 0 0]);
end
hold off

% --- Tukey–Kramer post-hoc via randomized-block ANOVA (Fly as random) ---
if ~isnan(pval) && pval < 0.05
    disp(['--- Tukey–Kramer post-hoc for ' nameData ' ---'])
    [~,~,statsA] = anovan(tbl.Value, {tbl.Condition, tbl.Fly}, ...
        'random', 2, ...
        'model', 'linear', ...
        'display', 'off', ...
        'varnames', {'Condition','Fly'});

    [c,~,~,gnames] = multcompare(statsA, ...
        'ctype','tukey-kramer', ...
        'dimension', 1, ...
        'display','off');

    posthocTbl = table( ...
        gnames(c(:,1)), gnames(c(:,2)), ...
        c(:,3), c(:,4), c(:,5), c(:,6), ...
        'VariableNames', {'Cond1','Cond2','LowerCI','Diff','UpperCI','pValue'});
    disp(posthocTbl)
end

% Save zoomed plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_fit';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');
disp('Complete.')

%% Plot error vs turn response (r/l combined)

condColors = lines(nObj);  % default colormap
figure; set(gcf, 'Position', [100 100 300 600]);
hold on;

for o = 1:nObj
    % Extract data for object condition o across all flies
    thisData = fixationEvO_rl(:, o, :);  % (posBins × 1 × nFlies)
    thisData = squeeze(thisData);     % (posBins × nFlies)

    % Compute mean and SEM across flies
    meanEv = mean(thisData, 2, 'omitnan');
    semEv = std(thisData, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(thisData), 2));

    % Only plot positions with sufficient data
    validPts = any(~isnan(thisData), 2);

    % Plot mean
    hLine = plot(posBins(validPts), meanEv(validPts), ...
                          'Color', condColors(o,:), 'LineWidth', settings.lwAvg);
    lineHandles(o) = hLine(1);  % ensure only one handle is stored

    % Plot SEM patch (excluded from legend)
    semPatch = patch([posBins(validPts)'; flipud(posBins(validPts)')], ...
                     [meanEv(validPts) - semEv(validPts); flipud(meanEv(validPts) + semEv(validPts))], ...
                     condColors(o,:), 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none');
    uistack(semPatch, 'bottom');  % push patch behind line
end

xline(0); yline(0);
xlabel('Object Error (deg)');
ylabel('Angular Velocity (deg/s)');
xlim([-60 60]);
ylim([-200 200]);

% Use only line handles for legend
legend(lineHandles, trialTypes, 'Location', 'southeast');

% Save combined plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_combined_rl';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox,'f');

% Save zoomed combined plot
ylim([-130 130]); xlim([-40 40])
% Save zoomed plot as PNG and SVG
cd(folder.summary)
plotname = 'error_v_turn_combined_rl_20x';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');
cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], folder.dropbox,'f');
disp('Complete.')

%% end
disp('ALL ANALYSES COMPLETE.')
end

