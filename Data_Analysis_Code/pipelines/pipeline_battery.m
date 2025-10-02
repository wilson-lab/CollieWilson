% pipeline_battery
%
% Pipeline Function
% Pulls all processed files from ALL flies in a given oscillatory experiment
% and performs necessary analyses and plots accordingly. Can be used for both
% behavior-only experiments and ephys experiments.
%
% INPUTS
% exptFolder - overarching experiment folder
% trialTypes - list of trial conditions (e.g., different target shapes/sizes)
%
% The function pools data across all flies, performs analysis on spike rate,
% target motion, cross-correlations, pursuit performance, and tracking indices.
% Results are saved as plots, including raster, summary figures, and comparisons
% across trial types.
%
% 07/07/2025 - MC adapted from oscillatory pipeline
%
function pipeline_battery(exptFolder,trialTypes)
%% initialize
disp('STARTING ANALYSES FOR BATTERY EXPERIMENT...')
close all

% load in processing settings
settings = processSettings();

% set filename info and create necessary directories
filebase = strrep(exptFolder,' ','_');
% generate folder structures as needed
folder = generateFolders(exptFolder);

% number of target conditions used by this battery
nCond = length(trialTypes);

%% set plotting variables
switch exptFolder
    case 'AOTU019 Battery'
        sr_limit = [0 100];
    case 'AOTU025 Battery'
        sr_limit = [0 100];
end
ps_limit = [-150 150];
fwd_limit = [0 12];
ang_limit = [-250 250];
sid_limit = [-7 7];

cm = [0.7 0.7 0.7 ; flip(colormap(jet(nCond-1)))]; %color map
close all

%% load in and pool all trials from each experiment folder
disp('Loading in datasets...')
% pull file info
cd(folder.int)
allFiles = dir('*int.mat');
nFlies = length(allFiles);

% initialize data storage arrays
nFliesThresh = nFlies; %number of flies used in behavior analysis
nt_t = 0; %counter
normSweep = [];

for nt = 1:nFlies
    % load this trial
    disp(['Loading fly ' num2str(nt) '/' num2str(nFlies)])
    thisTrial = allFiles(nt).name;
    thisFly = thisTrial(6:16);
    cd(folder.int); load(thisTrial)

    % determine if this fly exhibitted sufficient walking behavior
    flyRunTime(nt,1) = (sum(int_forward>settings.runThreshE,'all')/length(int_time))*60;

    % for each trial type
    for s = 1:nCond
        % pull this fly data for ONLY this trial type
        thisForward = int_forward(:,:,s);
        thisSideway = int_sideway(:,:,s);
        thisAngular = int_angular(:,:,s);
        thisPanelPs = int_panelps(:,:,s);
        thisSpikeRt = int_spikert(:,:,s);
        thisTime = int_time;

        % analyze visual tuning for this trial type
        [thisSweep,thisMean_all] = osc_v_output(thisPanelPs,thisForward,thisSpikeRt,-1,thisTime,settings); %all
        % fetch size of arrays
        checkSweep(1) = size(thisMean_all,1);
        thisLen = size(thisMean_all,1);

        if nt == 1
            normSweep(s) = thisLen;
            sweepPos{s} = thisSweep; % store sweep positions
        end

        % Ensure consistent length: trim or pad if necessary
        if thisLen >= normSweep(s)
            thisMean_all = thisMean_all(1:normSweep(s), :); % trim to expected size
        else
            padAmount = normSweep(s) - thisLen;
            thisMean_all = [thisMean_all; nan(padAmount, size(thisMean_all,2))]; % pad with NaNs
        end

        % store visual tuning analysis from this trial type
        srvpos_all{nt,s} = thisMean_all;

        % only include flies that showed sufficient running in behavior analyses
        if flyRunTime(nt,1)>settings.minRunTime
            if s==1
                nt_t = nt_t+1; %update counter
            end
            % analyze pursuit performance for this trial type
            [fidelity,vigor,w_sr,w_time] = pursuit_performance(thisPanelPs,thisAngular,thisSpikeRt,int_time);
            allFidelity{nt_t,s} = fidelity;
            allVigor{nt_t,s} = vigor;
            allWinSR{nt_t,s} = w_sr;

        else
            if s==1
                disp([thisFly ' omitted from behavior analyses.'])
                nFliesThresh = nFliesThresh-1;
            end
        end
    end
end

%% By trial: plot spike rate vs target motion (ALL behavior bin only, using patch for SEM)
disp('By trial type: analyzing spike rate vs target motion (mean ± SEM)...')

% create one figure for all conditions (subplots = 1 per condition)
figure; set(gcf,'Position',[100 100 1500 900])
tiledlayout(2,nCond/2,"TileSpacing","compact")

% initialize panel range
panel_pos = sweepPos{1};
idx_min = find(panel_pos == min(panel_pos), 1, 'first');
idx_max = find(panel_pos == max(panel_pos), 1, 'last');
idx_range = idx_min:idx_max;
idx_mid = round(mean(idx_range));

% initialize storage arrays
store_maxFR = [];
store_peakpos = [];
store_meanFR = [];

% for each condition
for s = 1:nCond
    % only analyze 'all' trial set
    thisTrialset = cat(2,srvpos_all{:,s}); % time x trials
    thisSweepMean = mean(thisTrialset,2,'omitnan');
    thisSweepSEM = std(thisTrialset,0,2,'omitnan') ./ sqrt(nFlies);
    thisPanelps = panel_pos(1:length(thisSweepMean));
    % initialize time vector
    t = int_time(1:length(thisSweepMean));
    cT = mean(t);

    % store max firing rate during ipsi sweep
    [maxFR, idx_peak] = max(thisTrialset(idx_mid:idx_max,:));
    store_maxFR(s,:) = maxFR;
    store_peakpos(s,:) = thisPanelps(idx_peak+idx_min);
    % store mean firing rate during ipsi sweep
    store_meanFR(s,:) = mean(thisTrialset(idx_mid:idx_max,:));

    % plot
    nexttile
    fs = gca; yyaxis right
    plot(t, thisPanelps, 'k', 'LineWidth', 2)
    fs.YAxis(1).Color = 'k'; fs.YAxis(2).Visible = 'off';
    yline(0,'Color','k'); xline(cT,'Color','k');

    yyaxis left
    % Plot SEM as patch
    xvals = t(:);
    y_upper = thisSweepMean + thisSweepSEM;
    y_lower = thisSweepMean - thisSweepSEM;
    patch([xvals; flipud(xvals)], [y_lower; flipud(y_upper)], 'k', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none', 'FaceColor', 'k'); hold on

    % Plot mean line
    plot(t, thisSweepMean, 'Color', settings.spkColor, 'LineWidth', 1.5)
    axis tight; ylim(sr_limit);
    ylabel(['All ' settings.spkLabel]); xlabel('Time (s)')
    title(trialTypes{s}, 'Interpreter', 'none')
end
sgtitle([strrep(filebase,'_',' ') ' ALL trials'])

% save plots
cd(folder.summary)
plotname = 'srvpos_all';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');

cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox,'f');


% zoom in to ipsiversive sweep only
for s = 1:nCond
    nexttile(s)
    xlim(sort([t(idx_min), t(idx_max)]))
end
% save plots
cd(folder.summary)
plotname = 'srvpos_all_ipsionly';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');

cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox,'f');


%% generate summary plot
figure; set(gcf,'Position',[100 100 1500 500])
tiledlayout(1,4,"TileSpacing","compact")

% Plot 1: Max firing rate per fly (connected points + median)
nexttile; hold on

% Each fly across conditions, skip last
for f = 1:nFlies
    plot(1:nCond-1, store_meanFR(1:nCond-1,f), '.-', 'Color', [0.8 0.8 0.8]);
end

% Median dash per condition, skip last
for s = 1:nCond-1
    yvals = store_meanFR(s,:);
    ymed  = median(yvals,'omitnan');
    plot([s-0.2 s+0.2], [ymed ymed], '-', 'LineWidth', 2, 'Color', cm(s,:));
end

xlim([0 nCond+1]); ylim([0 120])
xticks(1:nCond); xticklabels(trialTypes)
ylabel(['Mean ' settings.spkLabel])
title('Mean firing rate per fly')


% Plot 2 (middle): Mean ± SEM of max firing rate across flies
nexttile; hold on
meanFR = mean(store_meanFR, 2, 'omitnan');         % mean across flies
semFR = std(store_meanFR, 0, 2, 'omitnan') ./ sqrt(nFlies);  % SEM
% Plot bars
b = bar(1:nCond, meanFR, 'FaceColor', 'flat'); 
for s = 1:nCond
    b.CData(s,:) = cm(s,:); % set each bar color from colormap
end

% Add errorbars (line only, no caps, no markers)
er = errorbar(1:nCond, meanFR, semFR, ...
    'k', 'LineStyle', 'none', 'CapSize', 0, 'LineWidth', 1);

xlim([0 nCond+1]); ylim([0 200])
xticks(1:nCond); xticklabels(trialTypes)
ylabel(['Max ' settings.spkLabel])
title('Mean ± SEM firing rate')

% Plot 3: SEM of max FR per condition (dash line)
nexttile; hold on
% Reuse semFR from earlier, or recompute
semFR = std(store_maxFR, 0, 2, 'omitnan') ./ sqrt(nFlies);
% Colored markers for each condition
for s = 1:nCond
    plot([s-0.2, s+0.2], [semFR(s), semFR(s)], '-', ...
        'LineWidth', 2, 'Color', cm(s,:))
end
xlim([0 nCond+1])
ylim([0 20])
xticks(1:nCond); xticklabels(trialTypes)
ylabel(['SEM of ' settings.spkLabel])
title('SEM of max firing rate per condition')

% Plot 4: Panel position at max FR
nexttile; hold on
jitter = 0;
for s = 1:nCond
    tt_j = s + (-jitter + 2*jitter*rand(1, nFlies));
    yvals = store_peakpos(s,:);

    plot(tt_j, yvals, '.', 'Color', [0.8 0.8 0.8]);

    % Median dash
    ymed = median(yvals, 'omitnan');
    plot([s-0.2 s+0.2], [ymed ymed], '-', 'LineWidth', 2, 'Color', cm(s,:));
end
xlim([0 nCond+1])
xticks(1:nCond); xticklabels(trialTypes)
ylabel('Panel pos @ max FR')
title('Preferred position per fly')

% Save plots
cd(folder.summary)
plotname = 'srvpos_all_summary';
saveas(gcf, [plotname '.png']);
copyfile([plotname '.png'], folder.dropbox,'f');

cd(folder.vector)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg']);
copyfile([plotname '.svg'], folder.dropbox,'f');

disp('Complete.')


%% run repeated measures ANOVA
% Prepare data
data = store_maxFR(1:end-1, :)';  % [nFlies x (nCond - 1)]
nFlies = size(data, 1);
nCondUsed = size(data, 2);
condLabels = trialTypes(1:end-1);

% Stack into long format
flyID = repmat((1:nFlies)', nCondUsed, 1);
condition = repmat(condLabels(:)', nFlies, 1);
condition = condition(:);  % long vector
response = data(:);        % long vector

% Build table
data_tbl = table(flyID, categorical(condition), response, ...
    'VariableNames', {'Fly', 'Condition', 'MaxFR'});

% Fit LME model
lme = fitlme(data_tbl, 'MaxFR ~ Condition + (1|Fly)');

% ANOVA with Satterthwaite approximation
anova_tbl = anova(lme, 'DFMethod', 'satterthwaite');

% Convert if needed and save ANOVA
if ~istable(anova_tbl)
    anova_tbl = dataset2table(anova_tbl);
end

% Save ANOVA
writetable(anova_tbl, fullfile(folder.summary, 'anova_maxFR.csv'));

% === POST-HOC pairwise comparison between conditions ===
levels = categories(data_tbl.Condition);
n_comparisons = nchoosek(numel(levels), 2);

posthoc_array = perform_bonferroni_posthoc(lme, levels, n_comparisons, 'Condition');

% Convert to table and save
posthoc_tbl = cell2table(posthoc_array, 'VariableNames', ...
    {'Condition1', 'Condition2', 'Difference', 'SE', 'TStatistic', 'BonferroniPValue'});

writetable(posthoc_tbl, fullfile(folder.summary, 'posthoc_maxFR.csv'));

disp('LME ANOVA and post-hoc results saved.');

%% end
disp('ALL ANALYSES COMPLETE.')
end

