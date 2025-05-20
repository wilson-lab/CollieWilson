% analyze_jump_performance
%
% Analyzes jump performance for multiple genotypes and gains, focusing on the 
% time to reach object positions within specified thresholds. The function 
% generates plots for each genotype across gain conditions and calculates 
% statistical metrics.
%
% INPUTS:
% - kirJumpObj, wtJumpObj, naJumpObj : Cell arrays of object positions for 
%                                        each genotype.
% - kirJumpZeroT, wtJumpZeroT, naJumpZeroT : Cell arrays of times to zero 
%                                              for each genotype.
% - settings : Structure containing plotting settings (e.g., geneColor, 
%               pursuitGain, lwAvg, semAlpha, trialColor).
% - folder : Directory where results or plots should be saved.
%
% OUTPUTS:
% - kirZeroTAboveThresh, wtZeroTAboveThresh, naZeroTAboveThresh : Zero-time 
%   points between the minthreshold and maxthreshold for each genotype.
%
% PROCESS:
% The function initializes thresholds for object positions and retrieves gain 
% conditions. It then processes and analyzes jump data for each genotype, 
% creating scatter plots for object positions and the corresponding zero times 
% while highlighting values within specified thresholds. After processing, 
% kernel density estimates are plotted for the zero times across all genotypes 
% for each gain condition. Finally, it performs ANOVA on the results and 
% displays p-values for genotype and gain factors.
%
function [kirZeroTAboveThresh, wtZeroTAboveThresh, naZeroTAboveThresh] = analyze_jump_performance(kirJumpObj, wtJumpObj, naJumpObj, kirJumpZeroT, wtJumpZeroT, naJumpZeroT, settings, folder)
%% Initialize

% Set min and max thresholds
minthreshold = 50;
maxthreshold = 100;
% Fetch gain
nGain = size(kirJumpObj, 1);

%% Fetch data and analyze
% Initialize figure
figure; 
set(gcf, 'Position', [100 100 1200 900]); % Adjust figure size
tiledlayout(4, nGain, 'TileSpacing', 'compact'); % 3 rows for genotypes, nGain columns for each gain

% Initialize data storage arrays
kirAvgT = nan(1, nGain);
wtAvgT = nan(1, nGain);
naAvgT = nan(1, nGain);
kirZeroT = cell(1, nGain);
wtZeroT = cell(1, nGain);
naZeroT = cell(1, nGain);
kirObjPos = cell(1, nGain);
wtObjPos = cell(1, nGain);
naObjPos = cell(1, nGain);

% Loop over each genotype
for g = 1:3
    % Fetch data based on genotype
    switch g
        case 1
            jumpObj = kirJumpObj;
            jumpZeroT = kirJumpZeroT;
            avgT = kirAvgT;
            zeroT = kirZeroT;
            objPos = kirObjPos;
        case 2
            jumpObj = wtJumpObj;
            jumpZeroT = wtJumpZeroT;
            avgT = wtAvgT;
            zeroT = wtZeroT;
            objPos = wtObjPos;
        case 3
            jumpObj = naJumpObj;
            jumpZeroT = naJumpZeroT;
            avgT = naAvgT;
            zeroT = naZeroT;
            objPos = naObjPos;
    end

    % Loop over each gain condition
    for c = 1:nGain
        nexttile; hold on % Create a new tile for each gain within the genotype

        % Collect data for the current gain condition
        allObj = [];
        allZeroT = [];
        for i = 1:size(jumpObj, 2)
            allObj = [allObj, jumpObj{c, i}];
            allZeroT = [allZeroT, jumpZeroT{c, i}];
        end

        % Remove NaN values
        valid_idx = ~isnan(allObj) & ~isnan(allZeroT);
        allObj = allObj(valid_idx);
        allZeroT = allZeroT(valid_idx);

        % Filter data based on minthreshold and maxthreshold
        within_thresh_idx = allObj > minthreshold & allObj <= maxthreshold;

        % Scatter points outside the threshold in settings.trialColor
        scatter(allObj(~within_thresh_idx), allZeroT(~within_thresh_idx), '.', 'MarkerEdgeColor', settings.trialColor);

        % Scatter points within the threshold in black
        scatter(allObj(within_thresh_idx), allZeroT(within_thresh_idx), '.', 'k');

        % Store zeroT and object position points between minthreshold and maxthreshold
        zeroT{c} = allZeroT(within_thresh_idx);
        objPos{c} = allObj(within_thresh_idx);

        % Calculate the average correction time for object positions within the thresholds
        if sum(within_thresh_idx) > 0
            avgT(c) = mean(allZeroT(within_thresh_idx));
        else
            avgT(c) = NaN; % If no data within thresholds, set to NaN
        end

        % Add the horizontal line (yline) from 50 to 150 in the color of the genotype
        if ~isnan(avgT(c))
            plot([minthreshold maxthreshold], [avgT(c) avgT(c)], 'Color', settings.geneColor{g}, 'LineWidth', 1.5);
        end

        % Set plot limits
        xlim([0 200]); ylim([0 5]);

        % Set labels and titles
        if c == 1
            ylabel({'Time to Object', 'Cross Midline(s)'});
        end
        if g == 1
            title([num2str(settings.pursuitGain(c)) 'X']);
        end
    end
    
    % Store the average correction times, zeroT, and objPos for each genotype
    switch g
        case 1
            kirAvgT = avgT;
            kirZeroT = zeroT;
            kirObjPos = objPos;
        case 2
            wtAvgT = avgT;
            wtZeroT = zeroT;
            wtObjPos = objPos;
        case 3
            naAvgT = avgT;
            naZeroT = zeroT;
            naObjPos = objPos;
    end
end

% Plot kernel density estimates for all genotypes together for each gain
for c = 1:nGain
    nexttile; hold on
    
    % Kernel density estimation for each genotype's zeroTAboveThresh
    if ~isempty(kirZeroT{c})
        [fKir, xiKir] = ksdensity(kirZeroT{c});
        plot(xiKir, fKir, 'Color', settings.geneColor{1}, 'LineWidth', 2, 'DisplayName', 'KIR');
    end
    if ~isempty(wtZeroT{c})
        [fWt, xiWt] = ksdensity(wtZeroT{c});
        plot(xiWt, fWt, 'Color', settings.geneColor{2}, 'LineWidth', 2, 'DisplayName', 'WT');
    end
    if ~isempty(naZeroT{c})
        [fNa, xiNa] = ksdensity(naZeroT{c});
        plot(xiNa, fNa, 'Color', settings.geneColor{3}, 'LineWidth', 2, 'DisplayName', 'NA');
    end
    
    % Set plot limits and labels
    xlim([0 5]); ylim([0 1.2])
    ylabel('Density');
    xlabel('Time to Object (s)');
    
    % Add a legend
    if c == 1
        legend show;
    end

    % Set title for each gain
    title([num2str(settings.pursuitGain(c)) 'X']);
end


%% Statistical test
% Run ANCOVA
% pvalues = run_ancova_on_data(kirObjPos, wtObjPos, naObjPos, kirZeroT, wtZeroT, naZeroT, settings);
% pval_text = {['p(obj) = ' num2str(pvalues.ObjPos)], ['p(gene) = ' num2str(pvalues.Genotype)], ['p(gain) = ' num2str(pvalues.Gain)]};
% text(0.95, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% Run ANOVA
[pSeparate,~] = run_anova_on_jumps(kirZeroT, wtZeroT, naZeroT, 'JumpCorrection', folder);
pval_text = {['p(gene) = ' num2str(pSeparate.Genotype)], ['p(gain) = ' num2str(pSeparate.Gain)]};
text(0.95, 0.95, pval_text, 'Units', 'normalized', 'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

end
