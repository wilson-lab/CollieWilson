% LC10a Heatmap Analysis
%
% Description:
% This script analyzes the spatial distribution of LC10 neuron dendrites 
% and their synaptic output strength onto downstream neurons (e.g., A02, DNa02) 
% in the Drosophila visual system. The analysis includes loading heatmap 
% data, processing receptive field (RF) distributions, plotting both raw 
% and weighted RFs, and saving the results. Additionally, the script generates 
% a weighted input diagram for the synaptic connections to DNa02, visualized 
% as a graph with line thickness proportional to synapse count.
%
% Main Steps:
% 1. Load and process heatmap data for LC10a neurons.
% 2. Estimate and visualize the receptive fields for all LC10a neurons.
% 3. Analyze and plot the receptive fields for each neuron individually.
% 4. Calculate and plot distributions of LC10a inputs to downstream A02 neurons.
% 5. Generate a weighted input diagram for DNa02, showing synaptic connections 
%    with line thickness proportional to synapse count and vertically arranged 
%    input/output labels.
% 6. Save the analyzed data and plots.
%
% Outputs:
% - Plots of receptive fields for LC10a neurons and their A02 inputs.
% - Weighted input diagram for DNa02.
% - Data tables containing adjusted receptive fields.
% - .mat files storing the processed results for further modeling.
%
% Created: 06/24/24 MC
% Updated: 10/08/25 MC male and female

%% Initialize
clear
close all

% Get current file path and extract directories
thisFile = matlab.desktop.editor.getActiveFilename;
[mainPath,~,~] = fileparts(thisFile); % Main directory path
%flyPath = [mainPath '/AOTU_FemaleCNS'];
flyPath = [mainPath '/AOTU_MaleCNS'];

cd(fullfile(flyPath, 'heatmaps'))
heatmapFiles = dir('*lc10a_heatmap.png'); % Colorful heatmap files
greyscaleFiles = dir('*lc10a_bw_heatmap.png'); % Greyscale heatmap files for analysis
nNOI = size(greyscaleFiles,1); % Number of neurons of interest (NOI)

% Set parameters for smoothing and retina boundaries
gwin = 30; % Gaussian window size
cntr_boundary = -15; % Retina contralateral boundary (degrees)
ipsi_boundary = 165; % Retina ipsilateral boundary (degrees)
rTicks = cntr_boundary:15:ipsi_boundary; % Tick positions for plot (every 15 degrees)

%% Estimate receptive field (RF) position using all LC10a neurons
disp('Estimating A-P position...')
close all

% --- Load heatmap
cd(fullfile(flyPath, 'heatmaps'))
[ao_all, ~, alpha] = imread("all_lc10a.png");

% --- Generate smoothed A-P distribution
[ap_mean, ~] = heatmap_average(ao_all, alpha, gwin);

% --- Convert to binary presence/absence
minWT = 0.01;
lc_binary = double(ap_mean >= minWT);

% --- Determine retina position range
ss = find(ischange(lc_binary));
rRange = ss(1):ss(2);
rPos = linspace(cntr_boundary, ipsi_boundary, numel(rRange));

% --- Extend positions beyond calculated range
spacing = rPos(2) - rPos(1);
ext_start = rPos(1) - (ss(1)-1:-1:1) * spacing;
ext_stop  = rPos(end) + (1:(numel(lc_binary)-ss(2))) * spacing;
rPos_ext  = [ext_start, rPos, ext_stop];

% --- Plot
figure('Position', [100 100 450 750]);
tiledlayout(3,1,"TileSpacing","compact")

nexttile([2 1]); imshow(ao_all)
nexttile; hold on
plot(rPos_ext, lc_binary, 'k', 'LineWidth', 1.5)
plot(rPos_ext, ap_mean, 'r', 'LineWidth', 1.5)
xticks(rTicks); xtickangle(90)
xlabel('Retina Position (deg)')
axis tight
sgtitle('All LC10a')

% --- Save figure
cd(fullfile(flyPath, 'rf'))
plotname = 'rf_all_LC10a';
saveas(gcf, [plotname '.png'])
set(gcf, 'renderer', 'Painters')
saveas(gcf, [plotname '.svg'])

%% Prepare output weighting
% Read table
cd(fullfile(flyPath, 'data'))
dn_type_avg_DNa02 = readtable('dn_type_avg_DNa02.csv');

% Extract relevant columns
dn_subset = dn_type_avg_DNa02(:, {'type', 'DNa02'});

% Add normalized weight column (range 0–1)
dn_subset.norm_DNa02 = dn_subset.DNa02 ./ max(dn_subset.DNa02, [], 'omitnan');

%% Estimate distribution for each LC10a neuron separately
disp('Estimating LC10a distributions...')
close all

% Expect dn_subset to exist with columns: {'type','DNa02','norm_DNa02'}
types_in_table = string(dn_subset.type);

all_means           = zeros(nNOI, numel(rRange));
all_means_weighted  = zeros(nNOI, numel(rRange));

outdir = fullfile(flyPath, 'rf');
if ~exist(outdir, 'dir'); mkdir(outdir); end

for i = 1:nNOI
    gname = greyscaleFiles(i).name;
    cname = heatmapFiles(i).name;

    % Parse cell type for this neuron
    this_type = string(extractBefore(gname, '_lc10a'));

    % Load images
    [ao_grey, ~, alpha] = imread(fullfile(flyPath, 'heatmaps', gname));
    ao_color            = imread(fullfile(flyPath, 'heatmaps', cname));

    % Process heatmap -> A-P mean (suppress dv_density)
    [ap_mean, ~] = heatmap_average(ao_grey, alpha, gwin);

    % Plot (kept same)
    figure('Position', [100 100 450 750]);
    tiledlayout(3,1,"TileSpacing","compact")
    nexttile([2 1]); imshow(ao_color)
    nexttile; hold on
    x = rPos_ext(1:numel(ap_mean));          % ensure matching lengths
    plot(x, ap_mean, 'r', 'LineWidth', 1.5)
    xticks(rTicks); xtickangle(90); xlabel('Retina Position (deg)'); axis tight

    % Store unweighted mean over requested range
    idx = rRange(rRange <= numel(ap_mean));  % safe indexing into ap_mean
    all_means(i,:) = ap_mean(idx);

    % ===== Look up norm_DNa02 for this cell type and compute weighted mean =====
    hit = find(types_in_table == this_type, 1, 'first');
    if ~isempty(hit) && ~isnan(dn_subset.norm_DNa02(hit))
        w = dn_subset.norm_DNa02(hit);            % scalar in [0,1]
    else
        w = NaN;                                   % or set to 0 if you prefer: w = 0;
        % warning('No matching type or NaN norm_DNa02 for "%s".', this_type);
    end
    all_means_weighted(i,:) = w .* ap_mean(idx);   % weighted synapse density

    % Save
    plotname = "rf_lc10a_mean_" + this_type;
    saveas(gcf, fullfile(outdir, plotname + ".png"))
    set(gcf, 'renderer', 'Painters')
    saveas(gcf, fullfile(outdir, plotname + ".svg"))
end

disp('Complete.')

%% ===== 1) SAVE WEIGHTED MEANS AS CSV (one column per cell) =====

% Build per-cell names (type_id if available; fallback to type_i)
cell_names = strings(nNOI,1);
for i = 1:nNOI
    gname = greyscaleFiles(i).name;
    this_type = string(extractBefore(gname, '_'));
    cid = extractBetween(gname, this_type + "_", "_lc10");   % may be empty
    if ~isempty(cid)
        cell_names(i) = this_type + "_" + string(cid);
    else
        cell_names(i) = this_type + "_" + string(i);
    end
end

% X-axis positions (use rPos at the bins in rRange)
x_positions = rPos;               % column vector
T = table(x_positions', 'VariableNames', {'rPos_deg'});

% Append each cell's weighted RF as a column
for i = 1:nNOI
    cname = matlab.lang.makeValidName(cell_names(i));
    T.(cname) = all_means_weighted(i,:)';   % column
end

% Write CSV
csv_out = fullfile(outdir, 'weighted_rf_by_cell.csv');
writetable(T, csv_out);
fprintf('Saved weighted RFs to: %s\n', csv_out);

%% ===== 2) PLOT ALL WEIGHTED SYNAPSE DENSITIES TOGETHER (Spectral gradient) =====

figure('Color','w','Position',[100 100 700 420]); hold on

% --- Generate a Spectral-like continuous gradient for nNOI cells
% Define a rough spectral base palette (ColorBrewer-like)
spectral_base = [ ...
    158, 1, 66; ...
    213, 62, 79; ...
    244, 109, 67; ...
    253, 174, 97; ...
    254, 224, 139; ...
    230, 245, 152; ...
    171, 221, 164; ...
    102, 194, 165; ...
    50, 136, 189; ...
    94, 79, 162] / 255;

% Interpolate across that palette for however many cells you have
nColors = size(spectral_base,1);
x_base = linspace(1, nColors, nColors);
xq = linspace(1, nColors, nNOI);
spectral_colors = interp1(x_base, spectral_base, xq, 'linear');

% --- Plot each weighted mean with its own color
for i = 1:nNOI
    plot(x_positions, all_means_weighted(i,:), ...
        'Color', spectral_colors(i,:), 'LineWidth', 1.6);
end

xlabel('Retina Position (deg)');
ylabel('Weighted synapse density (a.u.)');
title('Weighted LC10a synapse densities by cell');
legend(cell_names, 'Interpreter','none', 'Location','eastoutside');
axis tight; box on; grid on

% --- Save outputs
saveas(gcf, fullfile(outdir, 'weighted_rf_all.png'));
set(gcf, 'Renderer', 'painters');
saveas(gcf, fullfile(outdir, 'weighted_rf_all.svg'));
