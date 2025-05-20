% ANALYZE_IINJ Analyzes neural and behavioral responses to current injection.
% 
% This function processes time series data of neural activity and behavior 
% following intracellular current injection. It normalizes data, identifies 
% trials where the fly exhibits a turn, computes the frequency of turning, 
% and extracts the maximum spike rate and angular velocity during the pulse period.
%
% INPUTS:
%   - Iinject: Injected current data (time x trials)
%   - SpikeRt: Spike rate data (time x trials)
%   - Forward: Forward velocity data (time x trials)
%   - Angular: Angular velocity data (time x trials)
%   - Sideway: Sideways velocity data (time x trials)
%   - buffDur: Duration before pulse onset (s)
%   - pulseDur: Duration of the pulse (s)
%   - ttime: Time array (s)
%
% OUTPUTS:
%   - turnFrequency: Fraction of trials where the fly turns
%   - turnStats: Struct containing max spike rate and max angular velocity
%
% CREATED: 03/04/2025 MC
%
function [turnFrequency, turnStats] = analyze_iinj(Iinject, SpikeRt, Forward, Angular, Sideway, buffDur, pulseDur, ttime)
%% Initialize
preIdx = fetchTimeIdx(ttime, buffDur); % Index for pre-pulse window
halfpulseIdx = fetchTimeIdx(ttime, pulseDur/2); % Index for halfpulse window
pulseIdx = fetchTimeIdx(ttime, pulseDur); % Index for pulse window
endIdx = size(Angular, 1); % Last index for time series

% Create time array matching filtered data length
dataTime = ttime(1:size(Iinject,1)) * 1000; % Convert to milliseconds

%% Compute baseline (average from 0 to preIdx)
baselineIinject = mean(Iinject(1:preIdx, :), 1, 'omitnan');
baselineSpikeRt = mean(SpikeRt(1:preIdx, :), 1, 'omitnan');
baselineForward = mean(Forward(1:preIdx, :), 1, 'omitnan');
baselineAngular = mean(Angular(1:preIdx, :), 1, 'omitnan');
baselineSideway = mean(Sideway(1:preIdx, :), 1, 'omitnan');

% Normalize data by subtracting baseline
normIinject = Iinject - baselineIinject;
normSpikeRt = SpikeRt - baselineSpikeRt;
normForward = Forward - baselineForward;
normAngular = Angular - baselineAngular;
normSideway = Sideway - baselineSideway;

%% Select trials where the average rotational velocity during the pulse is positive
pulseWindow = preIdx:endIdx;
turnTrials = mean(normAngular(pulseWindow, :), 1, 'omitnan') > 0;
noTurnTrials = mean(normAngular(pulseWindow, :), 1, 'omitnan') <= 0;

% Compute fraction of turning trials
turnFrequency = sum(turnTrials) / (sum(turnTrials) + sum(noTurnTrials));

% Create filtered data variables
filtIinject = normIinject;
filtSpikeRt = normSpikeRt;
filtForward = normForward;
filtAngular = normAngular;
filtSideway = normSideway;

% Keep only valid trials
filtIinject(:, ~turnTrials) = NaN;
filtSpikeRt(:, ~turnTrials) = NaN;
filtForward(:, ~turnTrials) = NaN;
filtAngular(:, ~turnTrials) = NaN;
filtSideway(:, ~turnTrials) = NaN;

%% Extract change values during pulse period
% Set analysis window
testWindow = preIdx+halfpulseIdx:preIdx+pulseIdx;

% Extract from ALL data
meanSpikertResponse = mean(normSpikeRt(testWindow,:), 2, 'omitnan');
meanAngularResponse = mean(normAngular(testWindow,:), 2, 'omitnan');

turnStats.spikert = mean(meanSpikertResponse);
turnStats.angular = mean(meanAngularResponse);

% Extract from ONLY turn data
meanSpikertResponseTurn = mean(filtSpikeRt(testWindow,:), 2, 'omitnan');
meanAngularResponseTurn = mean(filtAngular(testWindow,:), 2, 'omitnan');

turnStats.spikertTurnOnly = mean(meanSpikertResponseTurn);
turnStats.angularTurnOnly = mean(meanAngularResponseTurn);

%% Optional Plot
plotData = false; % Set to true to generate plots
if plotData
    figure;
    tiledlayout(2,3);
    
    % Plot raw normalized data
    nexttile; plot(dataTime, normSpikeRt, 'LineWidth', 0.5); hold on;
    plot(dataTime, mean(normSpikeRt,2,'omitnan'),'k', 'LineWidth', 1.5);
    xline(dataTime(preIdx),'k', 'LineWidth', 1.5); yline(0,'k', 'LineWidth', 1.5);
    title('Spike Rate'); xlabel('Time (ms)'); ylabel('Spike Rate'); axis tight; grid on; hold off;
    
    nexttile; plot(dataTime, normForward, 'LineWidth', 0.5); hold on;
    plot(dataTime, mean(normForward,2,'omitnan'),'k', 'LineWidth', 1.5);
    xline(dataTime(preIdx), 'k', 'LineWidth', 1.5); yline(0, 'k', 'LineWidth', 1.5);
    title('Forward Velocity'); xlabel('Time (ms)'); ylabel('Forward Velocity'); axis tight; grid on; hold off;
    
    nexttile; plot(dataTime, normAngular, 'LineWidth', 0.5); hold on;
    plot(dataTime, mean(normAngular,2,'omitnan'),'k', 'LineWidth', 1.5);
    xline(dataTime(preIdx), 'k', 'LineWidth', 1.5); yline(0, 'k', 'LineWidth', 1.5);
    title('Angular Velocity'); xlabel('Time (ms)'); ylabel('Angular Velocity'); axis tight; ylim([-400 400]); grid on; hold off;
    
    % Plot filtered data (turning trials only)
    nexttile; plot(dataTime, filtSpikeRt, 'LineWidth', 0.5); hold on;
    plot(dataTime, mean(filtSpikeRt,2,'omitnan'),'k', 'LineWidth', 1.5);
    xline(dataTime(preIdx), 'k', 'LineWidth', 1.5); yline(0, 'k', 'LineWidth', 1.5);
    title('Filtered Spike Rate'); xlabel('Time (ms)'); ylabel('Spike Rate'); axis tight; grid on; hold off;
    
    nexttile; plot(dataTime, filtForward, 'LineWidth', 0.5); hold on;
    plot(dataTime, mean(filtForward,2,'omitnan'),'k', 'LineWidth', 1.5);
    xline(dataTime(preIdx), 'k', 'LineWidth', 1.5); yline(0, 'k', 'LineWidth', 1.5);
    title('Filtered Forward Velocity'); xlabel('Time (ms)'); ylabel('Forward Velocity'); axis tight; grid on; hold off;
    
    nexttile; plot(dataTime, filtAngular, 'LineWidth', 0.5); hold on;
    plot(dataTime, mean(filtAngular,2,'omitnan'),'k', 'LineWidth', 1.5);
    xline(dataTime(preIdx), 'k', 'LineWidth', 1.5); yline(0, 'k', 'LineWidth', 1.5);
    title('Filtered Angular Velocity'); xlabel('Time (ms)'); ylabel('Angular Velocity'); axis tight; ylim([-400 400]); grid on; hold off;
end

end
