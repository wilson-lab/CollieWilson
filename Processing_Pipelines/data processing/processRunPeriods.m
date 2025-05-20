% processRunPeriods
% This function processes running periods in an experiment by merging short non-running gaps
% between consecutive running periods. It assumes that the input `run_idx` is a logical matrix
% indicating running (true) and non-running (false) periods. Short gaps between runs are merged
% if the gap duration is less than a specified maximum.
%
% INPUTS:
%   run_idx     - 3D logical array (time x trial x condition) indicating running periods
%   expttime    - Time vector (seconds)
%   max_gap_dur - Maximum allowable gap duration (seconds) to merge short gaps between runs
%
% OUTPUTS:
%   runIdxProcessed - 3D logical array (same size as run_idx) with merged running periods
%
% CREATED: [Date] MC
%
function runIdxProcessed = processRunPeriods(run_idx, expttime, max_gap_dur)
    % Convert max gap duration to indices
    dt = expttime(2) - expttime(1); % Time step
    max_gap_idx = max_gap_dur / dt; % Maximum gap duration in indices

    % Initialize processed run_idx
    runIdxProcessed = run_idx;

    % Loop through trials and conditions (assuming 3D run_idx: time x trial x condition)
    for trial = 1:size(run_idx, 2)
        for cond = 1:size(run_idx, 3)
            % Find start and end of running periods
            run_diff = diff([0; run_idx(:, trial, cond); 0]); % Edge detection
            run_starts = find(run_diff == 1); % Start of running
            run_ends = find(run_diff == -1) - 1; % End of running

            % Merge short non-running gaps between running periods
            for i = 1:(length(run_starts) - 1)
                % If the gap between the end of one run and the start of the next run is shorter than max_gap_idx
                if run_starts(i + 1) - run_ends(i) - 1 <= max_gap_idx
                    % Merge by extending the current run period to the end of the next run period
                    run_ends(i) = run_ends(i + 1); % Merge the run ends
                    run_starts(i + 1) = NaN; % Mark for removal
                    run_ends(i + 1) = NaN; % Mark for removal
                end
            end

            % Remove marked invalid periods (NaN entries)
            run_starts = run_starts(~isnan(run_starts));
            run_ends = run_ends(~isnan(run_ends));

            % Reset run_idx for this trial and condition
            runIdxProcessed(:, trial, cond) = false;
            for i = 1:length(run_starts)
                runIdxProcessed(run_starts(i):run_ends(i), trial, cond) = true;
            end
        end
    end
end
