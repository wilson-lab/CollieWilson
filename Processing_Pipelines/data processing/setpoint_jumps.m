% SETPOINT_JUMPS - This function analyzes the time required to correct setpoint deviations (jumps)
% based on object position during pursuit behavior. It computes the average time to cross zero after 
% small and large jumps, using a pre-defined threshold to separate the two types of jumps.

% INPUTS:
%   panelps   - 3D array of panel position data, where each slice represents a different condition, 
%               and each column represents a trial.
%   barjump   - 3D array of bar jump triggers used to detect jumps.
%   ttime     - Time vector corresponding to the temporal resolution of the data.

% OUTPUTS:
%   smallJumps - Average time (in seconds) to cross zero for small jumps (below threshold) for each condition.
%   largeJumps - Average time (in seconds) to cross zero for large jumps (above threshold) for each condition.

% Created: 09/03/24 by MC
% Updated: N/A

% The function:
% - Initializes parameters such as the pre/post-jump windows, threshold for jump magnitude, 
%   and sets up outputs for small and large jumps.
% - Processes each condition by identifying bar jumps, computing the position at the time of the jump, 
%   and determining how long it takes to cross zero after the jump.
% - Returns the average time to correct small (< angle threshold) and large (> angle threshold) 
%   jumps for each condition.
%
function [smallJumps, largeJumps] = setpoint_jumps(panelps, barjump, ttime)
%% initialize
% dataset info
nCond = size(panelps, 3);

% load processing settings
settings = processSettings();

% set how much time pre/post jump to analyze
preWin = 0.5; % s
pstWin = 5; % s
preIdx = fetchTimeIdx(ttime, preWin); % idx
pstIdx = fetchTimeIdx(ttime, pstWin); % idx

% set min change to be considered a jump
minPosChange = 10;

angle = 35;

% initialize outputs
smallJumps = NaN(1, nCond);
largeJumps = NaN(1, nCond);

% for each condition
for c = 1:nCond
    % Initialize arrays to store object positions and time to cross zero
    obj_at_jump = [];
    zero_cross_times = [];

    % estimate HD bias
    biasHD = mean(panelps(:,:,c), 'all', 'omitnan');

    % fetch data
    thispanelps = reshape(panelps(:,:,c) - biasHD, [], 1);
    thisjumptrg = reshape(barjump(:,:,c), [], 1);

    % find indices where jumps occur
    jump_right = find(diff(thisjumptrg) == 1);
    jump_left = find(diff(thisjumptrg) == 2);

    % Process jumps
    process_jumps = @(jump_indices) arrayfun(@(idx) process_single_jump(idx, thispanelps, preIdx, pstIdx, minPosChange, ttime, preWin), jump_indices, 'UniformOutput', false);
    
    % Process right and left jumps
    right_jump_data = process_jumps(jump_right);
    left_jump_data = process_jumps(jump_left);

    % Combine right and left data
    right_jump_data = [right_jump_data{:}];
    left_jump_data = [left_jump_data{:}];
    
    % Extract object positions and zero-crossing times
    obj_at_jump = [right_jump_data.obj_at_jump, left_jump_data.obj_at_jump];
    zero_cross_times = [right_jump_data.zero_cross_times, left_jump_data.zero_cross_times];

    % Calculate averages for small and large jumps
    small_jump_idx = obj_at_jump < angle;
    large_jump_idx = obj_at_jump > angle;

    if any(small_jump_idx)
        smallJumps(c) = mean(zero_cross_times(small_jump_idx), 'omitnan');
    end
    if any(large_jump_idx)
        largeJumps(c) = mean(zero_cross_times(large_jump_idx), 'omitnan');
    end
end

end

function result = process_single_jump(idx, panelps, preIdx, pstIdx, minPosChange, ttime, preWin)
    % Initialize result structure
    result.obj_at_jump = NaN;
    result.zero_cross_times = NaN;

    % Find the first point where the panel position changes significantly
    jump_detected_idx = find(abs(diff(panelps(idx:idx+pstIdx))) > minPosChange, 1, 'first');
    if isempty(jump_detected_idx)
        return; % No significant jump detected
    end
    jump_actual_idx = idx + jump_detected_idx;

    % Extract data before and after the jump
    thiswin_panelps = panelps(jump_actual_idx-preIdx:jump_actual_idx+pstIdx);

    % Omit if NaNs are present before the jump
    if any(isnan(thiswin_panelps(1:preIdx)))
        return; % Skip this trial
    end

    % Find zero-crossing time after the jump
    sign_change_window = diff(sign(thiswin_panelps(preIdx+1:end)));
    zero_cross_idx = find(sign_change_window ~= 0, 1, 'first');
    if isempty(zero_cross_idx)
        return; % No zero-crossing found
    end

    % Omit if NaNs are present between the jump and the zero-crossing
    if any(isnan(thiswin_panelps(preIdx+1:preIdx+zero_cross_idx)))
        return; % Skip this trial
    end

    % Store the object position at the jump and time to cross zero
    result.obj_at_jump = abs(thiswin_panelps(preIdx));  % Absolute object position at jump
    result.zero_cross_times = ttime(preIdx + zero_cross_idx) - preWin;  % Time to cross zero
end
