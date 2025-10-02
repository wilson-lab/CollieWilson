function [objAtJump, corrTimes] = setpoint_jumps2(panelps, barjump, ttime)
% SETPOINT_JUMP3 - For each detected jump, return the object position at the
% jump and the time to correct (first zero-crossing) without binning/averaging.
%
% INPUTS:
%   panelps  : [time x trials x conditions] panel position
%   barjump  : [time x trials x conditions] jump trigger (diff==1 right, diff==2 left)
%   ttime    : [time x 1] time vector
%   doPlot   : (optional) logical, if true plots successful vs unsuccessful traces
%
% OUTPUTS:
%   objAtJump{c} : vector of |position at jump| for every jump in condition c
%   corrTimes{c} : vector of correction times (s) to first zero-crossing for condition c
%
% Updated: 08/25/2025 - MC

doPlot = false;

%% setup
nCond = size(panelps, 3);
preWin = 0.5;                      % seconds before jump to include
pstWin = 10;                        % seconds after jump to search
dt = median(diff(ttime), "omitnan");
preIdx = max(1, round(preWin / dt));  % samples
pstIdx = max(1, round(pstWin / dt));  % samples
winLen = preIdx + pstIdx + 1;
tAxis = (-preIdx:pstIdx) * dt;        % time aligned to jump onset (0 at jump)

minPosChange = 10; % minimum step in panel position to consider a jump (deg)

objAtJump = cell(1, nCond);
corrTimes = cell(1, nCond);

% For optional plotting: collect aligned traces per condition
succTraces = cell(1, nCond);  % [winLen x nSucc]
failTraces = cell(1, nCond);  % [winLen x nFail)

%% per condition
for c = 1:nCond
    % bias-correct panel positions and vectorize time x trials
    biasHD = mean(panelps(:,:,c), 'all', 'omitnan');
    thispanelps = reshape(panelps(:,:,c) - biasHD, [], 1);
    thisjumptrg = reshape(barjump(:,:,c), [], 1);

    % jump indices (transition in trigger)
    jr = find(diff(thisjumptrg) == 1);  % right jumps
    jl = find(diff(thisjumptrg) == 2);  % left jumps
    jall = sort([jr; jl]);

    pos_list  = [];
    time_list = [];

    % trace collectors
    succ_mat = [];  % columns are traces
    fail_mat = [];

    L = numel(thispanelps);

    for k = 1:numel(jall)
        idx = jall(k);

        % quick bounds check for initial search region
        if idx + pstIdx > L, continue; end

        % find the *actual* jump onset: first significant step after idx
        dseg = abs(diff(thispanelps(idx:min(idx+pstIdx, L))));
        jrel = find(dseg > minPosChange, 1, 'first');
        if isempty(jrel)
            % no clear step—cannot even define a proper window → count as fail if we can extract a window
            wStart = idx - preIdx;
            wEnd   = idx + pstIdx;
            if wStart >= 1 && wEnd <= L
                seg0 = thispanelps(wStart:wEnd);
                if numel(seg0) == winLen
                    fail_mat = [fail_mat, seg0]; %#ok<AGROW>
                end
            end
            continue;
        end

        jOn = idx + jrel;  % absolute index of jump onset

        % window around the actual jump onset
        wStart = jOn - preIdx;
        wEnd   = jOn + pstIdx;
        if wStart < 1 || wEnd > L
            % out of bounds – cannot analyze; try to store as fail if a partial window exists (skip to keep lengths consistent)
            continue;
        end

        seg = thispanelps(wStart:wEnd);
        if numel(seg) ~= winLen
            continue; % safety
        end

        % require non-NaN pre-jump segment
        if any(isnan(seg(1:preIdx)))
            fail_mat = [fail_mat, seg]; %#ok<AGROW>
            continue;
        end

        % zero-crossing after the jump (relative to preIdx)
        sc = diff(sign(seg(preIdx+1:end)));
        zrel = find(sc ~= 0, 1, 'first'); % first sign change
        if isempty(zrel)
            % no zero crossing in window
            fail_mat = [fail_mat, seg]; %#ok<AGROW>
            continue;
        end

        % ensure no NaNs between jump and zero-crossing
        if any(isnan(seg(preIdx+1:preIdx+zrel)))
            fail_mat = [fail_mat, seg]; %#ok<AGROW>
            continue;
        end

        % SUCCESS: record |position at jump| and correction time
        pos_at_jump = abs(seg(preIdx));   % absolute position at jump
        t_correct   = zrel * dt;          % seconds to zero-cross from jump

        pos_list(end+1,1)  = pos_at_jump; %#ok<AGROW>
        time_list(end+1,1) = t_correct;   %#ok<AGROW>
        succ_mat = [succ_mat, seg];       %#ok<AGROW>
    end

    objAtJump{c} = pos_list;
    corrTimes{c} = time_list;

    % store traces for optional plotting
    succTraces{c} = succ_mat;
    failTraces{c} = fail_mat;
end

%% optional plotting
if doPlot
    % simple color palette per condition
    cmap = lines(max(1, nCond));

    tl = tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

    % Successful corrections
    ax1 = nexttile(tl,1); hold(ax1, 'on');
    for c = 1:nCond
        if ~isempty(succTraces{c})
            plot(ax1, tAxis, succTraces{c}, 'Color', [cmap(c,:) 0.25]); % light lines
            % mean trace
            mtrace = mean(succTraces{c}, 2, 'omitnan');
            plot(ax1, tAxis, mtrace, 'LineWidth', 2, 'Color', cmap(c,:));
        end
    end
    xline(ax1, 0, ':', 'LineWidth', 1);
    yline(ax1, 0, ':', 'LineWidth', 1);
    xlabel(ax1, 'Time from jump onset (s)');
    ylabel(ax1, 'Position (deg, bias-corrected)');
    title(ax1, 'Successful corrections');

    % Unsuccessful corrections
    ax2 = nexttile(tl,2); hold(ax2, 'on');
    for c = 1:nCond
        if ~isempty(failTraces{c})
            plot(ax2, tAxis, failTraces{c}, 'Color', [cmap(c,:) 0.35]);
            % mean (for context)
            mtrace = mean(failTraces{c}, 2, 'omitnan');
            plot(ax2, tAxis, mtrace, 'LineWidth', 2, 'Color', cmap(c,:));
        end
    end
    xline(ax2, 0, ':', 'LineWidth', 1);
    yline(ax2, 0, ':', 'LineWidth', 1);
    xlabel(ax2, 'Time from jump onset (s)');
    ylabel(ax2, 'Position (deg, bias-corrected)');
    title(ax2, 'Unsuccessful corrections');

    title(tl, 'Setpoint jump responses (aligned to jump onset)');
end
end
