function [mc, ran, p_anova] = rm_tukeykramer(thisData, trialTypes)
% CREATED: 09/08/2025 - MC
% RM-ANOVA (within-subjects) + Tukey–Kramer posthocs.
% thisData: [nFlies x nCond]; rows = flies, cols = conditions.

    arguments
        thisData double
        trialTypes = []
    end

    % --- Checks
    [nFlies, nCond] = size(thisData);
    if nCond < 2
        error('Need at least 2 conditions (columns).');
    end

    % --- Labels
    if isempty(trialTypes)
        trialTypes = "C" + string(1:nCond);
    else
        trialTypes = string(trialTypes(:).');
        if numel(trialTypes) ~= nCond
            error('Length of trialTypes must equal number of conditions.');
        end
    end
    varNames = matlab.lang.makeUniqueStrings( ...
                  matlab.lang.makeValidName("C_" + trialTypes));

    % --- Complete cases (fitrm requires no NaNs per row)
    validFly = all(isfinite(thisData), 2);
    if sum(validFly) < 2
        error('Not enough valid flies after removing NaNs (need >= 2).');
    end
    Y = thisData(validFly, :);

    % --- Wide table & WithinDesign (nCond x 1)
    Twide  = array2table(Y, 'VariableNames', cellstr(varNames));
    CondCats = categorical((1:nCond).', 1:nCond, trialTypes);
    Within = table(CondCats, 'VariableNames', {'Condition'});

    % --- Repeated-measures model & ANOVA
    rm  = fitrm(Twide, sprintf('%s-%s ~ 1', varNames(1), varNames(end)), ...
                'WithinDesign', Within);
    ran = ranova(rm, 'WithinModel', 'Condition');

    % --- Robust extraction of omnibus p for Condition
    p_anova = get_condition_pvalue(ran);

    % --- Tukey–Kramer posthocs
    mc = multcompare(rm, 'Condition', 'ComparisonType','tukey-kramer');

    % Normalize column names for readability across MATLAB versions
    vn = mc.Properties.VariableNames;
    if ismember('Condition_1', vn) && ismember('Condition_2', vn)
        mc.Properties.VariableNames(strcmp(vn,'Condition_1')) = {'Cond1'};
        mc.Properties.VariableNames(strcmp(vn,'Condition_2')) = {'Cond2'};
    elseif ismember('A', vn) && ismember('B', vn)
        mc.Properties.VariableNames(strcmp(vn,'A')) = {'Cond1'};
        mc.Properties.VariableNames(strcmp(vn,'B')) = {'Cond2'};
    end
    if iscategorical(mc.Cond1), mc.Cond1 = string(mc.Cond1); end
    if iscategorical(mc.Cond2), mc.Cond2 = string(mc.Cond2); end
end

% ===== helper: pick Condition row & best p-value column across versions =====
function p = get_condition_pvalue(ran)
    p = NaN;

    % Find the 'Condition' row
    row = false(height(ran),1);

    % 1) Some versions have a 'Term' variable
    if ismember('Term', ran.Properties.VariableNames)
        if iscategorical(ran.Term)
            row = ran.Term == 'Condition';
        else
            row = strcmp(string(ran.Term), 'Condition');
        end
    end

    % 2) Otherwise, use row names like 'Condition' or 'Condition:...' etc.
    if ~any(row) && ~isempty(ran.Properties.RowNames)
        rn = string(ran.Properties.RowNames);
        row = rn == "Condition" | startsWith(rn, "Condition");
    end

    % 3) As a final fallback, look for a variable literally named 'Condition'
    if ~any(row) && ismember('Condition', ran.Properties.VariableNames)
        v = ran.Condition;
        if iscategorical(v), v = string(v); end
        row = strcmp(string(v), 'Condition');
    end

    if ~any(row)
        warning('Could not locate Condition row in ranova table.');
        return
    end

    % Prefer GG-corrected p if present; otherwise try HF, LB, then plain p
    pcols = {'GG_pValue','pValueGG','HF_pValue','pValueHF','LB_pValue','pValueLB','pValue'};
    for k = 1:numel(pcols)
        if ismember(pcols{k}, ran.Properties.VariableNames)
            p = ran{row, pcols{k}};
            if iscell(p), p = cell2mat(p); end
            return
        }
    end

    % As absolute last resort, try to find any column with 'pValue' in its name
    anyp = contains(ran.Properties.VariableNames, 'pValue', 'IgnoreCase', true);
    if any(anyp)
        p = ran{row, find(anyp,1)};
        if iscell(p), p = cell2mat(p); end
    end
end
