function run_RL_anova_repeated(avg_srRL_Rightward, avg_srRL_Leftward, sweepPosR, sweepPosL)

[nFlies, nSweeps] = size(avg_srRL_Rightward);

% Find shared sweep positions between directions
[sharedSweepPos, idxR, idxL] = intersect(sweepPosR, sweepPosL);

% Restrict to sweep positions > 0 (ipsilateral only)
keepIdx = sharedSweepPos > 0;
sharedSweepPos = sharedSweepPos(keepIdx);
idxR = idxR(keepIdx);
idxL = idxL(keepIdx);

% Preallocate arrays
flyID = [];
sweepPos = [];
direction = {};
response = [];

% Loop over flies and selected sweeps
for nt = 1:nFlies
    for i = 1:length(sharedSweepPos)
        sR = idxR(i);
        sL = idxL(i);

        r_val = avg_srRL_Rightward(nt, sR);
        l_val = avg_srRL_Leftward(nt, sL);

        % Rightward
        flyID(end+1,1) = nt;
        sweepPos(end+1,1) = sharedSweepPos(i);
        direction{end+1,1} = 'Rightward';
        response(end+1,1) = r_val;

        % Leftward
        flyID(end+1,1) = nt;
        sweepPos(end+1,1) = sharedSweepPos(i);
        direction{end+1,1} = 'Leftward';
        response(end+1,1) = l_val;
    end
end

% Build table
T = table;
T.Fly = categorical(flyID);
T.SweepPos = categorical(sweepPos);  % Treat sweep positions as categories
T.Direction = categorical(direction);
T.Response = response;

% Fit linear mixed-effects model
lme = fitlme(T, 'Response ~ Direction * SweepPos + (1|Fly)');

% Display ANOVA table
anovaTbl = anova(lme);
disp('Repeated Measures ANOVA (SweepPos > 0 only):');
disp(anovaTbl);

% ---- Sanity Plot ----
figure('Name', 'R vs L per Sweep Position (>0 only)', 'Position', [100 100 600 400]); hold on;

% Extract sweep labels in numeric form
sweepLevels = categories(T.SweepPos);
sweepPosVals = str2double(sweepLevels);
nLevels = numel(sweepLevels);

% Build per-fly response matrices
R_mat = NaN(nFlies, nLevels);
L_mat = NaN(nFlies, nLevels);

for i = 1:nFlies
    for j = 1:nLevels
        idx_R = T.Fly == categorical(i) & T.SweepPos == sweepLevels{j} & T.Direction == "Rightward";
        idx_L = T.Fly == categorical(i) & T.SweepPos == sweepLevels{j} & T.Direction == "Leftward";

        if any(idx_R)
            R_mat(i, j) = T.Response(find(idx_R, 1));
        end
        if any(idx_L)
            L_mat(i, j) = T.Response(find(idx_L, 1));
        end
    end
end

% Plot per-fly traces
for i = 1:nFlies
    plot(sweepPosVals, R_mat(i,:), '-', 'Color', [0.8 0.8 1]);
    plot(sweepPosVals, L_mat(i,:), '-', 'Color', [1 0.8 0.8]);
end

% Plot means ± SEM
mean_R = mean(R_mat, 1, 'omitnan');
mean_L = mean(L_mat, 1, 'omitnan');
sem_R = std(R_mat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(R_mat), 1));
sem_L = std(L_mat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(L_mat), 1));

plot(sweepPosVals, mean_R, '-', 'Color', [0 0 1], 'LineWidth', 2, 'DisplayName', 'Rightward');
plot(sweepPosVals, mean_L, '-', 'Color', [1 0 0], 'LineWidth', 2, 'DisplayName', 'Leftward');

xlabel('Sweep Position (deg)');
ylabel('Response');
title('Mean ± SEM: Rightward vs Leftward (SweepPos > 0)');
grid on;

end
