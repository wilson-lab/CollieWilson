function [mdl, stats, tbl] = run_dir_sweep_lme(avg_srRL_Rightward, avg_srRL_Leftward)
% run_dir_sweep_lme
% Repeated-measures LME on firing rate with factors Direction, SweepPos, and their interaction.
% Aligns matrices as Rightward(:,1:end-1) vs Leftward(:,2:end) and includes Animal as a random effect.
%
% INPUTS:
%   avg_srRL_Rightward : [nAnimals x nCols] firing rates for rightward motion
%   avg_srRL_Leftward  : [nAnimals x nCols] firing rates for leftward motion
%
% OUTPUTS:
%   mdl   : fitted LinearMixedModel
%   stats : ANOVA table from mdl (Satterthwaite df)
%   tbl   : analysis table passed to fitlme (after NaN removal)
%
% CREATED: 09/29/2025 - MC

    % Align columns per spec
    R = avg_srRL_Rightward(:, 1:end-1);
    L = avg_srRL_Leftward(:,  2:end);

    % Basic checks
    if any(size(R) ~= size(L))
        error('Aligned matrices must match in size after trimming: size(R)=%s, size(L)=%s', ...
              mat2str(size(R)), mat2str(size(L)));
    end

    [nAnimals, nPos] = size(R);

    % Vectorize responses
    Resp = [R(:); L(:)];

    % Factors
    Direction = categorical( [repmat({'Rightward'}, nAnimals*nPos, 1); ...
                               repmat({'Leftward'},  nAnimals*nPos, 1)] );

    % Treat sweep position as categorical factor (bins)
    SweepPos = categorical( repmat((1:nPos)', 2*nAnimals, 1) );

    % Animal IDs
    Animal = categorical( repelem((1:nAnimals)', 2*nPos, 1) );

    % Build table and drop missing rows
    tbl = table(Resp, Direction, SweepPos, Animal, ...
        'VariableNames', {'Response','Direction','SweepPos','Animal'});

    valid = ~isnan(tbl.Response);
    if ~all(valid)
        tbl = tbl(valid, :);
        % (Optional) warn if many NaNs were removed:
        % warning('Removed %d rows with NaN responses.', sum(~valid));
    end

    % Fit mixed model: random intercept per Animal
    mdl = fitlme(tbl, 'Response ~ Direction*SweepPos + (1|Animal)');

    % ANOVA-style table with Satterthwaite df
    stats = anova(mdl, 'DFMethod', 'satterthwaite');
end
