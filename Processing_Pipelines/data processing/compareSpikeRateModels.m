function r2_values = compareSpikeRateModels(spikert, panelps, forward)
% CREATED: 03/30/2025 - MC
%
% Computes R² values for three models predicting spike rate:
%   1. Additive: panelps + forward
%   2. Multiplicative: panelps .* forward
%   3. Full: panelps + forward + (panelps .* forward)

    % Vectorize input matrices
    spike_flat = spikert(:);
    panelps_flat = panelps(:);
    forward_flat = forward(:);

    % Remove NaNs
    valid_idx = ~isnan(spike_flat) & ~isnan(panelps_flat) & ~isnan(forward_flat);
    spike_flat = spike_flat(valid_idx);
    panelps_flat = panelps_flat(valid_idx);
    forward_flat = forward_flat(valid_idx);

    % Define predictors
    X_add = [panelps_flat, forward_flat];                   % Additive model
    X_mult = panelps_flat .* forward_flat;                  % Multiplicative model
    X_full = [panelps_flat, forward_flat, X_mult];          % Full model

    % Fit linear models
    mdl_add = fitlm(X_add, spike_flat);
    mdl_mult = fitlm(X_mult, spike_flat);
    mdl_full = fitlm(X_full, spike_flat);

    % Store R² values in output struct
    r2_values.additive = mdl_add.Rsquared.Ordinary;
    r2_values.multiplicative = mdl_mult.Rsquared.Ordinary;
    r2_values.full = mdl_full.Rsquared.Ordinary;
end
