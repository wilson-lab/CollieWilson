function run_lme_and_tests(D, labelStr)
    fprintf('\n==== %s ====\n', labelStr);

    % Mixed model: Value ~ Visual + Arousal + (1|Animal)
    % (no interaction term; that cell is structurally missing in your design)
    lme = fitlme(D, 'Value ~ Visual + Arousal + (1|Animal)', ...
                    'DummyVarCoding','effects');

    disp('--- Fixed-effects ANOVA ---');
    disp(anova(lme,'DFMethod','Satterthwaite'));

    % Planned contrasts
    % With effects coding, the main-effect coefficients correspond to
    % mean-centered contrasts. But given the design:
    %  - Arousal effect is identified from Dark only (Dark+P1 vs Dark−P1).
    %  - Visual effect is identified under +P1 (Visual+P1 vs Dark+P1).
    %
    % We'll construct explicit contrasts on the fitted marginal means
    % by re-fitting with 'reference' coding for easier coef interpretation.
    lme_ref = fitlme(D, 'Value ~ Visual + Arousal + (1|Animal)', ...
                        'DummyVarCoding','reference');

    % Coefficient order (reference coding): (Intercept), Visual_Visual, Arousal_P1
    betaNames = lme_ref.CoefficientNames;
    ix_visual  = find(strcmp(betaNames,'Visual_Visual'));
    ix_arousal = find(strcmp(betaNames,'Arousal_P1'));

    % Contrast #1: Arousal effect in Dark (Dark+P1 - Dark−P1)
    % Under reference coding, Arousal_P1 represents (P1 - NoP1) at the reference Visual level (Dark).
    C1 = zeros(1, numel(betaNames)); C1(ix_arousal) = 1;
    [p1,F1,df1] = coefTest(lme_ref, C1);
    fprintf('Arousal effect in Dark (Dark+P1 vs Dark−P1): F(1,%g) = %.3f, p = %.4g\n', df1, F1, p1);

    % Contrast #2: Visual effect under +P1 (Visual+P1 - Dark+P1)
    % Visual_Visual represents (Visual - Dark) at the reference Arousal level (NoP1).
    % But we only have +P1 data for Visual, so we test the Visual main effect directly,
    % which in this design corresponds to (Visual+P1 - Dark+P1).
    C2 = zeros(1, numel(betaNames)); C2(ix_visual) = 1;
    [p2,F2,df2] = coefTest(lme_ref, C2);
    fprintf('Visual effect under +P1 (Visual+P1 vs Dark+P1): F(1,%g) = %.3f, p = %.4g\n', df2, F2, p2);

    % Optional: show fixed effects table
    disp('--- Fixed effects (reference coding) ---');
    disp(lme_ref.Coefficients);
end
