% run_ancova_on_data
% This function performs an ANCOVA (Analysis of Covariance) to evaluate the relationship 
% between object position (covariate), zero-time (dependent variable), genotype, and gain condition. 
% It tests for main effects of object position, genotype, and gain, as well as their interactions.
%
% INPUTS:
%   kirObjPos  - Object positions for KIR genotype (cell array by gain condition)
%   wtObjPos   - Object positions for WT genotype (cell array by gain condition)
%   naObjPos   - Object positions for NA genotype (cell array by gain condition)
%   kirZeroT   - Zero-time points for KIR genotype (cell array by gain condition)
%   wtZeroT    - Zero-time points for WT genotype (cell array by gain condition)
%   naZeroT    - Zero-time points for NA genotype (cell array by gain condition)
%   settings   - Structure containing settings, including gain values for each condition
%
% OUTPUTS:
%   pvalues - Structure containing p-values for the following effects:
%             - ObjPos: P-value for object position (covariate)
%             - Genotype: P-value for overall genotype effect
%             - Gain: P-value for overall gain effect
%             - ObjPos_Genotype: P-value for interaction between object position and genotype
%             - ObjPos_Gain: P-value for interaction between object position and gain
%
% CREATED: [Date] MC
%
function pvalues = run_ancova_on_data(kirObjPos, wtObjPos, naObjPos, kirZeroT, wtZeroT, naZeroT, settings)
% Initialize storage for all data
allObjPos = [];
allZeroT = [];
allGenotype = [];
allGain = [];

% Concatenate data for KIR
for c = 1:length(kirObjPos)
    allObjPos = [allObjPos; kirObjPos{c}'];
    allZeroT = [allZeroT; kirZeroT{c}'];
    allGenotype = [allGenotype; repmat({'KIR'}, length(kirObjPos{c}), 1)];
    allGain = [allGain; repmat(settings.pursuitGain(c), length(kirObjPos{c}), 1)];
end

% Concatenate data for WT
for c = 1:length(wtObjPos)
    allObjPos = [allObjPos; wtObjPos{c}'];
    allZeroT = [allZeroT; wtZeroT{c}'];
    allGenotype = [allGenotype; repmat({'WT'}, length(wtObjPos{c}), 1)];
    allGain = [allGain; repmat(settings.pursuitGain(c), length(wtObjPos{c}), 1)];
end

% Concatenate data for NA
for c = 1:length(naObjPos)
    allObjPos = [allObjPos; naObjPos{c}'];
    allZeroT = [allZeroT; naZeroT{c}'];
    allGenotype = [allGenotype; repmat({'NA'}, length(naObjPos{c}), 1)];
    allGain = [allGain; repmat(settings.pursuitGain(c), length(naObjPos{c}), 1)];
end

% Convert allGenotype and allGain to categorical for ANCOVA
allGenotype = categorical(allGenotype);
allGain = categorical(allGain);

% Perform ANCOVA
tbl = table(allObjPos, allZeroT, allGenotype, allGain, ...
    'VariableNames', {'ObjPos', 'ZeroT', 'Genotype', 'Gain'});

% Fit the linear model
ancovaModel = fitlm(tbl, 'ZeroT ~ ObjPos + Genotype + Gain + ObjPos:Genotype + ObjPos:Gain');

% Perform ANOVA on the fitted model and extract p-values for the overall terms
anovaResults = anova(ancovaModel, 'summary');

% Extract overall p-values for the main effects
pvalues.ObjPos = anovaResults.pValue(2);         % P-value for ObjPos (covariate)
pvalues.Genotype = anovaResults.pValue(3);       % P-value for overall Genotype effect
pvalues.Gain = anovaResults.pValue(4);           % P-value for overall Gain effect
pvalues.ObjPos_Genotype = anovaResults.pValue(5); % P-value for ObjPos:Genotype interaction
pvalues.ObjPos_Gain = anovaResults.pValue(6);     % P-value for ObjPos:Gain interaction

% Display overall ANOVA table (including p-values for main effects)
disp(anovaResults);

% % Plot interaction effects for Genotype and Gain
% figure;
% subplot(1, 2, 1);
% plotInteraction(ancovaModel, 'ObjPos', 'Genotype', 'predictions');
% title('Interaction between ObjPos and Genotype');
%
% subplot(1, 2, 2);
% plotInteraction(ancovaModel, 'ObjPos', 'Gain', 'predictions');
% title('Interaction between ObjPos and Gain');
end
