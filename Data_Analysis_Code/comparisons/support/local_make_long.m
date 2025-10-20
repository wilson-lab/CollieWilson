function L = local_make_long(T)
    animals = string(T.Properties.RowNames);
    vals = table2array(T);  % (nAnimals x 3)
    condNames = {'Dark_P1','Dark_noP1','Visual_P1'};
    % repeat animal names for 3 conditions
    Animal = categorical(repmat(animals, 3, 1));
    % stack values
    Value = vals(:,1);
    Value = [Value; vals(:,2); vals(:,3)];
    % condition labels aligned to stacked values
    Condition = categorical( ...
        [repmat(string(condNames{1}), numel(animals), 1); ...
         repmat(string(condNames{2}), numel(animals), 1); ...
         repmat(string(condNames{3}), numel(animals), 1)], ...
        string(condNames));
    % derive factorial coding
    Visual  = zeros(size(Condition));  % 0=Dark
    Arousal = zeros(size(Condition));  % 0=NoP1
    Visual(Condition=="Visual_P1") = 1;
    Arousal(Condition=="Dark_P1" | Condition=="Visual_P1") = 1;

    L = table(Animal, Value, Condition, ...
              categorical(Visual,[0 1],{'Dark','Visual'}), ...
              categorical(Arousal,[0 1],{'NoP1','P1'}), ...
              'VariableNames',{'Animal','Value','Condition','Visual','Arousal'});
    % drop missing rows (animals lacking a condition)
    L = rmmissing(L,'DataVariables','Value');
end