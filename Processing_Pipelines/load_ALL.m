% Master script for processing and reprocessing electrophysiology and 
% behavioral genetics experiments. The script handles oscillatory, motion 
% pulse, background, current injection, KIR pursuit, and depolarization 
% experiments. Users can specify whether to clear or retain previously 
% processed data for each experiment type.
%
%% Analyze oscillatory experiments
userClear = 'n'; % 'y' to clear or 'n' to keep previously preprocessed data
load_targetosc(userClear) % Load and process oscillatory target experiments

%% Analyze motion pulse experiments
clear; close all

% Motion pulse experiments to process
mopRuns = {'19_motion';'19_none';'25_motion';'25_none';'19_stationary';'25_stationary';'a02_motion'};
mopRuns = {'25_stationary'};

userClear = 'n'; % 'y' to clear or 'n' to keep previously preprocessed data
for mp = 1:length(mopRuns)
    thisRun = mopRuns{mp};
    disp(['ANALYZING: ' thisRun]) % Display current run
    load_targetpulse(thisRun, userClear); % Process current motion pulse run
end

%% Analyze background experiments
clear; close all

% Background experiments to process
bckRuns = {'19'; '25'};
userClear = 'n'; % 'y' to clear or 'n' to keep previously preprocessed data
for bg = 1:length(bckRuns)
    thisRun = bckRuns{bg};
    disp(['ANALYZING: ' thisRun]) % Display current run
    load_p1background(thisRun, userClear); % Process current background run
end

%% Analyze current injection experiments
clear; close all

% Process current injection hold experiments
% userClear = 'n'; % 'y' to clear or 'n' to keep previously preprocessed data
% load_iinj_hold(userClear); 

% Process acute current injection experiments
userClear = 'n'; % 'y' to clear or 'n' to keep previously preprocessed data
iinjRuns = {'19';'25'}; % Selected runs
for i = 1:length(iinjRuns)
    thisRun = iinjRuns{i};
    disp(['ANALYZING: ' thisRun]) % Display current run
    load_iinj_acute(thisRun,userClear);
end

%% KIR open loop pursuit behavior experiments
clear; close all

% Process KIR validation and open-loop pursuit experiments
userClear = 'n'; % 'y' to clear or 'n' to keep previously preprocessed data
%load_kir_validation(); % Process KIR validation
load_kir_openloop(userClear); % Process KIR open-loop experiments

%% KIR closed loop pursuit behavior experiments
clear; close all

% Process KIR validation and open-loop pursuit experiments
userClear = 'n'; % 'y' to clear or 'n' to keep previously preprocessed data
%load_kir_validation(); % Process KIR validation
load_kir_closedloop(userClear); % Process KIR closed-loop experiments

%% KIR menotaxis behavior experiments
clear; close all

% Process KIR menotaxis behavior experiments
userClear = 'n'; % 'y' to clear or 'n' to keep previously preprocessed data
load_kir_menotaxis(userClear)

%% DNa02 depolarization experiments
clear; close all

% Process DNa02 depolarization experiments
userClear = 'n'; % 'y' to clear or 'n' to keep previously preprocessed data
load_iinj_dna02(userClear)
