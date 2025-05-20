% pipeline_iinj_acute
%
% Pipeline Function
% Pulls all processed files from ALL flies in a given experiment, performs
% necessary analyses and plots accordingly. Can be used for both behavior
% only experiments and ephys experiments.
%
% INPUTS
% exptFolder - overarching experiment folder
%
% Created: 07/10/2023 by MC (adapted from motion pulse battery)
% Updated: 07/08/2024 by MC (cleaned up and simplified)
% Updated: 03/10/2025 by MC (adjusted for P1, depol only batteries)
%
function pipeline_iinj_acute(exptFolder)
%% Initialize Analysis
disp('STARTING ANALYSES FOR POOLED IINJ STIM...')

% Set filename and create necessary directories
filebase = strrep(exptFolder, ' ', '_');
folder = generateFolders(exptFolder);

% Load analysis settings
settings = processSettings();

% Define spike rate limits
sr_lim = [0 80];

% Set axis limits for trial-based and averaged data
% f = forward movement, a = angular movement, s = sideways movement
f_lim = [0 6];  
a_lim = [-75 75]; 
s_lim = [-0.8 0.8];

%% Load and Pool All Pulse Trials from Each Experiment
disp('Loading in datasets...')

% Navigate to the directory containing processed pulse trials
cd(folder.int)
pulseFiles = dir('*int.mat'); % Identify all pulse trial files
nFlies = length(pulseFiles);  % Total number of flies in the dataset

% Initialize data storage arrays
nt_P1on = 0; % Counter for flies with P1 on
nt_P1off = 0; % Counter for flies with P1 off

% Separate storage arrays for P1 on and P1 off conditions
i_v_spk_dep_P1on = []; i_v_spk_dep_P1off = [];
i_v_fwd_dep_P1on = []; i_v_fwd_dep_P1off = [];
i_v_ang_dep_P1on = []; i_v_ang_dep_P1off = [];
i_v_sid_dep_P1on = []; i_v_sid_dep_P1off = [];

i_v_spk_depH_P1on = []; i_v_fwd_depH_P1on = [];
i_v_ang_depH_P1on = []; i_v_sid_depH_P1on = [];
i_v_spk_depB_P1on = []; i_v_fwd_depB_P1on = [];
i_v_ang_depB_P1on = []; i_v_sid_depB_P1on = [];
i_v_spk_depA_P1on = []; i_v_fwd_depA_P1on = [];
i_v_ang_depA_P1on = []; i_v_sid_depA_P1on = [];

i_v_spk_hyp_P1on = []; i_v_spk_hyp_P1off = [];
i_v_fwd_hyp_P1on = []; i_v_fwd_hyp_P1off = [];
i_v_ang_hyp_P1on = []; i_v_ang_hyp_P1off = [];
i_v_sid_hyp_P1on = []; i_v_sid_hyp_P1off = [];

i_v_spk_ctr_P1on = []; i_v_spk_ctr_P1off = [];
i_v_fwd_ctr_P1on = []; i_v_fwd_ctr_P1off = [];
i_v_ang_ctr_P1on = []; i_v_ang_ctr_P1off = [];
i_v_sid_ctr_P1on = []; i_v_sid_ctr_P1off = [];

turnFreq_dep_P1on = []; turnFreq_dep_P1off = [];
turnFreq_hyp_P1on = []; turnFreq_hyp_P1off = [];
turnFreq_ctr_P1on = []; turnFreq_ctr_P1off = [];

spikertChange_dep_P1on = []; spikertChange_dep_P1off = [];
spikertChange_hyp_P1on = []; spikertChange_hyp_P1off = [];
spikertChange_ctr_P1on = []; spikertChange_ctr_P1off = [];

angularChange_dep_P1on = []; angularChange_dep_P1off = [];
angularChange_hyp_P1on = []; angularChange_hyp_P1off = [];
angularChange_ctr_P1on = []; angularChange_ctr_P1off = [];

% Loop through each fly's dataset
for nt = 1:nFlies
    % Load trial data
    thisTrial = pulseFiles(nt).name;
    disp(['Analyzing fly ' num2str(nt) '/' num2str(nFlies)])
    load(thisTrial)

    % Compute total time spent walking (converted to minutes)
    flyRunTime(nt,1) = (sum(int_forward > settings.runThreshE, 'all') / length(int_time)) * 60;

    % Determine the number of z-dimensions
    numStimConditions = size(int_iinject, 3);

    % Only include flies that meet the minimum running time threshold
    if flyRunTime(nt,1) > settings.minRunTime
        % Analyze behavior based on the number of stimulus conditions
        if numStimConditions == 2
            stimSelections = [1, 2]; % P1 off and P1 on available
        else
            stimSelections = 1; % Only P1 on available
        end

        for stimSelect = stimSelections
            % Extract relevant data for the selected condition
            thisInj = int_iinject(:,:,stimSelect);
            thisSpikert = int_spikert(:,:,stimSelect);
            thisForward = int_forward(:,:,stimSelect);
            thisAngular = int_angular(:,:,stimSelect);
            thisSideway = int_sideway(:,:,stimSelect);

            % Analyze current injection versus behavior for ONLY running timepoints
            [depolMean, hypolMean, controlMean, responseStats] = iinj_v_behavior(thisInj, thisSpikert, thisForward, thisAngular, thisSideway, int_time, settings.runThreshE);
            
            % Assign to appropriate storage arrays
            if numStimConditions == 2 && stimSelect == 1
                % P1 off
                nt_P1off = nt_P1off + 1;
                i_v_spk_dep_P1off(:,nt_P1off,1) = depolMean.spikert;
                i_v_fwd_dep_P1off(:,nt_P1off,1) = depolMean.forward;
                i_v_ang_dep_P1off(:,nt_P1off,1) = depolMean.angular;
                i_v_sid_dep_P1off(:,nt_P1off,1) = depolMean.sideway;
                i_v_spk_dep_P1off(:,nt_P1off,2) = depolMean.spikert_ipsi;
                i_v_fwd_dep_P1off(:,nt_P1off,2) = depolMean.forward_ipsi;
                i_v_ang_dep_P1off(:,nt_P1off,2) = depolMean.angular_ipsi;
                i_v_sid_dep_P1off(:,nt_P1off,2) = depolMean.sideway_ipsi;
                i_v_spk_dep_P1off(:,nt_P1off,3) = depolMean.spikert_contra;
                i_v_fwd_dep_P1off(:,nt_P1off,3) = depolMean.forward_contra;
                i_v_ang_dep_P1off(:,nt_P1off,3) = depolMean.angular_contra;
                i_v_sid_dep_P1off(:,nt_P1off,3) = depolMean.sideway_contra;
                
                i_v_spk_hyp_P1off(:,nt_P1off) = hypolMean.spikert;
                i_v_fwd_hyp_P1off(:,nt_P1off) = hypolMean.forward;
                i_v_ang_hyp_P1off(:,nt_P1off) = hypolMean.angular;
                i_v_sid_hyp_P1off(:,nt_P1off) = hypolMean.sideway;
                
                i_v_spk_ctr_P1off(:,nt_P1off) = controlMean.spikert;
                i_v_fwd_ctr_P1off(:,nt_P1off) = controlMean.forward;
                i_v_ang_ctr_P1off(:,nt_P1off) = controlMean.angular;
                i_v_sid_ctr_P1off(:,nt_P1off) = controlMean.sideway;

                turnFreq_dep_P1off(:,nt_P1off) = responseStats.dpFreq;
                turnFreq_hyp_P1off(:,nt_P1off) = responseStats.hpFreq;
                turnFreq_ctr_P1off(:,nt_P1off) = responseStats.ctrFreq;

                spikertChange_dep_P1off(:,nt_P1off) = responseStats.dpSpikert;
                spikertChange_hyp_P1off(:,nt_P1off) = responseStats.hpSpikert;
                spikertChange_ctr_P1off(:,nt_P1off) = responseStats.ctrSpikert;

                angularChange_dep_P1off(:,nt_P1off) = responseStats.dpAngular;
                angularChange_hyp_P1off(:,nt_P1off) = responseStats.hpAngular;
                angularChange_ctr_P1off(:,nt_P1off) = responseStats.ctrAngular;
            else
                % P1 on
                nt_P1on = nt_P1on + 1;
                i_v_spk_dep_P1on(:,nt_P1on,1) = depolMean.spikert;
                i_v_fwd_dep_P1on(:,nt_P1on,1) = depolMean.forward;
                i_v_ang_dep_P1on(:,nt_P1on,1) = depolMean.angular;
                i_v_sid_dep_P1on(:,nt_P1on,1) = depolMean.sideway;
                i_v_spk_dep_P1on(:,nt_P1on,2) = depolMean.spikert_ipsi;
                i_v_fwd_dep_P1on(:,nt_P1on,2) = depolMean.forward_ipsi;
                i_v_ang_dep_P1on(:,nt_P1on,2) = depolMean.angular_ipsi;
                i_v_sid_dep_P1on(:,nt_P1on,2) = depolMean.sideway_ipsi;
                i_v_spk_dep_P1on(:,nt_P1on,3) = depolMean.spikert_contra;
                i_v_fwd_dep_P1on(:,nt_P1on,3) = depolMean.forward_contra;
                i_v_ang_dep_P1on(:,nt_P1on,3) = depolMean.angular_contra;
                i_v_sid_dep_P1on(:,nt_P1on,3) = depolMean.sideway_contra;

                i_v_spk_depH_P1on(:,nt_P1on,1) = depolMean.spikertH;
                i_v_fwd_depH_P1on(:,nt_P1on,1) = depolMean.forwardH;
                i_v_ang_depH_P1on(:,nt_P1on,1) = depolMean.angularH;
                i_v_sid_depH_P1on(:,nt_P1on,1) = depolMean.sidewayH;
                i_v_spk_depB_P1on(:,nt_P1on,1) = depolMean.spikertB;
                i_v_fwd_depB_P1on(:,nt_P1on,1) = depolMean.forwardB;
                i_v_ang_depB_P1on(:,nt_P1on,1) = depolMean.angularB;
                i_v_sid_depB_P1on(:,nt_P1on,1) = depolMean.sidewayB;
                i_v_spk_depA_P1on(:,nt_P1on,1) = depolMean.spikertA;
                i_v_fwd_depA_P1on(:,nt_P1on,1) = depolMean.forwardA;
                i_v_ang_depA_P1on(:,nt_P1on,1) = depolMean.angularA;
                i_v_sid_depA_P1on(:,nt_P1on,1) = depolMean.sidewayA;
                
                i_v_spk_hyp_P1on(:,nt_P1on) = hypolMean.spikert;
                i_v_fwd_hyp_P1on(:,nt_P1on) = hypolMean.forward;
                i_v_ang_hyp_P1on(:,nt_P1on) = hypolMean.angular;
                i_v_sid_hyp_P1on(:,nt_P1on) = hypolMean.sideway;
                
                i_v_spk_ctr_P1on(:,nt_P1on) = controlMean.spikert;
                i_v_fwd_ctr_P1on(:,nt_P1on) = controlMean.forward;
                i_v_ang_ctr_P1on(:,nt_P1on) = controlMean.angular;
                i_v_sid_ctr_P1on(:,nt_P1on) = controlMean.sideway;

                turnFreq_dep_P1on(:,nt_P1on) = responseStats.dpFreq;
                turnFreq_hyp_P1on(:,nt_P1on) = responseStats.hpFreq;
                turnFreq_ctr_P1on(:,nt_P1on) = responseStats.ctrFreq;

                spikertChange_dep_P1on(:,nt_P1on) = responseStats.dpSpikert;
                spikertChange_hyp_P1on(:,nt_P1on) = responseStats.hpSpikert;
                spikertChange_ctr_P1on(:,nt_P1on) = responseStats.ctrSpikert;

                angularChange_dep_P1on(:,nt_P1on) = responseStats.dpAngular;
                angularChange_hyp_P1on(:,nt_P1on) = responseStats.hpAngular;
                angularChange_ctr_P1on(:,nt_P1on) = responseStats.ctrAngular;
            end
        end
    else
        % Exclude flies that did not meet the running threshold
        disp(['Fly ' num2str(nt) ' omitted due to poor behavior.'])
    end
end

% Convert time to milliseconds
msec_time = int_time .* 1000;
i_v_time = msec_time(1:size(i_v_spk_dep_P1on,1));

% Identify pulse onset times
pulse_idx = find(ischange(depolMean.iinject(:,1)));
pulse_time = msec_time(pulse_idx);

% Adjust angular and sideways movement to correct for turn bias ONLY for the first z-dimension
i_v_ang_dep_P1on(:,:,1) = i_v_ang_dep_P1on(:,:,1) - mean(i_v_ang_dep_P1on(1:pulse_idx(1),:,1), 1);
i_v_sid_dep_P1on(:,:,1) = i_v_sid_dep_P1on(:,:,1) - mean(i_v_sid_dep_P1on(1:pulse_idx(1),:,1), 1);
i_v_ang_hyp_P1on(:,:,1) = i_v_ang_hyp_P1on(:,:,1) - mean(i_v_ang_hyp_P1on(1:pulse_idx(1),:,1), 1);
i_v_sid_hyp_P1on(:,:,1) = i_v_sid_hyp_P1on(:,:,1) - mean(i_v_sid_hyp_P1on(1:pulse_idx(1),:,1), 1);
i_v_ang_ctr_P1on(:,:,1) = i_v_ang_ctr_P1on(:,:,1) - mean(i_v_ang_ctr_P1on(1:pulse_idx(1),:,1), 1);
i_v_sid_ctr_P1on(:,:,1) = i_v_sid_ctr_P1on(:,:,1) - mean(i_v_sid_ctr_P1on(1:pulse_idx(1),:,1), 1);

i_v_ang_dep_P1off(:,:,1) = i_v_ang_dep_P1off(:,:,1) - mean(i_v_ang_dep_P1off(1:pulse_idx(1),:,1), 1);
i_v_sid_dep_P1off(:,:,1) = i_v_sid_dep_P1off(:,:,1) - mean(i_v_sid_dep_P1off(1:pulse_idx(1),:,1), 1);
i_v_ang_hyp_P1off(:,:,1) = i_v_ang_hyp_P1off(:,:,1) - mean(i_v_ang_hyp_P1off(1:pulse_idx(1),:,1), 1);
i_v_sid_hyp_P1off(:,:,1) = i_v_sid_hyp_P1off(:,:,1) - mean(i_v_sid_hyp_P1off(1:pulse_idx(1),:,1), 1);
i_v_ang_ctr_P1off(:,:,1) = i_v_ang_ctr_P1off(:,:,1) - mean(i_v_ang_ctr_P1off(1:pulse_idx(1),:,1), 1);
i_v_sid_ctr_P1off(:,:,1) = i_v_sid_ctr_P1off(:,:,1) - mean(i_v_sid_ctr_P1off(1:pulse_idx(1),:,1), 1);

disp('All data loaded in.')

% Post-process confirm that FR was modulated sufficiently
disp('Running post-processing to confirm FR control meets criteria...')

% Remove animals where spikertChange(1,nt_t) < 10
valid_depol_p1on = spikertChange_dep_P1on(1,:) >= 10;
valid_hypol_p1on = spikertChange_hyp_P1on(1,:) <= -10;

valid_depol_p1off = spikertChange_dep_P1off(1,:) >= 10;
valid_hypol_p1off = spikertChange_hyp_P1off(1,:) <= -10;

% Apply mask to all stored arrays
i_v_spk_dep_P1on = i_v_spk_dep_P1on(:,valid_depol_p1on,:);
i_v_fwd_dep_P1on = i_v_fwd_dep_P1on(:,valid_depol_p1on,:);
i_v_ang_dep_P1on = i_v_ang_dep_P1on(:,valid_depol_p1on,:);
i_v_sid_dep_P1on = i_v_sid_dep_P1on(:,valid_depol_p1on,:);
i_v_spk_hyp_P1on = i_v_spk_hyp_P1on(:,valid_hypol_p1on,:);
i_v_fwd_hyp_P1on = i_v_fwd_hyp_P1on(:,valid_hypol_p1on,:);
i_v_ang_hyp_P1on = i_v_ang_hyp_P1on(:,valid_hypol_p1on,:);
i_v_sid_hyp_P1on = i_v_sid_hyp_P1on(:,valid_hypol_p1on,:);
i_v_spk_ctr_P1on = i_v_spk_ctr_P1on(:,valid_depol_p1on,:);
i_v_fwd_ctr_P1on = i_v_fwd_ctr_P1on(:,valid_depol_p1on,:);
i_v_ang_ctr_P1on = i_v_ang_ctr_P1on(:,valid_depol_p1on,:);
i_v_sid_ctr_P1on = i_v_sid_ctr_P1on(:,valid_depol_p1on,:);

turnFreq_dep_P1on = turnFreq_dep_P1on(valid_depol_p1on);
turnFreq_hyp_P1on = turnFreq_hyp_P1on(valid_hypol_p1on);
turnFreq_ctr_P1on = turnFreq_ctr_P1on(valid_depol_p1on);
spikertChange_dep_P1on = spikertChange_dep_P1on(valid_depol_p1on);
spikertChange_hyp_P1on = spikertChange_hyp_P1on(valid_hypol_p1on);
spikertChange_ctr_P1on = spikertChange_ctr_P1on(valid_depol_p1on);
angularChange_dep_P1on = angularChange_dep_P1on(valid_depol_p1on);
angularChange_hyp_P1on = angularChange_hyp_P1on(valid_hypol_p1on);
angularChange_ctr_P1on = angularChange_ctr_P1on(valid_depol_p1on);

i_v_spk_dep_P1off = i_v_spk_dep_P1off(:,valid_depol_p1off,:);
i_v_fwd_dep_P1off = i_v_fwd_dep_P1off(:,valid_depol_p1off,:);
i_v_ang_dep_P1off = i_v_ang_dep_P1off(:,valid_depol_p1off,:);
i_v_sid_dep_P1off = i_v_sid_dep_P1off(:,valid_depol_p1off,:);
i_v_spk_hyp_P1off = i_v_spk_hyp_P1off(:,valid_hypol_p1off,:);
i_v_fwd_hyp_P1off = i_v_fwd_hyp_P1off(:,valid_hypol_p1off,:);
i_v_ang_hyp_P1off = i_v_ang_hyp_P1off(:,valid_hypol_p1off,:);
i_v_sid_hyp_P1off = i_v_sid_hyp_P1off(:,valid_hypol_p1off,:);
i_v_spk_ctr_P1off = i_v_spk_ctr_P1off(:,valid_depol_p1off,:);
i_v_fwd_ctr_P1off = i_v_fwd_ctr_P1off(:,valid_depol_p1off,:);
i_v_ang_ctr_P1off = i_v_ang_ctr_P1off(:,valid_depol_p1off,:);
i_v_sid_ctr_P1off = i_v_sid_ctr_P1off(:,valid_depol_p1off,:);

turnFreq_dep_P1off = turnFreq_dep_P1off(valid_depol_p1off);
turnFreq_hyp_P1off = turnFreq_hyp_P1off(valid_hypol_p1off);
turnFreq_ctr_P1off = turnFreq_ctr_P1off(valid_depol_p1off);
spikertChange_dep_P1off = spikertChange_dep_P1off(valid_depol_p1off);
spikertChange_hyp_P1off = spikertChange_hyp_P1off(valid_hypol_p1off);
spikertChange_ctr_P1off = spikertChange_ctr_P1off(valid_depol_p1off);
angularChange_dep_P1off = angularChange_dep_P1off(valid_depol_p1off);
angularChange_hyp_P1off = angularChange_hyp_P1off(valid_hypol_p1off);
angularChange_ctr_P1off = angularChange_ctr_P1off(valid_depol_p1off);

n_P1on = nt_P1on- sum(~valid_depol_p1on);
n_P1off = nt_P1off - sum(~valid_hypol_p1off);

%% Plot average behavioral responses for each condition WITHOUT P1

% Mean Calculation (handling NaNs)
mean_i_v_spk_dep_P1off = mean(i_v_spk_dep_P1off(:,:,1), 2, 'omitnan');
mean_i_v_fwd_dep_P1off = mean(i_v_fwd_dep_P1off(:,:,1), 2, 'omitnan');
mean_i_v_ang_dep_P1off = mean(i_v_ang_dep_P1off(:,:,1), 2, 'omitnan');
mean_i_v_sid_dep_P1off = mean(i_v_sid_dep_P1off(:,:,1), 2, 'omitnan');
mean_i_v_spk_hyp_P1off = mean(i_v_spk_hyp_P1off(:,:,1), 2, 'omitnan');
mean_i_v_fwd_hyp_P1off = mean(i_v_fwd_hyp_P1off(:,:,1), 2, 'omitnan');
mean_i_v_ang_hyp_P1off = mean(i_v_ang_hyp_P1off(:,:,1), 2, 'omitnan');
mean_i_v_sid_hyp_P1off = mean(i_v_sid_hyp_P1off(:,:,1), 2, 'omitnan');
mean_i_v_spk_ctr_P1off = mean(i_v_spk_ctr_P1off(:,:,1), 2, 'omitnan');
mean_i_v_fwd_ctr_P1off = mean(i_v_fwd_ctr_P1off(:,:,1), 2, 'omitnan');
mean_i_v_ang_ctr_P1off = mean(i_v_ang_ctr_P1off(:,:,1), 2, 'omitnan');
mean_i_v_sid_ctr_P1off = mean(i_v_sid_ctr_P1off(:,:,1), 2, 'omitnan');

% SEM Calculation (handling NaNs)
sem_i_v_spk_dep_P1off = std(i_v_spk_dep_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_fwd_dep_P1off = std(i_v_fwd_dep_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_ang_dep_P1off = std(i_v_ang_dep_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_sid_dep_P1off = std(i_v_sid_dep_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_spk_hyp_P1off = std(i_v_spk_hyp_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_fwd_hyp_P1off = std(i_v_fwd_hyp_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_ang_hyp_P1off = std(i_v_ang_hyp_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_sid_hyp_P1off = std(i_v_sid_hyp_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_spk_ctr_P1off = std(i_v_spk_ctr_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_fwd_ctr_P1off = std(i_v_fwd_ctr_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_ang_ctr_P1off = std(i_v_ang_ctr_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_sid_ctr_P1off = std(i_v_sid_ctr_P1off(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);

figure; set(gcf,'Position',[100 100 1000 900]);
tiledlayout(4,3); % 4 rows (variables), 3 columns (conditions)

% Define labels
row_labels = {'Firing rate (spikes/s)', 'Forward velocity (mm/s)', 'Angular velocity (deg/s)', 'Sideways velocity (mm/s)'};
col_labels = {'Depolarize', 'Hyperpolarize', 'Control'};

mean_colors = {settings.spkColor,settings.velColor{1},settings.velColor{2},settings.velColor{3}};
ylim_values = {sr_lim, f_lim, a_lim, s_lim}; % Store limits for easy access

mean_data = {
    mean_i_v_spk_dep_P1off, mean_i_v_spk_hyp_P1off, mean_i_v_spk_ctr_P1off;
    mean_i_v_fwd_dep_P1off, mean_i_v_fwd_hyp_P1off, mean_i_v_fwd_ctr_P1off;
    mean_i_v_ang_dep_P1off, mean_i_v_ang_hyp_P1off, mean_i_v_ang_ctr_P1off;
    mean_i_v_sid_dep_P1off, mean_i_v_sid_hyp_P1off, mean_i_v_sid_ctr_P1off;
};

sem_data = {
    sem_i_v_spk_dep_P1off, sem_i_v_spk_hyp_P1off, sem_i_v_spk_ctr_P1off;
    sem_i_v_fwd_dep_P1off, sem_i_v_fwd_hyp_P1off, sem_i_v_fwd_ctr_P1off;
    sem_i_v_ang_dep_P1off, sem_i_v_ang_hyp_P1off, sem_i_v_ang_ctr_P1off;
    sem_i_v_sid_dep_P1off, sem_i_v_sid_hyp_P1off, sem_i_v_sid_ctr_P1off;
};

% Plotting
for r = 1:4
    for c = 1:3
        nexttile;
        
        % Plot mean and SEM
        hold on;
        fill([i_v_time, fliplr(i_v_time)], ...
             [mean_data{r, c} - sem_data{r, c}; flipud(mean_data{r, c} + sem_data{r, c})], ...
             'k', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none'); % SEM shading
        plot(i_v_time, mean_data{r, c}, 'Color',mean_colors{r}, 'LineWidth', 1.5); % Mean line
        hold off;

        % Labels and Formatting
        if r == 1
            title(col_labels{c}); % Column Titles
        end
        if c == 1
            ylabel(row_labels{r}); % Row Labels
        end
        xlabel('Time');
        xlim([min(i_v_time) max(i_v_time)]);
        ylim(ylim_values{r}); % Apply y-axis limits
        yline(0,'k')
        xline(pulse_time(1),'k')
        xline(pulse_time(2),'k')
    end
end

% Title for entire figure
sgtitle([strrep(filebase, '_', ' ') ' Thresholded, P1- (n = ' num2str(n_P1off) ')']);
% Save plot in multiple formats
cd(folder.summary);
plotname = ['run_iinj_v_behavior_none' '.png'];
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector);
set(gcf, 'renderer', 'Painters');
plotname = ['run_iinj_v_behavior_none' '.svg'];
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');


%% Plot average behavioral responses for each condition WITH P1

% Plot ALL behavior data
% Mean Calculation (handling NaNs)
mean_i_v_spk_dep_P1on = mean(i_v_spk_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_fwd_dep_P1on = mean(i_v_fwd_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_ang_dep_P1on = mean(i_v_ang_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_sid_dep_P1on = mean(i_v_sid_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_spk_hyp_P1on = mean(i_v_spk_hyp_P1on(:,:,1), 2, 'omitnan');
mean_i_v_fwd_hyp_P1on = mean(i_v_fwd_hyp_P1on(:,:,1), 2, 'omitnan');
mean_i_v_ang_hyp_P1on = mean(i_v_ang_hyp_P1on(:,:,1), 2, 'omitnan');
mean_i_v_sid_hyp_P1on = mean(i_v_sid_hyp_P1on(:,:,1), 2, 'omitnan');
mean_i_v_spk_ctr_P1on = mean(i_v_spk_ctr_P1on(:,:,1), 2, 'omitnan');
mean_i_v_fwd_ctr_P1on = mean(i_v_fwd_ctr_P1on(:,:,1), 2, 'omitnan');
mean_i_v_ang_ctr_P1on = mean(i_v_ang_ctr_P1on(:,:,1), 2, 'omitnan');
mean_i_v_sid_ctr_P1on = mean(i_v_sid_ctr_P1on(:,:,1), 2, 'omitnan');

% SEM Calculation (handling NaNs)
sem_i_v_spk_dep_P1on = std(i_v_spk_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_fwd_dep_P1on = std(i_v_fwd_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_ang_dep_P1on = std(i_v_ang_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_sid_dep_P1on = std(i_v_sid_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_spk_hyp_P1on = std(i_v_spk_hyp_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_fwd_hyp_P1on = std(i_v_fwd_hyp_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_ang_hyp_P1on = std(i_v_ang_hyp_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_sid_hyp_P1on = std(i_v_sid_hyp_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1off);
sem_i_v_spk_ctr_P1on = std(i_v_spk_ctr_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_fwd_ctr_P1on = std(i_v_fwd_ctr_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_ang_ctr_P1on = std(i_v_ang_ctr_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_sid_ctr_P1on = std(i_v_sid_ctr_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);

figure; set(gcf,'Position',[100 100 1000 900]);
tiledlayout(4,3); % 4 rows (variables), 3 columns (conditions)

% Define labels
row_labels = {'Firing rate (spikes/s)', 'Forward velocity (mm/s)', 'Angular velocity (deg/s)', 'Sideways velocity (mm/s)'};
col_labels = {'Depolarize', 'Hyperpolarize', 'Control'};

mean_colors = {settings.spkColor,settings.velColor{1},settings.velColor{2},settings.velColor{3}};
ylim_values = {sr_lim, f_lim, a_lim, s_lim}; % Store limits for easy access

mean_data = {
    mean_i_v_spk_dep_P1on, mean_i_v_spk_hyp_P1on, mean_i_v_spk_ctr_P1on;
    mean_i_v_fwd_dep_P1on, mean_i_v_fwd_hyp_P1on, mean_i_v_fwd_ctr_P1on;
    mean_i_v_ang_dep_P1on, mean_i_v_ang_hyp_P1on, mean_i_v_ang_ctr_P1on;
    mean_i_v_sid_dep_P1on, mean_i_v_sid_hyp_P1on, mean_i_v_sid_ctr_P1on;
};

sem_data = {
    sem_i_v_spk_dep_P1on, sem_i_v_spk_hyp_P1on, sem_i_v_spk_ctr_P1on;
    sem_i_v_fwd_dep_P1on, sem_i_v_fwd_hyp_P1on, sem_i_v_fwd_ctr_P1on;
    sem_i_v_ang_dep_P1on, sem_i_v_ang_hyp_P1on, sem_i_v_ang_ctr_P1on;
    sem_i_v_sid_dep_P1on, sem_i_v_sid_hyp_P1on, sem_i_v_sid_ctr_P1on;
};

% Plotting
for r = 1:4
    for c = 1:3
        nexttile;
        
        % Plot mean and SEM
        hold on;
        fill([i_v_time, fliplr(i_v_time)], ...
             [mean_data{r, c} - sem_data{r, c}; flipud(mean_data{r, c} + sem_data{r, c})], ...
             'k', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none'); % SEM shading
        plot(i_v_time, mean_data{r, c}, 'Color',mean_colors{r}, 'LineWidth', 1.5); % Mean line
        hold off;

        % Labels and Formatting
        if r == 1
            title(col_labels{c}); % Column Titles
        end
        if c == 1
            ylabel(row_labels{r}); % Row Labels
        end
        xlabel('Time');
        xlim([min(i_v_time) max(i_v_time)]);
        ylim(ylim_values{r}); % Apply y-axis limits
        yline(0,'k')
        xline(pulse_time(1),'k')
        xline(pulse_time(2),'k')
    end
end

% Title for entire figure
sgtitle([strrep(filebase, '_', ' ') ' Thresholded, P1+ (n = ' num2str(n_P1on) '/' num2str(n_P1off) ')']);
% Save plot in multiple formats
cd(folder.summary);
plotname = ['run_iinj_v_behavior_p1' '.png'];
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector);
set(gcf, 'renderer', 'Painters');
plotname = ['run_iinj_v_behavior_p1' '.svg'];
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

%% Run linear mixed-effects ANOVA on Δ forward and angular velocity (P1on)

% Set time windows (baseline: 0–500ms, response: 1000–1500ms)
baseline_idx = i_v_time >= 0 & i_v_time < 500;
response_idx = i_v_time >= 1000 & i_v_time < 1500;

% Initialize storage
deltaFwd = [];
deltaAng = [];
flyID = [];
condition = {};

% Helper to extract Δ and label
get_deltas = @(data) mean(data(response_idx,:),1,'omitnan') - mean(data(baseline_idx,:),1,'omitnan');

% ----- DEPOL -----
if exist('i_v_fwd_dep_P1on','var')
    n_depol = size(i_v_fwd_dep_P1on, 2);
    deltaFwd = [deltaFwd; get_deltas(i_v_fwd_dep_P1on(:,:,1))'];
    deltaAng = [deltaAng; get_deltas(i_v_ang_dep_P1on(:,:,1))'];
    flyID    = [flyID; (1:n_depol)'];
    condition = [condition; repmat({'Depol'}, n_depol, 1)];
end

% ----- HYPERPOL -----
if exist('i_v_fwd_hyp_P1on','var')
    n_hyper = size(i_v_fwd_hyp_P1on, 2);
    deltaFwd = [deltaFwd; get_deltas(i_v_fwd_hyp_P1on)'];
    deltaAng = [deltaAng; get_deltas(i_v_ang_hyp_P1on)'];
    flyID    = [flyID; (1:n_hyper)'];
    condition = [condition; repmat({'Hyperpol'}, n_hyper, 1)];
end

% ----- CONTROL -----
if exist('i_v_fwd_ctr_P1on','var')
    n_ctrl = size(i_v_fwd_ctr_P1on, 2);
    deltaFwd = [deltaFwd; get_deltas(i_v_fwd_ctr_P1on)'];
    deltaAng = [deltaAng; get_deltas(i_v_ang_ctr_P1on)'];
    flyID    = [flyID; (1:n_ctrl)'];
    condition = [condition; repmat({'Control'}, n_ctrl, 1)];
end

% Build table
T = table(deltaFwd, deltaAng, categorical(condition), flyID, ...
    'VariableNames', {'DeltaFwd', 'DeltaAng', 'Condition', 'Fly'});

% Fit linear mixed-effects models
lme_fwd = fitlme(T, 'DeltaFwd ~ Condition + (1|Fly)');
lme_ang = fitlme(T, 'DeltaAng ~ Condition + (1|Fly)');

% Display results
disp('LME on Δ Forward Velocity:');
disp(anova(lme_fwd));

disp('LME on Δ Angular Velocity:');
disp(anova(lme_ang));

% Get fixed effects and their covariance
fe = fixedEffects(lme_fwd);
covFE = lme_fwd.CoefficientCovariance;
names = lme_fwd.CoefficientNames;
levels = categories(T.Condition);

% Get actual names in order: base + contrast terms
% Intercept is baseline (e.g., Control), contrasts are relative to that
baseline = levels{1}; % Usually the first alphabetical, like Control
all_contrasts = [];

% Build pairwise contrast table manually
for i = 1:length(levels)
    for j = i+1:length(levels)
        g1 = levels{i}; g2 = levels{j};

        % Contrast vector
        c = zeros(size(fe));
        % Group 1 is baseline
        if ~strcmp(g1, baseline)
            c(strcmp(names, ['Condition_' g1])) = -1;
        end
        if ~strcmp(g2, baseline)
            c(strcmp(names, ['Condition_' g2])) = 1;
        end

        % Compute stats
        est = c' * fe;
        se = sqrt(c' * covFE * c);
        tval = est / se;
        df = lme_fwd.DFE;
        pval = 2 * tcdf(-abs(tval), df);

        all_contrasts = [all_contrasts; {g1, g2, est, se, tval, df, pval}];
    end
end

% Turn into table
posthoc_fwd = cell2table(all_contrasts, ...
    'VariableNames', {'Group1','Group2','Estimate','SE','tStat','DF','pValue'});

% Tukey correction using multcompare-like approach
k = height(posthoc_fwd);
posthoc_fwd.pAdj_Tukey = min(1, posthoc_fwd.pValue * k); % Conservative Tukey-style correction

disp('Posthoc (Tukey-corrected) Δ Forward Velocity');
disp(posthoc_fwd);

% Repeat for angular
fe = fixedEffects(lme_ang);
covFE = lme_ang.CoefficientCovariance;
names = lme_ang.CoefficientNames;
all_contrasts = [];

for i = 1:length(levels)
    for j = i+1:length(levels)
        g1 = levels{i}; g2 = levels{j};

        c = zeros(size(fe));
        if ~strcmp(g1, baseline)
            c(strcmp(names, ['Condition_' g1])) = -1;
        end
        if ~strcmp(g2, baseline)
            c(strcmp(names, ['Condition_' g2])) = 1;
        end

        est = c' * fe;
        se = sqrt(c' * covFE * c);
        tval = est / se;
        df = lme_ang.DFE;
        pval = 2 * tcdf(-abs(tval), df);

        all_contrasts = [all_contrasts; {g1, g2, est, se, tval, df, pval}];
    end
end

posthoc_ang = cell2table(all_contrasts, ...
    'VariableNames', {'Group1','Group2','Estimate','SE','tStat','DF','pValue'});
posthoc_ang.pAdj_Tukey = min(1, posthoc_ang.pValue * k);

disp('Posthoc (Tukey-corrected) Δ Angular Velocity');
disp(posthoc_ang);


%% Plot average behavioral responses for each condition WITH P1, binned by pre-pulse turning

% Mean Calculation (handling NaNs)
mean_i_v_spk_dep_P1on = mean(i_v_spk_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_fwd_dep_P1on = mean(i_v_fwd_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_ang_dep_P1on = mean(i_v_ang_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_sid_dep_P1on = mean(i_v_sid_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_spk_dep_P1on_ipsi = mean(i_v_spk_dep_P1on(:,:,2), 2, 'omitnan');
mean_i_v_fwd_dep_P1on_ipsi = mean(i_v_fwd_dep_P1on(:,:,2), 2, 'omitnan');
mean_i_v_ang_dep_P1on_ipsi = mean(i_v_ang_dep_P1on(:,:,2), 2, 'omitnan');
mean_i_v_sid_dep_P1on_ipsi = mean(i_v_sid_dep_P1on(:,:,2), 2, 'omitnan');
mean_i_v_spk_dep_P1on_ctra = mean(i_v_spk_dep_P1on(:,:,3), 2, 'omitnan');
mean_i_v_fwd_dep_P1on_ctra = mean(i_v_fwd_dep_P1on(:,:,3), 2, 'omitnan');
mean_i_v_ang_dep_P1on_ctra = mean(i_v_ang_dep_P1on(:,:,3), 2, 'omitnan');
mean_i_v_sid_dep_P1on_ctra = mean(i_v_sid_dep_P1on(:,:,3), 2, 'omitnan');

% SEM Calculation (handling NaNs)
sem_i_v_spk_dep_P1on = std(i_v_spk_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_fwd_dep_P1on = std(i_v_fwd_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_ang_dep_P1on = std(i_v_ang_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_sid_dep_P1on = std(i_v_sid_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_spk_dep_P1on_ipsi = std(i_v_spk_dep_P1on(:,:,2), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_fwd_dep_P1on_ipsi = std(i_v_fwd_dep_P1on(:,:,2), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_ang_dep_P1on_ipsi = std(i_v_ang_dep_P1on(:,:,2), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_sid_dep_P1on_ipsi = std(i_v_sid_dep_P1on(:,:,2), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_spk_dep_P1on_ctra = std(i_v_spk_dep_P1on(:,:,3), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_fwd_dep_P1on_ctra = std(i_v_fwd_dep_P1on(:,:,3), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_ang_dep_P1on_ctra = std(i_v_ang_dep_P1on(:,:,3), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_sid_dep_P1on_ctra = std(i_v_sid_dep_P1on(:,:,3), 0, 2, 'omitnan') ./ sqrt(n_P1on);

figure; set(gcf,'Position',[100 100 1000 900]);
tiledlayout(4,3); % 4 rows (variables), 3 columns (conditions)

% Define labels
row_labels = {'Firing rate (spikes/s)', 'Forward velocity (mm/s)', 'Angular velocity (deg/s)', 'Sideways velocity (mm/s)'};
col_labels = {'All', 'Ipsi Turn Prior', 'Contra Turn Prior'};

mean_colors = {settings.spkColor,settings.velColor{1},settings.velColor{2},settings.velColor{3}};
ylim_values = {sr_lim, f_lim, [-150 150], [-2 2]}; % Store limits for easy access

mean_data = {
    mean_i_v_spk_dep_P1on, mean_i_v_spk_dep_P1on_ipsi, mean_i_v_spk_dep_P1on_ctra;
    mean_i_v_fwd_dep_P1on, mean_i_v_fwd_dep_P1on_ipsi, mean_i_v_fwd_dep_P1on_ctra;
    mean_i_v_ang_dep_P1on, mean_i_v_ang_dep_P1on_ipsi, mean_i_v_ang_dep_P1on_ctra;
    mean_i_v_sid_dep_P1on, mean_i_v_sid_dep_P1on_ipsi, mean_i_v_sid_dep_P1on_ctra;
};

sem_data = {
    sem_i_v_spk_dep_P1on, sem_i_v_spk_dep_P1on_ipsi, sem_i_v_spk_dep_P1on_ctra;
    sem_i_v_fwd_dep_P1on, sem_i_v_fwd_dep_P1on_ipsi, sem_i_v_fwd_dep_P1on_ctra;
    sem_i_v_ang_dep_P1on, sem_i_v_ang_dep_P1on_ipsi, sem_i_v_ang_dep_P1on_ctra;
    sem_i_v_sid_dep_P1on, sem_i_v_sid_dep_P1on_ipsi, sem_i_v_sid_dep_P1on_ctra;
};

% Plotting
for r = 1:4
    for c = 1:3
        nexttile;
        
        % Plot mean and SEM
        hold on;
        fill([i_v_time, fliplr(i_v_time)], ...
             [mean_data{r, c} - sem_data{r, c}; flipud(mean_data{r, c} + sem_data{r, c})], ...
             'k', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none'); % SEM shading
        plot(i_v_time, mean_data{r, c}, 'Color',mean_colors{r}, 'LineWidth', 1.5); % Mean line
        hold off;

        % Labels and Formatting
        if r == 1
            title(col_labels{c}); % Column Titles
        end
        if c == 1
            ylabel(row_labels{r}); % Row Labels
        end
        xlabel('Time');
        xlim([min(i_v_time) max(i_v_time)]);
        ylim(ylim_values{r}); % Apply y-axis limits
        yline(0,'k')
        xline(pulse_time(1),'k')
        xline(pulse_time(2),'k')
    end
end

% Title for entire figure
sgtitle([strrep(filebase, '_', ' ') ' Thresholded, P1+ (n = ' num2str(n_P1on) ')']);
% Save plot in multiple formats
cd(folder.summary);
plotname = ['run_iinj_v_behavior_p1_prebined' '.png'];
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector);
set(gcf, 'renderer', 'Painters');
plotname = ['run_iinj_v_behavior_p1_prebinned' '.svg'];
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');


%% Plot average behavioral responses for each condition WITH P1, binned by forward velocity
% Mean Calculation (handling NaNs)
mean_i_v_spk_dep_P1on = mean(i_v_spk_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_fwd_dep_P1on = mean(i_v_fwd_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_ang_dep_P1on = mean(i_v_ang_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_sid_dep_P1on = mean(i_v_sid_dep_P1on(:,:,1), 2, 'omitnan');
mean_i_v_spk_dep_P1on_H = mean(i_v_spk_depH_P1on, 2, 'omitnan');
mean_i_v_fwd_dep_P1on_H = mean(i_v_fwd_depH_P1on, 2, 'omitnan');
mean_i_v_ang_dep_P1on_H = mean(i_v_ang_depH_P1on, 2, 'omitnan');
mean_i_v_sid_dep_P1on_H = mean(i_v_sid_depH_P1on, 2, 'omitnan');
mean_i_v_spk_dep_P1on_B = mean(i_v_spk_depB_P1on, 2, 'omitnan');
mean_i_v_fwd_dep_P1on_B = mean(i_v_fwd_depB_P1on, 2, 'omitnan');
mean_i_v_ang_dep_P1on_B = mean(i_v_ang_depB_P1on, 2, 'omitnan');
mean_i_v_sid_dep_P1on_B = mean(i_v_sid_depB_P1on, 2, 'omitnan');
mean_i_v_spk_dep_P1on_A = mean(i_v_spk_depA_P1on, 2, 'omitnan');
mean_i_v_fwd_dep_P1on_A = mean(i_v_fwd_depA_P1on, 2, 'omitnan');
mean_i_v_ang_dep_P1on_A = mean(i_v_ang_depA_P1on, 2, 'omitnan');
mean_i_v_sid_dep_P1on_A = mean(i_v_sid_depA_P1on, 2, 'omitnan');

% SEM Calculation (handling NaNs)
sem_i_v_spk_dep_P1on = std(i_v_spk_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_fwd_dep_P1on = std(i_v_fwd_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_ang_dep_P1on = std(i_v_ang_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_sid_dep_P1on = std(i_v_sid_dep_P1on(:,:,1), 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_spk_dep_P1on_H = std(i_v_spk_depH_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_fwd_dep_P1on_H = std(i_v_fwd_depH_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_ang_dep_P1on_H = std(i_v_ang_depH_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_sid_dep_P1on_H = std(i_v_sid_depH_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_spk_dep_P1on_B = std(i_v_spk_depB_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_fwd_dep_P1on_B = std(i_v_fwd_depB_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_ang_dep_P1on_B = std(i_v_ang_depB_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_sid_dep_P1on_B = std(i_v_sid_depB_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_spk_dep_P1on_A = std(i_v_spk_depA_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_fwd_dep_P1on_A = std(i_v_fwd_depA_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_ang_dep_P1on_A = std(i_v_ang_depA_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);
sem_i_v_sid_dep_P1on_A = std(i_v_sid_depA_P1on, 0, 2, 'omitnan') ./ sqrt(n_P1on);

figure; set(gcf,'Position',[100 100 1200 900]);
tiledlayout(4,4); % 4 rows (variables), 4 columns (conditions)

% Define labels
row_labels = {'Firing rate (spikes/s)', 'Forward velocity (mm/s)', 'Angular velocity (deg/s)', 'Sideways velocity (mm/s)'};
col_labels = {'0.3 Inclusion', '0.8 Inclusion', 'Fwd below Median', 'Fwd above Median'};

mean_colors = {settings.spkColor,settings.velColor{1},settings.velColor{2},settings.velColor{3}};
ylim_values = {sr_lim, f_lim, a_lim, s_lim}; % Store limits for easy access

mean_data = {
    mean_i_v_spk_dep_P1on, mean_i_v_spk_dep_P1on_H, mean_i_v_spk_dep_P1on_B, mean_i_v_spk_dep_P1on_A;
    mean_i_v_fwd_dep_P1on, mean_i_v_fwd_dep_P1on_H, mean_i_v_fwd_dep_P1on_B, mean_i_v_fwd_dep_P1on_A;
    mean_i_v_ang_dep_P1on, mean_i_v_ang_dep_P1on_H, mean_i_v_ang_dep_P1on_B, mean_i_v_ang_dep_P1on_A;
    mean_i_v_sid_dep_P1on, mean_i_v_sid_dep_P1on_H, mean_i_v_sid_dep_P1on_B, mean_i_v_sid_dep_P1on_A;
};

sem_data = {
    sem_i_v_spk_dep_P1on, sem_i_v_spk_dep_P1on_H, sem_i_v_spk_dep_P1on_B, sem_i_v_spk_dep_P1on_A;
    sem_i_v_fwd_dep_P1on, sem_i_v_fwd_dep_P1on_H, sem_i_v_fwd_dep_P1on_B, sem_i_v_fwd_dep_P1on_A;
    sem_i_v_ang_dep_P1on, sem_i_v_ang_dep_P1on_H, sem_i_v_ang_dep_P1on_B, sem_i_v_ang_dep_P1on_A;
    sem_i_v_sid_dep_P1on, sem_i_v_sid_dep_P1on_H, sem_i_v_sid_dep_P1on_B, sem_i_v_sid_dep_P1on_A;
};

% Plotting
for r = 1:4
    for c = 1:4
        nexttile;
        
        % Plot mean and SEM
        hold on;
        fill([i_v_time, fliplr(i_v_time)], ...
             [mean_data{r, c} - sem_data{r, c}; flipud(mean_data{r, c} + sem_data{r, c})], ...
             'k', 'FaceAlpha', settings.semAlpha, 'EdgeColor', 'none'); % SEM shading
        plot(i_v_time, mean_data{r, c}, 'Color',mean_colors{r}, 'LineWidth', 1.5); % Mean line
        hold off;

        % Labels and Formatting
        if r == 1
            title(col_labels{c}); % Column Titles
        end
        if c == 1
            ylabel(row_labels{r}); % Row Labels
        end
        xlabel('Time');
        xlim([min(i_v_time) max(i_v_time)]);
        ylim(ylim_values{r}); % Apply y-axis limits
        yline(0,'k')
        xline(pulse_time(1),'k')
        xline(pulse_time(2),'k')
    end
end

% Title for entire figure
sgtitle([strrep(filebase, '_', ' ') ' Thresholded, P1+ (n = ' num2str(n_P1on) ')']);
% Save plot in multiple formats
cd(folder.summary);
plotname = ['run_iinj_v_behavior_p1_fwdbined' '.png'];
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector);
set(gcf, 'renderer', 'Painters');
plotname = ['run_iinj_v_behavior_p1_fwdbinned' '.svg'];
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

%% Plot scatter between changes in FR v changes in turning WITHOUT P1

% Define conditions
spikertChange = {spikertChange_dep_P1off, spikertChange_hyp_P1off, spikertChange_ctr_P1off};
angularChange = {angularChange_dep_P1off, angularChange_hyp_P1off, angularChange_ctr_P1off};
conditions = {'Depolarize', 'Hyperpolarize', 'Control'};

% Define axis limits
spike_lim = [-100 100];
ang_lim = [-100 100];

% Initialize figure
figure; set(gcf, 'Position', [100 100 1500 600]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot for each condition
for i = 1:3
    nexttile;
    
    % Extract data and remove NaNs
    valid_idx = ~isnan(spikertChange{i}) & ~isnan(angularChange{i});
    x_data = spikertChange{i}(valid_idx);
    y_data = angularChange{i}(valid_idx);
    
    % Scatter plot
    scatter(x_data, y_data, 20, 'k', 'filled'); 
    hold on;

    % Labels and title
    ylabel('Change in Angular Velocity (deg/s)');
    xlabel('Change in Spike Rate');
    xlim(spike_lim);
    ylim(ang_lim);
    title(conditions{i});
    
    % Formatting
    axis square;
    grid on;
end

% Title for entire figure
sgtitle([strrep(filebase, '_', ' ') ' Thresholded, P1- (n = ' num2str(n_P1off) ')']);

% Save plot in multiple formats
cd(folder.summary);
plotname = 'run_iinj_v_behavior_scatter_none.png';
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector);
set(gcf, 'renderer', 'Painters');
plotname = 'run_iinj_v_behavior_scatter_none.svg';
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');


%% Plot scatter between changes in FR v changes in turning WITH P1

% Define conditions
spikertChange = {spikertChange_dep_P1on, spikertChange_hyp_P1on, spikertChange_ctr_P1on};
angularChange = {angularChange_dep_P1on, angularChange_hyp_P1on, angularChange_ctr_P1on};
conditions = {'Depolarize', 'Hyperpolarize', 'Control'};

% Define axis limits
spike_lim = [-100 100];
ang_lim = [-100 100];

% Initialize figure
figure; set(gcf, 'Position', [100 100 1500 600]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot for each condition
for i = 1:3
    nexttile;
    
    % Extract data and remove NaNs
    valid_idx = ~isnan(spikertChange{i}) & ~isnan(angularChange{i});
    x_data = spikertChange{i}(valid_idx);
    y_data = angularChange{i}(valid_idx);
    
    % Scatter plot
    scatter(x_data, y_data, 20, 'k', 'filled'); 
    hold on;

    % Labels and title
    ylabel('Change in Angular Velocity (deg/s)');
    xlabel('Change in Spike Rate');
    xlim(spike_lim);
    ylim(ang_lim);
    title(conditions{i});
    
    % Formatting
    axis square;
    grid on;
end

% Title for entire figure
sgtitle([strrep(filebase, '_', ' ') ' Thresholded, P1+ (n = ' num2str(n_P1on) '/' num2str(n_P1off) ')']);

% Save plot in multiple formats
cd(folder.summary);
plotname = 'run_iinj_v_behavior_scatter.png';
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector);
set(gcf, 'renderer', 'Painters');
plotname = 'run_iinj_v_behavior_scatter.svg';
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

%% Plot Turn Frequencies WITHOUT P1

% Define conditions and corresponding data
spikertChange = {spikertChange_dep_P1off, spikertChange_hyp_P1off, spikertChange_ctr_P1off};
turnFreq = {turnFreq_dep_P1off, turnFreq_hyp_P1off, turnFreq_ctr_P1off};
conditions = {'Depolarize', 'Hyperpolarize', 'Control'};

% Define axis limits
spike_lim = [-100 100];  % X-axis: Change in spike rate
turn_lim = [0 1];        % Y-axis: Turn frequency

% Initialize figure
figure; set(gcf, 'Position', [100 100 1500 420]);
tiledlayout(1,4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot for each condition
for i = 1:3
    nexttile;
    
    % Extract data and remove NaNs
    valid_idx = ~isnan(spikertChange{i}) & ~isnan(turnFreq{i});
    x_data = spikertChange{i}(valid_idx);
    y_data = turnFreq{i}(valid_idx);
    
    % Scatter plot
    scatter(x_data, y_data, 20, 'k', 'filled'); 
    hold on;
    
    % Labels and title
    ylabel('Turn Frequency');
    xlabel('Change in Spike Rate');
    xlim(spike_lim);
    ylim(turn_lim);
    yline(0.5, 'r', 'LineWidth', 1.5);
    title(conditions{i});
    
    % Formatting
    axis square;
    grid on;
end

% Plot all conditions together
nexttile;
hold on;
for i = 1:3
    % Extract valid data points
    valid_idx = ~isnan(turnFreq{i});
    y_data = turnFreq{i}(valid_idx);
    
    % Create x-axis positions for scatter points (Jitter to avoid overlap)
    x_data = i + (rand(size(y_data)) - 0.5) * 0.1; % Add jitter for visibility
    
    % Scatter plot
    scatter(x_data, y_data, 20, 'k', 'filled');
end
hold off;

% Formatting
ylabel('Turn Frequency');
xlim([0.5 3.5]);
ylim(turn_lim);
yline(0.5, 'r', 'LineWidth', 1.5);
grid on;
xticks(1:3);
xticklabels(conditions);
xtickangle(60);

% Title for entire figure
sgtitle([strrep(filebase, '_', ' ') ' Thresholded, P1- (n = ' num2str(n_P1off) ')']);

% Save plot in multiple formats
cd(folder.summary);
plotname = 'run_iinj_v_behavior_freq_none.png';
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector);
set(gcf, 'renderer', 'Painters');
plotname = 'run_iinj_v_behavior_freq_none.svg';
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

%% Plot Turn Frequencies WITH P1

% Define conditions and corresponding data
spikertChange = {spikertChange_dep_P1on, spikertChange_hyp_P1on, spikertChange_ctr_P1on};
turnFreq = {turnFreq_dep_P1on, turnFreq_hyp_P1on, turnFreq_ctr_P1on};
conditions = {'Depolarize', 'Hyperpolarize', 'Control'};

% Define axis limits
spike_lim = [-100 100];  % X-axis: Change in spike rate
turn_lim = [0 1];        % Y-axis: Turn frequency

% Initialize figure
figure; set(gcf, 'Position', [100 100 1500 420]);
tiledlayout(1,4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot for each condition
for i = 1:3
    nexttile;
    
    % Extract data and remove NaNs
    valid_idx = ~isnan(spikertChange{i}) & ~isnan(turnFreq{i});
    x_data = spikertChange{i}(valid_idx);
    y_data = turnFreq{i}(valid_idx);
    
    % Scatter plot
    scatter(x_data, y_data, 20, 'k', 'filled'); 
    hold on;
    
    % Labels and title
    ylabel('Turn Frequency');
    xlabel('Change in Spike Rate');
    xlim(spike_lim);
    ylim(turn_lim);
    yline(0.5, 'r', 'LineWidth', 1.5);
    title(conditions{i});
    
    % Formatting
    axis square;
    grid on;
end

% Plot all conditions together
nexttile;
hold on;
for i = 1:3
    % Extract valid data points
    valid_idx = ~isnan(turnFreq{i});
    y_data = turnFreq{i}(valid_idx);
    
    % Create x-axis positions for scatter points (Jitter to avoid overlap)
    x_data = i + (rand(size(y_data)) - 0.5) * 0.1; % Add jitter for visibility
    
    % Scatter plot
    scatter(x_data, y_data, 20, 'k', 'filled');
end
hold off;

% Formatting
ylabel('Turn Frequency');
xlim([0.5 3.5]);
ylim(turn_lim);
yline(0.5, 'r', 'LineWidth', 1.5);
grid on;
xticks(1:3);
xticklabels(conditions);
xtickangle(60);

% Title for entire figure
sgtitle([strrep(filebase, '_', ' ') ' Thresholded, P1+ (n = ' num2str(n_P1on) '/' num2str(n_P1off) ')']);

% Save plot in multiple formats
cd(folder.summary);
plotname = 'run_iinj_v_behavior_freq.png';
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

% Save vectorized plot
cd(folder.vector);
set(gcf, 'renderer', 'Painters');
plotname = 'run_iinj_v_behavior_freq.svg';
saveas(gcf, plotname);
copyfile(plotname, folder.dropbox, 'f');

%% end
disp('ALL ANALYSES COMPLETE.')

end

