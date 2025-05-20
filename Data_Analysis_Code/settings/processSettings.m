% processSettings
%
% Initializes and returns a structure containing shared settings across 
% various processing pipelines related to electrophysiology, behavior, 
% vision, genetics, and visualization. These settings are used to ensure 
% consistency and standardization in data processing and analysis.
%
% OUTPUT:
% - settings : Structure containing various settings for data processing,
%              including parameters for electrophysiology, behavioral 
%              analysis, visual processing, and plotting.
%
% SETTINGS INCLUDE:
% - Electrophysiology parameters such as liquid junction potential and 
%   current injection settings.
% - General behavior thresholds for running speed, minimum run time, and 
%   lag estimates for visual-motor behavior.
% - Genetic parameters, including possible genotypes and fixation times.
% - Menotaxis parameters, such as maximum stop time and window sizes for 
%   jump analysis.
% - Formatting and plotting specifications, including colors, line weights, 
%   and labels for different measurements (e.g., velocity and acceleration).
%
% The function establishes critical parameters used throughout the analysis 
% to ensure accurate and reproducible results across experiments.
%
function settings = processSettings()
%% general electrophysiology
% set liquid junction potential subtraction for voltage traces
% source: Gouwens, N. W. & Wilson, R.I. Signal propagation in Drosophila central neurons. J.Neurosci. 29, 6239-6249
settings.ljp = 13; %mV

% current injection
settings.acuteLabel = {'Depol';'Hyperol';'Ctr';'Depol +P1';'Hyperol -P1';'Ctr +P1'};
settings.holdLabel = {'No P1 No iHold';'No P1 with iHold';'P1 No iHold';'P1 with iHold'};

%% general vision
% binocular overlap region
settings.binoc = 15; %deg

%% general behavior
% min run threshold
settings.runThreshE = 2; %mm/s, ephys
settings.runThreshB = 5; %mm/s, behavior
% min time spent running
settings.minRunTime = 60; %s

% min peak prominence for xcorr
settings.minXCorrProm = 0.15; % r value, for fr v vel
settings.minXCorrPromVm = 0.075; % r value, for vm v vel

% start/stop transition time to exclude
settings.ssExclude = 0.250; %s
settings.ssThresh = 0.5; %m/s used to find when fwd changes

% lag estimate for visual-motor behavior
settings.visuomotorLag = 0.170; %s

% lag estimates for spikerate v each directional velocity
settings.fwdLag = 0; %s
settings.angLag = 0.200; %s
settings.sidLag = 0.100; %s

%% behavioral genetics
% set possible genotypes
settings.geneLabel = {'KIR';'WT';'NA'};

% for closed-loop experiments
% set min fixation time
settings.minFixationTime = 20; %s
% set gain options
settings.pursuitGain = [95, 110, 140, 170, 200]; %yaw gain
%settings.pursuitGain = [65, 80, 95, 110, 140, 170, 200]; %yaw gain
% set jump size
settings.pursuitJump = [60, -60]; %deg
settings.HDBins = 60; % set the number of heading bins

% setpoint bin size
settings.spBin = 5; %deg
% +/- how close the target can be to be considered at the setpoint HD
settings.spHD = 10; %deg

% set what is considered the front FOV (e.g., likely pursuit)
settings.fFOV = 50; %+/- deg

% set bootstrap parameters
settings.nBootstrap = 10000;
settings.pValSig = 0.05;

%% menotaxis
% max time spent stopped when menotaxing
settings.maxTimeStop = 4; %s
% set jump window size for pulling data before/after
settings.mtaxPreWin = 15; %sec, window before jump
settings.mtaxPstWin = 20; %sec, window after jump
settings.minMtaxCor = 2; %min number of jump corrections to be included in analysis

% NOTE: the following parameters are hardcoded in python socket-client code
% set voltages that were used to denote jump triggers
settings.voltages_inuse = 1:6;

% set jump sizes corresponding to said jump triggers
settings.jumps_inuse = [30, -30, 60, -60, 180, -180]; %deg

%% general formatting
settings.trialColor = [0.7 0.7 0.7];
settings.lwTri = 0.5; %lineweight
settings.lwAvg = 1.5; %lineweight
settings.semAlpha = 0.1; %sem opacity
settings.letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

%% plotting variables
% spikerate
settings.spkLabel = 'Firing Rate (spk/s)';
settings.spkColor = '#77AC30';

% velocity
settings.behaviorGroup = {'All'; 'Quiescent'; 'Pursuit'};
settings.velLabel = {'Forward Velocity (mm/s)'; 'Angular Velocity (deg/s)'; 'Sideways Velocity (mm/s)'}; %velocity names
settings.accLabel = {'Forward Acceleration (mm/s^2)'; 'Angular Acceleration (deg/s^2)'; 'Sideways Acceleration (mm/s^2)'}; %velocity names
settings.velColor = {'#D95319';'#0072BD';'#7E2F8E'}; %velocity colors

% motion pulse
settings.mopColor = {"#0072BD";"#7E2F8E"};

% background
settings.bckColor = {[0.6 0.6 0.6]; "#A2142F"};

% behavioral genetics
settings.geneColor = {'#77AC30',"#D95319","#EDB120"};

% menotaxis
settings.HDColor = "#77AC30";

end