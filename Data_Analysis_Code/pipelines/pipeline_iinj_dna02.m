% pipeline_iinj_dna02
%
% Pipeline Function
% Pulls all processed files from ALL flies in a given experiment, performs
% necessary analyses and plots accordingly. Can be used for both behavior
% only experiments and ephys experiments.
%
% INPUTS
% exptFolder - Experiment folder
% pulseSpeeds - Cell array of pulses presented (e.g., 0 25 75 dps)
%
% Created: 06/14/2023 by MC
% Updated: 07/05/2024 by MC (streamlined)
%
function pipeline_iinj_dna02(exptFolder)
disp('ANALYZING POOLED TARGET PULSE EXPT...')
%% initialize folders
% set filename info and create necessary directories
filebase = strrep(exptFolder,' ','_');
% generate file directories as needed
folder = generateFolders(exptFolder);
% load processing settings
settings = processSettings();

%% initialize experiment variables
% fetch number of sweeps and plot dimensions depending on experiment paradigm
nSweep = 9;
tileX = 3;
plotWidth = 800;

% set conditions
currCond = {'Baseline';'Depolarized'};
% initialize plotting variables
% spike rate
sr_limitA = [0 100]; %averages
sr_limitD = [-35 60]; %difference

% voltage
vm_diff = [-3 3]; %difference

% directional velocity
ang_limit = [-100 100]; %angular

% general
runSelect = [-1,0,settings.runThreshE];
n_limit = [-1.5 1.5]; %normalized average
n_limitD = [0 1.5]; %normalized difference
ds_limit = [-1 1]; %direction selectivity

% sweep colors
nSweep_ipsi = ceil(nSweep/2);
nSweep_cntr = floor(nSweep/2);
cm1 = flip(colormap(bone(nSweep_cntr+3)));
cm1(end,:) = []; cm1(1,:) = [];
cm2 = flip(colormap(hsv(nSweep_ipsi)));
color_rightward = [cm1(2:end,:) ; cm2];
color_leftward = [cm1 ; cm2(1:end-1,:)];
close all

%% load in and pool all trials from each experiment folder

disp('Loading in datasets...')
% find all files in this directory
% pull all file info
cd(folder.int)
allFiles = dir('*int.mat');
nFlies = length(allFiles);

% initialize data storage arrays
restPanelPs = cell(1,nFlies);
restForward = cell(1,nFlies);
restSpikeRt = cell(1,nFlies);
restVoltage = cell(1,nFlies);
depolPanelPs = cell(1,nFlies);
depolForward = cell(1,nFlies);
depolSpikeRt = cell(1,nFlies);
depolVoltage = cell(1,nFlies);

for nt = 1:nFlies
    % load in this fly
    cd(folder.int)
    thisTrial = allFiles(nt).name;
    flyShortNames{nt} = strrep(thisTrial(6:16),'_',' ');
    disp(['Loading: ' flyShortNames{nt}])
    load(thisTrial)

    % store trials at low voltage
    restPanelPs{nt} = int_panelps(:,:,1);
    restForward{nt} = int_forward(:,:,1);
    restSpikeRt{nt} = int_spikert(:,:,1);
    [mf_voltage] = spikeFilter(int_voltage(:,:,1),int_time); % median filter voltage signal to remove spikes
    restVoltage{nt} = mf_voltage;
    % store trials at high voltage
    depolPanelPs{nt} = int_panelps(:,:,2);
    depolForward{nt} = int_forward(:,:,2);
    depolSpikeRt{nt} = int_spikert(:,:,2);
    [mf_voltage] = spikeFilter(int_voltage(:,:,2),int_time); % median filter voltage signal to remove spikes
    depolVoltage{nt} = mf_voltage;

end

%% analyze each motion pulse vs spikerate response
disp('Analying motion pulse vs spikerate response...')
cd(folder.summary)

% for each current condition
for c = 1:2
    % select condition
    thisCond = currCond{c};

    % for each behavioral classification
    % run (1) all, (2) quiescent only, (3) running only
    for r = 1:3
        % pull run info
        thisName = settings.behaviorGroup{r};
        thisRun = runSelect(r);
        % initialize
        pulse_srR = [];
        pulse_srL = [];
        mean_pulse_srR = [];
        mean_pulse_srL = [];
        sem_pulse_srR = [];
        sem_pulse_srL = [];

        % for each fly
        for nt = 1:nFlies
            % select data
            switch c
                case 1
                    thisPanelPs = restPanelPs{nt};
                    thisForward = restForward{nt};
                    thisSpikert = restSpikeRt{nt};
                case 2
                    thisPanelPs = depolPanelPs{nt};
                    thisForward = depolForward{nt};
                    thisSpikert = depolSpikeRt{nt};
            end
            % determine relationship between pulse and spikerate
            [~, thisMean] = pulse_v_output(thisPanelPs,thisForward,thisSpikert,int_time,1,1,nSweep,thisRun);
            % if first run, store motion pulse positions
            if nt == 1
                % sweep positions
                pulse_posR(:,:) = thisMean.panelpsR;
                pulse_posL(:,:) = thisMean.panelpsL;
                pulseDur = size(pulse_posR,1);
                % sweep centers for each position
                sweepPosR = round(pulse_posR(round(pulseDur/2),1,:));
                sweepPosL = round(pulse_posL(round(pulseDur/2),1,:));
                % sweep indices
                sweepIdx = find(~isnan(pulse_posR(:,1)));
                sweepIdx2End = sweepIdx(1):length(pulse_posR(:,1));
            end
            % store spikerate averages for each fly
            pulse_srR(:,nt,:) = thisMean.varOutR;
            pulse_srL(:,nt,:) = thisMean.varOutL;
            pulse_srRL(:,nt,:) = thisMean.varOutR - flip(thisMean.varOutL,3);
            % for each sweep, store peak spikerate for each fly
            for s = 1:nSweep
                peak_srR(nt,s) = max(pulse_srR(sweepIdx,nt,s));
                peak_srL(nt,s) = max(pulse_srL(sweepIdx,nt,s));
                peak_srRL(nt,s) = max(pulse_srRL(sweepIdx,nt,s));
            end
        end
        % fetch adjusted n
        nFlies_adj = sum(~isnan(pulse_srRL(1,:,1)));

        % calculate spikerate means and sem
        mean_pulse_srR(:,:) = mean(pulse_srR,2,'omitnan');
        mean_pulse_srL(:,:) = mean(pulse_srL,2,'omitnan');
        mean_pulse_srRL(:,:) = mean(pulse_srRL,2,'omitnan');
        sem_pulse_srR(:,:) = std(pulse_srR,0,2,'omitnan')./sqrt(nFlies_adj);
        sem_pulse_srL(:,:) = std(pulse_srL,0,2,'omitnan')./sqrt(nFlies_adj);
        sem_pulse_srRL(:,:) = std(pulse_srRL,0,2,'omitnan')./sqrt(nFlies_adj);

        % plot averages for each speed
        % initialize
        figure; set(gcf,'Position',[100 100 plotWidth 800])
        tiledlayout(4,tileX,"TileSpacing","compact")
        % initialize remaining variables
        time_pulse = int_time(1:size(pulse_posR,1))*1000;
        time_limit = [0 max(time_pulse)];
        ps_max = ceil(max(pulse_posR,[],'all'));
        ps_limit = [-ps_max ps_max];

        % plot target position
        % plot rightward target
        nexttile; hold on
        for s = 1:nSweep
            plot(time_pulse,pulse_posR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        title('Rightward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(ps_limit); yline(0)
        % plot leftward target
        nexttile; hold on
        for s = nSweep:-1:1
            plot(time_pulse,pulse_posL(:,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
        end
        title('Leftward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(ps_limit); yline(0)

        % plot R-L
        nexttile; hold on
        for s = 1:nSweep
            plot(time_pulse,pulse_posR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        title('R-L'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(ps_limit); yline(0)

        % plot spikerate
        % plot rightward turning
        nexttile([3 1]); hold on
        for s = 1:nSweep
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srR(:,s)-sem_pulse_srR(:,s); flipud(mean_pulse_srR(:,s)+sem_pulse_srR(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = color_rightward(s,:);
            plot(time_pulse,mean_pulse_srR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        ylabel('Firing Rate (spikes/s)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(sr_limitA); yline(0)
        % plot leftward spikerate
        nexttile([3 1]); hold on
        for s = 1:nSweep
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srL(:,s)-sem_pulse_srL(:,s); flipud(mean_pulse_srL(:,s)+sem_pulse_srL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = color_leftward(s,:);
            plot(time_pulse,mean_pulse_srL(:,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
        end
        ylabel('Firing Rate (spikes/s)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(sr_limitA); yline(0)

        % plot R-L turning
        nexttile([3 1]); hold on
        for s = 1:nSweep
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srRL(:,s)-sem_pulse_srRL(:,s); flipud(mean_pulse_srRL(:,s)+sem_pulse_srRL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = color_rightward(s,:);
            plot(time_pulse,mean_pulse_srRL(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        ylabel('Expected Firing Rate Difference (spikes/s)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(sr_limitD); yline(0)

        sgtitle([strrep([filebase ' ' thisCond ' ' thisName],'_','/') ' (n = ' num2str(nFlies_adj) ')'])
        % save plot
        cd(folder.summary)
        plotname = strjoin({'pulse_v_fr', thisCond, 'dps', thisName},'_');
        saveas(gcf,[plotname '.png']);
        copyfile([plotname '.png'], folder.dropbox,'f');
        % save vectorized plot
        cd(folder.vector)
        set(gcf,'renderer','Painters')
        saveas(gcf, [plotname '.svg'])
        copyfile([plotname '.svg'], folder.dropbox,'f');

        % plot pulse speeds together
        % initialize
        figure; set(gcf,'Position',[100 100 1500 900])
        tiledlayout(6,nSweep,"TileSpacing","compact")
        % plot rightward target positions overlayed
        for s = 1:nSweep
            nexttile; hold on
            plot(time_pulse,pulse_posR(:,s),'Color',settings.mopColor{1},'LineWidth',settings.lwAvg)
            if s==1
                ylabel('Target Pos (deg)');
            end
            xlim(time_limit); ylim(ps_limit); yline(0)
        end
        % plot rightward firing rate averages overlayed
        for s = 1:nSweep
            nexttile([2 1]); hold on
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srR(:,s)-sem_pulse_srR(:,s); flipud(mean_pulse_srR(:,s)+sem_pulse_srR(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = settings.mopColor{1};
            plot(time_pulse,mean_pulse_srR(:,s),'Color',settings.mopColor{1},'LineWidth',settings.lwAvg)
            if s==1
                ylabel('Firing Rate (spikes/s)');
            end
            xlim(time_limit); ylim(sr_limitA); yline(0)
        end
        % plot leftward target positions overlayed
        for s = 1:nSweep
            nexttile; hold on
            plot(time_pulse,pulse_posL(:,s),'Color',settings.mopColor{1},'LineWidth',settings.lwAvg)

            if s==1
                ylabel('Target Pos (deg)');
            end
            xlim(time_limit); ylim(ps_limit); yline(0)
        end
        % plot leftwrad firing rate averages overlayed
        for s = 1:nSweep
            nexttile([2 1]); hold on
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srL(:,s)-sem_pulse_srL(:,s); flipud(mean_pulse_srL(:,s)+sem_pulse_srL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = settings.mopColor{1};
            plot(time_pulse,mean_pulse_srL(:,s),'Color',settings.mopColor{1},'LineWidth',settings.lwAvg)

            if s==1
                ylabel('Firing Rate (spikes/s)');
            end
            xlim(time_limit); ylim(sr_limitA); yline(0)
        end

        sgtitle(strrep([filebase ' ' thisCond ' ' thisName],'_','/'))
        % save plot
        cd(folder.summary)
        plotname = strjoin({'pulse_v_fr_overlay' ,thisCond,thisName},'_');
        saveas(gcf,[plotname '.png']);
        copyfile([plotname '.png'], folder.dropbox,'f');
        % save vectorized plot
        cd(folder.vector)
        set(gcf,'renderer','Painters')
        saveas(gcf, [plotname '.svg'])
        copyfile([plotname '.svg'], folder.dropbox,'f');
    end
end

%% analyze each motion pulse vs voltage response
disp('Analying motion pulse vs voltage response...')
cd(folder.summary)

% for each current condition
for c = 1:2
    % select condition
    thisCond = currCond{c};

    % for each behavioral classification
    % run (1) all, (2) quiescent only, (3) running only
    for r = 1:3
        % pull run info
        thisName = settings.behaviorGroup{r};
        thisRun = runSelect(r);
        % initialize
        pulse_srR = [];
        pulse_srL = [];
        mean_pulse_srR = [];
        mean_pulse_srL = [];
        sem_pulse_srR = [];
        sem_pulse_srL = [];

        % for each fly
        for nt = 1:nFlies
            % select data
            switch c
                case 1
                    thisPanelPs = restPanelPs{nt};
                    thisForward = restForward{nt};
                    thisSpikert = restVoltage{nt};
                case 2
                    thisPanelPs = depolPanelPs{nt};
                    thisForward = depolForward{nt};
                    thisSpikert = depolVoltage{nt};
            end
            % determine relationship between pulse and spikerate
            [~, thisMean] = pulse_v_output(thisPanelPs,thisForward,thisSpikert,int_time,1,1,nSweep,thisRun);
            % if first run, store motion pulse positions
            if nt == 1
                % sweep positions
                pulse_posR(:,:) = thisMean.panelpsR;
                pulse_posL(:,:) = thisMean.panelpsL;
                pulseDur = size(pulse_posR,1);
                % sweep centers for each position
                sweepPosR = round(pulse_posR(round(pulseDur/2),1,:));
                sweepPosL = round(pulse_posL(round(pulseDur/2),1,:));
                % sweep indices
                sweepIdx = find(~isnan(pulse_posR(:,1)));
                sweepIdx2End = sweepIdx(1):length(pulse_posR(:,1));
            end
            % store spikerate averages for each fly
            pulse_srR(:,nt,:) = thisMean.varOutR;
            pulse_srL(:,nt,:) = thisMean.varOutL;
            pulse_srRL(:,nt,:) = thisMean.varOutR - flip(thisMean.varOutL,3);
            % for each sweep, store peak spikerate for each fly
            for s = 1:nSweep
                peak_srR(nt,s) = max(pulse_srR(sweepIdx,nt,s));
                peak_srL(nt,s) = max(pulse_srL(sweepIdx,nt,s));
                peak_srRL(nt,s) = max(pulse_srRL(sweepIdx,nt,s));
            end
        end
        % fetch adjusted n
        nFlies_adj = sum(~isnan(pulse_srRL(1,:,1)));

        % calculate spikerate means and sem
        mean_pulse_srR(:,:) = mean(pulse_srR,2,'omitnan');
        mean_pulse_srL(:,:) = mean(pulse_srL,2,'omitnan');
        mean_pulse_srRL(:,:) = mean(pulse_srRL,2,'omitnan');
        sem_pulse_srR(:,:) = std(pulse_srR,0,2,'omitnan')./sqrt(nFlies_adj);
        sem_pulse_srL(:,:) = std(pulse_srL,0,2,'omitnan')./sqrt(nFlies_adj);
        sem_pulse_srRL(:,:) = std(pulse_srRL,0,2,'omitnan')./sqrt(nFlies_adj);

        % plot averages for each speed
        % set average baseline
        vm_baseline = mean(mean_pulse_srR(:,:),'all','omitnan');
        vm_limit = vm_baseline+vm_diff;
        % initialize
        figure; set(gcf,'Position',[100 100 plotWidth 800])
        tiledlayout(4,tileX,"TileSpacing","compact")
        % initialize remaining variables
        time_pulse = int_time(1:size(pulse_posR,1))*1000;
        time_limit = [0 max(time_pulse)];
        ps_max = ceil(max(pulse_posR,[],'all'));
        ps_limit = [-ps_max ps_max];

        % plot target position
        % plot rightward target
        nexttile; hold on
        for s = 1:nSweep
            plot(time_pulse,pulse_posR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        title('Rightward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(ps_limit); yline(0)
        % plot leftward target
        nexttile; hold on
        for s = nSweep:-1:1
            plot(time_pulse,pulse_posL(:,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
        end
        title('Leftward Sweeps'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(ps_limit); yline(0)

        % plot R-L
        nexttile; hold on
        for s = 1:nSweep
            plot(time_pulse,pulse_posR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        title('R-L'); ylabel('Target Pos (deg)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(ps_limit); yline(0)

        % plot spikerate
        % plot rightward turning
        nexttile([3 1]); hold on
        for s = 1:nSweep
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srR(:,s)-sem_pulse_srR(:,s); flipud(mean_pulse_srR(:,s)+sem_pulse_srR(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = color_rightward(s,:);
            plot(time_pulse,mean_pulse_srR(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        ylabel('Voltage (mV)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(vm_limit); yline(0)
        % plot leftward spikerate
        nexttile([3 1]); hold on
        for s = 1:nSweep
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srL(:,s)-sem_pulse_srL(:,s); flipud(mean_pulse_srL(:,s)+sem_pulse_srL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = color_leftward(s,:);
            plot(time_pulse,mean_pulse_srL(:,s),'Color',color_leftward(s,:),'LineWidth',settings.lwAvg)
        end
        ylabel('Voltage (mV)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(vm_limit); yline(0)

        % plot R-L turning
        nexttile([3 1]); hold on
        for s = 1:nSweep
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srRL(:,s)-sem_pulse_srRL(:,s); flipud(mean_pulse_srRL(:,s)+sem_pulse_srRL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = color_rightward(s,:);
            plot(time_pulse,mean_pulse_srRL(:,s),'Color',color_rightward(s,:),'LineWidth',settings.lwAvg)
        end
        ylabel('Expected Voltage Difference (mV)'); xlabel('Time (msec)')
        xlim(time_limit); ylim(vm_diff); yline(0)

        sgtitle([strrep([filebase ' ' thisCond ' ' thisName],'_','/') ' (n = ' num2str(nFlies_adj) ')'])
        % save plot
        cd(folder.summary)
        plotname = strjoin({'pulse_v_vm', thisCond, 'dps', thisName},'_');
        saveas(gcf,[plotname '.png']);
        copyfile([plotname '.png'], folder.dropbox,'f');
        % save vectorized plot
        cd(folder.vector)
        set(gcf,'renderer','Painters')
        saveas(gcf, [plotname '.svg'])
        copyfile([plotname '.svg'], folder.dropbox,'f');

        % plot pulse speeds together
        % initialize
        figure; set(gcf,'Position',[100 100 1500 600])
        tiledlayout(6,nSweep,"TileSpacing","compact")
        % plot rightward target positions overlayed
        for s = 1:nSweep
            nexttile; hold on
            plot(time_pulse,pulse_posR(:,s),'Color',settings.mopColor{1},'LineWidth',settings.lwAvg)
            if s==1
                ylabel('Target Pos (deg)');
            end
            xlim(time_limit); ylim(ps_limit); yline(0)
        end
        % plot rightward firing rate averages overlayed
        for s = 1:nSweep
            nexttile([2 1]); hold on
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srR(:,s)-sem_pulse_srR(:,s); flipud(mean_pulse_srR(:,s)+sem_pulse_srR(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = settings.mopColor{1};
            plot(time_pulse,mean_pulse_srR(:,s),'Color',settings.mopColor{1},'LineWidth',settings.lwAvg)
            if s==1
                ylabel('Firing Rate (spikes/s)');
            end
            xlim(time_limit); ylim(vm_limit); yline(0)
        end
        % plot leftward target positions overlayed
        for s = 1:nSweep
            nexttile; hold on
            plot(time_pulse,pulse_posL(:,s),'Color',settings.mopColor{1},'LineWidth',settings.lwAvg)

            if s==1
                ylabel('Target Pos (deg)');
            end
            xlim(time_limit); ylim(ps_limit); yline(0)
        end
        % plot leftwrad firing rate averages overlayed
        for s = 1:nSweep
            nexttile([2 1]); hold on
            sp = patch([time_pulse'; flipud(time_pulse')],[mean_pulse_srL(:,s)-sem_pulse_srL(:,s); flipud(mean_pulse_srL(:,s)+sem_pulse_srL(:,s))], 'r', 'FaceAlpha',settings.semAlpha, 'EdgeColor','none');
            sp.FaceColor = settings.mopColor{1};
            plot(time_pulse,mean_pulse_srL(:,s),'Color',settings.mopColor{1},'LineWidth',settings.lwAvg)

            if s==1
                ylabel('Firing Rate (spikes/s)');
            end
            xlim(time_limit); ylim(vm_limit); yline(0)
        end

        sgtitle(strrep([filebase ' ' thisCond ' ' thisName],'_','/'))
        % save plot
        cd(folder.summary)
        plotname = strjoin({'pulse_v_vm_overlay' ,thisCond,thisName},'_');
        saveas(gcf,[plotname '.png']);
        copyfile([plotname '.png'], folder.dropbox,'f');
        % save vectorized plot
        cd(folder.vector)
        set(gcf,'renderer','Painters')
        saveas(gcf, [plotname '.svg'])
        copyfile([plotname '.svg'], folder.dropbox,'f');
    end
end

%% end
disp('ALL ANALYSES COMPLETE.')
end

