% pipeline_kir_pursuit_pool
%
% Pipeline Function
% Pulls all processed files from ALL flies in a given experiment, performs
% necessary analyses and plots accordingly.
%
% INPUTS
% exptFolder- overarching experiment folder
%
% 06/07/2024 - MC adapted from kir pursuit pipeline
%

function pipeline_noise_pool(exptFolder,trialTypes)
%% initialize
disp('STARTING ANALYSES FOR POOLED KIR PURSUIT...')
close all

% set filename info and create necessary directories
filebase = strrep(exptFolder,' ','_');
mainfolder = ['E:\' exptFolder];
cd(mainfolder)
intFolder = [mainfolder '\interpolated']; %for saving interpolated data
if ~exist(intFolder, 'dir')
    mkdir(intFolder)
end
summaryFolder = [mainfolder '\summary']; %for saving plots
if ~exist(summaryFolder, 'dir')
    mkdir(summaryFolder)
end
vectorFolder = [summaryFolder '\vector']; %for saving plots
if ~exist(vectorFolder, 'dir')
    mkdir(vectorFolder)
end
dropboxFolder = ['C:\Users\wilson\Dropbox (HMS)\Data\' exptFolder '\summary']; %for saving data to dropbox
if ~exist(dropboxFolder, 'dir')
    mkdir(dropboxFolder)
end

cd(intFolder)
% find all files for each experiment set
intFiles = dir('*_int.mat');
nFlies = length(intFiles);

%% set key variables

genotypes = {'KIR';'WT';'NA'}; %genotype conditions
nGene = size(genotypes,1);
nTypes = length(trialTypes);

runThresh = 5; %mm/s, min fwd speed to be considered "running" bout
minRunTime = 6*60; %s, min time spend running to be considered "running well"

% set offset for overlapping data
o = 0.15;
offset = [-o; 0; o]; %used to offset overlapping plots

% set plot colors
sem_alpha = 0.2;
color_genotype = {'#77AC30',"#D95319","#EDB120"};
color_trials = [0.6 0.6 0.6];

% set plot limits
yHist = [0 0.2];
ySP = [0 0.25];

%% load in data
disp('Loading in and analyzing datasets...')
nKIR = 0; nWT = 0; nNA = 0;
for e = 1:nFlies
    disp(['Processing fly ' num2str(e) '/' num2str(nFlies) '...'])
    % load this trial
    cd(intFolder)
    thisTrial = intFiles(e).name;
    load(thisTrial)

    % find timepoints where the fly was running and for how long in total
    [run_forward,run_angular,run_panelps,~] = pursuitFinder(int_forward,int_angular,int_panelps,0,int_time,runThresh);
    int_time(end+1) = 60;
    thisRunTime = sum(int_time(sum(~isnan(run_forward))+1)); %sec
    thisRunSpeed = mean(run_forward,'all','omitnan'); %avg run speed
    thisTurnSpeed = mean(abs(run_angular),'all','omitnan'); %avg turn speed

    % measure distribution around setpoint (0)
    [hBins,hData,setpointProb] = setpoint_distribution(int_panelps);
    [hBins,hDataR,setpointProbR] = setpoint_distribution(run_panelps);

    % pool this trial data according to the assigned trial condition
    if contains(thisTrial,'KIR') %kir perturbation flies
        nKIR = nKIR+1;
        kirRunTime(:,nKIR) = reshape(thisRunTime,[],1);
        kirRunSpeed(:,nKIR) = reshape(thisRunSpeed,[],1);
        kirTurnSpeed(:,nKIR) = reshape(thisTurnSpeed,[],1);

        kirSetpoint(:,nKIR) = setpointProb;
        kirHist(:,:,nKIR) = hData;
        kirSetpointR(:,nKIR) = setpointProbR;
        kirHistR(:,:,nKIR) = hDataR;
    elseif contains(thisTrial,'WT') %wildtype control flies
        nWT = nWT+1;
        wtRunTime(:,nWT) = reshape(thisRunTime,[],1);
        wtRunSpeed(:,nWT) = reshape(thisRunSpeed,[],1);
        wtTurnSpeed(:,nWT) = reshape(thisTurnSpeed,[],1);

        wtSetpoint(:,nWT) = setpointProb;
        wtHist(:,:,nWT) = hData;
        wtSetpointR(:,nWT) = setpointProbR;
        wtHistR(:,:,nWT) = hDataR;
    elseif contains(thisTrial,'NA') %na control flies
        nNA=nNA+1;
        naRunTime(:,nNA) = reshape(thisRunTime,[],1);
        naRunSpeed(:,nNA) = reshape(thisRunSpeed,[],1);
        naTurnSpeed(:,nNA) = reshape(thisTurnSpeed,[],1);

        naSetpoint(:,nNA) = setpointProb;
        naHist(:,:,nNA) = hData;
        naSetpointR(:,nNA) = setpointProbR;
        naHistR(:,:,nNA) = hDataR;
    end
end
disp('Complete.')


%% compare heading distributions

% initialize
figure; set(gcf,'Position',[100 100 1500 600])
p = 1; %counter

% for each genotype
for g = 1%:nGene
    switch g
        case 1 %kir
            thisHist = kirHist;
            thisN = nKIR;
        case 2 %wt
            thisHist = wtHist;
            thisN = nWT;
        case 3 %na
            thisHist = naHist;
            thisN = nNA;
    end
    % calculate hist mean and sem
    meanHist = mean(thisHist,3);
    semHist = std(thisHist,[],3)./sqrt(thisN);
    
    for t = 1:nTypes
        % plot by genotypes
        subplot(nGene,nTypes,p); hold on
        plot(hBins,thisHist(:,t,:),'Color',color_trials)
        plot(hBins,meanHist(:,t),'Color',color_genotype{g},'LineWidth',1.5)
        xline(0); ylim(yHist)
        p = p+1; %update counter

        if t==1
            xlabel('Obj Pos'); ylabel(genotypes{g});
        end
    end

end

% save plot
sgtitle('Target Distributions All')
cd(summaryFolder)
plotname = 'hist_target';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], dropboxFolder,'f');
% save vectorized plot
cd(vectorFolder)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], dropboxFolder,'f');


% repeat with running data only
% initialize
figure; set(gcf,'Position',[100 100 1500 600])
p = 1; %counter

% for each genotype
for g = 1%:nGene
    switch g
        case 1 %kir
            thisHist = kirHistR;
            thisN = nKIR;
        case 2 %wt
            thisHist = wtHistR;
            thisN = nWT;
        case 3 %na
            thisHist = naHistR;
            thisN = nNA;
    end
    % calculate hist mean and sem
    meanHist = mean(thisHist,3);
    semHist = std(thisHist,[],3)./sqrt(thisN);
    
    for t = 1:nTypes
        % plot by genotypes
        subplot(nGene,nTypes,p); hold on
        plot(hBins,thisHist(:,t,:),'Color',color_trials)
        plot(hBins,meanHist(:,t),'Color',color_genotype{g},'LineWidth',1.5)
        xline(0); ylim(yHist)
        p = p+1; %update counter

        if t==1
            xlabel('Obj Pos'); ylabel(genotypes{g});
        end
    end

end

% save plot
sgtitle('Target Distributions Running')
cd(summaryFolder)
plotname = 'hist_target_r';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], dropboxFolder,'f');
% save vectorized plot
cd(vectorFolder)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], dropboxFolder,'f');

%% compare setpoint probabilities

% initialize
figure; set(gcf,'Position',[100 100 1500 400])

% for each genotype
for g = 1%:nGene
    switch g
        case 1 %kir
            thisSP = kirSetpoint;
            thisN = nKIR;
        case 2 %wt
            thisSP = wtSetpoint;
            thisN = nWT;
        case 3 %na
            thisSP = naSetpoint;
            thisN = nNA;
    end
    % calculate hist mean and sem
    meanSP = mean(thisSP,3);
    semSP = std(thisSP,[],3)./sqrt(thisN);
    
    % plot by genotypes
    subplot(1,nGene+1,g); hold on
    plot(trialTypes,thisSP,'Color',color_trials)
    errorbar(trialTypes,meanSP,semSP,'o','Color',color_genotype{g},'LineWidth',1.5)
    ylim(ySP); xlabel('Noise Amp'); xticks(trialTypes)
    if g==1
        ylabel('Setpoint Probability')
    end
end

% save plot
sgtitle('Setpoint Probabilities All')
cd(summaryFolder)
plotname = 'setpoint_target';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], dropboxFolder,'f');
% save vectorized plot
cd(vectorFolder)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], dropboxFolder,'f');


% repeat with running data only
% initialize
figure; set(gcf,'Position',[100 100 1500 400])

% for each genotype
for g = 1%:nGene
    switch g
        case 1 %kir
            thisSP = kirSetpointR;
            thisN = nKIR;
        case 2 %wt
            thisSP = wtSetpointR;
            thisN = nWT;
        case 3 %na
            thisSP = naSetpointR;
            thisN = nNA;
    end
    % calculate hist mean and sem
    meanSP = mean(thisSP,3);
    semSP = std(thisSP,[],3)./sqrt(thisN);
    
    % plot by genotypes
    subplot(1,nGene+1,g); hold on
    plot(trialTypes,thisSP,'Color',color_trials)
    errorbar(trialTypes,meanSP,semSP,'o','Color',color_genotype{g},'LineWidth',1.5)
    ylim(ySP); xlabel('Noise Amp'); xticks(trialTypes)
    if g==1
        ylabel('Setpoint Probability')
    end
end

% save plot
sgtitle('Setpoint Probabilities Running')
cd(summaryFolder)
plotname = 'setpoint_target_r';
saveas(gcf,[plotname '.png']);
copyfile([plotname '.png'], dropboxFolder,'f');
% save vectorized plot
cd(vectorFolder)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])
copyfile([plotname '.svg'], dropboxFolder,'f');

%% end
disp('ALL ANALYSES COMPLETE.')
end

