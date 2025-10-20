% setpoint_jumps
% analysis function for assessing steering performance during pursuit
% bar jumps. the setpoint is assumed to be 0 (fixation) for all analyses.
%
% INPUT
% panelps - panel position, yaw gain modified
% barjump - panel bar jump triggers
% forward - forward velocity of the fly
% ttime - time
% optPlot - 1/0 to plot or not
%
% OUTPUT
% jump_out - performance variables
%
% 08/21/24 - MC created
% 09/02/24 - MC updated inclusion criteria (must fixate prior)
% 09/03/24 - MC improved correction criteria, added post-jump stats
%
function [jump_out] = setpoint_jumps_old(panelps,barjump,forward,angular,ttime,optPlot)
%% initialize
% dataset info
nCond = size(panelps,3);
nTrials = size(panelps,2);

% load processing settings
settings = processSettings();

% set how much time pre/pst jump to analyze
preWin = 1; %s
pstWin = 10; %s
bufWin = 2; %s
preIdx = fetchTimeIdx(ttime,preWin); %idx
pstIdx = fetchTimeIdx(ttime,pstWin); %idx
bufIdx = fetchTimeIdx(ttime,bufWin); %idx
winSize = preIdx + pstIdx;

% set amount of pre jump window fly must have fixated for
preWinFixationMin = 0.75;
pstWinFixationMin = 0.5;
preWinFixationMinIdx = round(preIdx*preWinFixationMin);
pstWinFixationMinIdx = round(pstIdx*pstWinFixationMin);
pstJumpCorrectionMin = 10; %deg

% % adjust for slight offset between trigger and when panels actually jump
jOff = 0.035; %s
joIdx = fetchTimeIdx(ttime,jOff); %idx

%% find points where fly was fixating the target
fixationOutput = fixationFinder(panelps,forward,ttime,0);
fixIdx = fixationOutput.idx_run; %fetch indices when fly was actively fixating

%% fetch and analyze jumps during pursuit
% initialize data storage arrays
timeCorrect = [];
nCorrect = [];
nMissed = [];

trialsPanelps = {};
trialsAngular = {};
meanPanelps = [];
meanAngular = [];

meanIAE = [];
meanISE = [];
meanCSTD = [];
meanCVar = [];

meanFFT = [];
meanPSD = [];

% for each condition
for c = 1:nCond
    % initialize
    x = 1; %counter
    jumpPanelps = [];
    jumpAngular = [];
    jumpCorrect = [];

    iae = [];
    ise = [];
    cSTD = [];
    cVar = [];

    jFFT = [];
    jPSD = [];

    % estimate HD bias
    biasHD = mean(panelps(:,:,c),'all','omitnan');

    % for each trial
    for t = 1:nTrials
        % fetch data for this condition trial
        thisBarJump = barjump(:,t,c); %bar jump triggers
        thisPanelps = abs(panelps(:,t,c)); %panel position
        thisAngular = abs(angular(:,t,c)); %angular velocity
        thisFixIdx = fixIdx(:,t,c); %fixation tracking index

        % fictrac occassionally breaks, ensure this trial was run properly
        if sum(~isnan(thisPanelps))
            % find all jumps in this trial (left and right)
            idxJump = find(diff(thisBarJump)>0)+joIdx;
            
            % for each jump
            for j = 1:length(idxJump)
                % find this rightward jump
                idxPre = (idxJump(j)-preIdx+1):idxJump(j); %set pre jump indices
                idxPst = idxJump(j):(idxJump(j)+pstIdx); %set post jump indices
                idxPstB = (idxJump(j)+bufIdx):(idxJump(j)+bufIdx)+(pstIdx-bufIdx); %set post jump indices w/buffer
                idxFull = (idxJump(j)-preIdx+1):(idxJump(j)+pstIdx); %set full window indices

                % calculate heading frequency post jump
                freq_out = setpoint_freq(thisPanelps(idxPstB),ttime,0);
                fs = length(freq_out.fft);
                ps = length(freq_out.psd);

                % determine if fly was fixating prior to the bar jump
                preJumpFixate = thisFixIdx(idxPre); %store
                if sum(preJumpFixate)>preWinFixationMinIdx
                    % if fly fixated prior, store this jump data
                    jumpPanelps(1:winSize,x) = thisPanelps(idxFull);
                    jumpAngular(1:winSize,x) = thisAngular(idxFull);
                    
                    % determine if fly fixated post jump
                    pstJumpFixate = thisFixIdx(idxPstB);
                    if sum(pstJumpFixate)>pstWinFixationMinIdx
                        jumpCorrect(x) = findFirstCorrIndex(thisPanelps(idxPst),pstJumpCorrectionMin);
                        
                        if ~isnan(jumpCorrect(x))
                            % analyze behavior from correction to end of analysis window
                            pstCorIdx = jumpCorrect(x):idxFull(end);
                            pstJumpPanelps = thisPanelps(pstCorIdx);
                            pstJumpTime = ttime(1:length(pstCorIdx));
                            % calculate IAE (integral of the absolute error)
                            abs_error = pstJumpPanelps-biasHD;
                            iae(x) = trapz(pstJumpTime,abs_error);
                            % calculate ISE (integral of the squared error)
                            sqd_error = (pstJumpPanelps-biasHD).^2;
                            ise(x) = trapz(pstJumpTime,sqd_error);
                            % calculate circular standard deviation
                            cSTD(x) = circ_std(deg2rad(pstJumpPanelps),[],[],1);
                            % calculate circular variance
                            cVar(x) = 1-circ_var(deg2rad(pstJumpPanelps),[],[],1);
                            % store frequency
                            if x==1
                                jFFT(:,x) = freq_out.fft;
                                jPSD(:,x) = freq_out.psd;
                            else
                                jFFT(1:fs,x) = freq_out.fft;
                                jPSD(1:ps,x) = freq_out.psd;
                            end

                        else
                            iae(x) = nan;
                            ise(x) = nan;
                            cSTD(x) = nan;
                            cVar(x) = nan;

                            jFFT(:,x) = nan;
                            jPSD(:,x) = nan;
                        end
                    else
                        jumpCorrect(x) = nan;

                        iae(x) = nan;
                        ise(x) = nan;
                        cSTD(x) = nan;
                        cVar(x) = nan;

                        jFFT(:,x) = nan;
                        jPSD(:,x) = nan;
                    end
                    x = x+1; %update counter
                end
            end
        else %omit
        end
    end
    % fetch corrected trials
    corrIdx = ~isnan(jumpCorrect);

    % analyze trials from this conditionn
    timeCorrect(c) = mean(ttime(jumpCorrect(corrIdx)+1),'omitnan'); %correction time
    nCorrect(c) = sum(~isnan(jumpCorrect)); %number of corrected jumps
    nMissed(c) = sum(isnan(jumpCorrect)); %number of missed jumps

    trialsPanelps{c} = jumpPanelps(:,corrIdx); %trial panel data
    trialsAngular{c} = jumpAngular(:,corrIdx); %trial angular data
    if sum(corrIdx)
        meanPanelps(1:winSize,c) = mean(jumpPanelps(:,corrIdx),2,'omitnan'); %mean panel
        meanAngular(1:winSize,c) = mean(jumpAngular(:,corrIdx),2,'omitnan'); %mean angular

        meanIAE(c) = mean(iae,'omitnan');
        meanISE(c) = mean(ise,'omitnan');
        meanCSTD(c) = mean(cSTD,'omitnan');
        meanCVar(c) = mean(cVar,'omitnan');

        meanFFT(1:fs,c) = mean(jFFT,2,'omitnan');
        meanPSD(1:ps,c) = mean(jPSD,2,'omitnan');
    else
        meanPanelps(1:winSize,c) = nan;
        meanAngular(1:winSize,c) = nan;

        meanIAE(c) = nan;
        meanISE(c) = nan;
        meanCSTD(c) = nan;
        meanCVar(c) = nan;

        meanFFT(1:fs,c) = nan;
        meanPSD(1:ps,c) = nan;
    end
end

% fetch time, w/ 0 at time of jump
tjump = ttime(1:winSize)-preWin;

% store analyses for output
jump_out.correctTime = timeCorrect;
jump_out.nCorrect = nCorrect;
jump_out.nMissed = nMissed;

jump_out.trialPanelps = trialsPanelps;
jump_out.trialAngular = trialsAngular;
jump_out.meanPanelps = meanPanelps;
jump_out.meanAngular = meanAngular;

jump_out.iae = meanIAE;
jump_out.ise = meanISE;
jump_out.cstd = meanCSTD;
jump_out.cvar = meanCVar;

jump_out.fft = meanFFT;
jump_out.f_fft = freq_out.f_fft;
jump_out.psd = meanPSD;
jump_out.f_psd = freq_out.f_psd;

jump_out.tjump = tjump;

%% optional plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1500 600])
    tiledlayout(2,nCond,'TileSpacing','compact')

    for c = 1:nCond
        % fetch trial data
        thisTrialData = trialsPanelps{c};
        % plot correction panelps
        nexttile; hold on
        if ~isempty(thisTrialData)
            plot(tjump,thisTrialData,'Color',settings.trialColor,'LineWidth',settings.lwTri)
            plot(tjump,meanPanelps(:,c),'Color',settings.HDColor,'LineWidth',settings.lwAvg)
        end
        axis tight; ylim([0 200]); xline(0); xlabel('Time (s)'); ylabel('Object Pos (deg)')
        title([num2str(settings.pursuitGain(c)) 'X'])
    end

    for c = 1:nCond
        % fetch trial data
        thisTrialData = trialsAngular{c};
        % plot correction angular
        nexttile; hold on
        if ~isempty(thisTrialData)
            plot(tjump,thisTrialData,'Color',settings.trialColor,'LineWidth',settings.lwTri)
            plot(tjump,meanAngular(:,c),'Color',settings.velColor{2},'LineWidth',settings.lwAvg)
                end
        axis tight; ylim([0 300]); xline(0); xlabel('Time (s)'); ylabel('Angular Vel (deg/s)')
        title([num2str(settings.pursuitGain(c)) 'X'])
    end
end

end