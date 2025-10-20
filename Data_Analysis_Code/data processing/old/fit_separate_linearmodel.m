% fit_separate_linearmodel
% assess the visual and motor controbituions to changes in firing rate
% during pursuit by separately fitting a linear model to (1) firing vs
% behavior without a visual stimulus (pulse test) and (2) firing vs
% stimulus motion without behavior (quiescent)
%
% INPUT
% exptData - structure containing all necessary data from main experiment
% pulseData - structure containing data from initial opto pulse w/o stim
%
% OUTPUT
%
% CREATED: 06/12/2023 MC
%

function fit_separate_linearmodel(exptData,pulseData,expttime,highThresh)
%% pre-process pursuit data
% threshold for pursuit
[purForward,purAngular,~,purSpikeRt] = pursuitFinder(exptData.forward,exptData.angular,0,exptData.spikert,expttime,highThresh);
purPanelPs = exptData.panelps;
purPanelPs(isnan(purForward)) = nan;
% reshape and remove nans
purForward_r = reshape(purForward(~isnan(purForward)),[],1);
purAngular_r = reshape(purAngular(~isnan(purAngular)),[],1);
purSpikeRt_r = reshape(purSpikeRt(~isnan(purSpikeRt)),[],1);
purPanelPs_r = reshape(purPanelPs(~isnan(purPanelPs)),[],1);

%% pre-process visual data
% reshape visual data
visualSpikeRt_r = reshape(exptData.spikert,[],1);
visualForward_r = reshape(exptData.forward,[],1);
visualPanelPs_r = reshape(exptData.panelps,[],1);

% for visual motion, pull only quiescent behavior
visualSpikeRt_r_thresh = visualSpikeRt_r(visualForward_r<0.1 & visualForward_r>-0.1);
visualPanelPs_r_thresh = visualPanelPs_r(visualForward_r<0.1 & visualForward_r>-0.1);

%% pre-process pulse data
% reshape pulse data
pulseSpikeRt_r = reshape(pulseData.spikert,[],1);
pulseForward_r = reshape(pulseData.forward,[],1);
pulseAngular_r = reshape(pulseData.angular,[],1);

% for forward motion, pull only walking behavior
pulseForward_r_fwdthresh = pulseForward_r(pulseForward_r>0.1);
pulseSpikeRt_r_fwdthresh = pulseSpikeRt_r(pulseForward_r>0.1);
% for angular motion, pull only ipsilateral turning
pulseAngular_r_angthresh = pulseAngular_r(pulseAngular_r>5);
pulseSpikeRt_r_angthresh = pulseSpikeRt_r(pulseAngular_r>5);

%% estimate coefficient for each variable

% for visual
lnfitObj = mvregress(visualPanelPs_r_thresh,visualSpikeRt_r_thresh);
% for motor
lnfitFwd = polyfit(pulseForward_r_fwdthresh,pulseSpikeRt_r_fwdthresh,1);
lnfitAng = mvregress(pulseAngular_r_angthresh,pulseSpikeRt_r_angthresh);


%% generate fit estimates

fitObj = purPanelPs_r.*lnfitObj;
fitFwd = purForward_r.*lnfitFwd;
fitAng = purAngular_r.*lnfitAng;

%% determine correlation

rObj = corrcoef(fitObj,purSpikeRt_r);
rFwd = corrcoef(fitFwd,purSpikeRt_r);
rAng = corrcoef(fitAng,purSpikeRt_r);

end