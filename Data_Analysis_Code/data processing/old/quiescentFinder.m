% quiescentFinder.m
% Function used to pull only rest epochs by tresholding against forward
% movement
%
% INPUTS:
% allForward - forward velocity array containing all data
% allAngular - angular velocity array containing all data
% allSideway - sidways velocity array containing all data
% allSpikeRt - spike rate array containing all data
% expttime - time array for 1 trial
% moveThresh - minimum forward velocity for "moving" (e.g., 0.1)
% 
% INPUTS:
% restForward - forward velocity array containing only rest, others nan
% restAngular - angular velocity array containing only rest, others nan
% restSideway - sidways velocity array containing only rest, others nan
% restSpikeRt - spike rate array containing only rest, others nan
%
% CREATED: 06/16/2023 MC created from pursuitFinder
%

function [restForward,restAngular,restSideway,restSpikeRt] = quiescentFinder(allForward,allAngular,allSideway,allSpikeRt,expttime,moveThresh)
%% initialize

% copy data
restForward = allForward;
restAngular = allAngular;
restSideway = allSideway;
restSpikeRt = allSpikeRt;

% find all indices where fly was moving above threshold
movingIdx = find(allForward<-moveThresh | allForward>moveThresh);

%% pull quiescent data

restForward(movingIdx) = nan;
restAngular(movingIdx) = nan;
restSideway(movingIdx) = nan;
restSpikeRt(movingIdx) = nan;


end