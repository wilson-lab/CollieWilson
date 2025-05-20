% setpoint_polarhist
% This function generates a polar histogram of fly pursuit behavior based on 
% panel position (panelps). It outputs the vector strength and direction 
% for each condition, which represent the fly's directional preference during pursuit.
%
% INPUTS:
% panelps   - panel positions (degrees) representing fly's pursuit target positions
% forward   - forward velocities of the fly
% ttime     - time vector (in seconds) for the trial
% settings  - structure containing analysis settings (e.g., bin sizes, colors, gain values)
% optPlot   - 0 to omit plotting, 1 to display polar histogram plots
%
% OUTPUT:
% vector    - structure containing the following:
%             - hd: vector direction (in radians) for each condition
%             - strength: vector strength (magnitude) for each condition
%
% CREATED: [Date] MC
%
function vector = setpoint_polarhist(panelps,forward,ttime,settings,optPlot)
%% initialize
% fetch number of conditions
nCond = size(panelps,3);
% convert panelps to radians
panelps_rad = deg2rad(panelps);

%% restrict behavior space

% restrict to timepoints when the fly was running
[~,~,~,panelps_rad_r] = pursuitFinder(forward,0,0,panelps_rad,ttime,settings.runThreshB);

%% for each condition
% initialize
vector = [];
h_max = nan(1,nCond);
h_idx = nan(1,nCond);
figure; set(gcf,'Position',[100 100 1800 300])
tiledlayout(1,nCond,'TileSpacing','compact')

for c = 1:nCond
    thisPanelps = reshape(panelps_rad_r(:,:,c),[],1);
    thisGain = settings.pursuitGain(c);
    nexttile

    if any(~isnan(thisPanelps))
        h = polarhistogram(thisPanelps,settings.HDBins,'FaceColor',settings.HDColor,'FaceAlpha',0.5,'Normalization','probability');
        [h_max(c),h2_imax] = max(h.Values); %find most prominent heading
        hold on
        h_idx(c) = (h.BinEdges(h2_imax)+h.BinEdges(h2_imax+1))/2; %find prominent heading bin center
        polarplot([h_idx(c);h_idx(c)],[0;h_max(c)],'Color',"#7E2F8E",'LineWidth',2)
        title([num2str(thisGain) 'X'])
    else
        h_idx(c) = nan;
        h_max(c) = nan;
    end
end
% store HD vector direction and strength
vector.hd = h_idx; %rad, vector direction
vector.strength = h_max; %vector strength

%% optional plot
if ~optPlot
    close all
end
end