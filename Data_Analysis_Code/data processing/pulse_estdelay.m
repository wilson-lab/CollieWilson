% pulse_estdelay
% This function estimates the delay of the fly's turn response relative to a motion pulse 
% by analyzing when the fly accelerates in the direction of the pulse. It is best suited 
% for pulses in the front ipsilateral field (e.g., 20-60 degrees).
%
% INPUTS:
%   panelps   - Panel positions representing the motion pulses (degrees)
%   angular   - Angular velocity of the fly (degrees/second)
%   ttime     - Time vector (seconds)
%
% OUTPUTS:
%   delay_out - Delay time (in seconds) for the turn response relative to the pulse onset
%
% CREATED: [Date] MC
%
function delay_out = pulse_estdelay(panelps,angular,ttime)
%% initialize
nPulse = size(panelps,3);

% fetch data
pulsePanelps = reshape(panelps,[],nPulse);
pulseAngular = reshape(angular,[],nPulse);

% fetch time when pulses start
pStart = find(~isnan(pulsePanelps),1);

%% find when fly velocity changes following a motion pulse

% calculate derivative of velocity (acceleration)
diffAngular = diff(pulseAngular);

% find when change exceeds threshold
accThresh = 0.02;
peak_diffT = [];
for p = 1:nPulse
    delay_idx = find(diffAngular(pStart:end,p)>accThresh,1);
    peak_diffT(p) = ttime(delay_idx+pStart);
end

% (optional) plot
% clf(1);tiledlayout(4,1)
% ex = 6:7;
% nexttile; plot(ttime,pulsePanelps(:,ex)); xlim([0 max(ttime)])
% nexttile([2,1]); plot(ttime,pulseAngular(:,ex)); xlim([0 max(ttime)])
% xline(peak_diffT(ex))
% nexttile; plot(ttime(1:end-1),diffAngular(:,ex)); xlim([0 max(ttime)])

% calculate delay time relative to pulse onset time
delay_out = (peak_diffT - ttime(pStart))';

end