% convolveSpikeRate_input.m
%
% Analysis function that detects spikes in voltage data and applies a convolution 
% kernel to calculate the spike rate.
%
% INPUTS:
% - settings  : Struct containing electrophysiological setup settings, 
%               typically obtained from the ephysSettings() function.
% - exptData  : Struct containing processed experiment data, including scaled 
%               voltage and time arrays.
% - celltype  : Specifies the cell type, typically 'AOTU019' or similar, 
%               which can affect spike detection parameters.
% - ktype     : Type of convolution kernel to use, either 'gaussian' or 
%               'causal'.
%
% OUTPUTS:
% - pk_v      : Array of voltage values at detected spike peaks.
% - pk_r      : Binary raster array indicating spike occurrences (1 = spike, 
%               0 = no spike).
% - spikeRate : Array of the convolved spike rate, representing the frequency 
%               of spikes over time.
%
% PROCESS:
% The function begins by setting kernel parameters and fetching the relevant 
% data and sampling rate. It then establishes peak detection parameters based 
% on the specified cell type to optimize spike detection. After detecting 
% spikes, it creates a raster of spike occurrences and convolves this with 
% the specified kernel type to calculate the spike rate.
%
% The resulting spike rate is a smoothed representation of spike frequency 
% over time, which is useful for analyzing neuronal firing patterns.
%
% CREATED:  05/17/2023 by MC, created from main convolve function
% UPDATED:  05/29/2024 by MC, normalized area to 1
%           06/03/2024 by MC, fixed causal kernel
%           07/30/2024 by MC, made cell type specific, moved ksize to within function
%
function  [pk_v,pk_r,spikeRate] = convolveSpikeRate_input(settings,exptData,celltype,ktype)
%% set parameters
% set kernel size
ksize= 1000;

% fetch data and samplerate
spikeTrace = exptData.scaledVoltage;
t = exptData.t;
sr = settings.bob.sampRate;

% fetch peak detection parameters based on cell type
% AOTU019 spikes are exceptionally large compared to AOTU025/DNa02,
% therefore it is easier to have different parameters for the two groups
% for better spike detection
if contains(celltype,'AOTU019')
    minProm = 5;
    maxWidth = 0.5e-1 * sr;
    minDistance = 1.5e-3 * sr;
    minWidth = 0.5e-3 * sr;
else
    minProm = 2;
    maxWidth = 0.5e-1 * sr;
    minDistance = 1e-5 * sr;
    minWidth = 0.5e-3 * sr;
end

%% detect spikes
% find peaks
[pk_v,pk_i] = findpeaks(spikeTrace,...
    'MinPeakProminence', minProm,...
    'MinPeakWidth', minWidth,...
    'MaxPeakWidth', maxWidth,...
    'MinPeakDistance', minDistance);

% create spike raster
pk_r = zeros(length(t),1);
pk_r(pk_i) = 1;

% convolve to find spike rate
switch ktype
    case 'causal'
        % create kernel
        tau = ksize*2.5; %increase to increase smoothing
        tc = -(tau*10):(tau*10);
        tc_k = 0:(tau*10); %partial
        kernel = exp(-tc_k./tau);
        kernel = [zeros(1,tau*10) kernel]; %center to full time array
        kernel_norm = kernel./(trapz(tc,kernel)./sr);
        
        % convolve spikerate
        spikeRate = conv(pk_r,kernel_norm,'same');
        
    case 'gaussian'
        % create kernel
        sigma = ksize; %increase to increase smoothing
        tc = -(sigma*10):(sigma*10);
        kernel = gaussmf(tc,[sigma 0]);
        kernel_norm = kernel./(trapz(tc,kernel)./sr);

        % convolve spikerate
        spikeRate = conv(pk_r,kernel_norm,'same');
end

% (optional) plot to verify
if 0
    figure(10); clf;
    set(gcf,'Position',[100 100 1200 500])
    sr(1) = subplot(3,1,1:2);
    plot(t,spikeTrace,'k',t(pk_i),pk_v,'.','MarkerEdgeColor','#D95319')
    sr(2) = subplot(3,1,3);
    plot(t,spikeRate,'color','#D95319')
    ylabel('spikes/sec'); xlabel('time (sec.)')
    linkaxes(sr,'x'); axis tight
end

end

