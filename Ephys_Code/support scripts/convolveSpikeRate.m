% convolveSpikeRate.m
%
% Analysis function that detects spikes in voltage or current traces and 
% applies a convolution kernel to calculate the spike rate.
%
% INPUTS:
% - settings  : Struct containing electrophysiological setup settings, 
%               typically obtained from the ephysSettings() function.
% - exptData  : Struct containing processed experiment data, which includes 
%               scaled voltage or current traces and corresponding time data.
% - ktype     : Type of convolution kernel to use, either 'gaussian' or 
%               'causal'.
%
% OUTPUTS:
% - pk_v      : Array of voltage or current values at detected spike peaks.
% - pk_r      : Binary raster array indicating spike occurrences (1 = spike, 
%               0 = no spike).
% - spikeRate : Array of the convolved spike rate, representing the frequency 
%               of spikes over time.
%
% PROCESS:
% The function pulls the time and sample rate from the experimental data, 
% determining whether to process voltage or current traces based on the 
% available fields. It sets parameters for peak detection tailored to the 
% specific recording conditions. After detecting spikes, it creates a raster 
% of spike occurrences and convolves this with the specified kernel type to 
% calculate the spike rate.
%
% The resulting spike rate provides insight into the neuron's firing patterns 
% and how it responds to stimuli.
%
% CREATED:  12/14/2021 by MC
% UPDATED:  05/29/2024 by MC, normalized area to 1
%           06/03/2024 by MC, fixed causal kernel
%
function  [pk_v,pk_r,spikeRate] = convolveSpikeRate(settings,exptData,ktype)

%% pull trace and time

% pull time and sample rate
t = exptData.t;
sr = settings.bob.sampRate;

% pull voltage or current trace
if isfield(exptData,'scaledVoltage') %iclamp
    % set findpeak variables
    minProm = 5;
    maxWidth = 0.5e-1 * sr;
    minDistance = 1.5e-3 * sr;
    minWidth = 0.5e-3 * sr;
    
    ylbl = 'voltage (mV)';
    spikeTrace = exptData.scaledVoltage;
    inv = 1; %do not invert
elseif isfield(exptData,'scaledCurrent') %vclamp
    % set findpeak variables
    minProm = 6;
    maxWidth = .09e-2 * sr;
    minDistance = 1.5e-3 * sr;
    minWidth = 0.3e-3 * sr;
    
    ylbl = 'current (pA)';
    spikeTrace = exptData.scaledCurrent;
    inv = -1; %invert
end


%% detect spikes

% find peaks
[pk_v,pk_i] = findpeaks((spikeTrace .* inv),...
    'MinPeakProminence', minProm,...
    'MinPeakWidth', minWidth,...
    'MaxPeakWidth', maxWidth,...
    'MinPeakDistance', minDistance);
pk_v = pk_v * inv;

% create raster
pk_r = zeros(length(t),1);
pk_r(pk_i) = 1;


%% convolve to find spike rate

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
        sigma = 2000;
        %sigma = 4000; %increase to increase smoothing
        tc = -(sigma*10):(sigma*10);
        kernel = gaussmf(tc,[sigma 0]);
        kernel_norm = kernel./(trapz(tc,kernel)./sr);

        spikeRate = conv(pk_r,kernel_norm,'same');

end

% preview if desired
%plot(tc,kernel);hold on; plot(tc,kernel_norm);xline(0)

%% plot
if 1
    figure(10); clf;
    set(gcf,'Position',[100 100 1200 500])
    sr(1) = subplot(3,1,1:2);
    plot(t,spikeTrace,'k',t(pk_i),pk_v,'.','MarkerEdgeColor','#D95319')
    ylabel(ylbl)

    sr(2) = subplot(3,1,3);
    plot(t,spikeRate,'color','#D95319')
    ylabel('spikes/sec')
    xlabel('time (sec.)')

    linkaxes(sr,'x')
end
end

