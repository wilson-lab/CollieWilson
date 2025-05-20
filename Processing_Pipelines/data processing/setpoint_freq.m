% SETPOINT_FREQ - This function performs frequency analysis (FFT and PSD) during setpoint fixation.
% It computes the Fast Fourier Transform (FFT) and Power Spectral Density (PSD) for panel position data 
% to evaluate frequency components during pursuit behavior. The option to average over trials and apply
% median filtering is included.

% INPUTS:
%   panelps  - 3D array of panel position data (yaw gain-modified), with each slice representing 
%              a different condition and each column representing a trial.
%   ttime    - Time vector representing the sampling intervals.
%   optPlot  - Optional flag (1/0) for plotting the results (1 for yes, 0 for no).

% OUTPUTS:
%   freq_out - Structure containing frequency analysis results:
%              - freq_out.fft  : Averaged and filtered FFT for each condition.
%              - freq_out.f_fft: Frequency vector for FFT.
%              - freq_out.psd  : Averaged and filtered Power Spectral Density (PSD) for each condition.
%              - freq_out.f_psd: Frequency vector for PSD.

% Created: 09/03/24 by MC
% Updated: 09/04/24 by MC - Added PSD calculation.
% Updated: 09/25/24 by MC - Added trial averaging and median filtering option.
%
function [freq_out] = setpoint_freq(panelps,ttime,optPlot)
%% initialize
% dataset info
nCond = size(panelps,3);
nTrial = size(panelps,2);

% load processing settings
settings = processSettings();
Fs = fetchTimeIdx(ttime,1)-1 ; %sampling frequency

%% analyze setpoint
% initialize
storeFFT = [];
storePSD = [];

% for each condition
for c = 1:nCond
    tempFFT = [];
    tempPSD = [];
    
    % for each trial
    for t = 1:nTrial
        % fetch data
        thisPanelps = panelps(:,t,c);
        % remove nans - can either omit or replace with 0
        % thisPanelps(isnan(thisPanelps)) = [];
        thisPanelps(isnan(thisPanelps)) = 0;
        L = size(thisPanelps,1);

        % define Flat Top window
        window = flattopwin(L); % Create Flat Top window of length L

        % apply window to the signal
        windowedPanelps = thisPanelps .* window;

        % compute Fast Fourier Transform (FFT) estimate
        F = fft(windowedPanelps); % flat top window applied
        TSS = abs(F/L); % two-sided spectrum
        SSS = TSS(1:round(L/2+1)); % single-sided spectrum
        SSS(2:end-1) = 2*SSS(2:end-1); % correct the amplitude
        tempFFT(:,t) = SSS; % store individual trial FFT

        % compute Welch's Power Spectral Density (PSD) estimate
        [PW,f_psd] = pwelch(thisPanelps,window,[],[],Fs); % flat top window applied
        tempPSD(:,t) = PW; % store individual trial PSD
    end
    
    % average over trials
    avgFFT = mean(tempFFT, 2); 
    avgPSD = mean(tempPSD, 2);
    
    % (optional) apply median filtering to the averaged FFT and PSD
    avgFFT = medfilt1(avgFFT, 10);
    avgPSD = medfilt1(avgPSD, 10);

    % store filtered results
    storeFFT(:,c) = avgFFT; 
    storePSD(:,c) = avgPSD;
end

% create a frequency vector
f_fft = Fs*(0:round((L/2)))/L;

% store for output
freq_out.fft = storeFFT;
freq_out.f_fft = f_fft;

freq_out.psd = storePSD;
freq_out.f_psd = f_psd;

%% (optional) plot
if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1200 900])
    tiledlayout(nCond,2,'TileSpacing','compact') % 2 columns: FFT and PSD with filtered results
    
    for c = 1:nCond
        % Plot filtered FFT
        nexttile
        plot(f_fft,storeFFT(:,c))
        xlim([0 5]); xlabel('f (Hz)'); title(['FFT ' num2str(settings.pursuitGain(c)) 'X (Filtered)'])
        
        % Plot filtered PSD
        nexttile
        plot(f_psd,10*log10(storePSD(:,c)))
        xlim([0 5]); xlabel('f (Hz)'); title(['PSD ' num2str(settings.pursuitGain(c)) 'X (Filtered)'])
    end
end
end