% plotPatchCheck
%
% Plots the input resistance, spike rate, and baseline membrane potential 
% for a series of trials. This function visualizes the patch-clamp data 
% to assess the quality and stability of recordings across trials.
%
% INPUTS:
% - checkInputR : Array of input resistance values (in MOhm) for each trial.
% - checkVm     : Array of baseline membrane potentials (in mV) for each trial.
% - checkSR     : Array of spike rates (in spikes/sec) for each trial.
%
% OUTPUT:
% The function does not return any values, but it generates a figure with 
% three subplots displaying the respective data and saves the figure and 
% the data to files.
%
% The figure includes:
% - Subplot 1: Input resistance (IR) across trials.
% - Subplot 2: Spike rate (sp/s) across trials.
% - Subplot 3: Baseline membrane potential (Vm) across trials.
%
% The data is saved to 'patchcheck.mat' and the figure is saved as 
% 'patchcheck_scatterplot.png'.
%
% CREATED: [Date of creation] by [Your Name]
%
function  plotPatchCheck(checkInputR,checkVm,checkSR)
% pull trial N
trials = 1:size(checkInputR,2);

% plot input resistance and baseline membrane potential
figure(2); clf;
set(gcf,'Position',[100 100 700 400]);

s(1) = subplot(3,1,1);
scatter(trials,abs(checkInputR),'filled','MarkerFaceColor','#0072BD')
ylabel('IR (MOhm)')

s(2) = subplot(3,1,2);
scatter(trials,checkSR,'filled','MarkerFaceColor','#D95319')
ylabel('spikes/sec')

s(3) = subplot(3,1,3);
scatter(trials,checkVm,'filled','MarkerFaceColor','#A2142F')
ylabel('Vm (mV)')

sgtitle('patch status')
linkaxes(s,'x')


% save data
save('patchcheck.mat', 'checkInputR','checkVm','checkSR','-v7.3');
% save plot
saveas(gcf,'patchcheck_scatterplot.png');
end