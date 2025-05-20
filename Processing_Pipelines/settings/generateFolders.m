% generateFolders
%
% Creates a structured folder hierarchy for organizing files related to 
% experimental pipelines. This function ensures that necessary directories 
% for saving interpolated data, cross-correlation data, summary plots, 
% and Dropbox data backups are established. If a directory does not already 
% exist, it is created.
%
% INPUTS:
% - exptFolder : A string specifying the name of the experiment folder 
%                to be created under the main directory.
%
% OUTPUT:
% - folder : A structure containing paths to the created directories, 
%            including:
%   - main     : Path to the main experiment folder.
%   - int      : Path to the folder for saving interpolated data.
%   - xcorr    : Path to the folder for saving cross-correlation data.
%   - summary   : Path to the folder for saving summary plots.
%   - vector    : Path to the folder for saving vector plots.
%   - dropbox   : Path to the Dropbox folder for backing up data.
%   - compare   : Path to a predefined comparison folder for motion pulse data.
%
% The function changes the current directory to the main experiment folder 
% upon execution.
%
function folder = generateFolders(exptFolder)
folder.main = ['E:\' exptFolder];
cd(folder.main)
folder.int = [folder.main '\interpolated']; %for saving interpolated data
if ~exist(folder.int, 'dir')
    mkdir(folder.int)
end
folder.xcorr = [folder.int '\xcorr']; %for saving cross correlation data
if ~exist(folder.xcorr, 'dir')
    mkdir(folder.xcorr)
end
folder.accl = [folder.int '\accl']; %for saving acclimatization data
if ~exist(folder.accl, 'dir')
    mkdir(folder.accl)
end
folder.summary = [folder.main '\summary']; %for saving plots
if ~exist(folder.summary, 'dir')
    mkdir(folder.summary)
end
folder.vector = [folder.summary '\vector']; %for saving plots
if ~exist(folder.vector, 'dir')
    mkdir(folder.vector)
end
folder.dropbox = ['C:\Users\wilson\HMS Dropbox\Matt Collie\Data\' exptFolder]; %for saving data to dropbox
if ~exist(folder.dropbox, 'dir')
    mkdir(folder.dropbox)
end

folder.compare = 'E:\Compare Motion Pulse\data';
end