% startFicTrac.m
%
% Function for starting the FicTrac configuration and launching the 
% FicTrac system (and optionally Python scripts) for data acquisition 
% during experiments. This function checks if FicTrac is already 
% configured and, if not, sets it up using the appropriate configuration 
% files.
%
% INPUTS:
% - cellDirPath : Path to the directory where the configuration image 
%                 should be exported.
% - ftMode      : Specifies the mode for FicTrac operation:
%                 0 - Open-loop,
%                 1 - Closed-loop via frame rate,
%                 2 - Closed-loop via frame index,
%                 3 - Closed-loop via frame index with bar jumps,
%                 4 - Closed-loop via frame index with noise,
%                 5 - Leg camera mode.
%
% CREATED: 11/05/2021 by MC
% UPDATED: 04/20/2022 by MC, added open/closed loop alternatives.
% UPDATED: 10/30/2023 by MC, added leg camera option.
%
function startFicTrac(cellDirPath,ftMode)
%% configure fictrac
% set config file to be used
if ftMode==5
    thisConfig = 'config_leg.txt';
else
    system('Taskkill/IM cmd.exe'); %terminate any previous fictrac windows
    thisConfig = 'config_live.txt';
end

% only if it has not yet been configured for this experiment,
% as indicated by a moved configImg in the expt directory
if ~exist([cellDirPath '\fictrac-configImg.png'],'file')
    % start fictrac config
    fictracCfg = ['cd/ & cd /d C:\dev\fictrac\sample & ..\bin\Release\configGui.exe ' thisConfig ' &'];
    
    [~,~] = system(fictracCfg);
    check = input('[INPUT] Configure fictrac...');
    
    % move fictrac config image to data directory
    copyfile('C:\dev\fictrac\sample\fictrac-configImg.png', cellDirPath,'f');
else
    disp('Fictrac previously configured.')
end


%% start fictrac

% launch fictrac gui
system('Taskkill/IM cmd.exe'); %terminate any previous fictrac windows
disp('Launching fictrac...');
fictracStart = ['cd/ & cd /d C:\dev\fictrac\sample & ..\bin\Release\fictrac.exe ' thisConfig ' &'];
[~,~] = system(fictracStart);

%select python script to run
switch ftMode
    case 0
        disp('MODE SELECTION: open-loop');
        thisSocketClient = 'socket_client_360_ol'; %select open-loop
    case 1
        disp('MODE SELECTION: closed-loop - frame rate(4)');
        thisSocketClient = 'socket_client_360_cl_mode4'; %select closed-loop, frame rate
    case 2
        disp('MODE SELECTION: closed-loop - frame index(7)');
        thisSocketClient = 'socket_client_360_cl_mode7'; %select closed-loop, frame index
    case 3
        disp('MODE SELECTION: closed-loop - frame index(7) WITH bar jumps');
        thisSocketClient = 'socket_client_360_cl_mode7_jumps'; %select closed-loop, frame index
    case 4
        disp('MODE SELECTION: closed-loop - frame index(7) WITH noise');
        thisSocketClient = 'socket_client_360_cl_mode7_noise'; %select closed-loop, frame index
    case 5
        thisSocketClient = [];
end

if ~isempty(thisSocketClient)
    % start python script
    disp(['Launching ' thisSocketClient '...']);
    pythonStart = ['cd/ & cd /d C:\dev\fictrac\scripts & python ' thisSocketClient '.py &'];
    [~,~] = system(pythonStart);
else
    disp('No python script launched.')
end

disp('Fictrac ready.');


end

