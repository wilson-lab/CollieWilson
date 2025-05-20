% initUserDAQ.m
%
% Function to initialize user DAQ (i.e. DAQ that coordinates experiment and
%  collects appropriate voltage signals)
% Sets up channels and DAQ and channel settings. Does not set up any
%  experiment-specific output streams.
%
% INPUT:
%   settings - struct from ephysSettings() containing all DAQ settings
%   aIn - cell array of which analog input channels are actually used
%   aOut - cell arry of which analog output channels are actually used
%   dIn - cell array of which digital input channels are actually used
%   dOut - cell array of which digital output channels are actually used
%
% OUTPUT:
%   userDAQ - handle to user DAQ
%   aiCh - handle to analog input channels
%   aoCh - handle to analog output channels
%   diCh - handle to digital input channels
%   doCh - handle to digital output channels
%
% Created: 7/29/18
% Updated: 
%   8/2/18 - HHY
%   7/24/20 - HHY - updates to handling analog output
%

function [userDAQ, aiCh, aoCh, diCh, doCh] = initUserDAQ(...
    settings, aIn, aOut, dIn, dOut)
    
    disp('Initalizing DAQ');
    
    % start DAQ session
    daqreset;
    userDAQ = daq.createSession(settings.devVendor);
    
    % DAQ settings
    userDAQ.Rate = settings.bob.sampRate;
    
    % ADD CHANNELS
    
    % analog input channels
    if (~isempty(aIn)) % only if there are analog input channels
        aInInd = zeros(size(aIn));
        % map which analog input channels used to actual channel numbers
        for i = 1:length(aIn)
            aInInd(i) = settings.bob.aInChUsed(find(strcmpi(aIn{i},...
                settings.bob.aInChAssign)));
        end
        % actually add analog input channels
        aiCh = userDAQ.addAnalogInputChannel(settings.bob.devID, aInInd, ...
            settings.bob.aiMeasType);
        % set input type and range on channel
        for i = 1:length(aInInd)
            aiCh(i).TerminalConfig = settings.bob.aiInType;
            aiCh(i).Range = settings.bob.aiRange(i,:);
        end
    else
        % need to assign output arguments
        aiCh = [];
    end
    
    % analog output channels
    if (~isempty(aOut)) % only if there are analog output channels
        aOutInd = zeros(size(aOut));
        % map which analog input channels used to actual channel numbers
        for i = 1:length(aOut)
            aOutInd(i) = settings.bob.aOutChUsed(find(strcmpi(aOut{i},...
                settings.bob.aOutChAssign)));
        end
        % actually add analog output channels
        aoCh = userDAQ.addAnalogOutputChannel(settings.bob.devID, ...
            aOutInd, settings.bob.aoMeasType);
        % set input type and range on channels
        for i = 1:length(aOutInd)
            aoCh(i).TerminalConfig = settings.bob.aoOutType;
            aoCh(i).Range = settings.bob.aoRange(i,:);
        end
    else
        aoCh = []; % need to assign output arguments
    end
    
    % digital input channels
    if (~isempty(dIn)) % only if there are digital input channels
        % for all digital input channels
        for i = 1:length(dIn)
            % map channel name to index into dInChUsed
            dInInd = find(strcmpi(dIn{i}, settings.bob.dInChAssign));
            % convert dInChAssign matrix to channel IDs - element 1 is
            %  port, element 2 is line
            channelID = sprintf('Port%d/Line%d', ...
                settings.bob.dInChUsed(dInInd,1), ...
                settings.bob.dInChUsed(dInInd,2));
            % add digital input channel
            diCh(i) = userDAQ.addDigitalChannel(settings.bob.devID,...
                channelID, settings.bob.diType);
        end
    else
        diCh = []; % need to assign output arguments
    end
    
    % digital output channels
    if (~isempty(dOut)) % only if there are digital output channels
        % for all digital output channels
        for i = 1:length(dOut)
            % map channel name to index into dOutChUsed
            dOutInd = find(strcmpi(dOut{i}, settings.bob.dOutChAssign));
            % convert dOutChAssign matrix to channel IDs - element 1 is
            %  port, element 2 is line
            channelID = sprintf('Port%d/Line%d',...
                settings.bob.dOutChUsed(dOutInd,1),...
                settings.bob.dOutChUsed(dOutInd,2));
            % add digital output channel
            doCh(i) = userDAQ.addDigitalChannel(settings.bob.devID,...
                channelID, settings.bob.doType);
        end
    else
        doCh = []; % need to assign output arguments
    end
end