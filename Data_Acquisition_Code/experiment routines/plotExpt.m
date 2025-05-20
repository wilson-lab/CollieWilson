% plotExpt.m
%
% Function that generates a simple summary plot containing all relevant
% data acquired during a given trial. For electrophysiology (ephys) only, 
% plots only ephys data. For fictrac and ephys, plots both behavior and 
% ephys data.
%
% INPUTS:
%   exptData - processed data, can contain ephys, fictrac, and/or output.
%   exptMeta - processed electrophysiology metadata.
%
% OUTPUTS:
%
% Original: 04/05/2021 - MC
%           11/04/2021 - MC removed g3 and added g4.
%           02/06/2023 - MC made small modifications to display and output plots.
%
function [] = plotExpt(exptData,exptMeta)
    
    % reset counters
    n = 0;
    
    % set number of sublots based on number of expt variables
    checkEphys = contains(exptMeta.exptCond,'ephys','IgnoreCase',true);
    checkIInj = contains(exptMeta.exptCond,'inj','IgnoreCase',true);
    checkG4 = contains(exptMeta.exptCond,'g4','IgnoreCase',true);
    checkFicTrac = contains(exptMeta.exptCond,'fictrac','IgnoreCase',true);
    checkOpto = contains(exptMeta.exptCond,'stim','IgnoreCase',true); %not separate
    checkPython = contains(exptMeta.exptCond,'jump','IgnoreCase',true); %not separate

    s = checkEphys + checkIInj + checkG4 + (checkFicTrac*3);
    
    % set plot height based on number of expt variables
    if s>2
        plotHeight = 750;
    else
        plotHeight = 400;
    end
    
    % generate figure
    figure(1); clf;
    set(gcf,'Position',[100 100 1500 plotHeight])
    lw = 1;

    
    %% plot panel display
    if checkG4
        if std(exptData.g4displayXPos) %and only if not stationary
            n = n+1; % update counter
            
            % pull xpos data for pre-plot processing
            g4Pos_mod = exptData.g4displayXPos;
            
            % if function used center and remove noise
            if isfield(exptMeta,'func')
                % set midpoint (in front of fly) based on object size
                midPos = (88 - (exptMeta.objSize/2 - 1))/192 *360; %center position, in degrees
                % zero data across midpoing, with right + and left -
                g4Pos_mod = g4Pos_mod - midPos;

                % some experiments use a "hidden" position for when the target
                % is behind the fly in the empty column, and therefore not visible
                if contains(exptMeta.func,'pulse')
                    hiddenPos = (184 - (exptMeta.objSize/2))/192 * 360; %hidden position, in degrees
                    g4Pos_mod(exptData.g4displayXPos>hiddenPos) = nan;
                end

                % hide noise caused by data acquisition or motion across sides
                g4Pos_mod(abs(diff(g4Pos_mod))>2) = nan;
                g4Pos_mod(isoutlier(g4Pos_mod)) = nan;

                % else no function used, center but do not remove noise
            else
                g4Pos_mod = g4Pos_mod - 180;
            end


            % plot
            posLim = nanmax(abs(g4Pos_mod));
            ex(n) = subplot(s,1,n);
            p = plot(exptData.t, g4Pos_mod, 'Color','#77AC30','LineWidth',lw);
            %p.YData(abs(diff(p.YData))>180) = nan;
            y(n) = ylabel('Object (deg)');
            
            axis tight
            if posLim < 10 || isnan(posLim)
                ylim([-10 10])
            else
                ylim([-posLim posLim])
            end
            xlim([0 max(exptData.t)])
            yline(0,':','Color','k') %line at mideline
            
        end
    end
    
    
    %% plot current pulse
    if checkIInj
        n = n+1; % update counter
        
        ex(n) = subplot(s,1,n);
        plot(exptData.t, exptData.iInj, 'k')
        y(n) = ylabel('I-Inj (pA)');
        axis tight

    end


    %% plot ephys
    if checkEphys
        n = n+1; % advance counter
        
        switch exptMeta.mode
            % voltage clamp, scaled out is current
            case {'Track','V-Clamp'}
                ex(n) = subplot(s,1,n);
                plot(exptData.t, exptData.scaledCurrent, 'k')
                y(n) = ylabel('Current (pA)');
                axis tight

            % current clamp, scaled out is voltage
            case {'I=0','I-Clamp','I-Clamp Fast'}
                ex(n) = subplot(s,1,n);
                plot(exptData.t, exptData.scaledVoltage, 'k')
                y(n) = ylabel('Voltage (mV)');
                axis tight

        end
    end

    
    %% plot directional velocities (3)
    if checkFicTrac
        n = n+1; % update counter
        
        ex(n) = subplot(s,1,n);
        plot(exptData.tDS, exptData.angularVelocity, 'Color','#0072BD','LineWidth',lw)
        y(n) = ylabel('Angular (deg/s)');
        yline(0,':','Color','k')
        axis tight
        % if fly doesnt walk, set ylim above just baseline noise
        angMax = max(abs(exptData.angularVelocity));
        if angMax < 20
            ylim([-20 20])
        else
            ylim([-angMax angMax])
        end
        n = n+1; % update counter
        
        ex(n) = subplot(s,1,n);
        plot(exptData.tDS, exptData.sidewaysVelocity, 'Color','#7E2F8E','LineWidth',lw)
        y(n) = ylabel('Sideways (mm/s)');
        yline(0,':','Color','k')
        axis tight
        % if fly doesnt walk, set ylim above just baseline noise
        sidMax = max(abs(exptData.sidewaysVelocity));
        if sidMax < 1
            ylim([-1 1])
        else
            ylim([-sidMax sidMax])
        end
        n = n+1; % update counter
        
        ex(n) = subplot(s,1,n);
        plot(exptData.tDS, exptData.forwardVelocity, 'Color','#D95319','LineWidth',lw)
        yline(0,':','Color','k')
        y(n) = ylabel('Forward (mm/s)');
        axis tight
        % if fly doesnt walk, set ylim above just baseline noise
        fwdMax = exptData.forwardVelocity;
        if max(fwdMax) < 2
            ylim([min(exptData.forwardVelocity) 2])
        end
        
    end
    

    %% add opto stim line if needed
    if checkOpto
        % find both stim onsets and offsets based on output changes
        optoBuffer = [0; exptData.optoStim; 0]; %add buffer to catch stims at start/stop of trial
        optoOn = find(diff(optoBuffer)>0);
        optoOff = find(diff(optoBuffer)<0)-1;
        
        if ~isempty(optoOn)
            for ol=1:s
                subplot(s,1,ol)
                hold on
                xline(exptData.t(optoOn),'Color',"#77AC30",'LineWidth',4)
                xline(exptData.t(optoOff),'Color',"#A2142F",'LineWidth',4)
            end
        end
    end


    %% add python trigger if needed
    if checkPython
        try
            % find when panel xposition was jumped via python socket client
            jump_idx = find(diff(exptData.pythonJumpTrig)>0)+1;

            for ol=1
                subplot(s,1,ol)
                hold on
                xline(exptData.t(jump_idx),'Color',"#A2142F",'LineWidth',1)
            end
        end
    end


    %% add label
    % add at label at the end, regardless of the number of subplots
    xlabel('Time (s)')
    linkaxes(ex,'x');

    
end

