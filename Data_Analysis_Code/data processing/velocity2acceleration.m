% velocity2acceleration
% This function converts directional velocities (forward, angular, and sideways) 
% into their respective accelerations by calculating the derivative of each velocity
% with respect to time. It takes in arrays of forward velocity, angular velocity, 
% sideway velocity, and corresponding time points (ttime), then returns a structure 
% 'accel' containing the calculated forward, angular, and sideway accelerations.
%
% Inputs:
% - forward: a matrix where each column represents forward velocity over time for a trial.
% - angular: a matrix where each column represents angular velocity over time for a trial.
% - sideway: a matrix where each column represents sideways velocity over time for a trial.
% - ttime: a vector of time points corresponding to the velocities.
%
% Output:
% - accel: a structure containing the forward, angular, and sideway accelerations.
%
% Optional: If 'optPlot' is set to 1, the function will plot the forward and angular 
% velocities alongside their respective accelerations for a selected trial.
%
function [accel] = velocity2acceleration(forward,angular,sideway,ttime)
%% calculate acceleration as derivative of velocity

% determine change in time
dt = ttime(2) - ttime(1);

% initialize accel structure
accel = struct();

% calculate acceleration (dv/dt) only if input velocities are provided
if ~isempty(forward)
    accel.forward = diff(forward) ./ dt;
else
    accel.forward = []; % set as empty if no forward velocity provided
end

if ~isempty(angular)
    accel.angular = diff(angular) ./ dt;
else
    accel.angular = []; % set as empty if no angular velocity provided
end

if ~isempty(sideway)
    accel.sideway = diff(sideway) ./ dt;
else
    accel.sideway = []; % set as empty if no sideway velocity provided
end

% check (optional)
optPlot = 0;
if optPlot
    figure(1); clf(1)
    tSelect = 1; % trial to plot
    if ~isempty(forward)
        subplot(2,1,1)
        plot(forward(:,tSelect),'k'); hold on; plot(accel.forward(:,tSelect),'r')
    end
    if ~isempty(angular)
        subplot(2,1,2)
        plot(angular(:,tSelect),'k'); hold on; plot(accel.angular(:,tSelect),'r')
    end
end

end
