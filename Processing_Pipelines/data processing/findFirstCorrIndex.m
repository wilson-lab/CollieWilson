% findFirstCorrIndex
% This function identifies the first index in the panel positions array 
% where the value crosses below a specified minimum position (minPos) after a jump.
% It returns the index of this crossing, or NaN if no crossing occurs.
%
% INPUTS:
%   panelps - Array of panel positions (e.g., representing movement or object position)
%   minPos  - Minimum position value to check against (e.g., threshold value)
%
% OUTPUTS:
%   corrIdx - Index of the first crossing below minPos (NaN if no crossing is found)
%
% CREATED: [Date] MC
%
function corrIdx = findFirstCorrIndex(panelps,minPos)
% find first time panelps array crosses min position (e.g., 10) after a
% jump
corrIdx = find(panelps <= minPos, 1);
% check if no such element is found
if isempty(corrIdx)
    corrIdx = nan; % Return nan if not corrected
end
end