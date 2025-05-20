% findFirstCrossing
% This function identifies the first zero crossing in each trace of a given dataset.
% A zero crossing occurs where the sign of the values in the trace changes from positive
% to negative or vice versa. The function returns the index of the first zero crossing
% for each trace, or NaN if no crossings are detected.
%
% INPUTS:
%   traces - Matrix where each row represents a separate trace (e.g., time series data)
%
% OUTPUTS:
%   indices - Vector containing the indices of the first zero crossing for each trace
%             (NaN if no crossing is found)
%
% CREATED: [Date] MC
%
function indices = findFirstCrossing(traces)

% Initialize the output vector
indices = zeros(size(traces, 1), 1);

% Loop through each trace
for i = 1:size(traces, 1)
    trace = traces(i, :);
    % Find indices where the sign changes
    signChanges = find(diff(sign(trace)));
    % Determine the first zero crossing
    if ~isempty(signChanges)
        indices(i) = signChanges(1) + 1;
    else
        indices(i) = nan; % Return nan if there are no zero crossings
    end
end

end