function [x, y] = flattenJump(posCell, timeCell)
% posCell/timeCell: cell arrays (animals × conditions) with vectors per cell
    x = []; y = [];
    for ii = 1:numel(posCell)
        p = posCell{ii};
        t = timeCell{ii};
        if isempty(p) || isempty(t), continue; end
        p = p(:); t = t(:);
        n = min(numel(p), numel(t));  % keep matched pairs only
        p = p(1:n); t = t(1:n);
        m = isfinite(p) & isfinite(t);
        x = [x; p(m)];
        y = [y; t(m)];
    end
end