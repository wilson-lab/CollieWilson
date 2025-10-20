function p_adj = fdr_bh(p)
    [ps, ord] = sort(p(:));
    m = numel(ps);
    adj_sorted = (m ./ (1:m)') .* ps;           % (m/i) * p_(i)
    adj_sorted = flipud(cummin(flipud(adj_sorted))); % enforce monotonicity
    adj_sorted = min(adj_sorted, 1);
    p_adj = zeros(size(ps));
    p_adj(ord) = adj_sorted;
    p_adj = reshape(p_adj, size(p));
end