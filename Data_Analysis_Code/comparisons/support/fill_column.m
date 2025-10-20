function data_out = fill_column(data_out, all_names, condition_names, condition_vals, col_idx)
    for i = 1:numel(condition_names)
        idx = find(strcmp(all_names, condition_names{i}));
        if ~isempty(idx)
            data_out(idx, col_idx) = condition_vals(i);
        end
    end
end