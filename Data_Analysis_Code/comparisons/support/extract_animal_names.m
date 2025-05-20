function names = extract_animal_names(files, pattern)
    % Extract animal names based on a specified pattern
    names = cellfun(@(x) regexp(x, pattern, 'tokens', 'once'), {files.name}, 'UniformOutput', false);
    names = cellfun(@(x) x{1}, names, 'UniformOutput', false);
end