function B = parse_filters_band(s)
% PARSE_FILTERS_BAND  Parse the launch-file 'filters_band' string
%   ("8.0 10.0; 10.0 12.0; ...") into an [n_bands x 2] matrix.
%   Also accepts already-parsed yamlmatlab cell-of-[lo hi] pairs.
    if isnumeric(s)
        B = double(s);
        return;
    end
    if iscell(s)
        B = to_mat(s);
        return;
    end
    s = strtrim(char(s));
    parts = strsplit(s, ';');
    parts = parts(~cellfun(@isempty, strtrim(parts)));
    B = zeros(numel(parts), 2);
    for k = 1:numel(parts)
        nums = sscanf(parts{k}, '%f');
        if numel(nums) ~= 2
            error('parse_filters_band:bad_row', ...
                  'Band #%d does not have 2 numbers: "%s"', k, parts{k});
        end
        B(k, :) = nums(:).';
    end
end
