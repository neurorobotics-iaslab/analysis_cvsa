function M = to_mat(x)
% TO_MAT  Coerce yamlmatlab output (numeric / cell-of-rows / cell-of-cells)
%         into a 2-D numeric matrix.
    if isempty(x)
        M = [];
        return;
    end
    if isnumeric(x)
        M = double(x);
        return;
    end
    if iscell(x) && all(reshape(cellfun(@isscalar, x), [], 1))
        M = double(cell2mat(x));
        return;
    end
    if iscell(x)
        rows = cell(numel(x), 1);
        for i = 1:numel(x)
            rows{i} = to_vec(x{i});
        end
        M = vertcat(rows{:});
        return;
    end
    error('to_mat:unsupported', 'Cannot coerce class %s to matrix.', class(x));
end
