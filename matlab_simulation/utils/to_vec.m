function v = to_vec(x)
% TO_VEC  Coerce yamlmatlab output (scalar, numeric, or cell of scalars) into
%         a numeric row vector.
    if isempty(x)
        v = [];
    elseif isnumeric(x)
        v = double(x(:)).';
    elseif iscell(x)
        v = cellfun(@(e) double(e), x(:)).';
    else
        error('to_vec:unsupported', 'Cannot coerce class %s to vector.', class(x));
    end
end
