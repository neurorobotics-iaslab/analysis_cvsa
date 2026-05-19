function c = to_strcell(x)
% TO_STRCELL  Coerce yamlmatlab string lists into a cell array of char.
    if isempty(x)
        c = {};
    elseif ischar(x)
        c = {strtrim(x)};
    elseif iscell(x)
        c = cellfun(@(e) strtrim(char(e)), x(:), 'UniformOutput', false);
    elseif isstring(x)
        c = cellstr(strtrim(x(:)));
    else
        error('to_strcell:unsupported', 'Cannot coerce class %s to strcell.', class(x));
    end
end
