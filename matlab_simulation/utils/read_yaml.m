function y = read_yaml(path)
% READ_YAML  Load a YAML file and return it as a nested MATLAB struct/cell.
%   Handles two calling conventions (ReadYaml vs yaml.ReadYaml) to support
%   different versions of yamlmatlab. Raises an error with the download URL
%   if yamlmatlab is not on the path.
%   Used by: load_params_yaml, load_csp, load_slda.
%   Requires: yamlmatlab — https://github.com/jiri-cigler/yamlmatlab
    if exist('ReadYaml', 'file') == 2
        y = ReadYaml(path);
    elseif ~isempty(which('yaml.ReadYaml'))
        y = yaml.ReadYaml(path);
    else
        error('read_yaml:missing', ...
              'yamlmatlab not on path. Add it from https://github.com/jiri-cigler/yamlmatlab');
    end
end
