function y = read_yaml(path)
% READ_YAML  Thin wrapper around yamlmatlab's ReadYaml.
    if exist('ReadYaml', 'file') == 2
        y = ReadYaml(path);
    elseif ~isempty(which('yaml.ReadYaml'))
        y = yaml.ReadYaml(path);
    else
        error('read_yaml:missing', ...
              'yamlmatlab not on path. Add it from https://github.com/jiri-cigler/yamlmatlab');
    end
end
