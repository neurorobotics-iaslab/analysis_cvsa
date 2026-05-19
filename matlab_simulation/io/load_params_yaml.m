function [params, yaml_path] = load_params_yaml(gdf_path)
% LOAD_PARAMS_YAML  Locate and load the rosparam-dump YAML that the
%   bag_bci recorder writes next to every GDF. Same directory, same
%   basename, .yaml extension.
    [dir_, base, ~] = fileparts(gdf_path);
    yaml_path = fullfile(dir_, [base, '.yaml']);
    if ~exist(yaml_path, 'file')
        error('load_params_yaml:notfound', ...
              'Companion YAML not found:\n  %s', yaml_path);
    end
    log_step('load_params_yaml: loading "%s"', yaml_path);
    params = read_yaml(yaml_path);

    paradigm = params.integrator.paradigm;
    classes  = to_vec(params.integrator.classes);
    log_step('load_params_yaml: paradigm="%s", classes=[%s], samplerate=%d, framerate=%d', ...
             paradigm, num2str(classes), ...
             params.acquisition.samplerate, params.acquisition.framerate);
end
