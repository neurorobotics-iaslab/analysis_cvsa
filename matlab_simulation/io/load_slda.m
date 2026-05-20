function slda = load_slda(params, paradigm_key)
% LOAD_SLDA  Load the sLDA model. The companion YAML references the model
%           file by path under slda_node_{mi,cvsa}.path_slda_model; we
%           read that file (sLDACfg.params).
%   Returns a struct with:
%       .weights    [1 x n_features]
%       .intercept  scalar
%       .bands      [n_bands x 2]
%       .classes    [1 x n_classes] numeric
%       .n_components, .n_features
    node_name = ['slda_node_', paradigm_key];
    if ~isfield(params, node_name)
        error('load_slda:nonode', 'YAML has no "%s" section.', node_name);
    end
    model_path = params.(node_name).path_slda_model;
    if ~exist(model_path, 'file')
        error('load_slda:notfound', 'sLDA model not found:\n  %s', model_path);
    end
    log_step('load_slda[%s]: reading "%s"', paradigm_key, model_path);

    y = read_yaml(model_path);
    p = y.sLDACfg.params;

    W              = to_mat(p.slda_weights);
    slda.weights   = W(:).';
    slda.intercept = to_vec(p.slda_intercept);
    slda.bands     = to_mat(p.bands);

    cls = p.classes;
    if iscell(cls) && ~isempty(cls) && (isstring(cls{1}) || ischar(cls{1}))
        s = to_strcell(cls);
        slda.classes = cellfun(@str2double, s(:)).';
    else
        slda.classes = to_vec(cls);
    end
    slda.n_components = numel(to_vec(p.selected_components_indices));
    slda.n_features   = numel(slda.weights);

    % Safe parsing of Platt calibration parameters (default to standard sigmoid: a=1.0, b=0.0)
    if isfield(p, 'slda_calibrated_weights') && ~isempty(p.slda_calibrated_weights)
        slda.platt_a = to_vec(p.slda_calibrated_weights);
        slda.platt_a = slda.platt_a(1);
    else
        slda.platt_a = 1.0;
    end
    if isfield(p, 'slda_calibrated_intercept') && ~isempty(p.slda_calibrated_intercept)
        slda.platt_b = to_vec(p.slda_calibrated_intercept);
        slda.platt_b = slda.platt_b(1);
    else
        slda.platt_b = 0.0;
    end

    log_step('load_slda[%s]: %d features (%d comp x %d bands), classes=[%s], platt_a=%.3f, platt_b=%.3f', ...
             paradigm_key, slda.n_features, slda.n_components, ...
             size(slda.bands, 1), num2str(slda.classes), slda.platt_a, slda.platt_b);
end
