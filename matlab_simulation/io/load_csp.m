function csp = load_csp(params, paradigm_key)
% LOAD_CSP  Extract the CSP config embedded in the companion YAML under
%           processing_fbcsp_{mi,cvsa}.CspCfg.params for the given key
%           ('mi' or 'cvsa').
%   Returns a struct with:
%       .bands             [n_bands x 2]
%       .selected_channels {n_sel x 1} cellstr
%       .csp_matrices      {n_bands x 1} of [n_components x n_sel]
%       .n_bands, .n_components, .n_selected
    node_name = ['processing_fbcsp_', paradigm_key];
    if ~isfield(params, node_name)
        error('load_csp:nonode', 'YAML has no "%s" section.', node_name);
    end
    p = params.(node_name).CspCfg.params;

    csp.bands             = to_mat(p.bands);
    csp.selected_channels = to_strcell(p.selected_channels);

    raw_mats = p.csp_matrices;
    if ~iscell(raw_mats), raw_mats = {raw_mats}; end
    n = numel(raw_mats);
    csp.csp_matrices = cell(n, 1);
    for b = 1:n
        csp.csp_matrices{b} = to_mat(raw_mats{b});
    end
    csp.n_bands      = n;
    csp.n_components = size(csp.csp_matrices{1}, 1);
    csp.n_selected   = size(csp.csp_matrices{1}, 2);

    if size(csp.bands, 1) ~= n
        error('load_csp:bands_count', ...
              'CSP for %s: %d bands but %d csp_matrices.', paradigm_key, ...
              size(csp.bands, 1), n);
    end
    log_step('load_csp[%s]: %d bands, %d components, %d selected channels', ...
             paradigm_key, csp.n_bands, csp.n_components, csp.n_selected);
end
