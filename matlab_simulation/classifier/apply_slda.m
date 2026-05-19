function probs = apply_slda(features, slda, csp_bands)
% APPLY_SLDA  Apply log + sLDA decision function + sigmoid to every row of
%   `features`. Mirrors slda_bci/src/slda.py exactly.
%
%   features    [n_chunks x n_features]  raw mean-power features (NaN
%                                        rows before the ringbuffer
%                                        was full are passed through).
%   slda        struct from load_slda
%   csp_bands   [n_bands x 2] band order used at feature-extraction time.
%               If it differs from slda.bands, columns are reordered to
%               match the model's band order.
%
%   probs       [n_chunks x n_classes] sigmoid probabilities (binary
%               sLDA -> two columns that sum to 1). NaN rows pass through.

    n_chunks = size(features, 1);
    n_cls    = numel(slda.classes);
    probs    = NaN(n_chunks, n_cls);

    % Reorder feature columns from csp_bands order to slda.bands order
    perm = match_band_order(csp_bands, slda.bands);   % perm(i) = csp band that becomes slda band i
    n_comp_per_band = slda.n_components;
    col_idx = zeros(1, n_comp_per_band * numel(perm));
    for i = 1:numel(perm)
        col_idx((i-1)*n_comp_per_band + (1:n_comp_per_band)) = ...
            (perm(i)-1) * n_comp_per_band + (1:n_comp_per_band);
    end

    valid = ~any(isnan(features), 2);
    X     = features(valid, col_idx);
    Xlog  = log(X);

    score = Xlog * slda.weights(:) + slda.intercept;     % [n_valid x 1]
    p2    = 1 ./ (1 + exp(-score));                      % P(class 2)
    p1    = 1 - p2;                                      % P(class 1)

    probs(valid, :) = [p1, p2];
    log_step('apply_slda: %d valid frames classified, P(c2) mean=%.3f', ...
             sum(valid), mean(p2));
end


function perm = match_band_order(csp_bands, slda_bands)
% MATCH_BAND_ORDER  For each row in slda_bands, find the matching row in
%   csp_bands (within 1e-3 tol). Returns the perm vector.
    n = size(slda_bands, 1);
    perm = zeros(1, n);
    for i = 1:n
        d = max(abs(csp_bands - slda_bands(i, :)), [], 2);
        hit = find(d < 1e-3, 1);
        if isempty(hit)
            error('apply_slda:band', 'sLDA band [%g %g] not in CSP bands.', ...
                  slda_bands(i, 1), slda_bands(i, 2));
        end
        perm(i) = hit;
    end
end
