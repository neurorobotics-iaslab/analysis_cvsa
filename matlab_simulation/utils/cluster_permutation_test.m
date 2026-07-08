function result = cluster_permutation_test(diff_mat, time_axis, n_perm, t_thresh, alpha)
% CLUSTER_PERMUTATION_TEST  Nonparametric cluster-based permutation test
%   (Maris & Oostenveld 2007 style) for a paired difference time series, no
%   Statistics Toolbox required.
%
%   diff_mat  [n_rows x L] paired difference per row (trial, or per-subject
%             mean curve for a group-level test) x time frame. NaN allowed
%             (e.g. trials/subjects with fewer valid frames than others).
%   time_axis [1 x L] or [L x 1] time values in seconds, for reporting only.
%   n_perm    number of sign-flip permutations for the null distribution
%             (default 2000).
%   t_thresh  cluster-FORMING threshold on the per-frame t-like statistic
%             (default 2.0). This only decides candidate cluster boundaries
%             -- the actual corrected significance comes from the permuted
%             null distribution of the max cluster mass, not from this
%             threshold directly.
%   alpha     cluster-level significance threshold (default 0.05).
%
%   Method: at each frame, computes a per-frame t-like statistic (mean over
%   rows / SEM), thresholds it to find contiguous "candidate" clusters
%   (separately for the positive and negative direction), and sums the
%   statistic within each cluster ("cluster mass"). The null distribution is
%   built by randomly sign-flipping ENTIRE ROWS (not individual frames --
%   this preserves the temporal autocorrelation structure within a
%   trial/subject) and recomputing the largest |cluster mass| each time.
%   Each observed cluster's p-value is the fraction of permutations whose
%   largest |cluster mass| is >= that cluster's |mass| -- this controls the
%   family-wise error rate across the whole time course, unlike testing
%   every frame independently.
%
%   Shared by main_hybrid_advantage_probs.m (single-session, trial-level:
%   does P_fused(target)-P_MI(target) show a genuine cluster of influence
%   within the CVSA-influence window?) and main_group_analysis.m
%   (cross-subject, subject-level: using each subject's own mean diff curve
%   as one "row", does the effect generalise across the cohort?).
%
%   Returns a struct:
%     tstat     [1 x L] observed per-frame t-like statistic
%     clusters  struct array (one per candidate cluster): start_idx, end_idx,
%               start_t, end_t, mass, p, mean_delta, cohen_d, ci_lo, ci_hi --
%               the last four collapse each row's OWN mean over just this
%               cluster's time range into one scalar per row, then report:
%               mean_delta = mean of that scalar across rows; cohen_d =
%               paired Cohen's d on that same per-row scalar (row is the
%               unit -- trial for the single-session call, subject for the
%               group-level call); ci_lo/ci_hi = percentile bootstrap 95% CI
%               on mean_delta (2000 resamples of rows, with replacement).
%               This answers "how big is the effect, specifically where the
%               cluster says it's significant" -- a magnitude + uncertainty
%               complement to the cluster's p-value alone.
%     sig_mask  [1 x L] logical, true where the frame belongs to a cluster
%               with p < alpha
%     null_max  [1 x n_perm] max |cluster mass| under each permutation

    if nargin < 3 || isempty(n_perm),   n_perm = 2000; end
    if nargin < 4 || isempty(t_thresh), t_thresh = 2.0; end
    if nargin < 5 || isempty(alpha),    alpha = 0.05;  end

    [n_rows, L] = size(diff_mat);
    time_axis = time_axis(:)';

    obs_t = frame_tstat_local(diff_mat);
    [obs_ranges, obs_mass] = find_clusters_local(obs_t, t_thresh);

    null_max = zeros(1, n_perm);
    for p = 1:n_perm
        signs = sign(rand(n_rows, 1) - 0.5);
        signs(signs == 0) = 1;
        permuted = diff_mat .* signs;
        pt = frame_tstat_local(permuted);
        [~, pmass] = find_clusters_local(pt, t_thresh);
        if isempty(pmass)
            null_max(p) = 0;
        else
            null_max(p) = max(abs(pmass));
        end
    end

    n_obs = numel(obs_mass);
    clusters = struct('start_idx', {}, 'end_idx', {}, 'start_t', {}, 'end_t', {}, 'mass', {}, 'p', {}, ...
                       'mean_delta', {}, 'cohen_d', {}, 'ci_lo', {}, 'ci_hi', {});
    sig_mask = false(1, L);
    for k = 1:n_obs
        idx_range = obs_ranges{k};
        m = obs_mass(k);
        p_val = (sum(null_max >= abs(m)) + 1) / (n_perm + 1);
        clusters(k).start_idx = idx_range(1);
        clusters(k).end_idx   = idx_range(end);
        clusters(k).start_t   = time_axis(idx_range(1));
        clusters(k).end_t     = time_axis(idx_range(end));
        clusters(k).mass      = m;
        clusters(k).p         = p_val;
        if p_val < alpha
            sig_mask(idx_range) = true;
        end

        % ── Magnitude + uncertainty, localised to THIS cluster's own time
        %    range only (not the whole diff_mat) -- one scalar per row via
        %    that row's own mean over idx_range, then row is the unit. ────
        row_scalar = mean(diff_mat(:, idx_range), 2, 'omitnan');
        clusters(k).mean_delta = mean(row_scalar, 'omitnan');
        clusters(k).cohen_d    = local_cohen_d(row_scalar);
        [clusters(k).ci_lo, clusters(k).ci_hi] = local_bootstrap_ci(row_scalar, 2000);
    end

    result = struct();
    result.tstat    = obs_t;
    result.clusters = clusters;
    result.sig_mask = sig_mask;
    result.null_max = null_max;
end

function t = frame_tstat_local(mat)
    n = sum(~isnan(mat), 1);
    m = mean(mat, 1, 'omitnan');
    s = std(mat, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
    t = m ./ s;
    t(n < 2 | s == 0) = 0;   % undefined at this frame -- no evidence, not a cluster seed
end

function [ranges, mass] = find_clusters_local(tvec, thresh)
    L = numel(tvec);
    ranges = {};
    mass = [];
    supra = zeros(1, L);
    supra(tvec >= thresh)  = 1;
    supra(tvec <= -thresh) = -1;
    k = 1;
    while k <= L
        if supra(k) == 0
            k = k + 1;
            continue;
        end
        s = supra(k);
        k2 = k;
        while k2 <= L && supra(k2) == s
            k2 = k2 + 1;
        end
        idx_range = k:(k2-1);
        ranges{end+1} = idx_range; %#ok<AGROW>
        mass(end+1) = sum(tvec(idx_range)); %#ok<AGROW>
        k = k2;
    end
end

function d = local_cohen_d(x)
% LOCAL_COHEN_D  Paired/one-sample Cohen's d = mean(x)/std(x), row is the unit.
    x = x(~isnan(x));
    if numel(x) < 2, d = NaN; return; end
    s = std(x);
    if s == 0, d = NaN; return; end
    d = mean(x) / s;
end

function [ci_lo, ci_hi] = local_bootstrap_ci(x, n_boot)
% LOCAL_BOOTSTRAP_CI  Percentile bootstrap 95% CI on mean(x), resampling
%   rows with replacement -- no Statistics Toolbox required.
    x = x(~isnan(x));
    n = numel(x);
    if n < 2, ci_lo = NaN; ci_hi = NaN; return; end
    boot_means = zeros(n_boot, 1);
    for b = 1:n_boot
        idx = randi(n, n, 1);
        boot_means(b) = mean(x(idx));
    end
    boot_sorted = sort(boot_means);
    lo_i = max(1, round(0.025 * n_boot));
    hi_i = min(n_boot, round(0.975 * n_boot));
    ci_lo = boot_sorted(lo_i);
    ci_hi = boot_sorted(hi_i);
end
