%% MAIN_GROUP_ANALYSIS  Multi-subject aggregation of saved analysis .mat files.
%
%   Recursively scans a root folder (e.g. recordings/) for the .mat summaries
%   already saved by other scripts in this package:
%     <root>/<subject>/.../analysis_results/session_overview/session_summary.mat
%       (from main_session_overview)
%     <root>/<subject>/.../analysis_results/hybrid_advantage_integ/counterfactual_summary.mat
%       (from main_hybrid_advantage_integ)
%     <root>/<subject>/.../analysis_results/hybrid_advantage_probs/hybrid_advantage_probs_summary.mat
%       (from main_hybrid_advantage_probs)
%     <root>/<subject>/.../analysis_results/topo_erders/topo_erders_summary.mat
%       (from topo_erders.m, in analysis_gdf/)
%     <root>/<subject>/.../analysis_results/roc_analysis/roc_summary.mat
%       (from main_roc_analysis)
%
%   Subject ID is inferred from the folder structure: the first path
%   component directly under the selected root folder.
%
%   Multiple .mat files for the same subject (e.g. several evaluation
%   sessions) are pooled: hit/miss/to counts are summed, accuracy is
%   recomputed from the pooled counts; counterfactual accuracy is pooled
%   weighted by n_trials; fusion-mechanism and ERD/ERS metrics are averaged
%   (frame counts summed); ROC is built by RE-POOLING that subject's raw
%   per-frame (score, label) pairs across all their sessions (valid, since
%   the same deployed classifier/calibration produced every one of a
%   subject's frames) and recomputing one ROC curve from the pooled set.
%
%   TWO accuracy conventions are computed and kept side by side EVERYWHERE
%   (per-subject arrays, console tables, group-level tests, figures), always
%   clearly labelled so the two are never confused:
%     *_acc      = n_hit / n_total             (TIMEOUT counted as a failure)
%     *_acc_dec  = n_hit / (n_hit + n_miss)     ("decided trials only": how
%                  good is the classifier/integrator when it actually
%                  reaches a decision, independent of how often it doesn't)
%     *_fp_rate  = n_miss / n_total             ("false positive": the WRONG
%                  class's threshold was reached with confidence -- distinct
%                  from TIMEOUT, where neither threshold was reached)
%     *_to_rate  = n_to / n_total                (already existed)
%   Both real-session and counterfactual (simulated) data get this same
%   4-way split, plus a trial-pooled GRAND row/value.
%
%   Nine figures:
%     Fig 1  — per-subject real-session accuracy (TIMEOUT = fail) + TIMEOUT
%              rate per paradigm (MI/CVSA/Hybrid), grouped bars, plus a
%              final "GRAND" group built by POOLING every subject's trials
%              together (not the mean of the per-subject bars) -- the
%              population-level estimate, which can differ from mean-of-
%              means when subjects contribute unequal trial counts.
%     Fig 1b — per-subject real-session TIMES (TTH / T-miss / T-timeout,
%              one panel per paradigm), same per-subject + trial-pooled
%              GRAND group convention as Fig 1.
%     Fig 1c — real-session performance BREAKDOWN: decided-trial accuracy
%              (TIMEOUT excluded) on top, false-positive rate + timeout rate
%              below, per paradigm, per subject + GRAND. Gold star / red
%              triangle mark the BEST/WORST subject FOR THAT PARADIGM,
%              ranked by decided-trial accuracy (ties shown together).
%     Fig 2  — per-subject counterfactual (SIMULATED) accuracy (TIMEOUT =
%              fail; Hybrid/MI-only/CVSA-only) and Hybrid-MI / Hybrid-CVSA
%              accuracy deltas, with a group-level sign-flip permutation
%              test across SUBJECTS (not trials). SIMULATED means: the SAME
%              integrator/buffer/thresholds are driven counterfactually by
%              the Bayesian-fused, raw MI-only, or raw CVSA-only signal on
%              the SAME Hybrid-session trials -- isolating the fusion
%              algorithm's own effect from any real-session differences in
%              trial count, timing, or thresholds (this is NOT real data).
%     Fig 2b — simulated (counterfactual) performance BREAKDOWN: same layout
%              as Fig 1c, but for the 3 counterfactual streams instead of
%              the 3 real paradigms; best/worst marked PER STREAM.
%     Fig 3  — per-subject CVSA-fusion advantage (mean P_fused-P_MI) and
%              rescue-vs-cost frame-level effect, with the same group-level
%              sign-flip test across subjects — tests whether the fusion
%              mechanism described per-file by main_hybrid_advantage_probs
%              is consistent across the whole cohort, not just one subject.
%     Fig 3b — (console only, no new figure) group-level cluster-based
%              permutation test on P_fused-P_MI (see cluster_time_grid/
%              cluster_diff_subj), plus a Stouffer meta-analytic combination
%              of each subject's OWN within-session sign-flip test p-value
%              (from main_hybrid_advantage_probs.m) -- the sign-flip test on
%              fus_adv_subj above has a permutation null of only 2^n_subj
%              sign patterns, so with few subjects it structurally cannot
%              reach p<0.05 regardless of effect size; the meta-analysis
%              sidesteps that floor by combining each subject's own
%              well-powered (many-trial) within-session evidence instead of
%              first collapsing each subject to one mean.
%     Fig 4  — per-subject ERD/ERS class discrimination and CSP-weight
%              correlation (Pearson r), per task, with a group-level
%              sign-flip test on r across subjects — tests whether the
%              neurophysiological grounding shown per-file by topo_erders.m
%              holds across the whole cohort.
%     Fig 5  — classifier-level ROC (pre-integrator, see main_roc_analysis.m),
%              one panel per paradigm present: each subject's own re-pooled
%              ROC curve (thin) plus the CROSS-SUBJECT macro-average curve
%              (thick, ± SEM band) -- per-subject curves are interpolated
%              onto a common FPR grid and averaged, since pooling raw scores
%              directly across subjects would mix different classifiers'
%              probability scales. A group-level sign-flip test on
%              (subject AUC - 0.5) reports whether the classifier is
%              consistently above chance across the cohort.
%     Fig 6  — cross-subject correlates of real Hybrid performance: does a
%              subject's counterfactual advantage / fusion-mechanism strength
%              / neurophysiological (ERD-CSP) grounding / CVSA-only "quality"
%              predict how well they actually do in real sessions? Pearson r
%              (permutation p) across subjects, one scatter panel per
%              candidate predictor. The 4th panel (CVSA-only decided
%              accuracy vs Hybrid-MI counterfactual advantage) directly
%              supports the "CVSA helps even when imperfect" thesis: if the
%              advantage holds even for subjects with weak CVSA-only
%              accuracy, the benefit isn't contingent on CVSA being good.
%
%   Statistics: besides the per-figure group-level sign-flip tests above, a
%   subject-level Friedman omnibus test (rank-based, no Statistics Toolbox)
%   plus pairwise sign-flip tests (Hybrid-MI, Hybrid-CVSA, MI-CVSA) are run
%   on REAL-session accuracy, in BOTH conventions (TIMEOUT=fail and decided-
%   trials-only) — the behavioural claim that matters most for the paper,
%   complementing the counterfactual-only tests that already existed. The
%   counterfactual pairwise deltas (Hybrid-MI-only, Hybrid-CVSA-only) get
%   the same decided-trials-only twin. Best/worst subject rankings are
%   printed per paradigm (real) and per stream (simulated), ties included.
%
%   Console: per-subject tables (both accuracy conventions + false-positive
%   rate, timeout rate, TTH/T-miss/T-timeout) + the trial-pooled GRAND row +
%   group-level paired statistics (both conventions) + per-paradigm/per-
%   stream subject rankings + cross-subject performance correlates.
%
%   Saves group_summary.mat (per-subject + grand-average aggregated table)
%   and all SVGs under <root>/group_analysis/.

function main_group_analysis(root_dir, subjects_filter, show_figures)
%   Callable as a function:
%     main_group_analysis()                                  % default recordings root, ALL subjects
%     main_group_analysis(root_dir)                           % given root, ALL subjects
%     main_group_analysis(root_dir, {'a1','a2','a3'})        % only these subjects
%     main_group_analysis(root_dir, {'a1','a2','a3'}, true)  % also show figures on screen
%     main_group_analysis([], {'a1','a2'})                    % default root, subset of subjects
%   No GUI folder picker: root_dir defaults to DEFAULT_ROOT below (edit it
%   directly, or pass a path) rather than prompting. subjects_filter is an
%   explicit opt-in filter (cellstr of subject IDs); omit it (or pass {}) to
%   use every subject found under the root.

DEFAULT_ROOT = '/home/paolo/bci_vr_ws/recordings';

if nargin < 3, show_figures = false; end
if nargin < 2, subjects_filter = {}; end
if nargin < 1 || isempty(root_dir), root_dir = DEFAULT_ROOT; end
SHOW_FIGURES = show_figures;
if SHOW_FIGURES, fig_vis = 'on'; else, fig_vis = 'off'; end

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir,'utils'));

session_files = dir(fullfile(root_dir, '**', 'session_overview', 'session_summary.mat'));
cf_files      = dir(fullfile(root_dir, '**', 'hybrid_advantage_integ', 'counterfactual_summary.mat'));
probs_files   = dir(fullfile(root_dir, '**', 'hybrid_advantage_probs', 'hybrid_advantage_probs_summary.mat'));
erd_files     = dir(fullfile(root_dir, '**', 'topo_erders', 'topo_erders_summary.mat'));
roc_files     = dir(fullfile(root_dir, '**', 'roc_analysis', 'roc_summary.mat'));

% ── Optional subject filter: only include the requested subjects ────────
if ~isempty(subjects_filter)
    session_files = filter_files_by_subject(session_files, root_dir, subjects_filter);
    cf_files      = filter_files_by_subject(cf_files,      root_dir, subjects_filter);
    probs_files   = filter_files_by_subject(probs_files,   root_dir, subjects_filter);
    erd_files     = filter_files_by_subject(erd_files,     root_dir, subjects_filter);
    roc_files     = filter_files_by_subject(roc_files,     root_dir, subjects_filter);
    fprintf('Subject filter requested: %s\n', strjoin(subjects_filter, ', '));
end

fprintf('Found %d session_summary.mat, %d counterfactual_summary.mat, %d hybrid_advantage_probs_summary.mat, %d topo_erders_summary.mat, %d roc_summary.mat under %s\n', ...
        numel(session_files), numel(cf_files), numel(probs_files), numel(erd_files), numel(roc_files), root_dir);
if isempty(session_files) && isempty(cf_files) && isempty(probs_files) && isempty(erd_files) && isempty(roc_files)
    error('main_group_analysis:nodata', 'No .mat summaries found. Run main_session_overview / main_hybrid_advantage_integ / main_hybrid_advantage_probs / topo_erders / main_roc_analysis first.');
end

% ── Flatten real-session per-file records, tagged with subject ─────────────
flat_real = struct('subject', {}, 'paradigm', {}, 'n_hit', {}, 'n_miss', {}, 'n_to', {}, ...
                    'tth_vals', {}, 't_miss_vals', {}, 'to_vals', {});
for i = 1:numel(session_files)
    fpath = fullfile(session_files(i).folder, session_files(i).name);
    subj  = subject_of(fpath, root_dir);
    s = load(fpath); sf = s.summary_file;
    for k = 1:numel(sf)
        idx = numel(flat_real) + 1;
        flat_real(idx).subject     = subj;
        flat_real(idx).paradigm    = sf(k).paradigm;
        flat_real(idx).n_hit       = sf(k).n_hit;
        flat_real(idx).n_miss      = sf(k).n_miss;
        flat_real(idx).n_to        = sf(k).n_to;
        flat_real(idx).tth_vals    = sf(k).tth_vals;
        flat_real(idx).t_miss_vals = sf(k).t_miss_vals;
        flat_real(idx).to_vals     = sf(k).to_vals;
    end
end

% ── Flatten counterfactual per-file records, tagged with subject ───────────
flat_cf = struct('subject', {}, 'acc', {}, 'n_trials', {}, 't_hit_mean', {}, ...
                  'n_hit', {}, 'n_miss', {}, 'n_to', {}, ...
                  'buf_adv_mi_win_mean', {}, 'buf_adv_mi_win_pval', {}, 'buf_adv_mi_win_n_trials', {});
for i = 1:numel(cf_files)
    fpath = fullfile(cf_files(i).folder, cf_files(i).name);
    subj  = subject_of(fpath, root_dir);
    s = load(fpath); cf = s.counterfactual;
    idx = numel(flat_cf) + 1;
    flat_cf(idx).subject    = subj;
    flat_cf(idx).acc        = cf.acc;          % [1x3] Hybrid, MI-only, CVSA-only
    flat_cf(idx).n_trials   = cf.n_trials;
    flat_cf(idx).t_hit_mean = cf.t_hit_mean;   % [1x3]
    flat_cf(idx).n_hit      = cf.n_hit;        % [1x3] Hybrid, MI-only, CVSA-only
    flat_cf(idx).n_miss     = cf.n_miss;       % [1x3]
    flat_cf(idx).n_to       = cf.n_to;         % [1x3]
    % buf_adv_mi_win_* only present in counterfactual_summary.mat regenerated
    % after the cvsa_influence-window restricted test was added -- guarded
    % with isfield so older summaries don't crash the loader, same
    % backward-compatibility convention as flat_probs.cluster_diff_mean below.
    if isfield(cf, 'buf_adv_mi_win_pval')
        flat_cf(idx).buf_adv_mi_win_mean     = cf.buf_adv_mi_win_mean;
        flat_cf(idx).buf_adv_mi_win_pval     = cf.buf_adv_mi_win_pval;
        flat_cf(idx).buf_adv_mi_win_n_trials = cf.buf_adv_mi_win_n_trials;
    end
end

% ── Flatten hybrid_advantage_probs per-file records, tagged with subject ───
% cluster_* fields are only present in hybrid_advantage_probs_summary.mat
% files regenerated after the cluster-based permutation test was added --
% guarded with isfield so older summaries (from before that change) don't
% crash the loader, they just don't contribute to the group cluster test.
flat_probs = struct('subject', {}, 'fus_adv_mean', {}, 'n_rescued_fr', {}, 'n_hurt_fr', {}, ...
                     'rescue_delta', {}, 'cost_delta', {}, ...
                     'cluster_time_axis', {}, 'cluster_diff_mean', {}, 'cluster_n_trials', {}, ...
                     'fus_adv_pval', {}, 'fus_adv_cohend', {}, 'fus_adv_n_trials', {});
for i = 1:numel(probs_files)
    fpath = fullfile(probs_files(i).folder, probs_files(i).name);
    subj  = subject_of(fpath, root_dir);
    s = load(fpath); ps = s.probs_summary;
    for k = 1:numel(ps)
        idx = numel(flat_probs) + 1;
        flat_probs(idx).subject      = subj;
        flat_probs(idx).fus_adv_mean = ps(k).fus_adv_mean;
        flat_probs(idx).n_rescued_fr = ps(k).n_rescued_fr;
        flat_probs(idx).n_hurt_fr    = ps(k).n_hurt_fr;
        flat_probs(idx).rescue_delta = ps(k).rescue_delta;
        flat_probs(idx).cost_delta   = ps(k).cost_delta;
        if isfield(ps, 'cluster_diff_mean')
            flat_probs(idx).cluster_time_axis = ps(k).cluster_time_axis;
            flat_probs(idx).cluster_diff_mean = ps(k).cluster_diff_mean;
            flat_probs(idx).cluster_n_trials  = ps(k).cluster_n_trials;
        end
        % fus_adv_pval/fus_adv_cohend are only present in summaries
        % regenerated after the per-subject sign-flip test was added --
        % same isfield-guarded backward compatibility as the cluster fields.
        if isfield(ps, 'fus_adv_pval')
            flat_probs(idx).fus_adv_pval    = ps(k).fus_adv_pval;
            flat_probs(idx).fus_adv_cohend  = ps(k).fus_adv_cohend;
            flat_probs(idx).fus_adv_n_trials = ps(k).n_def;
        end
    end
end

% ── Flatten topo_erders per-band records, tagged with subject ──────────────
flat_erd = struct('subject', {}, 'paradigm', {}, 'band_origin', {}, 'discrimination', {}, 'csp_r', {});
for i = 1:numel(erd_files)
    fpath = fullfile(erd_files(i).folder, erd_files(i).name);
    subj  = subject_of(fpath, root_dir);
    s = load(fpath); es = s.erd_summary;
    for k = 1:numel(es)
        idx = numel(flat_erd) + 1;
        flat_erd(idx).subject        = subj;
        flat_erd(idx).paradigm       = es(k).paradigm;
        flat_erd(idx).band_origin    = es(k).band_origin;
        flat_erd(idx).discrimination = es(k).discrimination;
        flat_erd(idx).csp_r          = es(k).csp_r;
    end
end

% ── Flatten ROC per-paradigm records, tagged with subject (raw pooled
%    scores/labels, so they can be RE-POOLED per subject below) ────────────
flat_roc = struct('subject', {}, 'paradigm', {}, 'scores', {}, 'labels', {});
for i = 1:numel(roc_files)
    fpath = fullfile(roc_files(i).folder, roc_files(i).name);
    subj  = subject_of(fpath, root_dir);
    s = load(fpath); rs = s.roc_summary;
    for k = 1:numel(rs.roc_by_group)
        idx = numel(flat_roc) + 1;
        flat_roc(idx).subject  = subj;
        flat_roc(idx).paradigm = rs.roc_by_group(k).name;
        flat_roc(idx).scores   = rs.roc_by_group(k).scores;
        flat_roc(idx).labels   = rs.roc_by_group(k).labels;
    end
end

subjects = unique([{flat_real.subject}, {flat_cf.subject}, {flat_probs.subject}, {flat_erd.subject}, {flat_roc.subject}]);
n_subj   = numel(subjects);
pars     = {'mi', 'cvsa', 'hybrid'};
par_labels = {'MI', 'CVSA', 'Hybrid'};
cf_labels  = {'Hybrid', 'MI-only', 'CVSA-only'};   % counterfactual stream order, matches counterfactual.acc [1x3]

% ── Per-subject real accuracy/timeout-rate/times per paradigm (pooled across
%    that subject's files) ─────────────────────────────────────────────────
% TWO accuracy conventions are computed, kept side by side and clearly
% labelled everywhere (console headers + figure axis labels), rather than
% picking one:
%   real_acc     = n_hit / n_total            (TIMEOUT counted as a failure)
%   real_acc_dec = n_hit / (n_hit + n_miss)   (TIMEOUT excluded -- "decided
%                  trials only": how good is the classifier/integrator when
%                  it actually reaches a decision)
%   real_fp_rate = n_miss / n_total           ("false positive": the WRONG
%                  class's threshold was reached with confidence -- distinct
%                  from TIMEOUT, where neither threshold was reached)
real_acc     = nan(n_subj, 3);
real_acc_dec = nan(n_subj, 3);
real_fp_rate = nan(n_subj, 3);
real_to_rate = nan(n_subj, 3);
real_n_hit   = zeros(n_subj, 3);
real_n_miss  = zeros(n_subj, 3);
real_to_n    = zeros(n_subj, 3);
real_n       = zeros(n_subj, 3);
real_tth     = nan(n_subj, 3);
real_t_miss  = nan(n_subj, 3);
real_to_time = nan(n_subj, 3);
for si = 1:n_subj
    for pi = 1:3
        mask = strcmp({flat_real.subject}, subjects{si}) & strcmp({flat_real.paradigm}, pars{pi});
        if ~any(mask), continue; end
        nh = sum([flat_real(mask).n_hit]);
        nm = sum([flat_real(mask).n_miss]);
        nt = sum([flat_real(mask).n_to]);
        ntot = nh + nm + nt;
        real_n_hit(si,pi)   = nh;
        real_n_miss(si,pi)  = nm;
        real_acc(si,pi)     = nh / max(1, ntot);
        if (nh+nm) > 0, real_acc_dec(si,pi) = nh / (nh+nm); end
        real_fp_rate(si,pi) = nm / max(1, ntot);
        real_to_rate(si,pi) = nt / max(1, ntot);
        real_to_n(si,pi)    = nt;
        real_n(si,pi)       = ntot;
        tth_all = [flat_real(mask).tth_vals];
        if ~isempty(tth_all), real_tth(si,pi) = mean(tth_all); end
        tmiss_all = [flat_real(mask).t_miss_vals];
        if ~isempty(tmiss_all), real_t_miss(si,pi) = mean(tmiss_all); end
        toT_all = [flat_real(mask).to_vals];
        if ~isempty(toT_all), real_to_time(si,pi) = mean(toT_all); end
    end
end

% ── Grand average: pooled across ALL subjects' trials (not the mean of the
%    per-subject means above) — the population-level estimate, which can
%    differ from mean-of-subject-means when subjects contribute unequal
%    trial counts. ────────────────────────────────────────────────────────
grand_acc     = nan(1, 3);
grand_acc_dec = nan(1, 3);
grand_fp_rate = nan(1, 3);
grand_to_rate = nan(1, 3);
grand_to_n    = zeros(1, 3);
grand_n       = zeros(1, 3);
grand_tth     = nan(1, 3);
grand_t_miss  = nan(1, 3);
grand_to_time = nan(1, 3);
for pi = 1:3
    mask = strcmp({flat_real.paradigm}, pars{pi});
    if ~any(mask), continue; end
    nh = sum([flat_real(mask).n_hit]);
    nm = sum([flat_real(mask).n_miss]);
    nt = sum([flat_real(mask).n_to]);
    ntot = nh + nm + nt;
    grand_acc(pi)     = nh / max(1, ntot);
    if (nh+nm) > 0, grand_acc_dec(pi) = nh / (nh+nm); end
    grand_fp_rate(pi) = nm / max(1, ntot);
    grand_to_rate(pi) = nt / max(1, ntot);
    grand_to_n(pi)    = nt;
    grand_n(pi)       = ntot;
    tth_all = [flat_real(mask).tth_vals];
    if ~isempty(tth_all), grand_tth(pi) = mean(tth_all); end
    tmiss_all = [flat_real(mask).t_miss_vals];
    if ~isempty(tmiss_all), grand_t_miss(pi) = mean(tmiss_all); end
    toT_all = [flat_real(mask).to_vals];
    if ~isempty(toT_all), grand_to_time(pi) = mean(toT_all); end
end

% ── Per-subject counterfactual accuracy (pooled from raw n_hit/n_miss/n_to
%    counts, summed across that subject's hybrid files) ────────────────────
% Same two-convention split as the real-session block above: cf_acc (TIMEOUT
% = fail) and cf_acc_dec (decided trials only) + cf_fp_rate + cf_to_rate.
cf_n_hit  = zeros(n_subj, 3);   % Hybrid, MI-only, CVSA-only
cf_n_miss = zeros(n_subj, 3);
cf_n_to   = zeros(n_subj, 3);
cf_n      = zeros(n_subj, 1);   % total counterfactual trials (same underlying trials for all 3 streams)
for si = 1:n_subj
    mask = strcmp({flat_cf.subject}, subjects{si});
    if ~any(mask), continue; end
    rows = flat_cf(mask);
    cf_n_hit(si,:)  = sum(vertcat(rows.n_hit),  1);
    cf_n_miss(si,:) = sum(vertcat(rows.n_miss), 1);
    cf_n_to(si,:)   = sum(vertcat(rows.n_to),   1);
    cf_n(si) = sum([rows.n_trials]);
end
cf_tot     = cf_n_hit + cf_n_miss + cf_n_to;         % [n_subj x 3], == cf_n repeated across columns
cf_acc     = cf_n_hit ./ max(1, cf_tot);             % Hybrid, MI-only, CVSA-only -- TIMEOUT counted as fail
cf_fp_rate = cf_n_miss ./ max(1, cf_tot);
cf_to_rate = cf_n_to   ./ max(1, cf_tot);
cf_acc_dec = nan(n_subj, 3);                         % decided trials only
dec_denom_cf = cf_n_hit + cf_n_miss;
has_dec_cf   = dec_denom_cf > 0;
cf_acc_dec(has_dec_cf) = cf_n_hit(has_dec_cf) ./ dec_denom_cf(has_dec_cf);

% ── GRAND (trial-pooled across ALL subjects) for the simulated/counterfactual
%    streams, same convention as the real-session GRAND above ─────────────
if isempty(flat_cf)
    grand_cf_n_hit  = zeros(1, 3);
    grand_cf_n_miss = zeros(1, 3);
    grand_cf_n_to   = zeros(1, 3);
else
    grand_cf_n_hit  = sum(vertcat(flat_cf.n_hit),  1);   % [1x3]
    grand_cf_n_miss = sum(vertcat(flat_cf.n_miss), 1);
    grand_cf_n_to   = sum(vertcat(flat_cf.n_to),   1);
end
grand_cf_tot    = grand_cf_n_hit + grand_cf_n_miss + grand_cf_n_to;
grand_cf_acc     = grand_cf_n_hit  ./ max(1, grand_cf_tot);
grand_cf_fp_rate = grand_cf_n_miss ./ max(1, grand_cf_tot);
grand_cf_to_rate = grand_cf_n_to   ./ max(1, grand_cf_tot);
grand_cf_acc_dec = nan(1, 3);
grand_dec_denom_cf = grand_cf_n_hit + grand_cf_n_miss;
has_dec_grand_cf   = grand_dec_denom_cf > 0;
grand_cf_acc_dec(has_dec_grand_cf) = grand_cf_n_hit(has_dec_grand_cf) ./ grand_dec_denom_cf(has_dec_grand_cf);

% ── Per-subject counterfactual (INTEGRATOR) CVSA-help significance: mean
%    buf_adv_mi_win (Hybrid-MI buffer advantage, restricted to the first
%    cvsa_influence seconds of CF -- see main_hybrid_advantage_integ.m)
%    across that subject's files, plus a Stouffer meta-analytic combination
%    of each file's OWN within-session sign-flip p-value. Same
%    floor-avoidance rationale as fus_adv_meta_p_subj further below: a
%    group-level sign-flip test on the per-subject MEAN has only 2^n_subj
%    achievable p-values (with a handful of subjects it structurally cannot
%    reach p<0.05), while combining each subject's own well-powered
%    (many-trial) within-session evidence has no such floor and remains
%    valid at any cohort size. ───────────────────────────────────────────────
buf_adv_mi_win_subj       = nan(n_subj, 1);   % per-subject mean delta (files averaged)
buf_adv_mi_win_meta_p_subj = nan(n_subj, 1);  % per-subject one-sided combined p ("CVSA helps")
buf_adv_mi_win_meta_n_subj = zeros(n_subj, 1);
has_bufwin = ~cellfun(@isempty, {flat_cf.buf_adv_mi_win_pval});
if any(has_bufwin)
    rows_bw_all = flat_cf(has_bufwin);
    for si = 1:n_subj
        mask = strcmp({rows_bw_all.subject}, subjects{si});
        if ~any(mask), continue; end
        rows_bw = rows_bw_all(mask);
        buf_adv_mi_win_subj(si) = mean([rows_bw.buf_adv_mi_win_mean], 'omitnan');
        p2 = [rows_bw.buf_adv_mi_win_pval];        % two-sided, per session
        d  = [rows_bw.buf_adv_mi_win_mean];        % session-level sign, for direction
        nt = [rows_bw.buf_adv_mi_win_n_trials];    % session trial counts (Stouffer weight)
        p1 = p2 / 2;
        p1(d <= 0) = 1 - p1(d <= 0);                % wrong-direction sessions -> weak one-sided p
        [~, buf_adv_mi_win_meta_p_subj(si)] = stouffer_combine_local(p1, sqrt(max(nt, 1)));
        buf_adv_mi_win_meta_n_subj(si) = sum(nt);
    end
end
[~, p_bufadvwin_meta_grp] = stouffer_combine_local(buf_adv_mi_win_meta_p_subj);   % equal weight per subject

% ── Per-subject CVSA-fusion mechanism (mean across files) ──────────────────
fus_adv_subj    = nan(n_subj, 1);
rescue_subj     = nan(n_subj, 1);
cost_subj       = nan(n_subj, 1);
net_rescue_subj = nan(n_subj, 1);   % n_rescued_fr - n_hurt_fr, summed across files
for si = 1:n_subj
    mask = strcmp({flat_probs.subject}, subjects{si});
    if ~any(mask), continue; end
    rows = flat_probs(mask);
    fus_adv_subj(si) = mean([rows.fus_adv_mean], 'omitnan');
    rescue_subj(si)  = mean([rows.rescue_delta], 'omitnan');
    cost_subj(si)    = mean([rows.cost_delta],   'omitnan');
    net_rescue_subj(si) = sum([rows.n_rescued_fr]) - sum([rows.n_hurt_fr]);
end

% ── Meta-analytic combination of the within-session sign-flip tests on
%    fus_adv (main_hybrid_advantage_probs.m's fus_adv_pval/fus_adv_mean),
%    combined via Stouffer's Z method into one group-level p-value. This
%    exists alongside the group sign-flip test above (on fus_adv_subj) for a
%    reason: that test's own permutation null only has 2^n_subj sign
%    patterns, so with a handful of subjects its smallest achievable
%    p-value is floored well above 0.05 no matter how strong the effect --
%    it cannot detect significance by construction, not because the effect
%    is weak. Combining each subject's OWN within-session p-value (already
%    well-powered, since it uses that subject's full trial count) sidesteps
%    that floor. Still valid and still useful at large n_subj -- it uses
%    the full continuous per-subject evidence rather than collapsing each
%    subject to a single mean first, so it stays a meaningful complement to
%    the group sign-flip test, not just a workaround for small cohorts.
%    Two-step: (1) within a subject, combine that subject's own session-
%    level p-values weighted by sqrt(n_trials) (a subject with more than one
%    hybrid recording gets one combined p across sessions); (2) across
%    subjects, combine with EQUAL weight per subject (matching every other
%    group-level test in this script: subject is the unit). Both steps
%    convert the two-sided per-session p to a one-sided p in the
%    hypothesised "CVSA helps" direction (fus_adv_mean>0) first -- Stouffer's
%    method (like Fisher's) is only valid for combining evidence pointing
%    the same direction; a session pointing the wrong way gets a large
%    (weak) one-sided p rather than being dropped. ──────────────────────────
fus_adv_meta_p_subj = nan(n_subj, 1);   % combined one-sided p per subject
fus_adv_meta_n_subj = zeros(n_subj, 1); % total trials contributing, per subject
has_meta = ~cellfun(@isempty, {flat_probs.fus_adv_pval});
if any(has_meta)
    rows_meta_all = flat_probs(has_meta);
    for si = 1:n_subj
        mask = strcmp({rows_meta_all.subject}, subjects{si});
        if ~any(mask), continue; end
        rows_m = rows_meta_all(mask);
        p2 = [rows_m.fus_adv_pval];        % two-sided, per session
        d  = [rows_m.fus_adv_mean];        % session-level sign, for direction
        nt = [rows_m.fus_adv_n_trials];    % session trial counts (Stouffer weight)
        p1 = p2 / 2;
        p1(d <= 0) = 1 - p1(d <= 0);        % wrong-direction sessions -> weak one-sided p
        [~, fus_adv_meta_p_subj(si)] = stouffer_combine_local(p1, sqrt(max(nt, 1)));
        fus_adv_meta_n_subj(si) = sum(nt);
    end
end
[~, p_fusadv_meta_grp] = stouffer_combine_local(fus_adv_meta_p_subj);   % equal weight per subject

% ── Per-subject mean CVSA-fusion-effect curve, for the group-level (cross-
%    subject) cluster-based permutation test below: average that subject's
%    per-file mean P_fused(target)-P_MI(target) curves (each already
%    averaged over trials) onto a common time grid via interp1 -- same
%    onto-a-common-grid pattern as the ROC cross-subject averaging further
%    down. Requires hybrid_advantage_probs_summary.mat regenerated after the
%    cluster test was added (older summaries are skipped, not an error). ──
has_cluster = ~cellfun(@isempty, {flat_probs.cluster_diff_mean});
CLUSTER_TIME_GRID = [];
cluster_diff_subj = nan(n_subj, 0);
n_cluster_subj = 0;
if any(has_cluster)
    rows_c = flat_probs(has_cluster);
    t_max = min(cellfun(@(v) v(end), {rows_c.cluster_time_axis}));
    N_CLUSTER_GRID = 60;
    CLUSTER_TIME_GRID = linspace(0, t_max, N_CLUSTER_GRID);
    cluster_diff_subj = nan(n_subj, N_CLUSTER_GRID);
    for si = 1:n_subj
        mask = strcmp({rows_c.subject}, subjects{si});
        if ~any(mask), continue; end
        rows_si = rows_c(mask);
        curves = nan(numel(rows_si), N_CLUSTER_GRID);
        for r = 1:numel(rows_si)
            curves(r,:) = interp1(rows_si(r).cluster_time_axis, rows_si(r).cluster_diff_mean, CLUSTER_TIME_GRID, 'linear');
        end
        cluster_diff_subj(si,:) = mean(curves, 1, 'omitnan');
    end
    n_cluster_subj = sum(any(~isnan(cluster_diff_subj), 2));
end

% ── Per-subject ERD/ERS discrimination + CSP-weight correlation, per task ──
erd_discrim_subj = nan(n_subj, 3);   % mean over that subject's bands, per paradigm
erd_r_subj       = nan(n_subj, 3);
for si = 1:n_subj
    for pi = 1:3
        mask = strcmp({flat_erd.subject}, subjects{si}) & strcmp({flat_erd.paradigm}, pars{pi});
        if ~any(mask), continue; end
        rows = flat_erd(mask);
        erd_discrim_subj(si,pi) = mean([rows.discrimination], 'omitnan');
        erd_r_subj(si,pi)       = mean([rows.csp_r], 'omitnan');
    end
end

% ── Per-subject ROC: RE-POOL that subject's raw (score, label) frames across
%    all their sessions (valid -- same deployed classifier/calibration
%    produced every frame) and recompute one ROC curve + AUC per paradigm.
%    Then macro-average those per-subject curves across the cohort by
%    interpolating each onto a common FPR grid -- pooling raw scores ACROSS
%    subjects/classifiers directly would not be valid (different
%    classifiers/calibrations), exactly as in main_roc_analysis.m's
%    "ALL paradigms" comparison figure. ──────────────────────────────────────
n_perm = 2000;   % permutations for every sign-flip test below (defined here so
                 % the ROC group-level test can use it before the later tests do)
ROC_FPR_GRID = linspace(0, 1, 101);
roc_auc_subj      = nan(n_subj, 3);
roc_tpr_subj_grid = nan(n_subj, 3, numel(ROC_FPR_GRID));
for si = 1:n_subj
    for pi = 1:3
        mask = strcmp({flat_roc.subject}, subjects{si}) & strcmp({flat_roc.paradigm}, pars{pi});
        if ~any(mask), continue; end
        rows = flat_roc(mask);
        scores_all = double(vertcat(rows.scores));
        labels_all = vertcat(rows.labels);
        [fpr_s, tpr_s, ~, auc_s] = compute_roc_curve(scores_all, labels_all);
        if isnan(auc_s), continue; end
        roc_auc_subj(si,pi) = auc_s;
        [u_fpr, ~, ic] = unique(fpr_s);
        u_tpr = accumarray(ic, tpr_s, [], @max);
        roc_tpr_subj_grid(si,pi,:) = interp1(u_fpr, u_tpr, ROC_FPR_GRID, 'linear');
    end
end

roc_mean_tpr = nan(3, numel(ROC_FPR_GRID));
roc_sem_tpr  = nan(3, numel(ROC_FPR_GRID));
p_auc_grp    = nan(1, 3);
for pi = 1:3
    curves = reshape(roc_tpr_subj_grid(:,pi,:), n_subj, numel(ROC_FPR_GRID));
    valid_rows = ~all(isnan(curves), 2);
    n_v = sum(valid_rows);
    if n_v > 0
        roc_mean_tpr(pi,:) = mean(curves(valid_rows,:), 1, 'omitnan');
        roc_sem_tpr(pi,:)  = std(curves(valid_rows,:), 0, 1, 'omitnan') / sqrt(max(1,n_v));
    end
    p_auc_grp(pi) = sign_flip_test_local(roc_auc_subj(:,pi) - 0.5, n_perm);
end

% ── Group-level real-session accuracy across paradigms (subject-level) ─────
% Both accuracy conventions are tested, clearly labelled -- see the comment
% above the real_acc/real_acc_dec computation. Pairwise sign-flip (subject =
% statistical unit) mirrors the counterfactual tests further below, but on
% the REAL behavioural outcome -- the claim that matters most for the paper
% is that Hybrid beats unimodal in actual sessions, not only in the
% counterfactual replay.
d_real_hyb_mi  = real_acc(:,3) - real_acc(:,1);   % Hybrid - MI
d_real_hyb_cvs = real_acc(:,3) - real_acc(:,2);   % Hybrid - CVSA
d_real_mi_cvs  = real_acc(:,1) - real_acc(:,2);   % MI - CVSA
p_real_hyb_mi  = sign_flip_test_local(d_real_hyb_mi,  n_perm);
p_real_hyb_cvs = sign_flip_test_local(d_real_hyb_cvs, n_perm);
p_real_mi_cvs  = sign_flip_test_local(d_real_mi_cvs,  n_perm);
p_friedman_real = friedman_test_local(real_acc, n_perm);
% Cohen's d (subject is the unit) -- the effect-size companion to the
% p-value above: how big is the advantage, not just whether it's non-zero.
d_real_hyb_mi_cohen  = cohen_d_local(d_real_hyb_mi);
d_real_hyb_cvs_cohen = cohen_d_local(d_real_hyb_cvs);
d_real_mi_cvs_cohen  = cohen_d_local(d_real_mi_cvs);

% Decided-trials-only (TIMEOUT excluded) twin of the block above.
d_real_hyb_mi_dec  = real_acc_dec(:,3) - real_acc_dec(:,1);
d_real_hyb_cvs_dec = real_acc_dec(:,3) - real_acc_dec(:,2);
d_real_mi_cvs_dec  = real_acc_dec(:,1) - real_acc_dec(:,2);
p_real_hyb_mi_dec  = sign_flip_test_local(d_real_hyb_mi_dec,  n_perm);
p_real_hyb_cvs_dec = sign_flip_test_local(d_real_hyb_cvs_dec, n_perm);
d_real_hyb_mi_dec_cohen  = cohen_d_local(d_real_hyb_mi_dec);
d_real_hyb_cvs_dec_cohen = cohen_d_local(d_real_hyb_cvs_dec);
d_real_mi_cvs_dec_cohen  = cohen_d_local(d_real_mi_cvs_dec);
p_real_mi_cvs_dec  = sign_flip_test_local(d_real_mi_cvs_dec,  n_perm);
p_friedman_real_dec = friedman_test_local(real_acc_dec, n_perm);

% ── Best / worst subject ranking, PER PARADIGM (real) and PER STREAM (sim) ─
% Ranked on DECIDED-TRIAL accuracy (TIMEOUT excluded) -- the "quality of
% decision" metric, not diluted by how often the system failed to decide at
% all. Ties are kept together (best/worst can be more than one subject).
% See best_worst_ties() at the bottom of this file.
idx_best_real  = cell(1,3);
idx_worst_real = cell(1,3);
for pi = 1:3
    [idx_best_real{pi}, idx_worst_real{pi}] = best_worst_ties(real_acc_dec(:,pi));
end
idx_best_cf  = cell(1,3);
idx_worst_cf = cell(1,3);
for s = 1:3
    [idx_best_cf{s}, idx_worst_cf{s}] = best_worst_ties(cf_acc_dec(:,s));
end

% ── Console: per-subject real-session accuracy, BOTH conventions ───────────
fprintf('\n══════════════════ Per-subject real-session accuracy (TIMEOUT counted as fail) ══════════════════\n');
fprintf('  %-10s  %8s  %8s  %8s   %8s  %8s  %8s\n', 'subject', 'MI%', 'CVSA%', 'Hybrid%', 'n_MI', 'n_CVSA', 'n_Hyb');
for si = 1:n_subj
    fprintf('  %-10s  %8s  %8s  %8s   %8d  %8d  %8d\n', subjects{si}, ...
            fmt_pct(real_acc(si,1)), fmt_pct(real_acc(si,2)), fmt_pct(real_acc(si,3)), ...
            real_n(si,1), real_n(si,2), real_n(si,3));
end
fprintf('  %-10s  %8s  %8s  %8s   %8d  %8d  %8d\n', 'GRAND', fmt_pct(grand_acc(1)), fmt_pct(grand_acc(2)), fmt_pct(grand_acc(3)), ...
        grand_n(1), grand_n(2), grand_n(3));

fprintf('\n══════════════════ Per-subject real-session accuracy (DECIDED TRIALS ONLY, TIMEOUT excluded) + false-positive rate ══════════════════\n');
fprintf('  %-10s  %8s  %8s  %8s   %8s  %8s  %8s\n', 'subject', 'MI%', 'CVSA%', 'Hybrid%', 'FP-MI%', 'FP-CVSA%', 'FP-Hyb%');
for si = 1:n_subj
    fprintf('  %-10s  %8s  %8s  %8s   %8s  %8s  %8s\n', subjects{si}, ...
            fmt_pct(real_acc_dec(si,1)), fmt_pct(real_acc_dec(si,2)), fmt_pct(real_acc_dec(si,3)), ...
            fmt_pct(real_fp_rate(si,1)), fmt_pct(real_fp_rate(si,2)), fmt_pct(real_fp_rate(si,3)));
end
fprintf('  %-10s  %8s  %8s  %8s   %8s  %8s  %8s\n', 'GRAND', ...
        fmt_pct(grand_acc_dec(1)), fmt_pct(grand_acc_dec(2)), fmt_pct(grand_acc_dec(3)), ...
        fmt_pct(grand_fp_rate(1)), fmt_pct(grand_fp_rate(2)), fmt_pct(grand_fp_rate(3)));
fprintf('  (false positive = MISS: the WRONG class''s threshold was reached with confidence, distinct from TIMEOUT)\n');

fprintf('\n══════════════════ Group-level real-session accuracy across paradigms (n=%d subjects) ══════════════════\n', n_subj);
fprintf('  -- TIMEOUT counted as fail --\n');
fprintf('  Friedman omnibus (MI vs CVSA vs Hybrid) : p=%.4f  %s\n', p_friedman_real, stars(p_friedman_real));
fprintf('  Hybrid - MI    : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', 100*mean(d_real_hyb_mi,'omitnan'),  p_real_hyb_mi,  stars(p_real_hyb_mi),  d_real_hyb_mi_cohen,  cohen_d_label(d_real_hyb_mi_cohen));
fprintf('  Hybrid - CVSA  : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', 100*mean(d_real_hyb_cvs,'omitnan'), p_real_hyb_cvs, stars(p_real_hyb_cvs), d_real_hyb_cvs_cohen, cohen_d_label(d_real_hyb_cvs_cohen));
fprintf('  MI - CVSA      : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', 100*mean(d_real_mi_cvs,'omitnan'),  p_real_mi_cvs,  stars(p_real_mi_cvs),  d_real_mi_cvs_cohen,  cohen_d_label(d_real_mi_cvs_cohen));
fprintf('  -- DECIDED TRIALS ONLY (TIMEOUT excluded) --\n');
fprintf('  Friedman omnibus (MI vs CVSA vs Hybrid) : p=%.4f  %s\n', p_friedman_real_dec, stars(p_friedman_real_dec));
fprintf('  Hybrid - MI    : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', 100*mean(d_real_hyb_mi_dec,'omitnan'),  p_real_hyb_mi_dec,  stars(p_real_hyb_mi_dec),  d_real_hyb_mi_dec_cohen,  cohen_d_label(d_real_hyb_mi_dec_cohen));
fprintf('  Hybrid - CVSA  : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', 100*mean(d_real_hyb_cvs_dec,'omitnan'), p_real_hyb_cvs_dec, stars(p_real_hyb_cvs_dec), d_real_hyb_cvs_dec_cohen, cohen_d_label(d_real_hyb_cvs_dec_cohen));
fprintf('  MI - CVSA      : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', 100*mean(d_real_mi_cvs_dec,'omitnan'),  p_real_mi_cvs_dec,  stars(p_real_mi_cvs_dec),  d_real_mi_cvs_dec_cohen,  cohen_d_label(d_real_mi_cvs_dec_cohen));
fprintf('  (pairwise: two-sided sign-flip permutation, %d perms, subject is the statistical unit; d = paired Cohen''s d, |d|>0.2 small/>0.5 medium/>0.8 large; Friedman is the omnibus check run before the pairwise tests)\n', n_perm);
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

fprintf('\n══════════════════ Subject ranking by real-session accuracy, per paradigm (decided trials only) ══════════════════\n');
for pi = 1:3
    print_ranking_table(sprintf('Real -- %s', par_labels{pi}), subjects, real_acc_dec(:,pi), idx_best_real{pi}, idx_worst_real{pi});
end
fprintf('═══════════════════════════════════════════════════════════════════════════════\n');

fprintf('\n══════════════════ Per-subject real-session TIMEOUT rate (n_timeout) ══════════════════\n');
fprintf('  %-10s  %14s  %14s  %14s\n', 'subject', 'MI%', 'CVSA%', 'Hybrid%');
for si = 1:n_subj
    fprintf('  %-10s  %8s (%3d)  %8s (%3d)  %8s (%3d)\n', subjects{si}, ...
            fmt_pct(real_to_rate(si,1)), real_to_n(si,1), ...
            fmt_pct(real_to_rate(si,2)), real_to_n(si,2), ...
            fmt_pct(real_to_rate(si,3)), real_to_n(si,3));
end
fprintf('  %-10s  %8s (%3d)  %8s (%3d)  %8s (%3d)\n', 'GRAND', ...
        fmt_pct(grand_to_rate(1)), grand_to_n(1), ...
        fmt_pct(grand_to_rate(2)), grand_to_n(2), ...
        fmt_pct(grand_to_rate(3)), grand_to_n(3));

fprintf('\n══════════════════ Per-subject real-session times (s) ══════════════════\n');
fprintf('  %-10s  %22s  %22s  %22s\n', 'subject', 'TTH (MI/CVSA/Hyb)', 'T-miss (MI/CVSA/Hyb)', 'T-timeout (MI/CVSA/Hyb)');
for si = 1:n_subj
    fprintf('  %-10s  %6.2f/%6.2f/%6.2f  %6.2f/%6.2f/%6.2f  %6.2f/%6.2f/%6.2f\n', subjects{si}, ...
            real_tth(si,1), real_tth(si,2), real_tth(si,3), ...
            real_t_miss(si,1), real_t_miss(si,2), real_t_miss(si,3), ...
            real_to_time(si,1), real_to_time(si,2), real_to_time(si,3));
end
fprintf('  %-10s  %6.2f/%6.2f/%6.2f  %6.2f/%6.2f/%6.2f  %6.2f/%6.2f/%6.2f\n', 'GRAND', ...
        grand_tth(1), grand_tth(2), grand_tth(3), ...
        grand_t_miss(1), grand_t_miss(2), grand_t_miss(3), ...
        grand_to_time(1), grand_to_time(2), grand_to_time(3));
fprintf('  (GRAND = trial-pooled across all subjects, NOT the mean of the per-subject means)\n');

fprintf('\n══════════════════ Per-subject counterfactual accuracy (TIMEOUT counted as fail) ══════════════════\n');
fprintf('  %-10s  %8s  %8s  %8s   %8s  %8s\n', 'subject', 'Hybrid%', 'MI-only%', 'CVSA-only%', 'n_trials', 'Hyb-MI');
for si = 1:n_subj
    if isnan(cf_acc(si,1)), continue; end
    fprintf('  %-10s  %8s  %8s  %8s   %8d  %+7.1f%%\n', subjects{si}, ...
            fmt_pct(cf_acc(si,1)), fmt_pct(cf_acc(si,2)), fmt_pct(cf_acc(si,3)), ...
            cf_n(si), 100*(cf_acc(si,1)-cf_acc(si,2)));
end
fprintf('  %-10s  %8s  %8s  %8s   %8d\n', 'GRAND', fmt_pct(grand_cf_acc(1)), fmt_pct(grand_cf_acc(2)), fmt_pct(grand_cf_acc(3)), sum(cf_n));

fprintf('\n══════════════════ Per-subject counterfactual accuracy (DECIDED TRIALS ONLY) + false-positive rate ══════════════════\n');
fprintf('  %-10s  %8s  %8s  %8s   %8s  %8s  %8s\n', 'subject', 'Hybrid%', 'MI-only%', 'CVSA-only%', 'FP-Hyb%', 'FP-MI%', 'FP-CVSA%');
for si = 1:n_subj
    if isnan(cf_acc_dec(si,1)) && isnan(cf_acc_dec(si,2)) && isnan(cf_acc_dec(si,3)), continue; end
    fprintf('  %-10s  %8s  %8s  %8s   %8s  %8s  %8s\n', subjects{si}, ...
            fmt_pct(cf_acc_dec(si,1)), fmt_pct(cf_acc_dec(si,2)), fmt_pct(cf_acc_dec(si,3)), ...
            fmt_pct(cf_fp_rate(si,1)), fmt_pct(cf_fp_rate(si,2)), fmt_pct(cf_fp_rate(si,3)));
end
fprintf('  %-10s  %8s  %8s  %8s   %8s  %8s  %8s\n', 'GRAND', ...
        fmt_pct(grand_cf_acc_dec(1)), fmt_pct(grand_cf_acc_dec(2)), fmt_pct(grand_cf_acc_dec(3)), ...
        fmt_pct(grand_cf_fp_rate(1)), fmt_pct(grand_cf_fp_rate(2)), fmt_pct(grand_cf_fp_rate(3)));

fprintf('\n══════════════════ Subject ranking by counterfactual accuracy, per stream (decided trials only) ══════════════════\n');
for s = 1:3
    print_ranking_table(sprintf('Simulated -- %s', cf_labels{s}), subjects, cf_acc_dec(:,s), idx_best_cf{s}, idx_worst_cf{s});
end
fprintf('═══════════════════════════════════════════════════════════════════════════════\n');

% ── Group-level paired statistics across SUBJECTS ───────────────────────────
d_mi_subj  = cf_acc(:,1) - cf_acc(:,2);
d_cvs_subj = cf_acc(:,1) - cf_acc(:,3);
p_mi_grp  = sign_flip_test_local(d_mi_subj,  n_perm);
p_cvs_grp = sign_flip_test_local(d_cvs_subj, n_perm);
d_mi_grp_cohen  = cohen_d_local(d_mi_subj);
d_cvs_grp_cohen = cohen_d_local(d_cvs_subj);
d_mi_subj_dec  = cf_acc_dec(:,1) - cf_acc_dec(:,2);
d_cvs_subj_dec = cf_acc_dec(:,1) - cf_acc_dec(:,3);
p_mi_grp_dec  = sign_flip_test_local(d_mi_subj_dec,  n_perm);
p_cvs_grp_dec = sign_flip_test_local(d_cvs_subj_dec, n_perm);
d_mi_grp_dec_cohen  = cohen_d_local(d_mi_subj_dec);
d_cvs_grp_dec_cohen = cohen_d_local(d_cvs_subj_dec);
fprintf('\n══════════════════ Group-level counterfactual advantage (n=%d subjects) ══════════════════\n', sum(~isnan(d_mi_subj)));
fprintf('  -- TIMEOUT counted as fail --\n');
fprintf('  Hybrid - MI-only   : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
        100*mean(d_mi_subj,'omitnan'),  p_mi_grp,  stars(p_mi_grp),  d_mi_grp_cohen,  cohen_d_label(d_mi_grp_cohen));
fprintf('  Hybrid - CVSA-only : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
        100*mean(d_cvs_subj,'omitnan'), p_cvs_grp, stars(p_cvs_grp), d_cvs_grp_cohen, cohen_d_label(d_cvs_grp_cohen));
fprintf('  -- DECIDED TRIALS ONLY (TIMEOUT excluded) --\n');
fprintf('  Hybrid - MI-only   : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
        100*mean(d_mi_subj_dec,'omitnan'),  p_mi_grp_dec,  stars(p_mi_grp_dec),  d_mi_grp_dec_cohen,  cohen_d_label(d_mi_grp_dec_cohen));
fprintf('  Hybrid - CVSA-only : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
        100*mean(d_cvs_subj_dec,'omitnan'), p_cvs_grp_dec, stars(p_cvs_grp_dec), d_cvs_grp_dec_cohen, cohen_d_label(d_cvs_grp_dec_cohen));
fprintf('  (two-sided sign-flip permutation, %d perms; subject is the statistical unit; d = paired Cohen''s d, |d|>0.2 small/>0.5 medium/>0.8 large)\n', n_perm);
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

% ── Bonus: cross-subject correlates of real Hybrid performance ─────────────
% Does each subject's simulated counterfactual advantage predict their
% actual real-session advantage? (an individual-level cross-check of
% main_validate_counterfactual's pooled real-vs-simulated comparison.) Does
% the fusion-mechanism strength / neurophysiological (ERD-CSP) grounding
% predict how well a subject actually performs in real Hybrid sessions? And,
% supporting the "CVSA helps even when imperfect" thesis: does a subject's
% Hybrid-vs-MI-only advantage hold up even when their CVSA-only accuracy is
% weak (decided trials only -- CVSA "quality" in isolation)?
[r_cf_vs_real,  p_cf_vs_real]  = pearson_perm_test_local(d_mi_subj,        d_real_hyb_mi,   n_perm);
[r_fus_vs_real, p_fus_vs_real] = pearson_perm_test_local(fus_adv_subj,     real_acc(:,3),   n_perm);
[r_erd_vs_real, p_erd_vs_real] = pearson_perm_test_local(erd_r_subj(:,3),  real_acc(:,3),   n_perm);
[r_cvsa_quality, p_cvsa_quality] = pearson_perm_test_local(cf_acc_dec(:,3), d_mi_subj_dec, n_perm);
fprintf('\n══════════════════ Cross-subject correlates of real Hybrid performance ══════════════════\n');
fprintf('  Counterfactual adv. (Hyb-MI) vs real adv. (Hyb-MI) : r=%+.2f  p=%.4f  %s\n', r_cf_vs_real,  p_cf_vs_real,  stars(p_cf_vs_real));
fprintf('  Fusion advantage (P_fused-P_MI) vs real Hybrid acc : r=%+.2f  p=%.4f  %s\n', r_fus_vs_real, p_fus_vs_real, stars(p_fus_vs_real));
fprintf('  ERD/ERS-CSP grounding (r) vs real Hybrid acc       : r=%+.2f  p=%.4f  %s\n', r_erd_vs_real, p_erd_vs_real, stars(p_erd_vs_real));
fprintf('  CVSA-only quality vs Hybrid-MI counterfactual adv. : r=%+.2f  p=%.4f  %s  (does the advantage survive weak CVSA?)\n', ...
        r_cvsa_quality, p_cvsa_quality, stars(p_cvsa_quality));
[~, idx_weakest_cvsa] = min(cf_acc_dec(:,3));
if ~isnan(cf_acc_dec(idx_weakest_cvsa,3))
    fprintf('  Weakest CVSA-only subject: %s (CVSA-only=%s decided-acc)  ->  Hybrid-MI advantage=%+.1f%%\n', ...
            subjects{idx_weakest_cvsa}, fmt_pct(cf_acc_dec(idx_weakest_cvsa,3)), 100*d_mi_subj_dec(idx_weakest_cvsa));
end
fprintf('  (Pearson r across subjects, two-sided permutation p, %d perms)\n', n_perm);
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

% ── Group-level CVSA-fusion mechanism across SUBJECTS ───────────────────────
p_fusadv_grp = sign_flip_test_local(fus_adv_subj, n_perm);
d_resc_cost  = rescue_subj - cost_subj;
p_resc_grp   = sign_flip_test_local(d_resc_cost, n_perm);
d_fusadv_cohen = cohen_d_local(fus_adv_subj);
d_resc_cohen   = cohen_d_local(d_resc_cost);
fprintf('\n══════════════════ Group-level CVSA-fusion mechanism (n=%d subjects) ══════════════════\n', sum(~isnan(fus_adv_subj)));
fprintf('  %-10s  %10s  %10s  %10s  %12s\n', 'subject', 'fus_adv', 'rescue_d', 'cost_d', 'net_resc(fr)');
for si = 1:n_subj
    if isnan(fus_adv_subj(si)), continue; end
    fprintf('  %-10s  %+9.3f  %+9.3f  %+9.3f  %+11d\n', subjects{si}, ...
            fus_adv_subj(si), rescue_subj(si), cost_subj(si), net_rescue_subj(si));
end
fprintf('  mean fusion advantage (P_fused-P_MI) : %+.3f  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(fus_adv_subj, 'omitnan'), p_fusadv_grp, stars(p_fusadv_grp), d_fusadv_cohen, cohen_d_label(d_fusadv_cohen));
fprintf('  mean (rescue - cost) delta           : %+.3f  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(d_resc_cost, 'omitnan'), p_resc_grp, stars(p_resc_grp), d_resc_cohen, cohen_d_label(d_resc_cohen));
fprintf('  (two-sided sign-flip permutation, %d perms; subject is the statistical unit; d = Cohen''s d, |d|>0.2 small/>0.5 medium/>0.8 large)\n', n_perm);
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

% ── Meta-analytic (Stouffer) combination of per-subject within-session
%    sign-flip tests -- see the comment above fus_adv_meta_p_subj for why
%    this exists alongside the sign-flip test on fus_adv_subj just above. ──
if any(~isnan(fus_adv_meta_p_subj))
    fprintf('\n═══════════ Group-level CVSA-fusion meta-analysis (Stouffer, n=%d subjects) ═══════════\n', sum(~isnan(fus_adv_meta_p_subj)));
    fprintf('  %-10s  %14s  %10s\n', 'subject', 'within-subj p', 'n_trials');
    for si = 1:n_subj
        if isnan(fus_adv_meta_p_subj(si)), continue; end
        fprintf('  %-10s  %14.4f  %10d\n', subjects{si}, fus_adv_meta_p_subj(si), fus_adv_meta_n_subj(si));
    end
    fprintf('  combined one-sided p (does CVSA-fusion help, H1: mean advantage > 0) : p=%.4f  %s\n', ...
            p_fusadv_meta_grp, stars(p_fusadv_meta_grp));
    fprintf('  (each subject''s own within-session sign-flip test -- see main_hybrid_advantage_probs.m --\n');
    fprintf('   combined via Stouffer''s Z method, subject weighted equally; this test has no small-n floor,\n');
    fprintf('   unlike the sign-flip test on fus_adv_subj above whose null only has 2^n_subj sign patterns)\n');
    fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');
end

% ── Meta-analytic (Stouffer) combination for the INTEGRATOR (counterfactual
%    buffer) CVSA-help test -- same rationale/method as the P_fused-P_MI
%    meta-analysis just above, but on main_hybrid_advantage_integ.m's
%    buf_adv_mi_win_pval (Hybrid vs MI-only, cvsa_influence-window
%    restricted). Directly answers, per subject and on average: "was the
%    fused signal reliably closer to threshold than MI-only while CVSA was
%    still actively weighted?" ─────────────────────────────────────────────
if any(~isnan(buf_adv_mi_win_meta_p_subj))
    fprintf('\n═══════ Group-level CVSA-help (integrator) meta-analysis (Stouffer, n=%d subjects) ═══════\n', sum(~isnan(buf_adv_mi_win_meta_p_subj)));
    fprintf('  %-10s  %14s  %10s  %10s\n', 'subject', 'within-subj p', 'n_trials', 'mean delta');
    for si = 1:n_subj
        if isnan(buf_adv_mi_win_meta_p_subj(si)), continue; end
        fprintf('  %-10s  %14.4f  %10d  %+10.4f\n', subjects{si}, buf_adv_mi_win_meta_p_subj(si), ...
                buf_adv_mi_win_meta_n_subj(si), buf_adv_mi_win_subj(si));
    end
    fprintf('  combined one-sided p (does CVSA help, H1: Hybrid > MI-only in the cvsa_influence window) : p=%.4f  %s\n', ...
            p_bufadvwin_meta_grp, stars(p_bufadvwin_meta_grp));
    fprintf('  (each subject''s own within-session sign-flip test -- see main_hybrid_advantage_integ.m --\n');
    fprintf('   combined via Stouffer''s Z method, subject weighted equally; this test has no small-n floor)\n');
    fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');
end

% ── Group-level (cross-subject) cluster-based permutation test: does the
%    CVSA-fusion effect on P_fused(target)-P_MI(target) show a genuine,
%    temporally-localised period of influence that GENERALISES across the
%    cohort? Rows = subjects (each subject's own mean curve), the
%    permutation unit -- the statistically meaningful unit once N>1 subject
%    exists, same logic as main_hybrid_advantage_probs.m but one level up.
if ~isempty(CLUSTER_TIME_GRID)
    cluster_res_grp = cluster_permutation_test(cluster_diff_subj, CLUSTER_TIME_GRID, n_perm);
    fprintf('\n══════════════════ Group-level CVSA-fusion cluster test (n=%d subjects) ══════════════════\n', n_cluster_subj);
    if isempty(cluster_res_grp.clusters)
        fprintf('  no candidate cluster found (|t| never reached the cluster-forming threshold)\n');
    else
        for k = 1:numel(cluster_res_grp.clusters)
            cl = cluster_res_grp.clusters(k);
            sig_str = ''; if cl.p < 0.05, sig_str = '  <-- SIGNIFICANT'; end
            fprintf('  cluster %d: [%.2f - %.2f]s  mass=%+.1f  p=%.4f%s\n', k, cl.start_t, cl.end_t, cl.mass, cl.p, sig_str);
            fprintf('    within this cluster only: mean delta=%+.3f  95%% CI=[%+.3f, %+.3f]  d=%+.2f\n', ...
                    cl.mean_delta, cl.ci_lo, cl.ci_hi, cl.cohen_d);
        end
    end
    fprintf('  (cluster-forming threshold |t|>=2.0; cluster p-value = permutation on max |cluster mass|, per-subject sign-flip, %d perms; mean delta/CI/d = magnitude localised to that cluster''s own time range, subject is the unit)\n', n_perm);
    fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');
else
    cluster_res_grp = [];
    fprintf('\n[main_group_analysis] No subject has cluster-test data yet (regenerate hybrid_advantage_probs_summary.mat with the updated main_hybrid_advantage_probs.m) -- skipping group cluster test.\n');
end

% ── Group-level ERD/ERS discrimination + CSP grounding across SUBJECTS ─────
p_r_grp = nan(1,3);
fprintf('\n══════════════════ Group-level ERD/ERS discrimination & CSP grounding (n=%d subjects) ══════════════════\n', n_subj);
fprintf('  %-10s', 'subject');
for pi = 1:3, fprintf('  %9s  %9s', [par_labels{pi} '-discr'], [par_labels{pi} '-r']); end
fprintf('\n');
for si = 1:n_subj
    fprintf('  %-10s', subjects{si});
    for pi = 1:3
        fprintf('  %9.2f  %+9.2f', erd_discrim_subj(si,pi), erd_r_subj(si,pi));
    end
    fprintf('\n');
end
for pi = 1:3
    p_r_grp(pi) = sign_flip_test_local(erd_r_subj(:,pi), n_perm);
    fprintf('  %-8s  mean discrimination=%.2f%%  mean CSP-weight r=%+.2f  p=%.4f  %s\n', ...
            par_labels{pi}, mean(erd_discrim_subj(:,pi),'omitnan'), mean(erd_r_subj(:,pi),'omitnan'), ...
            p_r_grp(pi), stars(p_r_grp(pi)));
end
fprintf('  (p tests whether the CSP-weight correlation r is consistently non-zero across subjects)\n');
fprintf('═══════════════════════════════════════════════════════════════════════════════════════════════════════\n');

% ── Group-level classifier ROC/AUC across SUBJECTS (pre-integrator) ────────
fprintf('\n══════════════════ Group-level classifier ROC/AUC, pre-integrator (n=%d subjects) ══════════════════\n', n_subj);
fprintf('  %-10s  %10s  %10s  %10s\n', 'subject', 'MI-AUC', 'CVSA-AUC', 'Hybrid-AUC');
for si = 1:n_subj
    if all(isnan(roc_auc_subj(si,:))), continue; end
    fprintf('  %-10s  %10.3f  %10.3f  %10.3f\n', subjects{si}, ...
            roc_auc_subj(si,1), roc_auc_subj(si,2), roc_auc_subj(si,3));
end
for pi = 1:3
    fprintf('  %-8s  mean AUC=%.3f  p(AUC>0.5)=%.4f  %s\n', ...
            par_labels{pi}, mean(roc_auc_subj(:,pi),'omitnan'), p_auc_grp(pi), stars(p_auc_grp(pi)));
end
fprintf('  (AUC is per-subject, re-pooled across that subject''s sessions; p = two-sided sign-flip permutation on AUC-0.5, subject is the statistical unit)\n');
fprintf('═══════════════════════════════════════════════════════════════════════════════════════════════════════\n');

% ── Colors (consistent with main_session_overview / main_hybrid_advantage_integ) ─
COL_mi     = [0.85 0.30 0.10];
COL_cvsa   = [0.10 0.60 0.30];
COL_hybrid = [0.18 0.45 0.75];
COL_mi_only   = [0.85 0.40 0.10];
COL_cvsa_only = [0.49 0.18 0.56];

out_dir = fullfile(root_dir, 'group_analysis');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% ── Figure 1: per-subject real accuracy + timeout rate per paradigm, each
%    panel also showing a trial-pooled GRAND AVG column at the far right ──
bw = 0.25;
cols3 = {COL_mi, COL_cvsa, COL_hybrid};

fig1 = figure('Name', 'Group Analysis — Real Session Accuracy & Timeout Rate', 'Color', 'w', ...
              'NumberTitle', 'off', 'Visible', fig_vis);
set(fig1, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

ax1 = subplot(1, 2, 1); hold(ax1, 'on');
plot_subject_grand_bars(ax1, real_acc, grand_acc, subjects, par_labels, cols3, bw);
yline(ax1, 50, 'k:', 'HandleVisibility', 'off');
ylabel(ax1, 'HIT rate (%)');
legend(ax1, 'Location', 'south', 'FontSize', 9);
title(ax1, sprintf('Real-session accuracy, TIMEOUT counted as fail (n=%d subjects)', n_subj), 'FontWeight', 'bold');
grid(ax1, 'on');
% (best/worst per-paradigm markers are on Fig 1c, which uses the decided-
% trial accuracy convention -- the more meaningful ranking metric.)

ax1b = subplot(1, 2, 2); hold(ax1b, 'on');
plot_subject_grand_bars(ax1b, real_to_rate, grand_to_rate, subjects, par_labels, cols3, bw);
ylabel(ax1b, 'TIMEOUT rate (%)');
legend(ax1b, 'Location', 'best', 'FontSize', 9);
title(ax1b, 'Real-session timeout rate', 'FontWeight', 'bold');
grid(ax1b, 'on');

sgtitle(fig1, sprintf(['Group Analysis — Real Session Accuracy & Timeout Rate (n=%d subjects)\n' ...
        '"GRAND" = trial-pooled across all subjects (not the mean of per-subject means)'], n_subj), ...
        'Interpreter', 'none');

saveas(fig1, fullfile(out_dir, 'group_real_session_overview.svg'), 'svg');
if ~SHOW_FIGURES, close(fig1); end

% ── Figure 1b: per-subject real-session TTH / T-miss / T-timeout, one panel
%    per paradigm, with the trial-pooled GRAND AVG as the last group ──────
fig1b = figure('Name', 'Group Analysis — Real Session Times', 'Color', 'w', ...
               'NumberTitle', 'off', 'Visible', fig_vis);
set(fig1b, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

COL_tth   = [0.18 0.45 0.75];
COL_tmiss = [0.80 0.20 0.20];
COL_tto   = [0.55 0.55 0.10];
time_cols = {COL_tth, COL_tmiss, COL_tto};
time_labels = {'TTH', 'T-miss', 'T-timeout'};

for pi = 1:3
    axp = subplot(1, 3, pi); hold(axp, 'on');
    metric_mat  = [real_tth(:,pi), real_t_miss(:,pi), real_to_time(:,pi)];
    grand_metric = [grand_tth(pi), grand_t_miss(pi), grand_to_time(pi)];
    plot_subject_grand_bars(axp, metric_mat, grand_metric, subjects, time_labels, time_cols, bw, false);
    ylabel(axp, 'time (s)');
    if pi == 1, legend(axp, 'Location', 'best', 'FontSize', 9); end
    title(axp, par_labels{pi}, 'FontWeight', 'bold');
    grid(axp, 'on');
end
sgtitle(fig1b, sprintf(['Group Analysis — Real Session Times (n=%d subjects)\n' ...
        '"GRAND" = trial-pooled across all subjects'], n_subj), 'Interpreter', 'none');

saveas(fig1b, fullfile(out_dir, 'group_real_times.svg'), 'svg');
if ~SHOW_FIGURES, close(fig1b); end

% ── Figure 1c: real-session performance breakdown -- decided-trial accuracy
%    (TIMEOUT excluded) on top, false-positive rate + timeout rate below,
%    per paradigm, per subject + GRAND. Gold star / red triangle mark the
%    best/worst subject FOR THAT PARADIGM (ties shown together). ──────────
fig1c = figure('Name', 'Group Analysis — Real Session Performance Breakdown', 'Color', 'w', ...
               'NumberTitle', 'off', 'Visible', fig_vis);
set(fig1c, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

axA1 = subplot(2, 2, [1 2]); hold(axA1, 'on');
plot_subject_grand_bars(axA1, real_acc_dec, grand_acc_dec, subjects, par_labels, cols3, bw);
yline(axA1, 50, 'k:', 'HandleVisibility', 'off');
ylabel(axA1, 'HIT rate, decided trials only (%)');
legend(axA1, 'Location', 'south', 'FontSize', 9);
title(axA1, sprintf('Accuracy on decided trials, TIMEOUT excluded (n=%d subjects)', n_subj), 'FontWeight', 'bold');
grid(axA1, 'on');
for pi = 1:3
    offset_pi = (pi - 2) * bw;
    for si = idx_best_real{pi}
        y = min(100*real_acc_dec(si,pi) + 6, 108);
        plot(axA1, si + offset_pi, y, 'p', 'MarkerSize', 12, ...
             'MarkerFaceColor', [0.85 0.65 0.10], 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    for si = idx_worst_real{pi}
        y = min(100*real_acc_dec(si,pi) + 6, 108);
        plot(axA1, si + offset_pi, y, 'v', 'MarkerSize', 9, ...
             'MarkerFaceColor', [0.75 0.15 0.15], 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
end

axA2 = subplot(2, 2, 3); hold(axA2, 'on');
plot_subject_grand_bars(axA2, real_fp_rate, grand_fp_rate, subjects, par_labels, cols3, bw);
ylabel(axA2, 'False positive rate (%)');
title(axA2, 'MISS: wrong-class threshold reached', 'FontWeight', 'bold', 'FontSize', 9);
grid(axA2, 'on');

axA3 = subplot(2, 2, 4); hold(axA3, 'on');
plot_subject_grand_bars(axA3, real_to_rate, grand_to_rate, subjects, par_labels, cols3, bw);
ylabel(axA3, 'Timeout rate (%)');
title(axA3, 'TIMEOUT: no threshold reached', 'FontWeight', 'bold', 'FontSize', 9);
grid(axA3, 'on');

sgtitle(fig1c, sprintf(['Group Analysis — Real Session Performance Breakdown (n=%d subjects)\n' ...
        'gold star = best subject per paradigm, red triangle = worst (ties shown together), ranked by decided-trial accuracy  |  "GRAND" = trial-pooled'], n_subj), ...
        'Interpreter', 'none');

saveas(fig1c, fullfile(out_dir, 'group_real_accuracy_breakdown.svg'), 'svg');
if ~SHOW_FIGURES, close(fig1c); end

% ── Figure 2: counterfactual accuracy + group-level advantage ─────────────
fig2 = figure('Name', 'Group Analysis — Counterfactual Advantage', 'Color', 'w', ...
              'NumberTitle', 'off', 'Visible', fig_vis);
set(fig2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

ax2 = subplot(1, 2, 1); hold(ax2, 'on');
cols_cf = {COL_hybrid, COL_mi_only, COL_cvsa_only};
for pi = 1:3
    x = (1:n_subj) + (pi-2)*bw;
    v = 100*cf_acc(:,pi);
    bar(ax2, x, v, bw*0.9, 'FaceColor', cols_cf{pi}, 'EdgeColor', 'k', 'LineWidth', 1, ...
        'DisplayName', cf_labels{pi});
    m = mean(v, 'omitnan');
    if ~isnan(m)
        plot(ax2, [0.5, n_subj+0.5], [m m], '--', 'Color', cols_cf{pi}*0.6, ...
             'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end
yline(ax2, 50, 'k:', 'HandleVisibility', 'off');
set(ax2, 'XTick', 1:n_subj, 'XTickLabel', subjects, 'YLim', [0, 110], 'XLim', [0.4, n_subj+0.6]);
ylabel(ax2, 'simulated HIT rate (%)');
legend(ax2, 'Location', 'south', 'FontSize', 9);
title(ax2, 'Counterfactual accuracy per subject', 'FontWeight', 'bold');
grid(ax2, 'on');

ax3 = subplot(1, 2, 2); hold(ax3, 'on');
x1 = (1:n_subj) - 0.12;
x2 = (1:n_subj) + 0.12;
scatter(ax3, x1, 100*d_mi_subj,  70, COL_mi_only,   'filled', 'MarkerEdgeColor', 'k');
scatter(ax3, x2, 100*d_cvs_subj, 70, COL_cvsa_only, 'filled', 'MarkerEdgeColor', 'k');
yline(ax3, 0, 'k-', 'HandleVisibility', 'off');
m_mi  = mean(d_mi_subj,  'omitnan'); se_mi  = std(d_mi_subj,  'omitnan') / sqrt(max(1,sum(~isnan(d_mi_subj))));
m_cvs = mean(d_cvs_subj, 'omitnan'); se_cvs = std(d_cvs_subj, 'omitnan') / sqrt(max(1,sum(~isnan(d_cvs_subj))));
errorbar(ax3, n_subj+0.7, 100*m_mi,  100*se_mi,  'o', 'Color', COL_mi_only,   'MarkerFaceColor', COL_mi_only,   'LineWidth', 1.5, 'CapSize', 6);
errorbar(ax3, n_subj+1.0, 100*m_cvs, 100*se_cvs, 'o', 'Color', COL_cvsa_only, 'MarkerFaceColor', COL_cvsa_only, 'LineWidth', 1.5, 'CapSize', 6);
text(ax3, n_subj+0.7, 100*m_mi  + 100*se_mi  + 3, sprintf('%s', stars(p_mi_grp)),  'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
text(ax3, n_subj+1.0, 100*m_cvs + 100*se_cvs + 3, sprintf('%s', stars(p_cvs_grp)), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
set(ax3, 'XTick', [1:n_subj, n_subj+0.85], 'XTickLabel', [subjects, {'group'}], 'XLim', [0.4, n_subj+1.4]);
ylabel(ax3, '\Delta accuracy (Hybrid - unimodal), %');
legend(ax3, {'Hybrid - MI-only', 'Hybrid - CVSA-only'}, 'Location', 'best', 'FontSize', 9);
title(ax3, sprintf('Per-subject advantage  |  group: Hyb-MI=%+.1f%% (%s, d=%+.2f), Hyb-CVSA=%+.1f%% (%s, d=%+.2f)', ...
      100*m_mi, stars(p_mi_grp), d_mi_grp_cohen, 100*m_cvs, stars(p_cvs_grp), d_cvs_grp_cohen), 'FontWeight', 'bold', 'FontSize', 9);
grid(ax3, 'on');

sgtitle(fig2, sprintf(['Group Analysis — Counterfactual Advantage (n=%d subjects)\n' ...
        'SIMULATED, not real: SAME integrator/buffer/thresholds driven counterfactually by Hybrid-fused, MI-only, or CVSA-only signal\n' ...
        'on the SAME trials -- isolates the fusion algorithm''s own effect from real-session differences in trial count/timing/thresholds'], ...
        n_subj), 'Interpreter', 'none');

saveas(fig2, fullfile(out_dir, 'group_counterfactual_advantage.svg'), 'svg');
if ~SHOW_FIGURES, close(fig2); end

% ── Figure 2b: simulated (counterfactual) performance breakdown -- same
%    layout as Fig 1c, but for the 3 counterfactual streams instead of the
%    3 real paradigms. Best/worst marked PER STREAM (ties shown together). ─
fig2b = figure('Name', 'Group Analysis — Simulated Performance Breakdown', 'Color', 'w', ...
               'NumberTitle', 'off', 'Visible', fig_vis);
set(fig2b, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

axB1 = subplot(2, 2, [1 2]); hold(axB1, 'on');
plot_subject_grand_bars(axB1, cf_acc_dec, grand_cf_acc_dec, subjects, cf_labels, cols_cf, bw);
yline(axB1, 50, 'k:', 'HandleVisibility', 'off');
ylabel(axB1, 'simulated HIT rate, decided trials only (%)');
legend(axB1, 'Location', 'south', 'FontSize', 9);
title(axB1, sprintf('Counterfactual accuracy on decided trials, TIMEOUT excluded (n=%d subjects)', n_subj), 'FontWeight', 'bold');
grid(axB1, 'on');
for s = 1:3
    offset_s = (s - 2) * bw;
    for si = idx_best_cf{s}
        y = min(100*cf_acc_dec(si,s) + 6, 108);
        plot(axB1, si + offset_s, y, 'p', 'MarkerSize', 12, ...
             'MarkerFaceColor', [0.85 0.65 0.10], 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    for si = idx_worst_cf{s}
        y = min(100*cf_acc_dec(si,s) + 6, 108);
        plot(axB1, si + offset_s, y, 'v', 'MarkerSize', 9, ...
             'MarkerFaceColor', [0.75 0.15 0.15], 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
end

axB2 = subplot(2, 2, 3); hold(axB2, 'on');
plot_subject_grand_bars(axB2, cf_fp_rate, grand_cf_fp_rate, subjects, cf_labels, cols_cf, bw);
ylabel(axB2, 'False positive rate (%)');
title(axB2, 'MISS: wrong-class threshold reached', 'FontWeight', 'bold', 'FontSize', 9);
grid(axB2, 'on');

axB3 = subplot(2, 2, 4); hold(axB3, 'on');
plot_subject_grand_bars(axB3, cf_to_rate, grand_cf_to_rate, subjects, cf_labels, cols_cf, bw);
ylabel(axB3, 'Timeout rate (%)');
title(axB3, 'TIMEOUT: no threshold reached', 'FontWeight', 'bold', 'FontSize', 9);
grid(axB3, 'on');

sgtitle(fig2b, sprintf(['Group Analysis — Simulated (Counterfactual) Performance Breakdown (n=%d subjects)\n' ...
        'SAME integrator/thresholds, driven counterfactually by each stream on the same Hybrid-session trials\n' ...
        'gold star = best subject per stream, red triangle = worst (ties together)  |  "GRAND" = trial-pooled'], n_subj), ...
        'Interpreter', 'none');

saveas(fig2b, fullfile(out_dir, 'group_simulated_accuracy_breakdown.svg'), 'svg');
if ~SHOW_FIGURES, close(fig2b); end

% ── Figure 3: CVSA-fusion mechanism across subjects ────────────────────────
COL_rescue = [0.20 0.60 0.20];
COL_cost   = [0.80 0.20 0.20];

fig3 = figure('Name', 'Group Analysis — CVSA-Fusion Mechanism', 'Color', 'w', ...
              'NumberTitle', 'off', 'Visible', fig_vis);
set(fig3, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

ax4 = subplot(1, 2, 1); hold(ax4, 'on');
bar(ax4, 1:n_subj, fus_adv_subj, 0.5, 'FaceColor', COL_hybrid, 'EdgeColor', 'k');
yline(ax4, 0, 'k-', 'HandleVisibility', 'off');
m_fa = mean(fus_adv_subj, 'omitnan'); se_fa = std(fus_adv_subj, 'omitnan') / sqrt(max(1,sum(~isnan(fus_adv_subj))));
errorbar(ax4, n_subj+0.8, m_fa, se_fa, 'o', 'Color', COL_hybrid, 'MarkerFaceColor', COL_hybrid, 'LineWidth', 1.5, 'CapSize', 6);
text(ax4, n_subj+0.8, m_fa + se_fa*sign(m_fa+eps) + 0.01, stars(p_fusadv_grp), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
set(ax4, 'XTick', [1:n_subj, n_subj+0.8], 'XTickLabel', [subjects, {'group'}], 'XLim', [0.4, n_subj+1.2]);
ylabel(ax4, 'mean(P_{fused} - P_{MI})  on target class');
title(ax4, sprintf('Per-subject CVSA-fusion advantage  |  group mean=%+.3f (%s, d=%+.2f)', m_fa, stars(p_fusadv_grp), d_fusadv_cohen), ...
      'FontWeight', 'bold', 'FontSize', 9);
grid(ax4, 'on');

ax5 = subplot(1, 2, 2); hold(ax5, 'on');
for si = 1:n_subj
    if isnan(rescue_subj(si)), continue; end
    plot(ax5, [si-0.15, si+0.15], [cost_subj(si), rescue_subj(si)], '-', ...
         'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
end
scatter(ax5, (1:n_subj)-0.15, cost_subj,   60, COL_cost,   'filled', 'MarkerEdgeColor', 'k');
scatter(ax5, (1:n_subj)+0.15, rescue_subj, 60, COL_rescue, 'filled', 'MarkerEdgeColor', 'k');
yline(ax5, 0, 'k-', 'HandleVisibility', 'off');
set(ax5, 'XTick', 1:n_subj, 'XTickLabel', subjects, 'XLim', [0.4, n_subj+0.6]);
ylabel(ax5, 'P_{fused} - P_{MI}  on target class');
legend(ax5, {'cost (MI ok, CVSA wrong)', 'rescue (MI wrong, CVSA ok)'}, 'Location', 'best', 'FontSize', 8);
title(ax5, sprintf('Rescue vs cost per subject  |  group (rescue-cost)=%+.3f (%s, d=%+.2f)', ...
      mean(d_resc_cost,'omitnan'), stars(p_resc_grp), d_resc_cohen), 'FontWeight', 'bold', 'FontSize', 9);
grid(ax5, 'on');

sgtitle(fig3, sprintf('Group Analysis — CVSA-Fusion Mechanism (n=%d subjects)', n_subj), 'Interpreter', 'none');

saveas(fig3, fullfile(out_dir, 'group_fusion_mechanism.svg'), 'svg');
if ~SHOW_FIGURES, close(fig3); end

% ── Figure 3b: cross-subject cluster-based permutation test on the CVSA-
%    fusion effect (P_fused(target)-P_MI(target)) -- does the temporally-
%    localised influence found per-session generalise across the cohort?
%    Skipped if no subject has cluster-test data yet. ─────────────────────
if ~isempty(CLUSTER_TIME_GRID)
    fig3b = figure('Name', 'Group Analysis — CVSA-Fusion Cluster Test', 'Color', 'w', ...
                   'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig3b, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    ax3b1 = subplot(1,2,1); hold(ax3b1, 'on');
    m_d  = mean(cluster_diff_subj, 1, 'omitnan');
    n_e  = sum(~isnan(cluster_diff_subj), 1);
    s_d  = std(cluster_diff_subj, 0, 1, 'omitnan') ./ sqrt(max(n_e,1));
    sig = cluster_res_grp.sig_mask;
    k = 1;
    while k <= numel(sig)
        if sig(k)
            k2 = k;
            while k2 <= numel(sig) && sig(k2), k2 = k2 + 1; end
            fill(ax3b1, [CLUSTER_TIME_GRID(k) CLUSTER_TIME_GRID(k2-1) CLUSTER_TIME_GRID(k2-1) CLUSTER_TIME_GRID(k)], ...
                 [min(m_d-s_d,[],'omitnan') min(m_d-s_d,[],'omitnan') max(m_d+s_d,[],'omitnan') max(m_d+s_d,[],'omitnan')], ...
                 [1.0 0.90 0.50], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            k = k2;
        else
            k = k + 1;
        end
    end
    fill(ax3b1, [CLUSTER_TIME_GRID, fliplr(CLUSTER_TIME_GRID)], [m_d-s_d, fliplr(m_d+s_d)], ...
         COL_hybrid, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax3b1, CLUSTER_TIME_GRID, m_d, '-', 'Color', COL_hybrid, 'LineWidth', 2, 'DisplayName', 'mean(P_{fused}-P_{MI})');
    yline(ax3b1, 0, 'k-', 'HandleVisibility', 'off');
    set(ax3b1, 'XLim', [0, CLUSTER_TIME_GRID(end)]);
    xlabel(ax3b1, 'time from CF onset (s)', 'FontSize', 9);
    ylabel(ax3b1, 'P_{fused}(target) - P_{MI}(target)', 'FontSize', 9);
    title(ax3b1, sprintf('Cross-subject mean \\pm SEM (n=%d subjects)\nshaded = significant cluster (p<0.05)', n_cluster_subj), ...
          'FontSize', 10, 'FontWeight', 'bold');
    legend(ax3b1, 'FontSize', 8, 'Location', 'best');
    grid(ax3b1, 'on');

    ax3b2 = subplot(1,2,2); hold(ax3b2, 'on');
    plot(ax3b2, CLUSTER_TIME_GRID, cluster_res_grp.tstat, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'DisplayName', 't-statistic');
    yline(ax3b2, 2.0, 'r--', 'LineWidth', 1, 'DisplayName', 'cluster-forming |t|=2.0');
    yline(ax3b2, -2.0, 'r--', 'HandleVisibility', 'off');
    yline(ax3b2, 0, 'k-', 'HandleVisibility', 'off');
    t_top = max([cluster_res_grp.tstat, 2.5], [], 'omitnan');
    for k = 1:numel(cluster_res_grp.clusters)
        cl = cluster_res_grp.clusters(k);
        if cl.p < 0.05
            plot(ax3b2, [cl.start_t, cl.end_t], [t_top, t_top]*1.05, '-', ...
                 'Color', [0.9 0.6 0.1], 'LineWidth', 4, 'HandleVisibility', 'off');
            text(ax3b2, mean([cl.start_t, cl.end_t]), t_top*1.15, sprintf('p=%.3f', cl.p), ...
                 'HorizontalAlignment', 'center', 'FontSize', 8);
        end
    end
    set(ax3b2, 'XLim', [0, CLUSTER_TIME_GRID(end)]);
    xlabel(ax3b2, 'time from CF onset (s)', 'FontSize', 9);
    ylabel(ax3b2, 't-statistic  (mean/SEM across subjects)', 'FontSize', 9);
    title(ax3b2, 'Cluster-forming statistic', 'FontSize', 10, 'FontWeight', 'bold');
    legend(ax3b2, 'Location', 'best', 'FontSize', 8);
    grid(ax3b2, 'on');

    sgtitle(fig3b, sprintf(['Group Analysis — CVSA-Fusion Cluster Test (n=%d subjects)\n' ...
            'each subject contributes their own mean P_{fused}-P_{MI} curve; cluster permutation controls for multiple comparisons across time'], ...
            n_cluster_subj), 'Interpreter', 'none');

    saveas(fig3b, fullfile(out_dir, 'group_fusion_cluster_test.svg'), 'svg');
    if ~SHOW_FIGURES, close(fig3b); end
end

% ── Figure 3c: CVSA-help significance from the INTEGRATOR counterfactual
%    (main_hybrid_advantage_integ.m's buf_adv_mi_win) -- per-subject +
%    group-level answer to "does CVSA help", isolated to the cvsa_influence
%    window. Complements Figure 3 above, which uses the raw sLDA-probability
%    fusion advantage from main_hybrid_advantage_probs.m instead of the
%    integrator buffer output. Skipped if no subject has this field yet
%    (counterfactual_summary.mat regenerated after this test was added). ──
if any(~isnan(buf_adv_mi_win_subj))
    fig3c = figure('Name', 'Group Analysis — CVSA-Help Significance (Integrator)', 'Color', 'w', ...
                   'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig3c, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    has_bw = ~isnan(buf_adv_mi_win_subj);

    % Panel 1: per-subject mean advantage, bar coloured by that SUBJECT'S
    % OWN significance (buf_adv_mi_win_meta_p_subj, well-powered -- combines
    % all of that subject's own trials/sessions), plus a GROUP bar (Stouffer
    % meta-combination across subjects, no small-n permutation floor).
    ax3c1 = subplot(1, 2, 1); hold(ax3c1, 'on');
    bar_cols = repmat([0.6 0.6 0.6], n_subj, 1);
    sig_subj = has_bw & buf_adv_mi_win_meta_p_subj < 0.05;
    bar_cols(sig_subj, :) = repmat(COL_hybrid, sum(sig_subj), 1);
    b3c1 = bar(ax3c1, 1:n_subj, buf_adv_mi_win_subj, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'k');
    b3c1.CData = bar_cols;
    yline(ax3c1, 0, 'k-', 'HandleVisibility', 'off');
    for si = 1:n_subj
        if isnan(buf_adv_mi_win_subj(si)), continue; end
        if buf_adv_mi_win_subj(si) >= 0, va = 'bottom'; y_off = 0.003; else, va = 'top'; y_off = -0.003; end
        text(ax3c1, si, buf_adv_mi_win_subj(si) + y_off, ...
             sprintf('%s\np=%.3f', stars(buf_adv_mi_win_meta_p_subj(si)), buf_adv_mi_win_meta_p_subj(si)), ...
             'HorizontalAlignment', 'center', 'FontSize', 7, 'VerticalAlignment', va);
    end
    m_bw  = mean(buf_adv_mi_win_subj, 'omitnan');
    se_bw = std(buf_adv_mi_win_subj, 'omitnan') / sqrt(max(1, sum(has_bw)));
    errorbar(ax3c1, n_subj+0.8, m_bw, se_bw, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'CapSize', 6);
    text(ax3c1, n_subj+0.8, m_bw + se_bw*sign(m_bw+eps) + 0.01, stars(p_bufadvwin_meta_grp), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    set(ax3c1, 'XTick', [1:n_subj, n_subj+0.8], 'XTickLabel', [subjects, {'GROUP'}], 'XLim', [0.4, n_subj+1.2]);
    ylabel(ax3c1, 'mean buffer(target) advantage  [Hybrid - MI-only]');
    title(ax3c1, sprintf('Per-subject CVSA-help (integrator, cvsa-influence window)\ncolor = significant for that subject (p<0.05)  |  group p=%.4f %s', ...
          p_bufadvwin_meta_grp, stars(p_bufadvwin_meta_grp)), 'FontWeight', 'bold', 'FontSize', 9, 'Interpreter', 'none');
    grid(ax3c1, 'on');

    % Panel 2: the p-value itself per subject, so "is it significant for me"
    % can be read directly without cross-referencing the console table.
    ax3c2 = subplot(1, 2, 2); hold(ax3c2, 'on');
    for si = 1:n_subj
        if isnan(buf_adv_mi_win_meta_p_subj(si)), continue; end
        if buf_adv_mi_win_meta_p_subj(si) < 0.05, col = COL_hybrid; else, col = [0.6 0.6 0.6]; end
        bar(ax3c2, si, buf_adv_mi_win_meta_p_subj(si), 0.5, 'FaceColor', col, 'EdgeColor', 'k');
    end
    yline(ax3c2, 0.05, 'r--', 'p=0.05', 'LabelHorizontalAlignment', 'left');
    set(ax3c2, 'XTick', 1:n_subj, 'XTickLabel', subjects, 'XLim', [0.4, n_subj+0.6]);
    ylabel(ax3c2, 'within-subject one-sided p-value ("CVSA helps")');
    title(ax3c2, 'Is CVSA help significant for THIS subject?', 'FontWeight', 'bold', 'FontSize', 9);
    grid(ax3c2, 'on');

    sgtitle(fig3c, sprintf(['Group Analysis — CVSA-Help Significance, Integrator Buffer (n=%d subjects)\n' ...
            'Hybrid vs MI-only, restricted to the first cvsa_influence seconds of CF -- isolates CVSA''s actual influence window'], ...
            sum(has_bw)), 'Interpreter', 'none');

    saveas(fig3c, fullfile(out_dir, 'group_cvsa_help_integrator.svg'), 'svg');
    if ~SHOW_FIGURES, close(fig3c); end
end

% ── Figure 4: ERD/ERS discrimination + CSP grounding across subjects ───────
fig4 = figure('Name', 'Group Analysis — ERD/ERS Across Subjects', 'Color', 'w', ...
              'NumberTitle', 'off', 'Visible', fig_vis);
set(fig4, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

ax6 = subplot(1, 2, 1); hold(ax6, 'on');
for pi = 1:3
    x = (1:n_subj) + (pi-2)*bw;
    v = erd_discrim_subj(:,pi);
    bar(ax6, x, v, bw*0.9, 'FaceColor', cols3{pi}, 'EdgeColor', 'k', 'LineWidth', 1, ...
        'DisplayName', par_labels{pi});
    m = mean(v, 'omitnan');
    if ~isnan(m)
        plot(ax6, [0.5, n_subj+0.5], [m m], '--', 'Color', cols3{pi}*0.6, ...
             'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end
set(ax6, 'XTick', 1:n_subj, 'XTickLabel', subjects, 'XLim', [0.4, n_subj+0.6]);
ylabel(ax6, 'mean |ERD/ERS class1 - class2|  (%)');
legend(ax6, 'Location', 'best', 'FontSize', 9);
title(ax6, 'ERD/ERS class discrimination per subject, per task', 'FontWeight', 'bold');
grid(ax6, 'on');

ax7 = subplot(1, 2, 2); hold(ax7, 'on');
for pi = 1:3
    x = (1:n_subj) + (pi-2)*bw;
    v = erd_r_subj(:,pi);
    bar(ax7, x, v, bw*0.9, 'FaceColor', cols3{pi}, 'EdgeColor', 'k', 'LineWidth', 1, ...
        'DisplayName', par_labels{pi});
    m = mean(v, 'omitnan'); se = std(v, 'omitnan') / sqrt(max(1,sum(~isnan(v))));
    if ~isnan(m)
        plot(ax7, [0.5, n_subj+0.5], [m m], '--', 'Color', cols3{pi}*0.6, ...
             'LineWidth', 1.5, 'HandleVisibility', 'off');
        text(ax7, n_subj+0.55+(pi-2)*0.15, m, stars(p_r_grp(pi)), 'FontSize', 9, ...
             'Color', cols3{pi}*0.6, 'FontWeight', 'bold');
    end
end
yline(ax7, 0, 'k-', 'HandleVisibility', 'off');
set(ax7, 'XTick', 1:n_subj, 'XTickLabel', subjects, 'XLim', [0.4, n_subj+0.6]);
ylabel(ax7, 'ERD/ERS discrimination vs CSP-weight  (Pearson r)');
legend(ax7, 'Location', 'best', 'FontSize', 9);
title(ax7, 'CSP-weight grounding per subject, per task (group sign-flip test, dashed = group mean)', ...
      'FontWeight', 'bold', 'FontSize', 9);
grid(ax7, 'on');

sgtitle(fig4, sprintf('Group Analysis — ERD/ERS Across Subjects (n=%d subjects)', n_subj), 'Interpreter', 'none');

saveas(fig4, fullfile(out_dir, 'group_erd_across_subjects.svg'), 'svg');
if ~SHOW_FIGURES, close(fig4); end

% ── Figure 5: classifier-level ROC across subjects, one panel per paradigm ─
paradigms_with_roc = pars(any(~isnan(roc_auc_subj), 1));
if ~isempty(paradigms_with_roc)
    fig5 = figure('Name', 'Group Analysis — ROC Across Subjects', 'Color', 'w', ...
                  'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig5, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    for pi = 1:3
        if all(isnan(roc_auc_subj(:,pi))), continue; end
        axp = subplot(1, 3, pi); hold(axp, 'on');
        plot(axp, [0 1], [0 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
        col_light = cols3{pi} + (1 - cols3{pi}) * 0.65;   % tinted toward white, no alpha needed
        for si = 1:n_subj
            curve = reshape(roc_tpr_subj_grid(si,pi,:), 1, numel(ROC_FPR_GRID));
            if all(isnan(curve)), continue; end
            plot(axp, ROC_FPR_GRID, curve, '-', 'Color', col_light, ...
                 'LineWidth', 1, 'HandleVisibility', 'off');
        end
        m  = roc_mean_tpr(pi,:);
        se = roc_sem_tpr(pi,:);
        fill(axp, [ROC_FPR_GRID, fliplr(ROC_FPR_GRID)], [m+se, fliplr(m-se)], cols3{pi}, ...
             'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(axp, ROC_FPR_GRID, m, '-', 'Color', cols3{pi}, 'LineWidth', 2.5, ...
             'DisplayName', 'cross-subject mean \pm SEM');
        xlabel(axp, 'False Positive Rate'); ylabel(axp, 'True Positive Rate');
        xlim(axp, [0 1]); ylim(axp, [0 1]); axis(axp, 'square'); grid(axp, 'on');
        n_v = sum(~isnan(roc_auc_subj(:,pi)));
        title(axp, sprintf('%s  (n=%d subjects)\nmean AUC=%.3f  p(AUC>0.5)=%.4f %s', ...
              par_labels{pi}, n_v, mean(roc_auc_subj(:,pi),'omitnan'), p_auc_grp(pi), stars(p_auc_grp(pi))), ...
              'FontWeight', 'bold', 'FontSize', 9);
        legend(axp, 'Location', 'southeast', 'FontSize', 8);
    end
    sgtitle(fig5, sprintf(['Group Analysis — Classifier ROC/AUC, pre-integrator (n=%d subjects)\n' ...
            'thin = per-subject re-pooled curve, thick = cross-subject macro-average (interpolated onto a common FPR grid)'], n_subj), ...
            'Interpreter', 'none');

    saveas(fig5, fullfile(out_dir, 'group_roc_analysis.svg'), 'svg');
    if ~SHOW_FIGURES, close(fig5); end
else
    fprintf('\n[main_group_analysis] No ROC data found (run main_roc_analysis first) — skipping Fig 5.\n');
end

% ── Figure 6: cross-subject correlates of real Hybrid performance (bonus) ──
% Four scatter panels: does a subject's counterfactual advantage / fusion
% mechanism / neurophysiological grounding / CVSA-only "quality" predict how
% well they actually do in real Hybrid sessions? The 4th panel directly
% supports the "CVSA helps even when imperfect" thesis: if the Hybrid-MI
% advantage stays positive even for subjects with weak CVSA-only accuracy,
% adding CVSA is not merely helping when it happens to be good.
fig6 = figure('Name', 'Group Analysis — Performance Correlates', 'Color', 'w', ...
              'NumberTitle', 'off', 'Visible', fig_vis);
set(fig6, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

ax_c1 = subplot(2, 2, 1);
scatter_with_fit(ax_c1, 100*d_mi_subj, 100*d_real_hyb_mi, subjects, r_cf_vs_real, p_cf_vs_real, ...
    '\Delta acc, counterfactual (Hybrid-MI, %)', '\Delta acc, real (Hybrid-MI, %)', ...
    'Counterfactual advantage predicts real advantage?', COL_hybrid);

ax_c2 = subplot(2, 2, 2);
scatter_with_fit(ax_c2, fus_adv_subj, 100*real_acc(:,3), subjects, r_fus_vs_real, p_fus_vs_real, ...
    'fusion advantage  mean(P_{fused}-P_{MI})', 'real Hybrid HIT rate (%)', ...
    'Fusion mechanism predicts real accuracy?', COL_hybrid);

ax_c3 = subplot(2, 2, 3);
scatter_with_fit(ax_c3, erd_r_subj(:,3), 100*real_acc(:,3), subjects, r_erd_vs_real, p_erd_vs_real, ...
    'Hybrid ERD/ERS vs CSP-weight (Pearson r)', 'real Hybrid HIT rate (%)', ...
    'Neurophysiological grounding predicts real accuracy?', COL_hybrid);

ax_c4 = subplot(2, 2, 4);
scatter_with_fit(ax_c4, 100*cf_acc_dec(:,3), 100*d_mi_subj_dec, subjects, r_cvsa_quality, p_cvsa_quality, ...
    'CVSA-only decided-trial accuracy (%)', '\Delta acc, counterfactual (Hybrid-MI, %)', ...
    'Does the Hybrid advantage survive weak CVSA?', COL_cvsa_only);
yline(ax_c4, 0, 'k-', 'HandleVisibility', 'off');

sgtitle(fig6, sprintf(['Group Analysis — Cross-Subject Correlates of Real Hybrid Performance (n=%d subjects)\n' ...
        'each point = one subject; permutation test on Pearson r'], n_subj), 'Interpreter', 'none');

saveas(fig6, fullfile(out_dir, 'group_performance_correlates.svg'), 'svg');
if ~SHOW_FIGURES, close(fig6); end

% ── Save group summary .mat ────────────────────────────────────────────────
group_summary = struct();
group_summary.subjects = subjects;
group_summary.real_acc     = real_acc;      % [n_subj x 3] MI, CVSA, Hybrid -- TIMEOUT counted as fail
group_summary.real_acc_dec = real_acc_dec;  % [n_subj x 3] decided trials only, TIMEOUT excluded
group_summary.real_fp_rate = real_fp_rate;  % [n_subj x 3] MISS / n_total ("false positive": wrong class confidently triggered)
group_summary.real_to_rate = real_to_rate;  % [n_subj x 3]
group_summary.real_n_hit   = real_n_hit;    % [n_subj x 3]
group_summary.real_n_miss  = real_n_miss;   % [n_subj x 3]
group_summary.real_to_n    = real_to_n;     % [n_subj x 3] raw timeout counts
group_summary.real_n       = real_n;
group_summary.real_tth     = real_tth;
group_summary.real_t_miss  = real_t_miss;
group_summary.real_to_time = real_to_time;
group_summary.grand_acc     = grand_acc;     % [1 x 3] trial-pooled across all subjects
group_summary.grand_acc_dec = grand_acc_dec;
group_summary.grand_fp_rate = grand_fp_rate;
group_summary.grand_to_rate = grand_to_rate;
group_summary.grand_to_n    = grand_to_n;    % [1 x 3] raw timeout counts, trial-pooled
group_summary.grand_n       = grand_n;
group_summary.grand_tth     = grand_tth;
group_summary.grand_t_miss  = grand_t_miss;
group_summary.grand_to_time = grand_to_time;
group_summary.cf_acc     = cf_acc;       % [n_subj x 3] Hybrid, MI-only, CVSA-only -- TIMEOUT counted as fail
group_summary.cf_acc_dec = cf_acc_dec;   % [n_subj x 3] decided trials only
group_summary.cf_fp_rate = cf_fp_rate;
group_summary.cf_to_rate = cf_to_rate;
group_summary.cf_n_hit   = cf_n_hit;
group_summary.cf_n_miss  = cf_n_miss;
group_summary.cf_n_to    = cf_n_to;
group_summary.cf_n       = cf_n;
group_summary.grand_cf_acc     = grand_cf_acc;      % [1x3] trial-pooled across all subjects
group_summary.grand_cf_acc_dec = grand_cf_acc_dec;
group_summary.grand_cf_fp_rate = grand_cf_fp_rate;
group_summary.grand_cf_to_rate = grand_cf_to_rate;
group_summary.p_mi_group  = p_mi_grp;
group_summary.p_cvsa_group = p_cvs_grp;
group_summary.d_mi_group_cohen  = d_mi_grp_cohen;      % paired Cohen's d, subject is the unit (companion to p_mi_group)
group_summary.d_cvs_group_cohen = d_cvs_grp_cohen;
group_summary.p_mi_group_dec  = p_mi_grp_dec;    % decided-trials-only twin
group_summary.p_cvsa_group_dec = p_cvs_grp_dec;
group_summary.d_mi_group_dec_cohen  = d_mi_grp_dec_cohen;
group_summary.d_cvs_group_dec_cohen = d_cvs_grp_dec_cohen;
group_summary.fus_adv_subj    = fus_adv_subj;
group_summary.rescue_subj     = rescue_subj;
group_summary.cost_subj       = cost_subj;
group_summary.net_rescue_subj = net_rescue_subj;
group_summary.p_fusadv_group  = p_fusadv_grp;
group_summary.p_rescue_group  = p_resc_grp;
group_summary.d_fusadv_group_cohen = d_fusadv_cohen;   % paired Cohen's d (companion to p_fusadv_group)
group_summary.d_rescue_group_cohen = d_resc_cohen;     % paired Cohen's d (companion to p_rescue_group)
group_summary.fus_adv_meta_p_subj  = fus_adv_meta_p_subj;   % per-subject combined within-session one-sided p
group_summary.fus_adv_meta_n_subj  = fus_adv_meta_n_subj;   % per-subject total trials contributing
group_summary.p_fusadv_meta_group  = p_fusadv_meta_grp;     % Stouffer-combined group-level one-sided p (no small-n floor)
group_summary.buf_adv_mi_win_subj        = buf_adv_mi_win_subj;         % per-subject mean Hybrid-MI buffer advantage, cvsa_influence window
group_summary.buf_adv_mi_win_meta_p_subj = buf_adv_mi_win_meta_p_subj;  % per-subject combined within-session one-sided p ("CVSA helps")
group_summary.buf_adv_mi_win_meta_n_subj = buf_adv_mi_win_meta_n_subj;  % per-subject total trials contributing
group_summary.p_bufadvwin_meta_group     = p_bufadvwin_meta_grp;        % Stouffer-combined group-level one-sided p (no small-n floor)
group_summary.cluster_time_grid = CLUSTER_TIME_GRID;   % [1 x N] common time grid, empty if no subject has cluster data
group_summary.cluster_diff_subj = cluster_diff_subj;   % [n_subj x N] per-subject mean P_fused-P_MI curve
if ~isempty(CLUSTER_TIME_GRID)
    group_summary.cluster_clusters_group = cluster_res_grp.clusters;  % struct array: start_t/end_t/mass/p
    group_summary.cluster_tstat_group    = cluster_res_grp.tstat;
    group_summary.cluster_sig_mask_group = cluster_res_grp.sig_mask;
end
group_summary.erd_discrim_subj = erd_discrim_subj;   % [n_subj x 3] MI, CVSA, Hybrid
group_summary.erd_r_subj       = erd_r_subj;
group_summary.p_erd_r_group    = p_r_grp;            % [1x3]
group_summary.roc_auc_subj     = roc_auc_subj;        % [n_subj x 3] MI, CVSA, Hybrid -- per-subject, re-pooled across sessions
group_summary.roc_fpr_grid     = ROC_FPR_GRID;        % [1 x 101] common FPR grid used for cross-subject averaging
group_summary.roc_mean_tpr     = roc_mean_tpr;        % [3 x 101] cross-subject mean TPR curve, per paradigm
group_summary.roc_sem_tpr      = roc_sem_tpr;         % [3 x 101] cross-subject SEM band
group_summary.p_auc_group      = p_auc_grp;           % [1x3] sign-flip test on (subject AUC - 0.5)
group_summary.p_real_hyb_mi    = p_real_hyb_mi;       % subject-level sign-flip on REAL-session accuracy deltas (TIMEOUT=fail)
group_summary.p_real_hyb_cvs   = p_real_hyb_cvs;
group_summary.p_real_mi_cvs    = p_real_mi_cvs;
group_summary.p_friedman_real  = p_friedman_real;     % omnibus rank-based test, MI vs CVSA vs Hybrid, real accuracy
group_summary.d_real_hyb_mi_cohen  = d_real_hyb_mi_cohen;   % paired Cohen's d (companion to p_real_hyb_mi etc.)
group_summary.d_real_hyb_cvs_cohen = d_real_hyb_cvs_cohen;
group_summary.d_real_mi_cvs_cohen  = d_real_mi_cvs_cohen;
group_summary.p_real_hyb_mi_dec   = p_real_hyb_mi_dec;    % decided-trials-only twin
group_summary.p_real_hyb_cvs_dec  = p_real_hyb_cvs_dec;
group_summary.p_real_mi_cvs_dec   = p_real_mi_cvs_dec;
group_summary.p_friedman_real_dec = p_friedman_real_dec;
group_summary.d_real_hyb_mi_dec_cohen  = d_real_hyb_mi_dec_cohen;
group_summary.d_real_hyb_cvs_dec_cohen = d_real_hyb_cvs_dec_cohen;
group_summary.d_real_mi_cvs_dec_cohen  = d_real_mi_cvs_dec_cohen;
group_summary.idx_best_real    = {idx_best_real};     % 1x3 cell (MI,CVSA,Hybrid), each a vector of subject indices (ties), by decided-trial accuracy
group_summary.idx_worst_real   = {idx_worst_real};
group_summary.idx_best_cf      = {idx_best_cf};       % 1x3 cell (Hybrid,MI-only,CVSA-only), by decided-trial accuracy
group_summary.idx_worst_cf     = {idx_worst_cf};
group_summary.r_cf_vs_real     = r_cf_vs_real;        % cross-subject Pearson r: counterfactual adv. vs real adv. (Hybrid-MI)
group_summary.p_cf_vs_real     = p_cf_vs_real;
group_summary.r_fus_vs_real    = r_fus_vs_real;       % cross-subject Pearson r: fusion advantage vs real Hybrid accuracy
group_summary.p_fus_vs_real    = p_fus_vs_real;
group_summary.r_erd_vs_real    = r_erd_vs_real;       % cross-subject Pearson r: ERD-CSP grounding (Hybrid) vs real Hybrid accuracy
group_summary.p_erd_vs_real    = p_erd_vs_real;
group_summary.r_cvsa_quality   = r_cvsa_quality;      % cross-subject Pearson r: CVSA-only decided accuracy vs Hybrid-MI counterfactual advantage
group_summary.p_cvsa_quality   = p_cvsa_quality;
save(fullfile(out_dir, 'group_summary.mat'), 'group_summary');
fprintf('\nSaved group_summary.mat and figures to %s\n', out_dir);

end % main_group_analysis

% ── Local helpers ────────────────────────────────────────────────────────────

function plot_subject_grand_bars(ax, val_mat, grand_vec, subjects, labels, cols, bw, is_pct)
% PLOT_SUBJECT_GRAND_BARS  Grouped bars (one group per column of val_mat)
%   across subjects, plus a final "GRAND" group built from grand_vec
%   (trial-pooled, NOT the mean of the per-subject bars) -- visually offset
%   by a vertical separator line and a thicker bar edge.
    if nargin < 8, is_pct = true; end
    n_subj   = size(val_mat, 1);
    n_series = size(val_mat, 2);
    scale = 1; if is_pct, scale = 100; end
    x_subj  = 1:n_subj;
    x_grand = n_subj + 1.2;
    for k = 1:n_series
        offset = (k - (n_series+1)/2) * bw;
        bar(ax, x_subj + offset, scale*val_mat(:,k), bw*0.9, ...
            'FaceColor', cols{k}, 'EdgeColor', 'k', 'LineWidth', 1, 'DisplayName', labels{k});
        gv = scale * grand_vec(k);
        if ~isnan(gv)
            bar(ax, x_grand + offset, gv, bw*0.9, ...
                'FaceColor', cols{k}, 'EdgeColor', 'k', 'LineWidth', 2.2, 'HandleVisibility', 'off');
        end
    end
    xline(ax, n_subj + 0.6, 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    set(ax, 'XTick', [x_subj, x_grand], 'XTickLabel', [subjects, {'GRAND'}], 'XLim', [0.4, x_grand+0.9]);
    if is_pct, set(ax, 'YLim', [0, 110]); end
end

function subj = subject_of(fpath, root_dir)
    rel = erase(fpath, [char(root_dir) filesep]);
    parts = strsplit(rel, filesep);
    subj = parts{1};
end

function files_out = filter_files_by_subject(files_in, root_dir, subjects_filter)
% FILTER_FILES_BY_SUBJECT  Keep only dir() entries whose subject (first path
%   component under root_dir) is in subjects_filter.
    if isempty(files_in), files_out = files_in; return; end
    keep = false(numel(files_in), 1);
    for i = 1:numel(files_in)
        fpath = fullfile(files_in(i).folder, files_in(i).name);
        keep(i) = ismember(subject_of(fpath, root_dir), subjects_filter);
    end
    files_out = files_in(keep);
end

function s = fmt_pct(v)
    if isnan(v), s = '--'; else, s = sprintf('%.1f%%', 100*v); end
end

function [idx_best, idx_worst] = best_worst_ties(v)
% BEST_WORST_TIES  Indices of the max/min value in the column vector v,
%   INCLUDING ties (returns every subject tied for best, every subject tied
%   for worst). Returns empty row vectors if fewer than 2 valid (non-NaN)
%   entries, or if every valid entry is tied at the same value (no
%   meaningful best/worst distinction).
    valid = find(~isnan(v));
    if numel(valid) < 2
        idx_best = []; idx_worst = []; return;
    end
    vv = v(valid);
    idx_best  = valid(vv == max(vv))';
    idx_worst = valid(vv == min(vv))';
    if isequal(sort(idx_best), sort(idx_worst))
        idx_best = []; idx_worst = [];
    end
end

function print_ranking_table(label, subjects, v, idx_best, idx_worst)
% PRINT_RANKING_TABLE  Console ranking (descending) of subjects by v, with
%   BEST/WORST tags for indices in idx_best/idx_worst (ties share the tag).
    fprintf('\n──── %s ────\n', label);
    valid = find(~isnan(v));
    if isempty(valid)
        fprintf('  (no subject has valid data)\n');
        return;
    end
    [~, ord] = sort(v(valid), 'descend');
    ranked = valid(ord);
    for ri = 1:numel(ranked)
        si = ranked(ri);
        tag = '';
        if any(idx_best == si),  tag = [tag '  <-- BEST'];  end
        if any(idx_worst == si), tag = [tag '  <-- WORST']; end
        fprintf('  %2d. %-10s  %8s%s\n', ri, subjects{si}, fmt_pct(v(si)), tag);
    end
end

function s = stars(p)
    if isnan(p),      s = '';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,              s = 'n.s.';
    end
end

function d = cohen_d_local(diffs)
% COHEN_D_LOCAL  Paired Cohen's d = mean(diff)/std(diff) across subjects --
%   for paired/one-sample data this is the standard effect size on the
%   differences, the natural companion to the sign-flip test's p-value
%   (subject is the unit, same convention as sign_flip_test_local).
    diffs = diffs(~isnan(diffs));
    if numel(diffs) < 2, d = NaN; return; end
    s = std(diffs);
    if s == 0, d = NaN; return; end
    d = mean(diffs) / s;
end

function s = cohen_d_label(d)
% COHEN_D_LABEL  Conventional |d| magnitude bins (Cohen, 1988).
    if isnan(d), s = ''; return; end
    ad = abs(d);
    if     ad < 0.2, s = 'negligible';
    elseif ad < 0.5, s = 'small';
    elseif ad < 0.8, s = 'medium';
    else,            s = 'large';
    end
end

function p = sign_flip_test_local(diffs, n_perm)
    diffs = diffs(~isnan(diffs));
    n = numel(diffs);
    if n < 2, p = NaN; return; end
    obs = abs(mean(diffs));
    cnt = 0;
    for i = 1:n_perm
        signs = sign(rand(n, 1) - 0.5);
        if abs(mean(diffs .* signs)) >= obs - 1e-12
            cnt = cnt + 1;
        end
    end
    p = (cnt + 1) / (n_perm + 1);
end

function [z_comb, p_comb] = stouffer_combine_local(p_one_sided, weights)
% STOUFFER_COMBINE_LOCAL  Weighted Stouffer's Z-score method: combines
%   independent one-sided p-values, all testing the SAME directional
%   hypothesis, into one combined z-score/p-value. No Statistics Toolbox
%   required (norminv/normcdf hand-rolled below via erfcinv/erfc, both base
%   MATLAB). weights default to equal (1 each) if omitted -- unweighted
%   Stouffer, used for the across-subject combination so every subject
%   counts equally, matching every other group-level test in this script.
%   Returns NaN/NaN if no valid (non-NaN, positive-weight) inputs remain.
    if nargin < 2 || isempty(weights), weights = ones(size(p_one_sided)); end
    valid = ~isnan(p_one_sided) & ~isnan(weights) & weights > 0;
    p_one_sided = p_one_sided(valid);
    weights = weights(valid);
    if isempty(p_one_sided), z_comb = NaN; p_comb = NaN; return; end
    z_i = norminv_local(1 - p_one_sided);
    z_comb = sum(weights .* z_i) / sqrt(sum(weights.^2));
    p_comb = 1 - normcdf_local(z_comb);
end

function z = norminv_local(p)
% NORMINV_LOCAL  Standard normal inverse CDF, via erfcinv (base MATLAB, no
%   Statistics Toolbox). Clamped away from 0/1 to avoid +-Inf.
    p = min(max(p, eps), 1 - eps);
    z = -sqrt(2) * erfcinv(2 * p);
end

function p = normcdf_local(z)
% NORMCDF_LOCAL  Standard normal CDF, via erfc (base MATLAB).
    p = 0.5 * erfc(-z / sqrt(2));
end

function p = friedman_test_local(mat, n_perm)
% FRIEDMAN_TEST_LOCAL  Rank-based omnibus test for >=2 repeated measures
%   (here: paradigms) across subjects (rows), via exact-style permutation
%   (shuffle each subject's own ranks independently) -- no Statistics
%   Toolbox required, consistent with the sign-flip tests elsewhere in this
%   script.
    valid = ~any(isnan(mat), 2);
    mat = mat(valid, :);
    n = size(mat, 1);
    k = size(mat, 2);
    if n < 2, p = NaN; return; end
    obs = friedman_stat_local(mat);
    cnt = 0;
    for i = 1:n_perm
        permuted = mat;
        for r = 1:n
            permuted(r, :) = permuted(r, randperm(k));
        end
        if friedman_stat_local(permuted) >= obs - 1e-12
            cnt = cnt + 1;
        end
    end
    p = (cnt + 1) / (n_perm + 1);
end

function stat = friedman_stat_local(mat)
    [n, k] = size(mat);
    ranks = zeros(n, k);
    for r = 1:n
        ranks(r, :) = tiedrank_local(mat(r, :));
    end
    R = sum(ranks, 1);
    stat = (12 / (n*k*(k+1))) * sum(R.^2) - 3*n*(k+1);
end

function r = tiedrank_local(v)
% TIEDRANK_LOCAL  Ranks a row vector, averaging ranks within tied groups
%   (equivalent to Statistics Toolbox's tiedrank, hand-rolled to avoid the
%   dependency).
    [~, idx] = sort(v);
    r = zeros(size(v));
    r(idx) = 1:numel(v);
    u = unique(v);
    for i = 1:numel(u)
        mask = (v == u(i));
        if sum(mask) > 1
            r(mask) = mean(r(mask));
        end
    end
end

function [r, p] = pearson_perm_test_local(x, y, n_perm)
% PEARSON_PERM_TEST_LOCAL  Pearson correlation across subjects with a
%   two-sided permutation p-value (shuffle y), no Statistics Toolbox
%   required.
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid); y = y(valid);
    n = numel(x);
    if n < 3, r = NaN; p = NaN; return; end
    r = corr_local(x, y);
    obs = abs(r);
    if isnan(obs), p = NaN; return; end   % zero-variance predictor/outcome -- no correlation is defined
    cnt = 0;
    for i = 1:n_perm
        yp = y(randperm(n));
        if abs(corr_local(x, yp)) >= obs - 1e-12
            cnt = cnt + 1;
        end
    end
    p = (cnt + 1) / (n_perm + 1);
end

function r = corr_local(x, y)
    x = x(:); y = y(:);
    xm = x - mean(x); ym = y - mean(y);
    denom = sqrt(sum(xm.^2) * sum(ym.^2));
    if denom == 0, r = NaN; return; end
    r = sum(xm .* ym) / denom;
end

function scatter_with_fit(ax, x, y, labels, r, p, xlab, ylab, ttl, col)
% SCATTER_WITH_FIT  One cross-subject scatter panel: points + linear fit +
%   subject-name labels + Pearson r / permutation p in the title.
    x = x(:); y = y(:);
    valid = ~isnan(x) & ~isnan(y);
    hold(ax, 'on');
    scatter(ax, x(valid), y(valid), 70, col, 'filled', 'MarkerEdgeColor', 'k');
    if sum(valid) >= 2 && range(x(valid)) > 0
        pf = polyfit(x(valid), y(valid), 1);
        xx = linspace(min(x(valid)), max(x(valid)), 50);
        plot(ax, xx, polyval(pf, xx), '--', 'Color', col*0.6, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    idxv = find(valid);
    for k = 1:numel(idxv)
        text(ax, x(idxv(k)), y(idxv(k)), ['  ' labels{idxv(k)}], 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
    end
    xlabel(ax, xlab); ylabel(ax, ylab);
    title(ax, sprintf('%s\nr=%+.2f, p=%.4f %s', ttl, r, p, stars(p)), 'FontWeight', 'bold', 'FontSize', 9);
    grid(ax, 'on');
end
