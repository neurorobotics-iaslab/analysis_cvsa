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

function group_summary = main_group_analysis(root_dir, subjects_filter, show_figures, group_tag)
%   Callable as a function:
%     main_group_analysis()                                  % default recordings root, ALL subjects
%     main_group_analysis(root_dir)                           % given root, ALL subjects
%     main_group_analysis(root_dir, {'a1','a2','a3'})        % only these subjects
%     main_group_analysis(root_dir, {'a1','a2','a3'}, true)  % also show figures on screen
%     main_group_analysis([], {'a1','a2'})                    % default root, subset of subjects
%     main_group_analysis(root_dir, {'a1','a2'}, false, 'well') % tag output files 'well'
%   No GUI folder picker: root_dir defaults to DEFAULT_ROOT below (edit it
%   directly, or pass a path) rather than prompting. subjects_filter is an
%   explicit opt-in filter (cellstr of subject IDs); omit it (or pass {}) to
%   use every subject found under the root.
%
%   group_tag (optional, default 'all') is stitched into every output
%   filename ('<num>_<tag>_<name>.svg', 'group_summary_<tag>.mat') so that
%   several calls with different subject subsets (e.g. a manually-assigned
%   'well'/'bad' performer split -- see batch_group_analysis.m, which drives
%   this) can coexist under the same <root>/group_analysis/ folder without
%   overwriting each other. This function itself only ever analyzes ONE
%   group per call; running several groups is orchestrated by the caller.
%   Returns the group_summary struct so a caller can feed it (e.g. for a
%   'well' and a 'bad' run) into well_vs_bad_comparison.m.

DEFAULT_ROOT = '/home/paolo/bci_vr_ws/recordings';

if nargin < 4 || isempty(group_tag), group_tag = 'all'; end
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
                    'tth_vals', {}, 't_miss_vals', {}, 'to_vals', {}, ...
                    'sigqc_n_dead', {}, 'sigqc_n_noisy', {});
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
        % Backward-compatible: older session_summary.mat files (saved before
        % the raw-signal-quality check existed) simply don't have these.
        if isfield(sf(k), 'sigqc_n_dead')
            flat_real(idx).sigqc_n_dead  = sf(k).sigqc_n_dead;
            flat_real(idx).sigqc_n_noisy = sf(k).sigqc_n_noisy;
        else
            flat_real(idx).sigqc_n_dead  = [];
            flat_real(idx).sigqc_n_noisy = [];
        end
    end
end

% ── Flatten counterfactual per-file records, tagged with subject ───────────
flat_cf = struct('subject', {}, 'acc', {}, 'n_trials', {}, 't_hit_mean', {}, ...
                  'n_hit', {}, 'n_miss', {}, 'n_to', {}, ...
                  't_miss_mean', {}, 't_to_mean', {}, ...
                  'buf_adv_mi_win_mean', {}, 'buf_adv_mi_win_pval', {}, 'buf_adv_mi_win_n_trials', {}, ...
                  'buf_adv_cvs_win_mean', {}, ...
                  'oc_hyb', {}, 'oc_mi', {}, 'oc_cvsa', {}, 'cvsa_inf_window', {}, ...
                  'frac_correct_raw', {}, 'frac_correct_int', {}, 'frac_correct_n_trials', {});
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
    % t_miss_mean/t_to_mean already existed in counterfactual_summary.mat before
    % this loader read them -- isfield-guarded only for symmetry with the newer
    % buf_adv_mi_win_* fields below, not because they're actually new.
    if isfield(cf, 't_miss_mean'), flat_cf(idx).t_miss_mean = cf.t_miss_mean; end % [1x3]
    if isfield(cf, 't_to_mean'),   flat_cf(idx).t_to_mean   = cf.t_to_mean;   end % [1x3]
    % buf_adv_mi_win_* only present in counterfactual_summary.mat regenerated
    % after the cvsa_influence-window restricted test was added -- guarded
    % with isfield so older summaries don't crash the loader, same
    % backward-compatibility convention as flat_probs.cluster_diff_mean below.
    if isfield(cf, 'buf_adv_mi_win_pval')
        flat_cf(idx).buf_adv_mi_win_mean     = cf.buf_adv_mi_win_mean;
        flat_cf(idx).buf_adv_mi_win_pval     = cf.buf_adv_mi_win_pval;
        flat_cf(idx).buf_adv_mi_win_n_trials = cf.buf_adv_mi_win_n_trials;
    end
    if isfield(cf, 'buf_adv_cvs_win_mean')
        flat_cf(idx).buf_adv_cvs_win_mean = cf.buf_adv_cvs_win_mean;
    end
    % oc_hyb/oc_mi/oc_cvsa (raw per-trial [outcome_code, t_event_s]) only
    % present in counterfactual_summary.mat regenerated after they were
    % added -- lets main_group_analysis.m split trials by WHEN they
    % resolved (within vs after cvsa_influence) without re-running the
    % pipeline.
    if isfield(cf, 'oc_hyb')
        flat_cf(idx).oc_hyb  = cf.oc_hyb;
        flat_cf(idx).oc_mi   = cf.oc_mi;
        flat_cf(idx).oc_cvsa = cf.oc_cvsa;
        flat_cf(idx).cvsa_inf_window = cf.cvsa_inf_window;
    end
    % frac_correct_raw/int ([1x3] Hybrid/MI-only/CVSA-only, fraction of CF
    % frames with target-class signal > 0.5, ALL trials, whole CF duration)
    % only present in counterfactual_summary.mat regenerated after this
    % "time in correct zone" metric was added.
    if isfield(cf, 'frac_correct_raw')
        flat_cf(idx).frac_correct_raw = cf.frac_correct_raw;
        flat_cf(idx).frac_correct_int = cf.frac_correct_int;
        flat_cf(idx).frac_correct_n_trials = cf.frac_correct_n_trials;
    end
end

% ── Flatten hybrid_advantage_probs per-file records, tagged with subject ───
% cluster_* fields are only present in hybrid_advantage_probs_summary.mat
% files regenerated after the cluster-based permutation test was added --
% guarded with isfield so older summaries (from before that change) don't
% crash the loader, they just don't contribute to the group cluster test.
flat_probs = struct('subject', {}, 'fus_adv_mean', {}, 'n_rescued_fr', {}, 'n_hurt_fr', {}, ...
                     'n_both_ok_fr', {}, 'n_both_bad_fr', {}, ...
                     'n_rescued_fr_all', {}, 'n_hurt_fr_all', {}, 'n_both_ok_fr_all', {}, 'n_both_bad_fr_all', {}, ...
                     'n_rescued_fr_post', {}, 'n_hurt_fr_post', {}, 'n_both_ok_fr_post', {}, 'n_both_bad_fr_post', {}, ...
                     'n_trials_post', {}, ...
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
        flat_probs(idx).n_both_ok_fr  = ps(k).n_both_ok_fr;
        flat_probs(idx).n_both_bad_fr = ps(k).n_both_bad_fr;
        flat_probs(idx).rescue_delta = ps(k).rescue_delta;
        flat_probs(idx).cost_delta   = ps(k).cost_delta;
        % _all/_post window variants only present in summaries regenerated
        % after the whole-trial/post-cvsa_influence breakdown was added --
        % same isfield-guarded backward compatibility as the cluster fields.
        if isfield(ps, 'n_rescued_fr_all')
            flat_probs(idx).n_rescued_fr_all  = ps(k).n_rescued_fr_all;
            flat_probs(idx).n_hurt_fr_all     = ps(k).n_hurt_fr_all;
            flat_probs(idx).n_both_ok_fr_all  = ps(k).n_both_ok_fr_all;
            flat_probs(idx).n_both_bad_fr_all = ps(k).n_both_bad_fr_all;
            flat_probs(idx).n_rescued_fr_post  = ps(k).n_rescued_fr_post;
            flat_probs(idx).n_hurt_fr_post     = ps(k).n_hurt_fr_post;
            flat_probs(idx).n_both_ok_fr_post  = ps(k).n_both_ok_fr_post;
            flat_probs(idx).n_both_bad_fr_post = ps(k).n_both_bad_fr_post;
            flat_probs(idx).n_trials_post      = ps(k).n_trials_post;
        end
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

% ── Flatten topo_erders per-band records, tagged with subject + SESSION
%    TYPE (evaluation/calibration) -- topo_erders.m is run once on the
%    evaluation/ folder and once on the sibling calibration/ folder (per
%    run_subject_analysis.m), so a subject can have topo_erders_summary.mat
%    from BOTH; without this tag they would get silently pooled together
%    into one average, mixing two different recording contexts (same bug
%    class already fixed in group_topo_erders.m). Everything that consumes
%    flat_erd below (Fig 12, r_erd_vs_real) is restricted to 'evaluation'
%    to match what "real Hybrid accuracy" itself is measured on; 'calibration'
%    entries are used separately, only by the calibration-baseline
%    correlate near the end of this function. ──────────────────────────────
flat_erd = struct('subject', {}, 'session', {}, 'paradigm', {}, 'band_origin', {}, 'discrimination', {}, 'csp_r', {});
for i = 1:numel(erd_files)
    fpath = fullfile(erd_files(i).folder, erd_files(i).name);
    subj  = subject_of(fpath, root_dir);
    s = load(fpath); es = s.erd_summary;
    for k = 1:numel(es)
        idx = numel(flat_erd) + 1;
        flat_erd(idx).subject        = subj;
        flat_erd(idx).session        = session_of(fpath);
        flat_erd(idx).paradigm       = es(k).paradigm;
        flat_erd(idx).band_origin    = es(k).band_origin;
        flat_erd(idx).discrimination = es(k).discrimination;
        flat_erd(idx).csp_r          = es(k).csp_r;
    end
end
flat_erd_eval = flat_erd(strcmp({flat_erd.session}, 'evaluation'));

% ── Flatten ROC per-paradigm records, tagged with subject (raw pooled
%    scores/labels, so they can be RE-POOLED per subject below) ────────────
flat_roc = struct('subject', {}, 'paradigm', {}, 'scores', {}, 'labels', {}, ...
                   'scores_win', {}, 'labels_win', {}, 'window_s', {});
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
        % scores_win/labels_win/window_s (frames restricted to the first
        % cvsa_influence seconds of CF -- the only window where a fused
        % classifier can differ from pure MI/CVSA) only present in
        % roc_summary.mat regenerated after this windowed AUC was added.
        if isfield(rs.roc_by_group(k), 'scores_win')
            flat_roc(idx).scores_win = rs.roc_by_group(k).scores_win;
            flat_roc(idx).labels_win = rs.roc_by_group(k).labels_win;
            flat_roc(idx).window_s   = rs.roc_by_group(k).window_s;
        end
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
real_T_all   = nan(n_subj, 3);   % mean time-to-ANY-outcome (HIT+MISS+TIMEOUT pooled) -- Wolpaw ITR's T
% Raw-signal-quality (dead/noisy channel during CF, see main_session_overview.m)
% pooled from raw counts (not averaged rates) -- same "pool then divide"
% convention as everything else in this file. NaN entries (no matched
% outcome for that trial) are excluded from both numerator and denominator.
sigqc_n_flagged = zeros(n_subj, 3);
sigqc_n_total   = zeros(n_subj, 3);
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
        all_t = [tth_all, tmiss_all, toT_all];
        if ~isempty(all_t), real_T_all(si,pi) = mean(all_t); end

        dead_all  = [flat_real(mask).sigqc_n_dead];
        noisy_all = [flat_real(mask).sigqc_n_noisy];
        valid_t   = ~isnan(dead_all);
        flagged_t = valid_t & ((dead_all > 0) | (noisy_all > 0));
        sigqc_n_total(si,pi)   = sum(valid_t);
        sigqc_n_flagged(si,pi) = sum(flagged_t);
    end
end
sigqc_rate = sigqc_n_flagged ./ max(1, sigqc_n_total);   % [n_subj x 3] MI, CVSA, Hybrid
grand_sigqc_n_flagged = sum(sigqc_n_flagged, 1);
grand_sigqc_n_total   = sum(sigqc_n_total, 1);
grand_sigqc_rate      = grand_sigqc_n_flagged ./ max(1, grand_sigqc_n_total);   % [1x3], trial-pooled across subjects
% Overall per-subject rate, pooled ACROSS paradigms too (by raw counts) --
% used for the well/bad comparison and the cross-subject correlate below,
% since "how clean was this subject's EEG overall" is a subject-level
% question, not tied to one paradigm's session.
sigqc_rate_overall = sum(sigqc_n_flagged, 2) ./ max(1, sum(sigqc_n_total, 2));   % [n_subj x 1]

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
grand_T_all   = nan(1, 3);   % mean time-to-ANY-outcome, trial-pooled -- Wolpaw ITR's T
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
    all_t = [tth_all, tmiss_all, toT_all];
    if ~isempty(all_t), grand_T_all(pi) = mean(all_t); end
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

% ── Per-subject counterfactual (SIMULATED) times, stream order matches
%    cf_labels = {'Hybrid','MI-only','CVSA-only'}. Each metric weighted by
%    ITS OWN outcome count (TTH by n_hit, T-miss by n_miss, T-timeout by
%    n_to) when a subject contributes more than one hybrid_advantage_integ
%    run -- exact pooling, not a plain average of session means. ───────────
cf_tth_subj   = nan(n_subj, 3);
cf_tmiss_subj = nan(n_subj, 3);
cf_tto_subj   = nan(n_subj, 3);
has_cftime = ~cellfun(@isempty, {flat_cf.t_miss_mean});
if any(has_cftime)
    rows_ct_all = flat_cf(has_cftime);
    for si = 1:n_subj
        mask = strcmp({rows_ct_all.subject}, subjects{si});
        if ~any(mask), continue; end
        rows_ct = rows_ct_all(mask);
        cf_tth_subj(si,:)   = weighted_mean_cols(vertcat(rows_ct.t_hit_mean),  vertcat(rows_ct.n_hit));
        cf_tmiss_subj(si,:) = weighted_mean_cols(vertcat(rows_ct.t_miss_mean), vertcat(rows_ct.n_miss));
        cf_tto_subj(si,:)   = weighted_mean_cols(vertcat(rows_ct.t_to_mean),   vertcat(rows_ct.n_to));
    end
end
% ── Per-subject "time in correct zone" (fraction of CF frames, ALL trials,
%    with target-class signal > 0.5 -- raw pre-integrator / int
%    post-integrator), weighted by each file's OWN trial count when a
%    subject contributes more than one hybrid_advantage_integ run. Skipped
%    if no subject's counterfactual_summary.mat has been regenerated with
%    these fields yet. ────────────────────────────────────────────────────
fraccorr_raw_subj = nan(n_subj, 3);
fraccorr_int_subj = nan(n_subj, 3);
has_fraccorr = ~cellfun(@isempty, {flat_cf.frac_correct_raw});
if any(has_fraccorr)
    rows_fc_all = flat_cf(has_fraccorr);
    for si = 1:n_subj
        mask = strcmp({rows_fc_all.subject}, subjects{si});
        if ~any(mask), continue; end
        rows_fc = rows_fc_all(mask);
        w = vertcat(rows_fc.frac_correct_n_trials) * [1 1 1];   % same weight, all 3 streams
        fraccorr_raw_subj(si,:) = weighted_mean_cols(vertcat(rows_fc.frac_correct_raw), w);
        fraccorr_int_subj(si,:) = weighted_mean_cols(vertcat(rows_fc.frac_correct_int), w);
    end
end

% GRAND: pool the (already-exact) per-subject means, weighted by each
% subject's own outcome counts (cf_n_hit/cf_n_miss/cf_n_to, computed above).
grand_cf_tth   = weighted_mean_cols(cf_tth_subj,   cf_n_hit);
grand_cf_tmiss = weighted_mean_cols(cf_tmiss_subj, cf_n_miss);
grand_cf_tto   = weighted_mean_cols(cf_tto_subj,   cf_n_to);

% ── Mean time-to-ANY-outcome (Wolpaw ITR's T): pool TTH/T-miss/T-timeout,
%    weighted by each category's own trial count -- a TIMEOUT-heavy stream
%    is penalised both by lower accuracy AND by a longer T here. ───────────
cf_T_all       = pool3_weighted(cf_tth_subj,   cf_n_hit,   cf_tmiss_subj, cf_n_miss, cf_tto_subj, cf_n_to);
grand_cf_T_all = pool3_weighted(grand_cf_tth,  grand_cf_n_hit, grand_cf_tmiss, grand_cf_n_miss, grand_cf_tto, grand_cf_n_to);

% ── Information Transfer Rate (ITR, Wolpaw et al.): bits/trial from
%    accuracy, bits/min by dividing by T above. N_CLS=2 for every
%    paradigm/stream in this dataset (every recorded class pair seen is
%    2-class) -- an assumption, not read per-subject, since neither
%    session_summary.mat nor counterfactual_summary.mat currently saves the
%    class count.
%
%    IMPORTANT: bits/trial is computed from DECIDED-TRIALS-ONLY accuracy
%    (acc_dec), scaled by the decided-trial rate (1 - timeout_rate) --
%    NOT from the raw TIMEOUT=fail accuracy fed directly into Wolpaw's
%    formula. Reason: for N_CLS=2, Wolpaw's bits/trial is SYMMETRIC around
%    50% -- an always-wrong classifier (P=0) scores the same 1 bit/trial as
%    an always-right one (P=1), since in principle you could just invert
%    its answers. That's fine for a genuine forced-choice classifier, but
%    real_acc/cf_acc fold TIMEOUT into "wrong" too -- so a stream that NEVER
%    reaches ANY decision (100% timeout, e.g. a weak CVSA-only session)
%    would score exactly like a stream that decides confidently every
%    trial and is invertibly wrong, i.e. 1 bit/trial, when it should score
%    0 (a system that never answers carries no information at all). Scaling
%    the decided-trials bits/trial by the decided-trial rate fixes this: a
%    100%-timeout stream now correctly gets 0 bits/trial regardless of what
%    its (undefined) decided-trial accuracy would have been. ───────────────
N_CLS = 2;
real_itr_bpt = nan(n_subj, 3); real_itr_bpm = nan(n_subj, 3);
cf_itr_bpt   = nan(n_subj, 3); cf_itr_bpm   = nan(n_subj, 3);
grand_real_itr_bpt = nan(1, 3); grand_real_itr_bpm = nan(1, 3);
grand_cf_itr_bpt   = nan(1, 3); grand_cf_itr_bpm   = nan(1, 3);
for pi = 1:3
    for si = 1:n_subj
        real_itr_bpt(si,pi) = decided_rate_itr_local(real_acc_dec(si,pi), real_to_rate(si,pi), N_CLS);
        if ~isnan(real_T_all(si,pi)) && real_T_all(si,pi) > 0
            real_itr_bpm(si,pi) = real_itr_bpt(si,pi) * 60 / real_T_all(si,pi);
        end
        cf_itr_bpt(si,pi) = decided_rate_itr_local(cf_acc_dec(si,pi), cf_to_rate(si,pi), N_CLS);
        if ~isnan(cf_T_all(si,pi)) && cf_T_all(si,pi) > 0
            cf_itr_bpm(si,pi) = cf_itr_bpt(si,pi) * 60 / cf_T_all(si,pi);
        end
    end
    grand_real_itr_bpt(pi) = decided_rate_itr_local(grand_acc_dec(pi), grand_to_rate(pi), N_CLS);
    if ~isnan(grand_T_all(pi)) && grand_T_all(pi) > 0
        grand_real_itr_bpm(pi) = grand_real_itr_bpt(pi) * 60 / grand_T_all(pi);
    end
    grand_cf_itr_bpt(pi) = decided_rate_itr_local(grand_cf_acc_dec(pi), grand_cf_to_rate(pi), N_CLS);
    if ~isnan(grand_cf_T_all(pi)) && grand_cf_T_all(pi) > 0
        grand_cf_itr_bpm(pi) = grand_cf_itr_bpt(pi) * 60 / grand_cf_T_all(pi);
    end
end

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
buf_adv_cvs_win_subj      = nan(n_subj, 1);   % same, Hybrid vs CVSA-only (no Stouffer meta-analysis needed, just the raw mean)
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
        if isfield(rows_bw, 'buf_adv_cvs_win_mean')
            buf_adv_cvs_win_subj(si) = mean([rows_bw.buf_adv_cvs_win_mean], 'omitnan');
        end
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

% ── Frame-level rescue/hurt/both-ok/both-bad BREAKDOWN, as fractions of
%    classified frames, over THREE window definitions -- "on average, in how
%    many samples does the CVSA help vs hurt vs not matter (both right/both
%    wrong)?" Per-subject fraction (that subject's own files summed first)
%    + GRAND (frame-pooled across the whole cohort, the population estimate,
%    matching every other GRAND in this script). Categories match
%    main_hybrid_advantage_probs.m's Panel(1,1): rescued/hurt/both-ok/both-bad. ─
has_rescue_all = ~cellfun(@isempty, {flat_probs.n_rescued_fr_all});
rescue_frac_inf  = nan(n_subj, 4);
rescue_frac_all  = nan(n_subj, 4);
rescue_frac_post = nan(n_subj, 4);
n_trials_post_subj  = nan(n_subj, 1);
n_trials_probs_subj = nan(n_subj, 1);
for si = 1:n_subj
    mask = strcmp({flat_probs.subject}, subjects{si});
    if ~any(mask), continue; end
    rows = flat_probs(mask);
    c_inf = [sum([rows.n_rescued_fr]), sum([rows.n_hurt_fr]), sum([rows.n_both_ok_fr]), sum([rows.n_both_bad_fr])];
    if sum(c_inf) > 0, rescue_frac_inf(si,:) = c_inf / sum(c_inf); end
    mask_all = mask & has_rescue_all;
    if any(mask_all)
        rows_a = flat_probs(mask_all);
        c_all  = [sum([rows_a.n_rescued_fr_all]),  sum([rows_a.n_hurt_fr_all]),  sum([rows_a.n_both_ok_fr_all]),  sum([rows_a.n_both_bad_fr_all])];
        c_post = [sum([rows_a.n_rescued_fr_post]), sum([rows_a.n_hurt_fr_post]), sum([rows_a.n_both_ok_fr_post]), sum([rows_a.n_both_bad_fr_post])];
        if sum(c_all)  > 0, rescue_frac_all(si,:)  = c_all  / sum(c_all);  end
        if sum(c_post) > 0, rescue_frac_post(si,:) = c_post / sum(c_post); end
        n_trials_post_subj(si)  = sum([rows_a.n_trials_post]);
        n_trials_probs_subj(si) = sum([rows_a.fus_adv_n_trials]);   % total trials, for the "X/Y trials contribute" caveat
    end
end
grand_rescue_frac_inf  = pooled_frac_4({flat_probs.n_rescued_fr},  {flat_probs.n_hurt_fr},  {flat_probs.n_both_ok_fr},  {flat_probs.n_both_bad_fr});
grand_rescue_frac_all  = pooled_frac_4({flat_probs.n_rescued_fr_all},  {flat_probs.n_hurt_fr_all},  {flat_probs.n_both_ok_fr_all},  {flat_probs.n_both_bad_fr_all});
grand_rescue_frac_post = pooled_frac_4({flat_probs.n_rescued_fr_post}, {flat_probs.n_hurt_fr_post}, {flat_probs.n_both_ok_fr_post}, {flat_probs.n_both_bad_fr_post});
grand_n_trials_post  = sum([flat_probs.n_trials_post],  'omitnan');
grand_n_trials_probs = sum([flat_probs.fus_adv_n_trials], 'omitnan');

fprintf('\n══════ Frame-level CVSA effect breakdown: rescued / hurt / both-correct / both-wrong (n=%d subjects) ══════\n', sum(any(~isnan(rescue_frac_inf),2)));
fprintf('  -- within cvsa_influence window (all trials) --\n');
for si = 1:n_subj
    if all(isnan(rescue_frac_inf(si,:))), continue; end
    fprintf('  %-10s  rescued=%s  hurt=%s  both-ok=%s  both-wrong=%s\n', subjects{si}, ...
            fmt_pct(rescue_frac_inf(si,1)), fmt_pct(rescue_frac_inf(si,2)), fmt_pct(rescue_frac_inf(si,3)), fmt_pct(rescue_frac_inf(si,4)));
end
fprintf('  %-10s  rescued=%s  hurt=%s  both-ok=%s  both-wrong=%s  (GRAND, frame-pooled)\n', 'GRAND', ...
        fmt_pct(grand_rescue_frac_inf(1)), fmt_pct(grand_rescue_frac_inf(2)), fmt_pct(grand_rescue_frac_inf(3)), fmt_pct(grand_rescue_frac_inf(4)));
fprintf('  -- whole trial (no time restriction) --\n');
fprintf('  %-10s  rescued=%s  hurt=%s  both-ok=%s  both-wrong=%s  (GRAND, frame-pooled)\n', 'GRAND', ...
        fmt_pct(grand_rescue_frac_all(1)), fmt_pct(grand_rescue_frac_all(2)), fmt_pct(grand_rescue_frac_all(3)), fmt_pct(grand_rescue_frac_all(4)));
fprintf('  -- from cvsa_influence ONWARD (alpha~0 -- %d/%d trials contribute) --\n', grand_n_trials_post, grand_n_trials_probs);
fprintf('  %-10s  rescued=%s  hurt=%s  both-ok=%s  both-wrong=%s  (GRAND, frame-pooled)\n', 'GRAND', ...
        fmt_pct(grand_rescue_frac_post(1)), fmt_pct(grand_rescue_frac_post(2)), fmt_pct(grand_rescue_frac_post(3)), fmt_pct(grand_rescue_frac_post(4)));
fprintf('  (per-subject breakdown for whole-trial/post-cvsa_influence omitted from console for brevity -- see the figure and group_summary.mat)\n');
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

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
        mask = strcmp({flat_erd_eval.subject}, subjects{si}) & strcmp({flat_erd_eval.paradigm}, pars{pi});
        if ~any(mask), continue; end
        rows = flat_erd_eval(mask);
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
% Brier score (calibration: are the probabilities themselves trustworthy,
% not just their ranking?) and d-prime (classic signal-detection
% sensitivity index at the operating point P=0.5, same convention as the
% frame-accuracy figures elsewhere) -- computed on the SAME re-pooled
% per-subject (score, label) frames as the ROC/AUC above, so all three
% metrics describe the identical pooled data.
brier_subj  = nan(n_subj, 3);
dprime_subj = nan(n_subj, 3);
% Reliability diagram data: 10 equal-width probability bins: for each
% (subject, paradigm, bin), the mean predicted P(target) in that bin
% ("confidence") and the observed fraction of frames where target was
% really the label ("accuracy") -- perfect calibration = points on the
% diagonal. Averaged across subjects per paradigm below, same convention
% as the ROC macro-average.
CALIB_NBINS = 10;
calib_conf_subj = nan(n_subj, 3, CALIB_NBINS);
calib_obs_subj  = nan(n_subj, 3, CALIB_NBINS);
% Windowed (cvsa_influence-restricted) counterparts of AUC/Brier/d-prime --
% same re-pooled per-subject frames as above, but only frames within the
% first cvsa_influence seconds of each trial's CF, the ONLY window where a
% fused/Hybrid classifier can differ from pure MI (bayesian_fuse.m forces
% alpha=0 -- fused==MI exactly -- from cvsa_influence onward). Pooling ALL
% CF frames dilutes any early advantage with the later identical portion;
% this isolates it. Skipped (stays NaN) if no subject's roc_summary.mat has
% been regenerated with scores_win/labels_win yet.
has_roc_win = ~cellfun(@isempty, {flat_roc.scores_win});
roc_auc_win_subj = nan(n_subj, 3);
roc_tpr_win_subj_grid = nan(n_subj, 3, numel(ROC_FPR_GRID));
brier_win_subj   = nan(n_subj, 3);
dprime_win_subj  = nan(n_subj, 3);
window_s_subj    = nan(n_subj, 3);
calib_conf_win_subj = nan(n_subj, 3, CALIB_NBINS);
calib_obs_win_subj  = nan(n_subj, 3, CALIB_NBINS);
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

        [brier_subj(si,pi), dprime_subj(si,pi)] = brier_dprime_local(scores_all, labels_all);

        lbl = logical(labels_all);
        bin_idx = min(max(1, ceil(scores_all * CALIB_NBINS)), CALIB_NBINS);
        for b = 1:CALIB_NBINS
            in_b = bin_idx == b;
            if any(in_b)
                calib_conf_subj(si,pi,b) = mean(scores_all(in_b), 'omitnan');
                calib_obs_subj(si,pi,b)  = mean(double(lbl(in_b)), 'omitnan');
            end
        end

        mask_w = mask & has_roc_win;
        if any(mask_w)
            rows_w = flat_roc(mask_w);
            scores_w = double(vertcat(rows_w.scores_win));
            labels_w = vertcat(rows_w.labels_win);
            [fpr_w, tpr_w, ~, auc_w] = compute_roc_curve(scores_w, labels_w);
            roc_auc_win_subj(si,pi) = auc_w;
            [u_fpr_w, ~, ic_w] = unique(fpr_w);
            u_tpr_w = accumarray(ic_w, tpr_w, [], @max);
            roc_tpr_win_subj_grid(si,pi,:) = interp1(u_fpr_w, u_tpr_w, ROC_FPR_GRID, 'linear');
            [brier_win_subj(si,pi), dprime_win_subj(si,pi)] = brier_dprime_local(scores_w, labels_w);
            window_s_subj(si,pi) = mean([rows_w.window_s], 'omitnan');

            lbl_w = logical(labels_w);
            bin_idx_w = min(max(1, ceil(scores_w * CALIB_NBINS)), CALIB_NBINS);
            for b = 1:CALIB_NBINS
                in_b = bin_idx_w == b;
                if any(in_b)
                    calib_conf_win_subj(si,pi,b) = mean(scores_w(in_b), 'omitnan');
                    calib_obs_win_subj(si,pi,b)  = mean(double(lbl_w(in_b)), 'omitnan');
                end
            end
        end
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

% Windowed (cvsa_influence-restricted) macro-average ROC curve, same method
% as the ALL-frames version above -- feeds Fig 14a.
roc_mean_tpr_win = nan(3, numel(ROC_FPR_GRID));
roc_sem_tpr_win  = nan(3, numel(ROC_FPR_GRID));
for pi = 1:3
    curves_w = reshape(roc_tpr_win_subj_grid(:,pi,:), n_subj, numel(ROC_FPR_GRID));
    valid_rows_w = ~all(isnan(curves_w), 2);
    n_vw = sum(valid_rows_w);
    if n_vw > 0
        roc_mean_tpr_win(pi,:) = mean(curves_w(valid_rows_w,:), 1, 'omitnan');
        roc_sem_tpr_win(pi,:)  = std(curves_w(valid_rows_w,:), 0, 1, 'omitnan') / sqrt(max(1,n_vw));
    end
end

% Pairwise AUC comparison (Hybrid vs MI, Hybrid vs CVSA) -- p_auc_grp above
% only tests each paradigm against chance (0.5) separately; this is the
% "does Hybrid actually discriminate BETTER than either unimodal stream"
% test, subject is the unit, mirroring every other pairwise test here.
% Columns follow pars={'mi','cvsa','hybrid'}: col1=MI, col2=CVSA, col3=Hybrid.
d_auc_hyb_mi  = roc_auc_subj(:,3) - roc_auc_subj(:,1);
d_auc_hyb_cvs = roc_auc_subj(:,3) - roc_auc_subj(:,2);
d_auc_mi_cvs  = roc_auc_subj(:,1) - roc_auc_subj(:,2);   % MI - CVSA, same convention as d_real_mi_cvs
p_auc_hyb_mi  = sign_flip_test_local(d_auc_hyb_mi,  n_perm);
p_auc_hyb_cvs = sign_flip_test_local(d_auc_hyb_cvs, n_perm);
p_auc_mi_cvs  = sign_flip_test_local(d_auc_mi_cvs,  n_perm);
d_auc_hyb_mi_cohen  = cohen_d_local(d_auc_hyb_mi);
d_auc_hyb_cvs_cohen = cohen_d_local(d_auc_hyb_cvs);
d_auc_mi_cvs_cohen  = cohen_d_local(d_auc_mi_cvs);
fprintf('\n══════════════════ ROC/AUC pairwise comparison: which classifier discriminates better? (n=%d subjects) ══════════════════\n', n_subj);
fprintf('  Hybrid - MI   AUC delta : %+.3f  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(d_auc_hyb_mi,'omitnan'),  p_auc_hyb_mi,  stars(p_auc_hyb_mi),  d_auc_hyb_mi_cohen,  cohen_d_label(d_auc_hyb_mi_cohen));
fprintf('  Hybrid - CVSA AUC delta : %+.3f  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(d_auc_hyb_cvs,'omitnan'), p_auc_hyb_cvs, stars(p_auc_hyb_cvs), d_auc_hyb_cvs_cohen, cohen_d_label(d_auc_hyb_cvs_cohen));
fprintf('  MI - CVSA     AUC delta : %+.3f  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(d_auc_mi_cvs,'omitnan'),  p_auc_mi_cvs,  stars(p_auc_mi_cvs),  d_auc_mi_cvs_cohen,  cohen_d_label(d_auc_mi_cvs_cohen));
fprintf('  (paired two-sided sign-flip permutation on per-subject AUC delta, %d perms; subject is the statistical unit; AUC itself: trapezoidal area under the non-parametric ROC on that subject''s pooled per-frame scores, equivalent to Mann-Whitney U/(n_pos*n_neg) -- see compute_roc_curve.m)\n', n_perm);
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

% ── Group-level pairwise comparisons on Brier score (calibration -- lower is
%    better) and d-prime (signal-detection sensitivity at P=0.5 -- higher is
%    better), same per-subject pooled scores as AUC above, same pairwise
%    sign-flip + Cohen's d pattern. ─────────────────────────────────────────
d_brier_hyb_mi  = brier_subj(:,3) - brier_subj(:,1);   % Hybrid - MI (lower Brier = better, so Hybrid-MI<0 favors Hybrid)
d_brier_hyb_cvs = brier_subj(:,3) - brier_subj(:,2);
d_brier_mi_cvs  = brier_subj(:,1) - brier_subj(:,2);
p_brier_hyb_mi  = sign_flip_test_local(d_brier_hyb_mi,  n_perm);
p_brier_hyb_cvs = sign_flip_test_local(d_brier_hyb_cvs, n_perm);
p_brier_mi_cvs  = sign_flip_test_local(d_brier_mi_cvs,  n_perm);
d_brier_hyb_mi_cohen  = cohen_d_local(d_brier_hyb_mi);
d_brier_hyb_cvs_cohen = cohen_d_local(d_brier_hyb_cvs);
d_brier_mi_cvs_cohen  = cohen_d_local(d_brier_mi_cvs);

d_dprime_hyb_mi  = dprime_subj(:,3) - dprime_subj(:,1);
d_dprime_hyb_cvs = dprime_subj(:,3) - dprime_subj(:,2);
d_dprime_mi_cvs  = dprime_subj(:,1) - dprime_subj(:,2);
p_dprime_hyb_mi  = sign_flip_test_local(d_dprime_hyb_mi,  n_perm);
p_dprime_hyb_cvs = sign_flip_test_local(d_dprime_hyb_cvs, n_perm);
p_dprime_mi_cvs  = sign_flip_test_local(d_dprime_mi_cvs,  n_perm);
d_dprime_hyb_mi_cohen  = cohen_d_local(d_dprime_hyb_mi);
d_dprime_hyb_cvs_cohen = cohen_d_local(d_dprime_hyb_cvs);
d_dprime_mi_cvs_cohen  = cohen_d_local(d_dprime_mi_cvs);

fprintf('\n══════════════════ Classifier calibration (Brier) & sensitivity (d-prime), n=%d subjects ══════════════════\n', n_subj);
fprintf('  %-8s  %10s  %10s\n', 'stream', 'Brier (lo better)', 'd-prime (hi better)');
for pi = 1:3
    fprintf('  %-8s  %18.3f  %20.2f\n', par_labels{pi}, mean(brier_subj(:,pi),'omitnan'), mean(dprime_subj(:,pi),'omitnan'));
end
fprintf('  Hybrid - MI   : Brier delta=%+.3f p=%.4f %s d=%+.2f (%s)  |  d-prime delta=%+.2f p=%.4f %s d=%+.2f (%s)\n', ...
        mean(d_brier_hyb_mi,'omitnan'), p_brier_hyb_mi, stars(p_brier_hyb_mi), d_brier_hyb_mi_cohen, cohen_d_label(d_brier_hyb_mi_cohen), ...
        mean(d_dprime_hyb_mi,'omitnan'), p_dprime_hyb_mi, stars(p_dprime_hyb_mi), d_dprime_hyb_mi_cohen, cohen_d_label(d_dprime_hyb_mi_cohen));
fprintf('  Hybrid - CVSA : Brier delta=%+.3f p=%.4f %s d=%+.2f (%s)  |  d-prime delta=%+.2f p=%.4f %s d=%+.2f (%s)\n', ...
        mean(d_brier_hyb_cvs,'omitnan'), p_brier_hyb_cvs, stars(p_brier_hyb_cvs), d_brier_hyb_cvs_cohen, cohen_d_label(d_brier_hyb_cvs_cohen), ...
        mean(d_dprime_hyb_cvs,'omitnan'), p_dprime_hyb_cvs, stars(p_dprime_hyb_cvs), d_dprime_hyb_cvs_cohen, cohen_d_label(d_dprime_hyb_cvs_cohen));
fprintf('  MI - CVSA     : Brier delta=%+.3f p=%.4f %s d=%+.2f (%s)  |  d-prime delta=%+.2f p=%.4f %s d=%+.2f (%s)\n', ...
        mean(d_brier_mi_cvs,'omitnan'), p_brier_mi_cvs, stars(p_brier_mi_cvs), d_brier_mi_cvs_cohen, cohen_d_label(d_brier_mi_cvs_cohen), ...
        mean(d_dprime_mi_cvs,'omitnan'), p_dprime_mi_cvs, stars(p_dprime_mi_cvs), d_dprime_mi_cvs_cohen, cohen_d_label(d_dprime_mi_cvs_cohen));
fprintf('  (Brier = mean((P(target)-label)^2), pooled per-subject frames, LOWER=better calibrated; d-prime = z(hit rate)-z(false alarm rate) at P=0.5, HIGHER=better separated; paired sign-flip, %d perms, subject=unit)\n', n_perm);
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

% ── WINDOWED (cvsa_influence-restricted) pairwise comparisons: AUC, Brier,
%    d-prime -- same 3 pairwise tests as the ALL-frames versions above, but
%    computed only on frames within the first cvsa_influence seconds of CF,
%    where a fused classifier can actually differ from pure MI/CVSA. Skipped
%    if no subject's roc_summary.mat has scores_win yet. ───────────────────
if any(has_roc_win)
    d_aucw_hyb_mi  = roc_auc_win_subj(:,3) - roc_auc_win_subj(:,1);
    d_aucw_hyb_cvs = roc_auc_win_subj(:,3) - roc_auc_win_subj(:,2);
    d_aucw_mi_cvs  = roc_auc_win_subj(:,1) - roc_auc_win_subj(:,2);
    p_aucw_hyb_mi  = sign_flip_test_local(d_aucw_hyb_mi,  n_perm);
    p_aucw_hyb_cvs = sign_flip_test_local(d_aucw_hyb_cvs, n_perm);
    p_aucw_mi_cvs  = sign_flip_test_local(d_aucw_mi_cvs,  n_perm);
    d_aucw_hyb_mi_cohen  = cohen_d_local(d_aucw_hyb_mi);
    d_aucw_hyb_cvs_cohen = cohen_d_local(d_aucw_hyb_cvs);
    d_aucw_mi_cvs_cohen  = cohen_d_local(d_aucw_mi_cvs);

    d_brierw_hyb_mi  = brier_win_subj(:,3) - brier_win_subj(:,1);
    d_brierw_hyb_cvs = brier_win_subj(:,3) - brier_win_subj(:,2);
    d_brierw_mi_cvs  = brier_win_subj(:,1) - brier_win_subj(:,2);
    p_brierw_hyb_mi  = sign_flip_test_local(d_brierw_hyb_mi,  n_perm);
    p_brierw_hyb_cvs = sign_flip_test_local(d_brierw_hyb_cvs, n_perm);
    p_brierw_mi_cvs  = sign_flip_test_local(d_brierw_mi_cvs,  n_perm);
    d_brierw_hyb_mi_cohen  = cohen_d_local(d_brierw_hyb_mi);
    d_brierw_hyb_cvs_cohen = cohen_d_local(d_brierw_hyb_cvs);
    d_brierw_mi_cvs_cohen  = cohen_d_local(d_brierw_mi_cvs);

    d_dprimew_hyb_mi  = dprime_win_subj(:,3) - dprime_win_subj(:,1);
    d_dprimew_hyb_cvs = dprime_win_subj(:,3) - dprime_win_subj(:,2);
    d_dprimew_mi_cvs  = dprime_win_subj(:,1) - dprime_win_subj(:,2);
    p_dprimew_hyb_mi  = sign_flip_test_local(d_dprimew_hyb_mi,  n_perm);
    p_dprimew_hyb_cvs = sign_flip_test_local(d_dprimew_hyb_cvs, n_perm);
    p_dprimew_mi_cvs  = sign_flip_test_local(d_dprimew_mi_cvs,  n_perm);
    d_dprimew_hyb_mi_cohen  = cohen_d_local(d_dprimew_hyb_mi);
    d_dprimew_hyb_cvs_cohen = cohen_d_local(d_dprimew_hyb_cvs);
    d_dprimew_mi_cvs_cohen  = cohen_d_local(d_dprimew_mi_cvs);

    mean_window_s = mean(window_s_subj(:), 'omitnan');
    fprintf('\n══════════════════ WINDOWED (first %.1fs of CF only) classifier comparison, n=%d subjects ══════════════════\n', mean_window_s, n_subj);
    fprintf('  %-8s  %10s  %18s  %20s\n', 'stream', 'AUC', 'Brier (lo better)', 'd-prime (hi better)');
    for pi = 1:3
        fprintf('  %-8s  %10.3f  %18.3f  %20.2f\n', par_labels{pi}, mean(roc_auc_win_subj(:,pi),'omitnan'), mean(brier_win_subj(:,pi),'omitnan'), mean(dprime_win_subj(:,pi),'omitnan'));
    end
    fprintf('  Hybrid - MI   : AUC delta=%+.3f p=%.4f %s d=%+.2f  |  Brier delta=%+.3f p=%.4f %s d=%+.2f  |  d-prime delta=%+.2f p=%.4f %s d=%+.2f\n', ...
            mean(d_aucw_hyb_mi,'omitnan'), p_aucw_hyb_mi, stars(p_aucw_hyb_mi), d_aucw_hyb_mi_cohen, ...
            mean(d_brierw_hyb_mi,'omitnan'), p_brierw_hyb_mi, stars(p_brierw_hyb_mi), d_brierw_hyb_mi_cohen, ...
            mean(d_dprimew_hyb_mi,'omitnan'), p_dprimew_hyb_mi, stars(p_dprimew_hyb_mi), d_dprimew_hyb_mi_cohen);
    fprintf('  Hybrid - CVSA : AUC delta=%+.3f p=%.4f %s d=%+.2f  |  Brier delta=%+.3f p=%.4f %s d=%+.2f  |  d-prime delta=%+.2f p=%.4f %s d=%+.2f\n', ...
            mean(d_aucw_hyb_cvs,'omitnan'), p_aucw_hyb_cvs, stars(p_aucw_hyb_cvs), d_aucw_hyb_cvs_cohen, ...
            mean(d_brierw_hyb_cvs,'omitnan'), p_brierw_hyb_cvs, stars(p_brierw_hyb_cvs), d_brierw_hyb_cvs_cohen, ...
            mean(d_dprimew_hyb_cvs,'omitnan'), p_dprimew_hyb_cvs, stars(p_dprimew_hyb_cvs), d_dprimew_hyb_cvs_cohen);
    fprintf('  MI - CVSA     : AUC delta=%+.3f p=%.4f %s d=%+.2f  |  Brier delta=%+.3f p=%.4f %s d=%+.2f  |  d-prime delta=%+.2f p=%.4f %s d=%+.2f\n', ...
            mean(d_aucw_mi_cvs,'omitnan'), p_aucw_mi_cvs, stars(p_aucw_mi_cvs), d_aucw_mi_cvs_cohen, ...
            mean(d_brierw_mi_cvs,'omitnan'), p_brierw_mi_cvs, stars(p_brierw_mi_cvs), d_brierw_mi_cvs_cohen, ...
            mean(d_dprimew_mi_cvs,'omitnan'), p_dprimew_mi_cvs, stars(p_dprimew_mi_cvs), d_dprimew_mi_cvs_cohen);
    fprintf('  (paired two-sided sign-flip permutation on per-subject delta, %d perms, subject=unit; frames restricted to time-since-CF-onset < cvsa_influence, ~%.1fs)\n', n_perm, mean_window_s);
    fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');
end

% ── Group-level significance: TTH and ITR(bits/min), both real-session and
%    counterfactual, same paired sign-flip + Cohen's d pattern as accuracy
%    below (subject is the unit). TTH delta is negative when Hybrid is
%    FASTER (time, not a rate), unlike every other delta in this file. ─────
d_tth_hyb_mi    = real_tth(:,3) - real_tth(:,1);   % Hybrid - MI (real)
d_tth_hyb_cvs   = real_tth(:,3) - real_tth(:,2);   % Hybrid - CVSA (real)
d_cftth_hyb_mi  = cf_tth_subj(:,1) - cf_tth_subj(:,2);   % Hybrid - MI-only (counterfactual)
d_cftth_hyb_cvs = cf_tth_subj(:,1) - cf_tth_subj(:,3);   % Hybrid - CVSA-only (counterfactual)
p_tth_hyb_mi    = sign_flip_test_local(d_tth_hyb_mi,    n_perm);
p_tth_hyb_cvs   = sign_flip_test_local(d_tth_hyb_cvs,   n_perm);
p_cftth_hyb_mi  = sign_flip_test_local(d_cftth_hyb_mi,  n_perm);
p_cftth_hyb_cvs = sign_flip_test_local(d_cftth_hyb_cvs, n_perm);
d_tth_hyb_mi_cohen    = cohen_d_local(d_tth_hyb_mi);
d_tth_hyb_cvs_cohen   = cohen_d_local(d_tth_hyb_cvs);
d_cftth_hyb_mi_cohen  = cohen_d_local(d_cftth_hyb_mi);
d_cftth_hyb_cvs_cohen = cohen_d_local(d_cftth_hyb_cvs);

d_itr_hyb_mi    = real_itr_bpm(:,3) - real_itr_bpm(:,1);   % Hybrid - MI (real), bits/min
d_itr_hyb_cvs   = real_itr_bpm(:,3) - real_itr_bpm(:,2);   % Hybrid - CVSA (real)
d_cfitr_hyb_mi  = cf_itr_bpm(:,1) - cf_itr_bpm(:,2);        % Hybrid - MI-only (counterfactual)
d_cfitr_hyb_cvs = cf_itr_bpm(:,1) - cf_itr_bpm(:,3);        % Hybrid - CVSA-only (counterfactual)
p_itr_hyb_mi    = sign_flip_test_local(d_itr_hyb_mi,    n_perm);
p_itr_hyb_cvs   = sign_flip_test_local(d_itr_hyb_cvs,   n_perm);
p_cfitr_hyb_mi  = sign_flip_test_local(d_cfitr_hyb_mi,  n_perm);
p_cfitr_hyb_cvs = sign_flip_test_local(d_cfitr_hyb_cvs, n_perm);
d_itr_hyb_mi_cohen    = cohen_d_local(d_itr_hyb_mi);
d_itr_hyb_cvs_cohen   = cohen_d_local(d_itr_hyb_cvs);
d_cfitr_hyb_mi_cohen  = cohen_d_local(d_cfitr_hyb_mi);
d_cfitr_hyb_cvs_cohen = cohen_d_local(d_cfitr_hyb_cvs);

fprintf('\n══════════════════ Group-level TTH & ITR (n=%d subjects) ══════════════════\n', n_subj);
fprintf('  -- Real-session TTH (s) --\n');
fprintf('  Hybrid - MI    : mean delta=%+.2fs  p=%.4f  %s  d=%+.2f (%s)  (negative = Hybrid faster)\n', ...
        mean(d_tth_hyb_mi,'omitnan'),  p_tth_hyb_mi,  stars(p_tth_hyb_mi),  d_tth_hyb_mi_cohen,  cohen_d_label(d_tth_hyb_mi_cohen));
fprintf('  Hybrid - CVSA  : mean delta=%+.2fs  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(d_tth_hyb_cvs,'omitnan'), p_tth_hyb_cvs, stars(p_tth_hyb_cvs), d_tth_hyb_cvs_cohen, cohen_d_label(d_tth_hyb_cvs_cohen));
fprintf('  -- Counterfactual (simulated) TTH (s) --\n');
fprintf('  Hybrid - MI-only   : mean delta=%+.2fs  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(d_cftth_hyb_mi,'omitnan'),  p_cftth_hyb_mi,  stars(p_cftth_hyb_mi),  d_cftth_hyb_mi_cohen,  cohen_d_label(d_cftth_hyb_mi_cohen));
fprintf('  Hybrid - CVSA-only : mean delta=%+.2fs  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(d_cftth_hyb_cvs,'omitnan'), p_cftth_hyb_cvs, stars(p_cftth_hyb_cvs), d_cftth_hyb_cvs_cohen, cohen_d_label(d_cftth_hyb_cvs_cohen));
fprintf('  -- Real-session ITR (bits/min) --\n');
fprintf('  Hybrid - MI    : mean delta=%+.2f  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(d_itr_hyb_mi,'omitnan'),  p_itr_hyb_mi,  stars(p_itr_hyb_mi),  d_itr_hyb_mi_cohen,  cohen_d_label(d_itr_hyb_mi_cohen));
fprintf('  Hybrid - CVSA  : mean delta=%+.2f  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(d_itr_hyb_cvs,'omitnan'), p_itr_hyb_cvs, stars(p_itr_hyb_cvs), d_itr_hyb_cvs_cohen, cohen_d_label(d_itr_hyb_cvs_cohen));
fprintf('  -- Counterfactual (simulated) ITR (bits/min) --\n');
fprintf('  Hybrid - MI-only   : mean delta=%+.2f  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(d_cfitr_hyb_mi,'omitnan'),  p_cfitr_hyb_mi,  stars(p_cfitr_hyb_mi),  d_cfitr_hyb_mi_cohen,  cohen_d_label(d_cfitr_hyb_mi_cohen));
fprintf('  Hybrid - CVSA-only : mean delta=%+.2f  p=%.4f  %s  d=%+.2f (%s)\n', ...
        mean(d_cfitr_hyb_cvs,'omitnan'), p_cfitr_hyb_cvs, stars(p_cfitr_hyb_cvs), d_cfitr_hyb_cvs_cohen, cohen_d_label(d_cfitr_hyb_cvs_cohen));
fprintf('  (two-sided sign-flip permutation, %d perms; subject is the statistical unit; N_CLS=%d assumed for ITR)\n', n_perm, N_CLS);
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

% ── Group-level "time in correct zone": Hybrid vs MI-only / Hybrid vs
%    CVSA-only, on the per-subject mean fraction of CF frames with the
%    target-class signal > 0.5, both raw (pre-integrator) and int
%    (post-integrator). Paired sign-flip test + Cohen's d, subject is the
%    unit (each subject's own per-trial paired test already lives in
%    main_hybrid_advantage_integ.m's console output). ─────────────────────
d_fcr_mi  = fraccorr_raw_subj(:,1) - fraccorr_raw_subj(:,2);   % Hybrid - MI-only, raw
d_fcr_cvs = fraccorr_raw_subj(:,1) - fraccorr_raw_subj(:,3);   % Hybrid - CVSA-only, raw
d_fci_mi  = fraccorr_int_subj(:,1) - fraccorr_int_subj(:,2);   % Hybrid - MI-only, integrated
d_fci_cvs = fraccorr_int_subj(:,1) - fraccorr_int_subj(:,3);   % Hybrid - CVSA-only, integrated
p_fcr_mi  = sign_flip_test_local(d_fcr_mi,  n_perm); p_fcr_cvs = sign_flip_test_local(d_fcr_cvs, n_perm);
p_fci_mi  = sign_flip_test_local(d_fci_mi,  n_perm); p_fci_cvs = sign_flip_test_local(d_fci_cvs, n_perm);
d_fcr_mi_cohen  = cohen_d_local(d_fcr_mi);  d_fcr_cvs_cohen = cohen_d_local(d_fcr_cvs);
d_fci_mi_cohen  = cohen_d_local(d_fci_mi);  d_fci_cvs_cohen = cohen_d_local(d_fci_cvs);
if any(has_fraccorr)
    fprintf('\n══════════════════ Group-level "time in correct zone" (n=%d subjects) ══════════════════\n', sum(has_fraccorr));
    fprintf('  -- raw (pre-integrator) --\n');
    fprintf('  Hybrid - MI-only   : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
            100*mean(d_fcr_mi,'omitnan'),  p_fcr_mi,  stars(p_fcr_mi),  d_fcr_mi_cohen,  cohen_d_label(d_fcr_mi_cohen));
    fprintf('  Hybrid - CVSA-only : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
            100*mean(d_fcr_cvs,'omitnan'), p_fcr_cvs, stars(p_fcr_cvs), d_fcr_cvs_cohen, cohen_d_label(d_fcr_cvs_cohen));
    fprintf('  -- integrated (post-integrator control signal) --\n');
    fprintf('  Hybrid - MI-only   : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
            100*mean(d_fci_mi,'omitnan'),  p_fci_mi,  stars(p_fci_mi),  d_fci_mi_cohen,  cohen_d_label(d_fci_mi_cohen));
    fprintf('  Hybrid - CVSA-only : mean delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
            100*mean(d_fci_cvs,'omitnan'), p_fci_cvs, stars(p_fci_cvs), d_fci_cvs_cohen, cohen_d_label(d_fci_cvs_cohen));
    fprintf('  (two-sided sign-flip permutation, %d perms; subject is the statistical unit; fraction of ALL-trial, whole-CF-duration frames with target-class signal > 0.5)\n', n_perm);
    fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');
end

% ── Early vs late HIT split (COUNTERFACTUAL only): would MI-only / CVSA-only
%    have kept up? Real MI-only/CVSA-only sessions don't share trial indices
%    with Hybrid (different recordings), so this cross-stream "would the
%    alternative also have hit THIS SAME trial" comparison is only
%    well-posed on the counterfactual replay, where all 3 streams run on the
%    identical trials. For each subject: split Hybrid's own HIT trials into
%    "early" (resolved within cvsa_influence, while CVSA still actively
%    weighted) vs "late" (resolved after, alpha~0 -- fused signal is
%    already ~pure MI by construction), then check whether MI-only /
%    CVSA-only ALSO hit on those exact same trial indices. A low keep-up
%    rate in the "early" group (but not "late", where Hybrid~MI already)
%    would mean CVSA specifically enabled the fast resolutions. ────────────
has_oc = ~cellfun(@isempty, {flat_cf.oc_hyb});
keepup_mi_early = nan(n_subj, 1); keepup_mi_late = nan(n_subj, 1);
keepup_cv_early = nan(n_subj, 1); keepup_cv_late = nan(n_subj, 1);
n_early_subj    = zeros(n_subj, 1); n_late_subj = zeros(n_subj, 1);
n_hit_subj      = zeros(n_subj, 1);   % Hybrid's own HIT trials (= n_early + n_late), the denominator below
n_total_subj    = zeros(n_subj, 1);   % ALL Hybrid trials (HIT+MISS+TIMEOUT) -- context only, no longer the denominator
pool_early_mi = []; pool_early_cv = []; pool_late_mi = []; pool_late_cv = [];

% ── "Cost" of fusion (complements the keep-up rate above): on Hybrid's own
%    MISS/TIMEOUT trials (i.e. Hybrid FAILED), would MI-only or CVSA-only
%    alone have HIT instead on that exact same trial? A non-zero rescue
%    rate means fusion actively cost a win that a single stream would have
%    gotten -- the flip side of "does the unimodal stream keep up".
%    MISS (wrong class confidently triggered) and TIMEOUT (neither
%    threshold reached) are kept separate throughout this file (different
%    failure modes), plus a pooled "ALL failures" version. ────────────────
rescue_mi_all = nan(n_subj, 1);  rescue_cv_all = nan(n_subj, 1);
rescue_mi_miss = nan(n_subj, 1); rescue_cv_miss = nan(n_subj, 1);
rescue_mi_to  = nan(n_subj, 1);  rescue_cv_to  = nan(n_subj, 1);
n_fail_subj = zeros(n_subj, 1); n_miss_subj = zeros(n_subj, 1); n_to_subj = zeros(n_subj, 1);
pool_fail_mi = []; pool_fail_cv = []; pool_miss_mi = []; pool_miss_cv = []; pool_to_mi = []; pool_to_cv = [];

if any(has_oc)
    rows_oc_all = flat_cf(has_oc);
    for si = 1:n_subj
        mask = strcmp({rows_oc_all.subject}, subjects{si});
        if ~any(mask), continue; end
        rows_oc = rows_oc_all(mask);
        oc_h = vertcat(rows_oc.oc_hyb);   % [n x 2]: outcome_code, t_event_s
        oc_m = vertcat(rows_oc.oc_mi);
        oc_c = vertcat(rows_oc.oc_cvsa);
        cinf = mean([rows_oc.cvsa_inf_window], 'omitnan');
        is_hit   = oc_h(:,1) == 1;
        is_early = is_hit & oc_h(:,2) <  cinf;
        is_late  = is_hit & oc_h(:,2) >= cinf;
        n_early_subj(si) = sum(is_early);
        n_late_subj(si)  = sum(is_late);
        n_hit_subj(si)   = sum(is_hit);
        n_total_subj(si) = size(oc_h, 1);
        if any(is_early)
            keepup_mi_early(si) = mean(oc_m(is_early,1) == 1);
            keepup_cv_early(si) = mean(oc_c(is_early,1) == 1);
        end
        if any(is_late)
            keepup_mi_late(si) = mean(oc_m(is_late,1) == 1);
            keepup_cv_late(si) = mean(oc_c(is_late,1) == 1);
        end
        pool_early_mi = [pool_early_mi; oc_m(is_early,1)==1]; %#ok<AGROW>
        pool_early_cv = [pool_early_cv; oc_c(is_early,1)==1]; %#ok<AGROW>
        pool_late_mi  = [pool_late_mi;  oc_m(is_late,1)==1];  %#ok<AGROW>
        pool_late_cv  = [pool_late_cv;  oc_c(is_late,1)==1];  %#ok<AGROW>

        is_fail = ~is_hit;              % Hybrid MISS or TIMEOUT
        is_miss = oc_h(:,1) == 2;
        is_to   = oc_h(:,1) == 3;
        n_fail_subj(si) = sum(is_fail);
        n_miss_subj(si) = sum(is_miss);
        n_to_subj(si)   = sum(is_to);
        if any(is_fail)
            rescue_mi_all(si) = mean(oc_m(is_fail,1) == 1);
            rescue_cv_all(si) = mean(oc_c(is_fail,1) == 1);
        end
        if any(is_miss)
            rescue_mi_miss(si) = mean(oc_m(is_miss,1) == 1);
            rescue_cv_miss(si) = mean(oc_c(is_miss,1) == 1);
        end
        if any(is_to)
            rescue_mi_to(si) = mean(oc_m(is_to,1) == 1);
            rescue_cv_to(si) = mean(oc_c(is_to,1) == 1);
        end
        pool_fail_mi = [pool_fail_mi; oc_m(is_fail,1)==1]; %#ok<AGROW>
        pool_fail_cv = [pool_fail_cv; oc_c(is_fail,1)==1]; %#ok<AGROW>
        pool_miss_mi = [pool_miss_mi; oc_m(is_miss,1)==1]; %#ok<AGROW>
        pool_miss_cv = [pool_miss_cv; oc_c(is_miss,1)==1]; %#ok<AGROW>
        pool_to_mi   = [pool_to_mi;   oc_m(is_to,1)==1];   %#ok<AGROW>
        pool_to_cv   = [pool_to_cv;   oc_c(is_to,1)==1];   %#ok<AGROW>
    end
end
grand_keepup_mi_early = mean(pool_early_mi, 'omitnan'); grand_keepup_mi_late = mean(pool_late_mi, 'omitnan');
grand_keepup_cv_early = mean(pool_early_cv, 'omitnan'); grand_keepup_cv_late = mean(pool_late_cv, 'omitnan');

% What fraction of Hybrid's own HIT trials (n_hit_subj = n_early+n_late,
% NOT all Hybrid trials incl. MISS/TIMEOUT) landed in each bucket -- needed
% to read n_early/n_late as a rate, not just a raw count (e.g. "31 early
% hits" means little without knowing whether that's 31/40 or 31/200 HITs).
% By construction pct_early_of_trials + pct_late_of_trials == 100% (every
% HIT trial is either early or late, no third bucket).
pct_early_of_trials = 100 * n_early_subj ./ max(n_hit_subj, 1);
pct_late_of_trials  = 100 * n_late_subj  ./ max(n_hit_subj, 1);
grand_n_total_oc    = sum(n_total_subj);
grand_n_hit_oc      = sum(n_hit_subj);
grand_pct_early = 100 * sum(n_early_subj) / max(grand_n_hit_oc, 1);
grand_pct_late  = 100 * sum(n_late_subj)  / max(grand_n_hit_oc, 1);

d_keepup_mi = keepup_mi_early - keepup_mi_late;
d_keepup_cv = keepup_cv_early - keepup_cv_late;
p_keepup_mi = sign_flip_test_local(d_keepup_mi, n_perm);
p_keepup_cv = sign_flip_test_local(d_keepup_cv, n_perm);
d_keepup_mi_cohen = cohen_d_local(d_keepup_mi);
d_keepup_cv_cohen = cohen_d_local(d_keepup_cv);

fprintf('\n══════ Early vs late Hybrid-HIT: would MI-only/CVSA-only have kept up? (n=%d subjects) ══════\n', sum(has_oc));
fprintf('  %-10s  %6s  %6s  %6s  %8s  %8s  %8s  %12s  %12s  %12s  %12s\n', 'subject', 'n_early', 'n_late', 'n_hit', 'n_total', ...
        'early%', 'late%', 'MI keep-early', 'MI keep-late', 'CV keep-early', 'CV keep-late');
for si = 1:n_subj
    if isnan(keepup_mi_early(si)) && isnan(keepup_mi_late(si)), continue; end
    fprintf('  %-10s  %6d  %6d  %6d  %8d  %7.1f%%  %7.1f%%  %12s  %12s  %12s  %12s\n', subjects{si}, ...
            n_early_subj(si), n_late_subj(si), n_hit_subj(si), n_total_subj(si), pct_early_of_trials(si), pct_late_of_trials(si), ...
            fmt_pct(keepup_mi_early(si)), fmt_pct(keepup_mi_late(si)), fmt_pct(keepup_cv_early(si)), fmt_pct(keepup_cv_late(si)));
end
fprintf('  %-10s  %6d  %6d  %6d  %8d  %7.1f%%  %7.1f%%  %12s  %12s  %12s  %12s  (GRAND, trial-pooled)\n', 'GRAND', ...
        numel(pool_early_mi), numel(pool_late_mi), grand_n_hit_oc, grand_n_total_oc, grand_pct_early, grand_pct_late, ...
        fmt_pct(grand_keepup_mi_early), fmt_pct(grand_keepup_mi_late), fmt_pct(grand_keepup_cv_early), fmt_pct(grand_keepup_cv_late));
fprintf('  ("early%%"/"late%%" = share of Hybrid''s own HIT trials (n_hit = n_early+n_late, NOT n_total) that resolved in that window --\n');
fprintf('   e.g. GRAND: of the trials Hybrid actually hit, %.1f%% resolved within cvsa_influence, %.1f%% after it;\n', grand_pct_early, grand_pct_late);
fprintf('   n_total is shown only for context -- how many trials (incl. MISS/TIMEOUT) n_hit came out of)\n');
fprintf('  MI-only  early-late delta : %+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
        100*mean(d_keepup_mi,'omitnan'), p_keepup_mi, stars(p_keepup_mi), d_keepup_mi_cohen, cohen_d_label(d_keepup_mi_cohen));
fprintf('  CVSA-only early-late delta : %+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
        100*mean(d_keepup_cv,'omitnan'), p_keepup_cv, stars(p_keepup_cv), d_keepup_cv_cohen, cohen_d_label(d_keepup_cv_cohen));
fprintf('  (COUNTERFACTUAL only -- same trials, all 3 streams; NOTE: even "late" HITs can still be partly enabled by\n');
fprintf('   early CVSA-driven push -- the leaky integrator has memory across the WHOLE trial, so a high MI keep-up\n');
fprintf('   rate late is not automatically guaranteed just because alpha~0 by the time the trial actually resolves)\n');
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

% ── GRAND (trial-pooled) rescue rates + group-level MI-vs-CVSA test ───────
grand_rescue_mi_all  = mean(pool_fail_mi, 'omitnan');  grand_rescue_cv_all  = mean(pool_fail_cv, 'omitnan');
grand_rescue_mi_miss = mean(pool_miss_mi, 'omitnan');  grand_rescue_cv_miss = mean(pool_miss_cv, 'omitnan');
grand_rescue_mi_to   = mean(pool_to_mi, 'omitnan');    grand_rescue_cv_to   = mean(pool_to_cv, 'omitnan');
grand_n_fail_oc = sum(n_fail_subj); grand_n_miss_oc = sum(n_miss_subj); grand_n_to_oc = sum(n_to_subj);

d_rescue_all  = rescue_mi_all  - rescue_cv_all;
d_rescue_miss = rescue_mi_miss - rescue_cv_miss;
d_rescue_to   = rescue_mi_to   - rescue_cv_to;
p_rescue_all  = sign_flip_test_local(d_rescue_all,  n_perm);
p_rescue_miss = sign_flip_test_local(d_rescue_miss, n_perm);
p_rescue_to   = sign_flip_test_local(d_rescue_to,   n_perm);
d_rescue_all_cohen  = cohen_d_local(d_rescue_all);
d_rescue_miss_cohen = cohen_d_local(d_rescue_miss);
d_rescue_to_cohen   = cohen_d_local(d_rescue_to);

fprintf('\n══════ Cost of fusion: on Hybrid''s own MISS/TIMEOUT trials, would MI-only/CVSA-only have HIT instead? (n=%d subjects) ══════\n', sum(has_oc));
fprintf('  %-10s  %6s  %6s  %6s  %13s  %13s  %13s  %13s  %13s  %13s\n', 'subject', 'n_fail', 'n_miss', 'n_to', ...
        'MI resc-ALL', 'CV resc-ALL', 'MI resc-MISS', 'CV resc-MISS', 'MI resc-TO', 'CV resc-TO');
for si = 1:n_subj
    if isnan(rescue_mi_all(si)) && isnan(rescue_cv_all(si)), continue; end
    fprintf('  %-10s  %6d  %6d  %6d  %13s  %13s  %13s  %13s  %13s  %13s\n', subjects{si}, ...
            n_fail_subj(si), n_miss_subj(si), n_to_subj(si), ...
            fmt_pct(rescue_mi_all(si)), fmt_pct(rescue_cv_all(si)), fmt_pct(rescue_mi_miss(si)), fmt_pct(rescue_cv_miss(si)), ...
            fmt_pct(rescue_mi_to(si)), fmt_pct(rescue_cv_to(si)));
end
fprintf('  %-10s  %6d  %6d  %6d  %13s  %13s  %13s  %13s  %13s  %13s  (GRAND, trial-pooled)\n', 'GRAND', ...
        grand_n_fail_oc, grand_n_miss_oc, grand_n_to_oc, ...
        fmt_pct(grand_rescue_mi_all), fmt_pct(grand_rescue_cv_all), fmt_pct(grand_rescue_mi_miss), fmt_pct(grand_rescue_cv_miss), ...
        fmt_pct(grand_rescue_mi_to), fmt_pct(grand_rescue_cv_to));
fprintf('  ("rescue rate" = %% of Hybrid''s failed trials [that failure type] where the unimodal stream ALONE would have HIT --\n');
fprintf('   a non-zero rescue rate means fusion actively cost a win the single stream would have gotten on that exact trial)\n');
fprintf('  MI-only vs CVSA-only, ALL failures  : delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
        100*mean(d_rescue_all,'omitnan'), p_rescue_all, stars(p_rescue_all), d_rescue_all_cohen, cohen_d_label(d_rescue_all_cohen));
fprintf('  MI-only vs CVSA-only, MISS only     : delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
        100*mean(d_rescue_miss,'omitnan'), p_rescue_miss, stars(p_rescue_miss), d_rescue_miss_cohen, cohen_d_label(d_rescue_miss_cohen));
fprintf('  MI-only vs CVSA-only, TIMEOUT only  : delta=%+.1f%%  p=%.4f  %s  d=%+.2f (%s)\n', ...
        100*mean(d_rescue_to,'omitnan'), p_rescue_to, stars(p_rescue_to), d_rescue_to_cohen, cohen_d_label(d_rescue_to_cohen));
fprintf('  (COUNTERFACTUAL only -- same trials, all 3 streams; MISS = wrong class confidently triggered, TIMEOUT = neither reached)\n');
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

% ── Plain (naive) group-level sign-flip test on buf_adv_mi_win_subj itself
%    -- the companion this was missing next to the Stouffer meta-analysis
%    above (fus_adv_subj/p_fusadv_grp already has this pairing; this one
%    didn't, which is why Fig 3c's GROUP bar showed only the meta p and
%    looked inconsistent with its own visible error bar). Subject is the
%    unit, same 2^n_subj floor as every other naive group test here -- at
%    n=14 the smallest achievable p is 1/2^14 ~ 6e-5, so THIS test, not the
%    meta-analysis, is the one limited by cohort size; the meta-analysis
%    exists precisely to sidestep that limitation, not to replace this
%    test -- both are shown so the two can't be confused for each other. ──
p_bufadvwin_naive_grp   = sign_flip_test_local(buf_adv_mi_win_subj, n_perm);
d_bufadvwin_naive_cohen = cohen_d_local(buf_adv_mi_win_subj);
fprintf('\n  Naive group-level sign-flip test on buf_adv_mi_win_subj (subject is the unit, %d subjects): p=%.4f  %s  d=%+.2f (%s)\n', ...
        sum(~isnan(buf_adv_mi_win_subj)), p_bufadvwin_naive_grp, stars(p_bufadvwin_naive_grp), d_bufadvwin_naive_cohen, cohen_d_label(d_bufadvwin_naive_cohen));
fprintf('  (this is the test limited by cohort size, NOT the Stouffer meta-analysis above -- the two are reported\n');
fprintf('   side by side in group_cvsa_help_integrator.svg precisely so they are never mistaken for one another)\n');

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
% Is a subject's real Hybrid accuracy predicted by how often their raw EEG
% had a dead/noisy channel during CF (main_session_overview.m Fig 11,
% pooled across paradigms/sessions)? Distinguishes "this subject is a poor
% BCI performer" from "this subject had bad electrode contact" -- expected
% direction is NEGATIVE (more flagged trials -> lower accuracy).
[r_sigqc_vs_real, p_sigqc_vs_real] = pearson_perm_test_local(sigqc_rate_overall, real_acc(:,3), n_perm);
fprintf('\n══════════════════ Cross-subject correlates of real Hybrid performance ══════════════════\n');
fprintf('  Counterfactual adv. (Hyb-MI) vs real adv. (Hyb-MI) : r=%+.2f  p=%.4f  %s\n', r_cf_vs_real,  p_cf_vs_real,  stars(p_cf_vs_real));
fprintf('  Fusion advantage (P_fused-P_MI) vs real Hybrid acc : r=%+.2f  p=%.4f  %s\n', r_fus_vs_real, p_fus_vs_real, stars(p_fus_vs_real));
fprintf('  ERD/ERS-CSP grounding (r) vs real Hybrid acc       : r=%+.2f  p=%.4f  %s\n', r_erd_vs_real, p_erd_vs_real, stars(p_erd_vs_real));
fprintf('  CVSA-only quality vs Hybrid-MI counterfactual adv. : r=%+.2f  p=%.4f  %s  (does the advantage survive weak CVSA?)\n', ...
        r_cvsa_quality, p_cvsa_quality, stars(p_cvsa_quality));
fprintf('  Signal-quality flagged-trial rate vs real Hybrid acc: r=%+.2f  p=%.4f  %s  (negative = dirtier signal, lower accuracy)\n', ...
        r_sigqc_vs_real, p_sigqc_vs_real, stars(p_sigqc_vs_real));
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
sig_pairs_real_acc = {{3,1,p_real_hyb_mi}, {3,2,p_real_hyb_cvs}, {1,2,p_real_mi_cvs}};   % cols: MI,CVSA,Hybrid
plot_subject_grand_bars(ax1, real_acc, grand_acc, subjects, par_labels, cols3, bw, true, false, sig_pairs_real_acc);
yline(ax1, 50, 'k:', 'HandleVisibility', 'off');
ylabel(ax1, 'accuracy (hit/(hit+miss+timeout)) (%)');
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

saveas(fig1, fullfile(out_dir, sprintf('01_%s_real_session_overview.svg', group_tag)), 'svg');
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
        '"GRAND" = trial-pooled across all subjects\n' ...
        'TTH is HIT-trials-only (speed when CORRECT, not "speed to any answer")\n' ...
        'TTH significance (sign-flip test, subject=unit -- full test in 06_%s_speed_itr.svg):\n' ...
        'Hybrid-MI %+.2fs (%s, d=%+.2f)   Hybrid-CVSA %+.2fs (%s, d=%+.2f)'], ...
        n_subj, group_tag, mean(d_tth_hyb_mi,'omitnan'), stars(p_tth_hyb_mi), d_tth_hyb_mi_cohen, ...
        mean(d_tth_hyb_cvs,'omitnan'), stars(p_tth_hyb_cvs), d_tth_hyb_cvs_cohen), 'Interpreter', 'none');

saveas(fig1b, fullfile(out_dir, sprintf('02_%s_real_times.svg', group_tag)), 'svg');
if ~SHOW_FIGURES, close(fig1b); end

% ── Figure 1c: real-session performance breakdown -- decided-trial accuracy
%    (TIMEOUT excluded) on top, false-positive rate + timeout rate below,
%    per paradigm, per subject + GRAND. Gold star / red triangle mark the
%    best/worst subject FOR THAT PARADIGM (ties shown together). ──────────
fig1c = figure('Name', 'Group Analysis — Real Session Performance Breakdown', 'Color', 'w', ...
               'NumberTitle', 'off', 'Visible', fig_vis);
set(fig1c, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

axA1 = subplot(2, 3, [1 2 3]); hold(axA1, 'on');
sig_pairs_real_acc_dec = {{3,1,p_real_hyb_mi_dec}, {3,2,p_real_hyb_cvs_dec}, {1,2,p_real_mi_cvs_dec}};   % cols: MI,CVSA,Hybrid
plot_subject_grand_bars(axA1, real_acc_dec, grand_acc_dec, subjects, par_labels, cols3, bw, true, false, sig_pairs_real_acc_dec);
yline(axA1, 50, 'k:', 'HandleVisibility', 'off');
ylabel(axA1, 'accuracy (hit/(hit+miss)) (%)');
legend(axA1, 'Location', 'south', 'FontSize', 9);
title(axA1, sprintf('Accuracy on decided trials, TIMEOUT excluded (n=%d subjects)', n_subj), 'FontWeight', 'bold');
grid(axA1, 'on');
% Best/worst markers sit at a fixed offset to the right of each condition's
% box (plot_subject_grand_bars no longer places subjects along x, so these
% can't be aligned to a specific subject's own jittered dot -- they mark
% the VALUE, stacked at that y height, next to the box they belong to).
for pi = 1:3
    for si = idx_best_real{pi}
        y = min(100*real_acc_dec(si,pi), 108);
        plot(axA1, pi + 0.42, y, 'p', 'MarkerSize', 12, ...
             'MarkerFaceColor', [0.85 0.65 0.10], 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    for si = idx_worst_real{pi}
        y = min(100*real_acc_dec(si,pi), 108);
        plot(axA1, pi + 0.42, y, 'v', 'MarkerSize', 9, ...
             'MarkerFaceColor', [0.75 0.15 0.15], 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
end

axA2 = subplot(2, 3, 4); hold(axA2, 'on');
plot_subject_grand_bars(axA2, real_fp_rate, grand_fp_rate, subjects, par_labels, cols3, bw);
ylabel(axA2, 'False positive rate (%)');
title(axA2, 'MISS: wrong-class threshold reached', 'FontWeight', 'bold', 'FontSize', 9);
grid(axA2, 'on');

axA3 = subplot(2, 3, 5); hold(axA3, 'on');
plot_subject_grand_bars(axA3, real_to_rate, grand_to_rate, subjects, par_labels, cols3, bw);
ylabel(axA3, 'Timeout rate (%)');
title(axA3, 'TIMEOUT: no threshold reached', 'FontWeight', 'bold', 'FontSize', 9);
grid(axA3, 'on');

axA4 = subplot(2, 3, 6);
plot_acc_vs_timeout_scatter(axA4, real_acc_dec, real_to_rate, cols3, par_labels);

sgtitle(fig1c, sprintf(['Group Analysis — Real Session Performance Breakdown (n=%d subjects)\n' ...
        'gold star = best subject per paradigm, red triangle = worst (ties shown together), ranked by decided-trial accuracy  |  "GRAND" = trial-pooled'], n_subj), ...
        'Interpreter', 'none');

saveas(fig1c, fullfile(out_dir, sprintf('03_%s_real_accuracy_breakdown.svg', group_tag)), 'svg');
if ~SHOW_FIGURES, close(fig1c); end

% ── Figure 2: counterfactual accuracy + group-level advantage ─────────────
fig2 = figure('Name', 'Group Analysis — Counterfactual Advantage', 'Color', 'w', ...
              'NumberTitle', 'off', 'Visible', fig_vis);
set(fig2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

ax2 = subplot(1, 2, 1); hold(ax2, 'on');
cols_cf = {COL_hybrid, COL_mi_only, COL_cvsa_only};
plot_subject_grand_bars(ax2, cf_acc, grand_cf_acc, subjects, cf_labels, cols_cf, bw);
yline(ax2, 50, 'k:', 'HandleVisibility', 'off');
ylabel(ax2, 'simulated accuracy (hit/(hit+miss+timeout)) (%)');
legend(ax2, 'Location', 'south', 'FontSize', 9);
title(ax2, 'Counterfactual accuracy per subject', 'FontWeight', 'bold');
grid(ax2, 'on');

ax3 = subplot(1, 2, 2); hold(ax3, 'on');
% Box = mean +/- SD across subjects, whiskers = true min/max, dots = each
% subject's own advantage -- directly shows whether the advantage is
% UNIFORM across the cohort (tight box) or HETEROGENEOUS (wide box/dots
% spanning positive to negative -- helps some subjects a lot, others not
% at all or even negatively).
grand_delta_cf = [grand_cf_acc(1)-grand_cf_acc(2), grand_cf_acc(1)-grand_cf_acc(3)];   % trial-pooled Hybrid-MI/Hybrid-CVSA
plot_delta_boxplot(ax3, [d_mi_subj, d_cvs_subj], grand_delta_cf, ...
                    {'Hybrid - MI-only', 'Hybrid - CVSA-only'}, {COL_mi_only, COL_cvsa_only}, 100, [p_mi_grp, p_cvs_grp]);
m_mi  = mean(d_mi_subj,  'omitnan');
m_cvs = mean(d_cvs_subj, 'omitnan');
ylabel(ax3, '\Delta accuracy (hit/(hit+miss+timeout)), Hybrid - unimodal, %');
legend(ax3, 'Location', 'best', 'FontSize', 9);
title(ax3, sprintf('Per-subject advantage (box=mean+-SD, dots=subjects)  |  group: Hyb-MI=%+.1f%% (%s, d=%+.2f), Hyb-CVSA=%+.1f%% (%s, d=%+.2f)', ...
      100*m_mi, stars(p_mi_grp), d_mi_grp_cohen, 100*m_cvs, stars(p_cvs_grp), d_cvs_grp_cohen), 'FontWeight', 'bold', 'FontSize', 9);
grid(ax3, 'on');

sgtitle(fig2, sprintf(['Group Analysis — Counterfactual Advantage (n=%d subjects)\n' ...
        'SIMULATED, not real: SAME integrator/buffer/thresholds driven counterfactually by Hybrid-fused, MI-only, or CVSA-only signal\n' ...
        'on the SAME trials -- isolates the fusion algorithm''s own effect from real-session differences in trial count/timing/thresholds'], ...
        n_subj), 'Interpreter', 'none');

saveas(fig2, fullfile(out_dir, sprintf('04_%s_counterfactual_advantage.svg', group_tag)), 'svg');
if ~SHOW_FIGURES, close(fig2); end

% ── Figure 2d (04b): "time in correct zone" -- fraction of EVERY trial's
%    own CF duration (5s window, ALL trials: hit/miss/timeout alike, no
%    outcome-based truncation) that the target-class control signal spent
%    on the correct side of 0.5. Complements Fig 2's binary hit/miss/timeout
%    view: a stream can lose on the binary outcome yet still spend most of
%    the trial leaning the right way (or vice versa). RAW = pre-integrator
%    classifier signal; INT = post-integrator control signal (what the VR
%    feedback actually shows). Skipped if no subject's counterfactual_
%    summary.mat has these fields yet (regenerate via main_hybrid_advantage_
%    integ.m). ──────────────────────────────────────────────────────────
if any(has_fraccorr)
    fig2d = figure('Name', 'Group Analysis — Time in Correct Zone', 'Color', 'w', ...
                   'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig2d, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    x1 = (1:n_subj) - 0.12; x2 = (1:n_subj) + 0.12;

    ax2d1 = subplot(2, 2, 1); hold(ax2d1, 'on');
    for pi = 1:3
        x = (1:n_subj) + (pi-2)*bw;
        v = 100*fraccorr_raw_subj(:,pi);
        bar(ax2d1, x, v, bw*0.9, 'FaceColor', cols_cf{pi}, 'EdgeColor', 'k', 'LineWidth', 1, 'DisplayName', cf_labels{pi});
        m = mean(v, 'omitnan');
        if ~isnan(m), plot(ax2d1, [0.5, n_subj+0.5], [m m], '--', 'Color', cols_cf{pi}*0.6, 'LineWidth', 1.5, 'HandleVisibility', 'off'); end
    end
    set(ax2d1, 'XTick', 1:n_subj, 'XTickLabel', subjects, 'XLim', [0.4, n_subj+0.6], 'YLim', [0 100]);
    ylabel(ax2d1, '% of CF frames with P(target)>0.5'); legend(ax2d1, 'Location', 'best', 'FontSize', 8);
    title(ax2d1, 'RAW (pre-integrator classifier signal), per subject', 'FontWeight', 'bold', 'FontSize', 9);
    grid(ax2d1, 'on');

    ax2d2 = subplot(2, 2, 2); hold(ax2d2, 'on');
    scatter(ax2d2, x1, 100*d_fcr_mi,  70, COL_mi_only,   'filled', 'MarkerEdgeColor', 'k');
    scatter(ax2d2, x2, 100*d_fcr_cvs, 70, COL_cvsa_only, 'filled', 'MarkerEdgeColor', 'k');
    yline(ax2d2, 0, 'k-', 'HandleVisibility', 'off');
    plot_group_marker(ax2d2, n_subj, 100*d_fcr_mi,  COL_mi_only,   p_fcr_mi,  0.7);
    plot_group_marker(ax2d2, n_subj, 100*d_fcr_cvs, COL_cvsa_only, p_fcr_cvs, 1.0);
    set(ax2d2, 'XTick', [1:n_subj, n_subj+0.85], 'XTickLabel', [subjects, {'group'}], 'XLim', [0.4, n_subj+1.4]);
    ylabel(ax2d2, '\Delta %% time in correct zone [Hybrid - unimodal]');
    legend(ax2d2, {'Hybrid - MI-only', 'Hybrid - CVSA-only'}, 'Location', 'best', 'FontSize', 8);
    title(ax2d2, sprintf('RAW: %+.1f%% (%s, d=%+.2f) vs MI-only, %+.1f%% (%s, d=%+.2f) vs CVSA-only', ...
          100*mean(d_fcr_mi,'omitnan'), stars(p_fcr_mi), d_fcr_mi_cohen, ...
          100*mean(d_fcr_cvs,'omitnan'), stars(p_fcr_cvs), d_fcr_cvs_cohen), 'FontWeight', 'bold', 'FontSize', 9);
    grid(ax2d2, 'on');

    ax2d3 = subplot(2, 2, 3); hold(ax2d3, 'on');
    for pi = 1:3
        x = (1:n_subj) + (pi-2)*bw;
        v = 100*fraccorr_int_subj(:,pi);
        bar(ax2d3, x, v, bw*0.9, 'FaceColor', cols_cf{pi}, 'EdgeColor', 'k', 'LineWidth', 1, 'DisplayName', cf_labels{pi});
        m = mean(v, 'omitnan');
        if ~isnan(m), plot(ax2d3, [0.5, n_subj+0.5], [m m], '--', 'Color', cols_cf{pi}*0.6, 'LineWidth', 1.5, 'HandleVisibility', 'off'); end
    end
    set(ax2d3, 'XTick', 1:n_subj, 'XTickLabel', subjects, 'XLim', [0.4, n_subj+0.6], 'YLim', [0 100]);
    ylabel(ax2d3, '% of CF frames with buffer(target)>0.5'); legend(ax2d3, 'Location', 'best', 'FontSize', 8);
    title(ax2d3, 'INTEGRATED (post-integrator control signal), per subject', 'FontWeight', 'bold', 'FontSize', 9);
    grid(ax2d3, 'on');

    ax2d4 = subplot(2, 2, 4); hold(ax2d4, 'on');
    scatter(ax2d4, x1, 100*d_fci_mi,  70, COL_mi_only,   'filled', 'MarkerEdgeColor', 'k');
    scatter(ax2d4, x2, 100*d_fci_cvs, 70, COL_cvsa_only, 'filled', 'MarkerEdgeColor', 'k');
    yline(ax2d4, 0, 'k-', 'HandleVisibility', 'off');
    plot_group_marker(ax2d4, n_subj, 100*d_fci_mi,  COL_mi_only,   p_fci_mi,  0.7);
    plot_group_marker(ax2d4, n_subj, 100*d_fci_cvs, COL_cvsa_only, p_fci_cvs, 1.0);
    set(ax2d4, 'XTick', [1:n_subj, n_subj+0.85], 'XTickLabel', [subjects, {'group'}], 'XLim', [0.4, n_subj+1.4]);
    ylabel(ax2d4, '\Delta %% time in correct zone [Hybrid - unimodal]');
    legend(ax2d4, {'Hybrid - MI-only', 'Hybrid - CVSA-only'}, 'Location', 'best', 'FontSize', 8);
    title(ax2d4, sprintf('INTEGRATED: %+.1f%% (%s, d=%+.2f) vs MI-only, %+.1f%% (%s, d=%+.2f) vs CVSA-only', ...
          100*mean(d_fci_mi,'omitnan'), stars(p_fci_mi), d_fci_mi_cohen, ...
          100*mean(d_fci_cvs,'omitnan'), stars(p_fci_cvs), d_fci_cvs_cohen), 'FontWeight', 'bold', 'FontSize', 9);
    grid(ax2d4, 'on');

    sgtitle(fig2d, sprintf(['Group Analysis — Time in Correct Zone (n=%d subjects)\n' ...
            '%% of EVERY trial''s own CF duration (ALL trials: hit+miss+timeout, no outcome-based truncation) ' ...
            'with the target-class signal on the correct side of 0.5'], sum(has_fraccorr)), 'Interpreter', 'none');

    saveas(fig2d, fullfile(out_dir, sprintf('04b_%s_time_in_correct_zone.svg', group_tag)), 'svg');
    if ~SHOW_FIGURES, close(fig2d); end
end

% ── Figure 2b: simulated (counterfactual) performance breakdown -- same
%    layout as Fig 1c, but for the 3 counterfactual streams instead of the
%    3 real paradigms. Best/worst marked PER STREAM (ties shown together). ─
fig2b = figure('Name', 'Group Analysis — Simulated Performance Breakdown', 'Color', 'w', ...
               'NumberTitle', 'off', 'Visible', fig_vis);
set(fig2b, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

axB1 = subplot(2, 3, [1 2 3]); hold(axB1, 'on');
sig_pairs_cf_acc_dec = {{1,2,p_mi_grp_dec}, {1,3,p_cvs_grp_dec}};   % cols: Hybrid,MI-only,CVSA-only (no precomputed MI-vs-CVSA test)
plot_subject_grand_bars(axB1, cf_acc_dec, grand_cf_acc_dec, subjects, cf_labels, cols_cf, bw, true, false, sig_pairs_cf_acc_dec);
yline(axB1, 50, 'k:', 'HandleVisibility', 'off');
ylabel(axB1, 'simulated accuracy (hit/(hit+miss)) (%)');
legend(axB1, 'Location', 'south', 'FontSize', 9);
title(axB1, sprintf('Counterfactual accuracy on decided trials, TIMEOUT excluded (n=%d subjects)', n_subj), 'FontWeight', 'bold');
grid(axB1, 'on');
% Best/worst markers sit at a fixed offset to the right of each stream's
% box (see the identical comment on Fig 1c's axA1 above).
for s = 1:3
    for si = idx_best_cf{s}
        y = min(100*cf_acc_dec(si,s), 108);
        plot(axB1, s + 0.42, y, 'p', 'MarkerSize', 12, ...
             'MarkerFaceColor', [0.85 0.65 0.10], 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    for si = idx_worst_cf{s}
        y = min(100*cf_acc_dec(si,s), 108);
        plot(axB1, s + 0.42, y, 'v', 'MarkerSize', 9, ...
             'MarkerFaceColor', [0.75 0.15 0.15], 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
end

axB2 = subplot(2, 3, 4); hold(axB2, 'on');
plot_subject_grand_bars(axB2, cf_fp_rate, grand_cf_fp_rate, subjects, cf_labels, cols_cf, bw);
ylabel(axB2, 'False positive rate (%)');
title(axB2, 'MISS: wrong-class threshold reached', 'FontWeight', 'bold', 'FontSize', 9);
grid(axB2, 'on');

axB3 = subplot(2, 3, 5); hold(axB3, 'on');
plot_subject_grand_bars(axB3, cf_to_rate, grand_cf_to_rate, subjects, cf_labels, cols_cf, bw);
ylabel(axB3, 'Timeout rate (%)');
title(axB3, 'TIMEOUT: no threshold reached', 'FontWeight', 'bold', 'FontSize', 9);
grid(axB3, 'on');

axB4 = subplot(2, 3, 6);
plot_acc_vs_timeout_scatter(axB4, cf_acc_dec, cf_to_rate, cols_cf, cf_labels);

sgtitle(fig2b, sprintf(['Group Analysis — Simulated (Counterfactual) Performance Breakdown (n=%d subjects)\n' ...
        'SAME integrator/thresholds, driven counterfactually by each stream on the same Hybrid-session trials\n' ...
        'gold star = best subject per stream, red triangle = worst (ties together)  |  "GRAND" = trial-pooled'], n_subj), ...
        'Interpreter', 'none');

saveas(fig2b, fullfile(out_dir, sprintf('05_%s_simulated_accuracy_breakdown.svg', group_tag)), 'svg');
if ~SHOW_FIGURES, close(fig2b); end

% ── Figure 2c: speed (TTH) and communication rate (ITR), real-session AND
%    counterfactual -- the direct answer to "is Hybrid not just more
%    accurate but also faster / a better bits/min channel?", which none of
%    the figures above test explicitly (they show mean+-SEM bars with no
%    paired significance, or don't cover ITR at all). Same per-subject-
%    delta-scatter + group-errorbar-with-stars layout as Fig 2's ax3. ─────
fig2c = figure('Name', 'Group Analysis — Speed (TTH) & Communication Rate (ITR)', 'Color', 'w', ...
               'NumberTitle', 'off', 'Visible', fig_vis);
set(fig2c, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Box = mean +/- SD across subjects, whiskers = true min/max, dots = each
% subject's own delta, diamond = trial-pooled GRAND delta -- same
% uniform-vs-heterogeneous reading as Fig 04's advantage panel.
ax2c1 = subplot(2, 2, 1);
grand_tth_mi  = grand_tth(3)-grand_tth(1); grand_tth_cvs = grand_tth(3)-grand_tth(2);
plot_delta_boxplot(ax2c1, [d_tth_hyb_mi, d_tth_hyb_cvs], [grand_tth_mi, grand_tth_cvs], ...
                    {'Hybrid - MI', 'Hybrid - CVSA'}, {COL_mi_only, COL_cvsa_only}, 1, [p_tth_hyb_mi, p_tth_hyb_cvs]);
ylabel(ax2c1, '\Delta TTH (s)  [Hybrid - unimodal]  (negative = Hybrid faster)');
legend(ax2c1, 'Location', 'best', 'FontSize', 8);
title(ax2c1, sprintf('REAL sessions: %+.2fs (%s, d=%+.2f) vs MI, %+.2fs (%s, d=%+.2f) vs CVSA', ...
      mean(d_tth_hyb_mi,'omitnan'), stars(p_tth_hyb_mi), d_tth_hyb_mi_cohen, ...
      mean(d_tth_hyb_cvs,'omitnan'), stars(p_tth_hyb_cvs), d_tth_hyb_cvs_cohen), 'FontWeight', 'bold', 'FontSize', 9);
grid(ax2c1, 'on');

ax2c2 = subplot(2, 2, 2);
grand_cftth_mi  = grand_cf_tth(1)-grand_cf_tth(2); grand_cftth_cvs = grand_cf_tth(1)-grand_cf_tth(3);
plot_delta_boxplot(ax2c2, [d_cftth_hyb_mi, d_cftth_hyb_cvs], [grand_cftth_mi, grand_cftth_cvs], ...
                    {'Hybrid - MI-only', 'Hybrid - CVSA-only'}, {COL_mi_only, COL_cvsa_only}, 1, [p_cftth_hyb_mi, p_cftth_hyb_cvs]);
ylabel(ax2c2, '\Delta TTH (s)  [Hybrid - unimodal]');
legend(ax2c2, 'Location', 'best', 'FontSize', 8);
title(ax2c2, sprintf('COUNTERFACTUAL: %+.2fs (%s, d=%+.2f) vs MI-only, %+.2fs (%s, d=%+.2f) vs CVSA-only', ...
      mean(d_cftth_hyb_mi,'omitnan'), stars(p_cftth_hyb_mi), d_cftth_hyb_mi_cohen, ...
      mean(d_cftth_hyb_cvs,'omitnan'), stars(p_cftth_hyb_cvs), d_cftth_hyb_cvs_cohen), 'FontWeight', 'bold', 'FontSize', 9);
grid(ax2c2, 'on');

ax2c3 = subplot(2, 2, 3);
grand_itr_mi  = grand_real_itr_bpm(3)-grand_real_itr_bpm(1); grand_itr_cvs = grand_real_itr_bpm(3)-grand_real_itr_bpm(2);
plot_delta_boxplot(ax2c3, [d_itr_hyb_mi, d_itr_hyb_cvs], [grand_itr_mi, grand_itr_cvs], ...
                    {'Hybrid - MI', 'Hybrid - CVSA'}, {COL_mi_only, COL_cvsa_only}, 1, [p_itr_hyb_mi, p_itr_hyb_cvs]);
ylabel(ax2c3, '\Delta ITR (bits/min)  [Hybrid - unimodal]');
legend(ax2c3, 'Location', 'best', 'FontSize', 8);
title(ax2c3, sprintf('REAL sessions: %+.2f (%s, d=%+.2f) vs MI, %+.2f (%s, d=%+.2f) vs CVSA  (N_{CLS}=%d)', ...
      mean(d_itr_hyb_mi,'omitnan'), stars(p_itr_hyb_mi), d_itr_hyb_mi_cohen, ...
      mean(d_itr_hyb_cvs,'omitnan'), stars(p_itr_hyb_cvs), d_itr_hyb_cvs_cohen, N_CLS), 'FontWeight', 'bold', 'FontSize', 9);
grid(ax2c3, 'on');

ax2c4 = subplot(2, 2, 4);
grand_cfitr_mi  = grand_cf_itr_bpm(1)-grand_cf_itr_bpm(2); grand_cfitr_cvs = grand_cf_itr_bpm(1)-grand_cf_itr_bpm(3);
plot_delta_boxplot(ax2c4, [d_cfitr_hyb_mi, d_cfitr_hyb_cvs], [grand_cfitr_mi, grand_cfitr_cvs], ...
                    {'Hybrid - MI-only', 'Hybrid - CVSA-only'}, {COL_mi_only, COL_cvsa_only}, 1, [p_cfitr_hyb_mi, p_cfitr_hyb_cvs]);
ylabel(ax2c4, '\Delta ITR (bits/min)  [Hybrid - unimodal]');
legend(ax2c4, 'Location', 'best', 'FontSize', 8);
title(ax2c4, sprintf('COUNTERFACTUAL: %+.2f (%s, d=%+.2f) vs MI-only, %+.2f (%s, d=%+.2f) vs CVSA-only', ...
      mean(d_cfitr_hyb_mi,'omitnan'), stars(p_cfitr_hyb_mi), d_cfitr_hyb_mi_cohen, ...
      mean(d_cfitr_hyb_cvs,'omitnan'), stars(p_cfitr_hyb_cvs), d_cfitr_hyb_cvs_cohen), 'FontWeight', 'bold', 'FontSize', 9);
grid(ax2c4, 'on');

sgtitle(fig2c, sprintf(['Group Analysis — Speed (TTH) & Communication Rate (ITR), n=%d subjects\n' ...
        'Left: REAL recorded sessions.  Right: COUNTERFACTUAL (same Hybrid trials, driven by each stream'' signal)\n' ...
        'ITR = Wolpaw bits/min, N_{CLS}=%d assumed\n' ...
        'Test: sign-flip permutation, subject=unit, + paired Cohen''s d\n' ...
        'TTH/ITR use HIT-trials-only speed -- "faster" cannot mean "rushing to a wrong answer" (see MISS rate in Fig 03/05)'], n_subj, N_CLS), 'Interpreter', 'tex');

saveas(fig2c, fullfile(out_dir, sprintf('06_%s_speed_itr.svg', group_tag)), 'svg');
if ~SHOW_FIGURES, close(fig2c); end

% ── Figure 2d: early vs late Hybrid-HIT -- would MI-only/CVSA-only have
%    kept up? COUNTERFACTUAL only (see comment above the computation block).
%    Slopegraph style mirrors Fig 3's rescue-vs-cost panel below. ──────────
if any(has_oc)
    fig2d = figure('Name', 'Group Analysis — Early vs Late Hybrid-HIT Keep-Up Rate', 'Color', 'w', ...
                   'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig2d, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    % x-tick labels carry each subject's OWN Hybrid early%/late% split
    % (out of Hybrid's own HIT trials only -- early%+late% == 100%) -- so
    % "58% keep-up" can be read against "how many HIT trials this bucket
    % actually represents". Kept as SHORT plain subject names for XTickLabel -- a wide multi-line
    % label per subject (tried previously) risks MATLAB's automatic tick-
    % label thinning once there are enough subjects that they'd overlap,
    % which silently drops a subset of labels instead of shrinking them.
    % The early%/late% info is now a small rotated text annotation above
    % each subject's own "early" point instead, immune to that thinning
    % since it's regular plotted text, not a tick label.

    ax2d1 = subplot(1, 2, 1); hold(ax2d1, 'on');
    for si = 1:n_subj
        if isnan(keepup_mi_early(si)) || isnan(keepup_mi_late(si)), continue; end
        plot(ax2d1, [si-0.15, si+0.15], 100*[keepup_mi_early(si), keepup_mi_late(si)], '-', ...
             'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    scatter(ax2d1, (1:n_subj)-0.15, 100*keepup_mi_early, 60, COL_mi_only, 'filled', 'MarkerEdgeColor', 'k');
    scatter(ax2d1, (1:n_subj)+0.15, 100*keepup_mi_late,  60, COL_mi_only, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    for si = 1:n_subj
        if isnan(pct_early_of_trials(si)), continue; end
        text(ax2d1, si-0.15, 100.5, sprintf('e%.0f/l%.0f%%', pct_early_of_trials(si), pct_late_of_trials(si)), ...
             'FontSize', 6, 'Rotation', 90, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
    set(ax2d1, 'XTick', 1:n_subj, 'XTickLabel', subjects, 'XLim', [0.4, n_subj+0.6], 'YLim', [0 115]);
    xlabel(ax2d1, 'text above bars: early%/late% = share of THIS subject''s Hybrid HIT trials that resolved early vs late');
    ylabel(ax2d1, 'MI-only HIT rate on Hybrid''s own HIT trials (%)');
    legend(ax2d1, {'early (within cvsa\_influence, filled)', 'late (after, hollow)'}, 'Location', 'best', 'FontSize', 8);
    title(ax2d1, sprintf('MI-only keep-up  |  GRAND: early=%s, late=%s  |  delta=%+.1f%% (%s, d=%+.2f)', ...
          fmt_pct(grand_keepup_mi_early), fmt_pct(grand_keepup_mi_late), ...
          100*mean(d_keepup_mi,'omitnan'), stars(p_keepup_mi), d_keepup_mi_cohen), 'FontWeight', 'bold', 'FontSize', 9);
    grid(ax2d1, 'on');

    ax2d2 = subplot(1, 2, 2); hold(ax2d2, 'on');
    for si = 1:n_subj
        if isnan(keepup_cv_early(si)) || isnan(keepup_cv_late(si)), continue; end
        plot(ax2d2, [si-0.15, si+0.15], 100*[keepup_cv_early(si), keepup_cv_late(si)], '-', ...
             'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    scatter(ax2d2, (1:n_subj)-0.15, 100*keepup_cv_early, 60, COL_cvsa_only, 'filled', 'MarkerEdgeColor', 'k');
    scatter(ax2d2, (1:n_subj)+0.15, 100*keepup_cv_late,  60, COL_cvsa_only, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    for si = 1:n_subj
        if isnan(pct_early_of_trials(si)), continue; end
        text(ax2d2, si-0.15, 100.5, sprintf('e%.0f/l%.0f%%', pct_early_of_trials(si), pct_late_of_trials(si)), ...
             'FontSize', 6, 'Rotation', 90, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
    set(ax2d2, 'XTick', 1:n_subj, 'XTickLabel', subjects, 'XLim', [0.4, n_subj+0.6], 'YLim', [0 115]);
    xlabel(ax2d2, 'text above bars: early%/late% = share of THIS subject''s Hybrid HIT trials that resolved early vs late');
    ylabel(ax2d2, 'CVSA-only HIT rate on Hybrid''s own HIT trials (%)');
    legend(ax2d2, {'early (within cvsa\_influence, filled)', 'late (after, hollow)'}, 'Location', 'best', 'FontSize', 8);
    title(ax2d2, sprintf('CVSA-only keep-up  |  GRAND: early=%s, late=%s  |  delta=%+.1f%% (%s, d=%+.2f)', ...
          fmt_pct(grand_keepup_cv_early), fmt_pct(grand_keepup_cv_late), ...
          100*mean(d_keepup_cv,'omitnan'), stars(p_keepup_cv), d_keepup_cv_cohen), 'FontWeight', 'bold', 'FontSize', 9);
    grid(ax2d2, 'on');

    sgtitle(fig2d, sprintf(['Group Analysis — Early vs Late Hybrid-HIT: would the unimodal stream have kept up? (n=%d subjects)\n' ...
            'COUNTERFACTUAL, same trials as Hybrid''s own HIT trials, split by whether Hybrid resolved within or after cvsa_influence\n' ...
            'GRAND: of Hybrid''s own %d HIT trials, %.1f%% resolved within cvsa_influence ("early") and %.1f%% after it ("late")\n' ...
            'NOTE: the leaky integrator has memory across the WHOLE trial, so even "late" HITs can still be partly enabled by an earlier CVSA push --\n' ...
            'a high MI keep-up rate in the "late" group is NOT automatically guaranteed just because alpha~0 by the time the trial resolves'], ...
            sum(has_oc), grand_n_hit_oc, grand_pct_early, grand_pct_late), 'Interpreter', 'none');

    saveas(fig2d, fullfile(out_dir, sprintf('07_%s_early_late_keepup.svg', group_tag)), 'svg');
    if ~SHOW_FIGURES, close(fig2d); end
end

% ── Figure 07b: cost of fusion -- on Hybrid's own MISS/TIMEOUT trials, would
%    MI-only or CVSA-only ALONE have hit instead? Complements Fig 07 (keep-
%    up = benefit side): a non-zero rescue rate here means fusion actively
%    cost a win the single stream would have gotten on that exact trial.
%    MISS/TIMEOUT kept separate (different failure modes) + a pooled "ALL
%    failures" panel, same paired MI-vs-CVSA sign-flip test + Cohen's d as
%    everywhere else (subject is the unit). ────────────────────────────────
if any(has_oc)
    fig2e = figure('Name', 'Group Analysis — Cost of Fusion (Rescue Rate)', 'Color', 'w', ...
                   'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig2e, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    resc_panels = {
        'ALL failures (MISS+TIMEOUT)', rescue_mi_all,  rescue_cv_all,  d_rescue_all,  p_rescue_all,  d_rescue_all_cohen,  grand_rescue_mi_all,  grand_rescue_cv_all;
        'MISS only',                   rescue_mi_miss, rescue_cv_miss, d_rescue_miss, p_rescue_miss, d_rescue_miss_cohen, grand_rescue_mi_miss, grand_rescue_cv_miss;
        'TIMEOUT only',                rescue_mi_to,   rescue_cv_to,   d_rescue_to,   p_rescue_to,   d_rescue_to_cohen,   grand_rescue_mi_to,   grand_rescue_cv_to;
    };

    for pi = 1:3
        ax = subplot(1, 3, pi); hold(ax, 'on');
        mi_v = 100 * resc_panels{pi,2}; cv_v = 100 * resc_panels{pi,3};
        d_v  = resc_panels{pi,4}; p_v = resc_panels{pi,5}; dcohen_v = resc_panels{pi,6};
        grand_mi_v = 100*resc_panels{pi,7}; grand_cv_v = 100*resc_panels{pi,8};

        jitter_mi = (rand(n_subj,1)-0.5)*0.15; jitter_cv = (rand(n_subj,1)-0.5)*0.15;
        scatter(ax, ones(n_subj,1)+jitter_mi, mi_v, 30, COL_mi_only, 'filled', 'MarkerFaceAlpha', 0.6);
        scatter(ax, 2*ones(n_subj,1)+jitter_cv, cv_v, 30, COL_cvsa_only, 'filled', 'MarkerFaceAlpha', 0.6);
        for si = 1:n_subj
            if isnan(mi_v(si)) || isnan(cv_v(si)), continue; end
            plot(ax, [1+jitter_mi(si), 2+jitter_cv(si)], [mi_v(si), cv_v(si)], '-', ...
                 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
        plot(ax, 1, grand_mi_v, 'd', 'MarkerFaceColor', COL_mi_only, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.5);
        plot(ax, 2, grand_cv_v, 'd', 'MarkerFaceColor', COL_cvsa_only, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.5);

        yv = [mi_v(:); cv_v(:)];
        yr_top = max([yv(~isnan(yv)); 5], [], 'omitnan') + 12;
        draw_sig_bracket_local(ax, 1, 2, yr_top, sprintf('%s (d=%+.2f)', stars(p_v), dcohen_v));

        set(ax, 'XTick', [1 2], 'XTickLabel', {'MI-only', 'CVSA-only'}, 'XLim', [0.5 2.5], 'YLim', [0, yr_top+8]);
        ylabel(ax, 'rescue rate: % of Hybrid failures that stream alone would have HIT');
        title(ax, sprintf('%s\nGRAND: MI=%s, CVSA=%s  |  delta=%+.1f%% (%s, d=%+.2f)', ...
              resc_panels{pi,1}, fmt_pct(resc_panels{pi,7}), fmt_pct(resc_panels{pi,8}), ...
              mean(d_v,'omitnan')*100, stars(p_v), dcohen_v), 'FontWeight', 'bold', 'FontSize', 9);
        grid(ax, 'on');
    end

    sgtitle(fig2e, sprintf(['Group Analysis — Cost of Fusion: on Hybrid''s own MISS/TIMEOUT trials, would a single stream have hit? (n=%d subjects)\n' ...
            'COUNTERFACTUAL, same trials as Hybrid''s own failures; filled dots = per-subject rescue rate, diamond = GRAND (trial-pooled)\n' ...
            'A non-zero rescue rate means fusion actively cost a win the unimodal stream would have gotten on that exact trial -- the flip side of Fig 07''s keep-up rate'], ...
            sum(has_oc)), 'Interpreter', 'none');

    saveas(fig2e, fullfile(out_dir, sprintf('07b_%s_fusion_cost.svg', group_tag)), 'svg');
    if ~SHOW_FIGURES, close(fig2e); end
end

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

saveas(fig3, fullfile(out_dir, sprintf('08_%s_fusion_mechanism.svg', group_tag)), 'svg');
if ~SHOW_FIGURES, close(fig3); end

% ── Figure 3e: frame-level rescue/hurt/both-ok/both-bad breakdown, 3 panels
%    (within cvsa_influence / whole trial / from cvsa_influence onward) --
%    "on average, in how many samples does the CVSA help vs hurt vs not
%    matter?" Skipped if no subject's hybrid_advantage_probs_summary.mat has
%    been regenerated with the _all/_post window variants yet. ────────────
if any(has_rescue_all)
    COL_both_ok  = [0.55 0.55 0.55];
    COL_both_bad = [0.15 0.15 0.15];
    rescue_cols = {COL_rescue, COL_cost, COL_both_ok, COL_both_bad};
    rescue_labels = {'rescued', 'hurt', 'both correct', 'both wrong'};

    fig3e = figure('Name', 'Group Analysis — Frame-Level Rescue/Hurt Breakdown', 'Color', 'w', ...
                   'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig3e, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    ax3e1 = subplot(3, 1, 1); hold(ax3e1, 'on');
    plot_subject_grand_bars(ax3e1, rescue_frac_inf, grand_rescue_frac_inf, subjects, rescue_labels, rescue_cols, bw, true, true);
    ylabel(ax3e1, '% of frames'); legend(ax3e1, 'Location', 'eastoutside', 'FontSize', 8);
    title(ax3e1, sprintf('Within cvsa\\_influence window (ALL trials)  |  GRAND: rescued=%s, hurt=%s', ...
          fmt_pct(grand_rescue_frac_inf(1)), fmt_pct(grand_rescue_frac_inf(2))), 'FontWeight', 'bold', 'FontSize', 9, 'Interpreter', 'none');
    grid(ax3e1, 'on');

    ax3e2 = subplot(3, 1, 2); hold(ax3e2, 'on');
    plot_subject_grand_bars(ax3e2, rescue_frac_all, grand_rescue_frac_all, subjects, rescue_labels, rescue_cols, bw, true, true);
    ylabel(ax3e2, '% of frames'); legend(ax3e2, 'Location', 'eastoutside', 'FontSize', 8);
    title(ax3e2, sprintf('WHOLE trial (no time restriction)  |  GRAND: rescued=%s, hurt=%s', ...
          fmt_pct(grand_rescue_frac_all(1)), fmt_pct(grand_rescue_frac_all(2))), 'FontWeight', 'bold', 'FontSize', 9, 'Interpreter', 'none');
    grid(ax3e2, 'on');

    ax3e3 = subplot(3, 1, 3); hold(ax3e3, 'on');
    plot_subject_grand_bars(ax3e3, rescue_frac_post, grand_rescue_frac_post, subjects, rescue_labels, rescue_cols, bw, true, true);
    ylabel(ax3e3, '% of frames'); legend(ax3e3, 'Location', 'eastoutside', 'FontSize', 8);
    title(ax3e3, sprintf('From cvsa\\_influence ONWARD, alpha~0  |  NOT all trials contribute (%d/%d overall)  |  GRAND: rescued=%s, hurt=%s', ...
          grand_n_trials_post, grand_n_trials_probs, fmt_pct(grand_rescue_frac_post(1)), fmt_pct(grand_rescue_frac_post(2))), ...
          'FontWeight', 'bold', 'FontSize', 9, 'Interpreter', 'none');
    grid(ax3e3, 'on');

    sgtitle(fig3e, sprintf(['Group Analysis — Frame-Level CVSA Effect Breakdown (n=%d subjects)\n' ...
            'rescued = MI wrong -> fused correct (CVSA helped)  |  hurt = MI correct -> fused wrong (CVSA hurt)\n' ...
            'per-subject %% of classified frames + trial-pooled GRAND; same 3 GDF hybrid files as elsewhere'], ...
            sum(any(~isnan(rescue_frac_inf),2))), 'Interpreter', 'none');

    saveas(fig3e, fullfile(out_dir, sprintf('09_%s_rescue_hurt_breakdown.svg', group_tag)), 'svg');
    if ~SHOW_FIGURES, close(fig3e); end
end

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

    saveas(fig3b, fullfile(out_dir, sprintf('10_%s_fusion_cluster_test.svg', group_tag)), 'svg');
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
    % TWO different p-values shown for the GROUP bar, explicitly labelled --
    % NOT the same test: "naive" is the plain sign-flip test on these 14
    % bars' own cross-subject variability (floor-limited: smallest possible
    % p is 1/2^n_subj here, ~6e-5 at n=14) -- what the error bar visually
    % represents. "meta" is the Stouffer combination of each subject's own
    % well-powered within-session test (hundreds of trials each), which has
    % no such floor. Showing only the meta star next to a bar with a
    % visible error bar (as this figure used to do) looks inconsistent --
    % e.g. "p=0.0000" beside a bar whose own spread looks unremarkable --
    % because the two numbers answer different questions.
    text(ax3c1, n_subj+0.8, m_bw + se_bw*sign(m_bw+eps) + 0.01, ...
         sprintf('naive: p=%.3f %s\nmeta: p=%.4f %s', p_bufadvwin_naive_grp, stars(p_bufadvwin_naive_grp), ...
                  p_bufadvwin_meta_grp, stars(p_bufadvwin_meta_grp)), ...
         'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold');
    set(ax3c1, 'XTick', [1:n_subj, n_subj+0.8], 'XTickLabel', [subjects, {'GROUP'}], 'XLim', [0.4, n_subj+1.2]);
    ylabel(ax3c1, 'mean buffer(target) advantage  [Hybrid - MI-only]');
    title(ax3c1, sprintf(['Per-subject CVSA-help (integrator, cvsa-influence window)\n' ...
          'bar colour = that SUBJECT''S OWN test significant (p<0.05)  |  GROUP bar height/error = simple mean+-SEM of these bars\n' ...
          '"naive" p = sign-flip test ON these bars (floor-limited at n=%d) -- "meta" p = Stouffer combination of each subject''s own well-powered test (no floor)'], ...
          n_subj), 'FontWeight', 'bold', 'FontSize', 9, 'Interpreter', 'none');
    grid(ax3c1, 'on');

    % Panel 2: each subject's OWN within-session p-value (from that
    % subject's own trials, well-powered) -- a DIFFERENT thing from Panel
    % 1's cross-subject GROUP test: this answers "was it significant for
    % just me", Panel 1's GROUP bar answers "is it significant across the
    % whole cohort". Bar height IS the p-value itself (not a delta/effect),
    % so lower bars = more significant, per subject.
    ax3c2 = subplot(1, 2, 2); hold(ax3c2, 'on');
    for si = 1:n_subj
        if isnan(buf_adv_mi_win_meta_p_subj(si)), continue; end
        if buf_adv_mi_win_meta_p_subj(si) < 0.05, col = COL_hybrid; else, col = [0.6 0.6 0.6]; end
        bar(ax3c2, si, buf_adv_mi_win_meta_p_subj(si), 0.5, 'FaceColor', col, 'EdgeColor', 'k');
    end
    yline(ax3c2, 0.05, 'r--', 'p=0.05', 'LabelHorizontalAlignment', 'left');
    set(ax3c2, 'XTick', 1:n_subj, 'XTickLabel', subjects, 'XLim', [0.4, n_subj+0.6]);
    ylabel(ax3c2, 'p-value (that subject''s own within-session test; lower = more significant)');
    title(ax3c2, ['Is CVSA help significant for THIS subject ALONE?' newline ...
          '(each subject''s own well-powered test on their own trials -- not the cross-subject GROUP test in Panel 1)'], ...
          'FontWeight', 'bold', 'FontSize', 9);
    grid(ax3c2, 'on');

    sgtitle(fig3c, sprintf(['Group Analysis — CVSA-Help Significance, Integrator Buffer (n=%d subjects)\n' ...
            'Hybrid vs MI-only, restricted to the first cvsa_influence seconds of CF -- isolates CVSA''s actual influence window'], ...
            sum(has_bw)), 'Interpreter', 'none');

    saveas(fig3c, fullfile(out_dir, sprintf('11_%s_cvsa_help_integrator.svg', group_tag)), 'svg');
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

saveas(fig4, fullfile(out_dir, sprintf('12_%s_erd_across_subjects.svg', group_tag)), 'svg');
if ~SHOW_FIGURES, close(fig4); end

% ── Figure 5: classifier-level ROC across subjects, one panel per paradigm ─
paradigms_with_roc = pars(any(~isnan(roc_auc_subj), 1));
if ~isempty(paradigms_with_roc)
    fig5 = figure('Name', 'Group Analysis — ROC Across Subjects', 'Color', 'w', ...
                  'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig5, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    % One distinct (solid-line-only) colour per subject via an evenly-spaced
    % hue sweep -- no per-subject dash/dot styles. The per-subject legend
    % now lives in its own dedicated tile OUTSIDE the ROC axes (tile 5,
    % below), so there is no in-plot clutter cost to a unique-per-subject
    % colour even when cohort size makes some hues only subtly different;
    % the legend (still readable, just off to the side) is what actually
    % identifies each curve.
    subj_roc_colors = hsv(n_subj);

    % Layout: 3 ROC panels + 1 pairwise-AUC panel + a 5th tile dedicated
    % ONLY to the per-subject legend -- tiledlayout's Layout.Tile placement
    % puts it fully OUTSIDE every plot axes (not overlapping the curves,
    % and not eating into panel 4), which a subplot()-based 'Location'
    % legend cannot do. Subject curves themselves never carry their own
    % per-axes legend entry (HandleVisibility off); their handles are
    % collected once, from the last populated ROC panel, and legended
    % explicitly into tile 5 below.
    tl5 = tiledlayout(fig5, 1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
    last_roc_pi = find(any(~isnan(roc_auc_subj), 1), 1, 'last');
    subj_legend_h = gobjects(0);
    subj_legend_labels = {};
    for pi = 1:3
        if all(isnan(roc_auc_subj(:,pi))), continue; end
        axp = nexttile(tl5, pi); hold(axp, 'on');
        plot(axp, [0 1], [0 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
        show_subj_legend = (pi == last_roc_pi);
        for si = 1:n_subj
            curve = reshape(roc_tpr_subj_grid(si,pi,:), 1, numel(ROC_FPR_GRID));
            if all(isnan(curve)), continue; end
            hplot = plot(axp, ROC_FPR_GRID, curve, 'Color', subj_roc_colors(si,:), ...
                         'LineWidth', 1.2, 'HandleVisibility', 'off');
            if show_subj_legend
                subj_legend_h(end+1)      = hplot;        %#ok<AGROW>
                subj_legend_labels{end+1} = subjects{si};  %#ok<AGROW>
            end
        end
        m  = roc_mean_tpr(pi,:);
        se = roc_sem_tpr(pi,:);
        fill(axp, [ROC_FPR_GRID, fliplr(ROC_FPR_GRID)], [m+se, fliplr(m-se)], cols3{pi}, ...
             'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(axp, ROC_FPR_GRID, m, '-', 'Color', cols3{pi}, 'LineWidth', 3.5, ...
             'DisplayName', 'cross-subject mean \pm SEM');
        xlabel(axp, 'False Positive Rate'); ylabel(axp, 'True Positive Rate');
        xlim(axp, [0 1]); ylim(axp, [0 1]); axis(axp, 'square'); grid(axp, 'on');
        n_v = sum(~isnan(roc_auc_subj(:,pi)));
        title(axp, sprintf('%s  (n=%d subjects)\nmean AUC=%.3f  (pairwise comparisons in panel 4)', ...
              par_labels{pi}, n_v, mean(roc_auc_subj(:,pi),'omitnan')), ...
              'FontWeight', 'bold', 'FontSize', 9);
        legend(axp, 'Location', 'southeast', 'FontSize', 8);   % just the mean +/- SEM entry
    end

    % 4th panel: direct, PAIRWISE AUC comparison across paradigms -- Hybrid
    % vs MI, Hybrid vs CVSA, and MI vs CVSA -- so "which one discriminates
    % best" is readable at a glance instead of comparing 3 separate panel
    % titles. No vs-chance test here (that answer is implicit: all 3 AUCs
    % are visibly >>0.5); the only questions this panel answers are the 3
    % pairwise ones. For Hybrid, the "classifier" being scored is the
    % cosine-annealed Bayesian-fused P(c) pooled across all CF frames (see
    % main_roc_analysis) -- a moving-weight blend of MI/CVSA, not a single
    % fixed classifier like the MI-only/CVSA-only panels -- so its AUC is a
    % genuine but time-averaged summary, called out in the title. All 3
    % brackets are the SAME test: paired two-sided sign-flip permutation on
    % the per-subject AUC delta, subject is the unit, + Cohen's d.
    ax_auc = nexttile(tl5, 4); hold(ax_auc, 'on');
    for pi = 1:3
        v = roc_auc_subj(:,pi);
        m = mean(v, 'omitnan'); se = std(v, 'omitnan') / sqrt(max(1, sum(~isnan(v))));
        bar(ax_auc, pi, m, 0.6, 'FaceColor', cols3{pi}, 'EdgeColor', 'k', 'LineWidth', 1);
        errorbar(ax_auc, pi, m, se, 'k', 'LineWidth', 1.2, 'CapSize', 6);
        text(ax_auc, pi, m + se + 0.015, sprintf('%.3f', m), ...
             'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
    end
    y_top = min(1, max(roc_auc_subj(:), [], 'omitnan') + 0.22);
    draw_sig_bracket_local(ax_auc, 1, 3, y_top,        sprintf('Hybrid vs MI: %s (d=%+.2f)', stars(p_auc_hyb_mi), d_auc_hyb_mi_cohen));
    draw_sig_bracket_local(ax_auc, 2, 3, y_top - 0.07, sprintf('Hybrid vs CVSA: %s (d=%+.2f)', stars(p_auc_hyb_cvs), d_auc_hyb_cvs_cohen));
    draw_sig_bracket_local(ax_auc, 1, 2, y_top - 0.14, sprintf('MI vs CVSA: %s (d=%+.2f)', stars(p_auc_mi_cvs), d_auc_mi_cvs_cohen));
    set(ax_auc, 'XTick', 1:3, 'XTickLabel', par_labels, 'YLim', [0.4, y_top+0.10]);
    ylabel(ax_auc, 'AUC');
    title(ax_auc, sprintf(['Pairwise AUC comparison (n=%d subjects)\n' ...
          'test: paired sign-flip permutation on per-subject AUC delta, subject=unit, + Cohen''s d\n' ...
          '(Hybrid AUC = time-averaged, annealed MI/CVSA blend, not a single fixed classifier)'], ...
          n_subj), ...
          'FontWeight', 'bold', 'FontSize', 8);
    grid(ax_auc, 'on');

    if ~isempty(subj_legend_h)
        lgd5 = legend(subj_legend_h, subj_legend_labels, 'FontSize', 7, 'NumColumns', 2);
        lgd5.Layout.Tile = 5;
        lgd5.Title.String = 'Subjects';
    end

    sgtitle(fig5, sprintf(['Group Analysis — Classifier ROC/AUC, pre-integrator, ALL CF frames (n=%d subjects)\n' ...
            'thin = per-subject re-pooled curve, thick = cross-subject macro-average\n' ...
            'Windowed (cvsa_influence-only) version: 14a_%s_roc_window.svg'], n_subj, group_tag), ...
            'Interpreter', 'none');

    saveas(fig5, fullfile(out_dir, sprintf('13a_%s_roc_analysis.svg', group_tag)), 'svg');
    if ~SHOW_FIGURES, close(fig5); end

    % ── Figure 5b (13b): classifier calibration (Brier + reliability diagram)
    %    and sensitivity (d-prime) -- complements AUC/ROC, which only tests
    %    whether scores RANK classes correctly, not whether the probability
    %    VALUES themselves are trustworthy. Calibration matters here because
    %    both the integrator (threshold crossing) and the Bayesian fusion
    %    (P_CVSA^alpha) use the actual probability value, not just its order. ─
    fig5b = figure('Name', 'Group Analysis — Calibration & Sensitivity', 'Color', 'w', ...
                   'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig5b, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    ax5b1 = subplot(1, 3, 1); hold(ax5b1, 'on');
    for pi = 1:3
        v = brier_subj(:,pi);
        m = mean(v, 'omitnan'); se = std(v, 'omitnan') / sqrt(max(1, sum(~isnan(v))));
        bar(ax5b1, pi, m, 0.6, 'FaceColor', cols3{pi}, 'EdgeColor', 'k', 'LineWidth', 1);
        errorbar(ax5b1, pi, m, se, 'k', 'LineWidth', 1.2, 'CapSize', 6);
        text(ax5b1, pi, m + se + 0.008, sprintf('%.3f', m), 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
    end
    yb_top = min(0.5, max(brier_subj(:), [], 'omitnan') + 0.10);
    draw_sig_bracket_local(ax5b1, 1, 3, yb_top,        sprintf('MI vs Hybrid: %s (d=%+.2f)', stars(p_brier_hyb_mi), d_brier_hyb_mi_cohen));
    draw_sig_bracket_local(ax5b1, 2, 3, yb_top - 0.035, sprintf('CVSA vs Hybrid: %s (d=%+.2f)', stars(p_brier_hyb_cvs), d_brier_hyb_cvs_cohen));
    draw_sig_bracket_local(ax5b1, 1, 2, yb_top - 0.07,  sprintf('MI vs CVSA: %s (d=%+.2f)', stars(p_brier_mi_cvs), d_brier_mi_cvs_cohen));
    set(ax5b1, 'XTick', 1:3, 'XTickLabel', par_labels, 'YLim', [0, yb_top+0.05]);
    ylabel(ax5b1, 'Brier score (LOWER = better calibrated)');
    title(ax5b1, sprintf('Calibration (n=%d subjects)\nmean((P(target)-label)^2), pooled per-subject frames', n_subj), 'FontWeight', 'bold', 'FontSize', 8);
    grid(ax5b1, 'on');

    ax5b2 = subplot(1, 3, 2); hold(ax5b2, 'on');
    for pi = 1:3
        v = dprime_subj(:,pi);
        m = mean(v, 'omitnan'); se = std(v, 'omitnan') / sqrt(max(1, sum(~isnan(v))));
        bar(ax5b2, pi, m, 0.6, 'FaceColor', cols3{pi}, 'EdgeColor', 'k', 'LineWidth', 1);
        errorbar(ax5b2, pi, m, se, 'k', 'LineWidth', 1.2, 'CapSize', 6);
        text(ax5b2, pi, m + se + 0.05, sprintf('%.2f', m), 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
    end
    yd_top = max(dprime_subj(:), [], 'omitnan') + 0.9;
    draw_sig_bracket_local(ax5b2, 1, 3, yd_top,       sprintf('MI vs Hybrid: %s (d=%+.2f)', stars(p_dprime_hyb_mi), d_dprime_hyb_mi_cohen));
    draw_sig_bracket_local(ax5b2, 2, 3, yd_top - 0.3, sprintf('CVSA vs Hybrid: %s (d=%+.2f)', stars(p_dprime_hyb_cvs), d_dprime_hyb_cvs_cohen));
    draw_sig_bracket_local(ax5b2, 1, 2, yd_top - 0.6, sprintf('MI vs CVSA: %s (d=%+.2f)', stars(p_dprime_mi_cvs), d_dprime_mi_cvs_cohen));
    set(ax5b2, 'XTick', 1:3, 'XTickLabel', par_labels, 'YLim', [min(0,min(dprime_subj(:),[],'omitnan')-0.2), yd_top+0.4]);
    yline(ax5b2, 0, 'k:', 'HandleVisibility', 'off');
    ylabel(ax5b2, 'd-prime (HIGHER = better separated)');
    title(ax5b2, 'Sensitivity index at P=0.5 (classic signal-detection d-prime)', 'FontWeight', 'bold', 'FontSize', 8);
    grid(ax5b2, 'on');

    ax5b3 = subplot(1, 3, 3); hold(ax5b3, 'on');
    plot(ax5b3, [0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', 'perfect calibration');
    for pi = 1:3
        conf_mean = squeeze(mean(calib_conf_subj(:,pi,:), 1, 'omitnan'))';
        obs_mean  = squeeze(mean(calib_obs_subj(:,pi,:),  1, 'omitnan'))';
        n_contrib = squeeze(sum(~isnan(calib_obs_subj(:,pi,:)), 1))';
        ok = n_contrib > 0;
        plot(ax5b3, conf_mean(ok), obs_mean(ok), '-o', 'Color', cols3{pi}, 'MarkerFaceColor', cols3{pi}, ...
             'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', par_labels{pi});
    end
    xlabel(ax5b3, 'mean predicted P(target) in bin ("confidence")');
    ylabel(ax5b3, 'observed fraction target=1 ("accuracy")');
    xlim(ax5b3, [0 1]); ylim(ax5b3, [0 1]); axis(ax5b3, 'square');
    legend(ax5b3, 'Location', 'best', 'FontSize', 8);
    title(ax5b3, sprintf('Reliability diagram (%d bins, cross-subject mean)\nabove diagonal = under-confident, below = over-confident', CALIB_NBINS), 'FontWeight', 'bold', 'FontSize', 8);
    grid(ax5b3, 'on');

    sgtitle(fig5b, sprintf(['Group Analysis — Classifier Calibration & Sensitivity, ALL CF frames (n=%d subjects)\n' ...
            'Complements 13a_%s_roc_analysis.svg: AUC tests RANKING only, not whether probability VALUES are trustworthy\n' ...
            'Both the integrator threshold and the Bayesian fusion P_CVSA^alpha use the actual value\n' ...
            'Windowed (cvsa_influence-only) version: 14b_%s_calibration_dprime_window.svg'], n_subj, group_tag, group_tag), 'Interpreter', 'none');

    saveas(fig5b, fullfile(out_dir, sprintf('13b_%s_calibration_dprime.svg', group_tag)), 'svg');
    if ~SHOW_FIGURES, close(fig5b); end

    % ── Figure 14a: WINDOWED classifier ROC/AUC -- exact mirror of Fig 13a
    %    (13/fig5 above), but restricted to the first cvsa_influence seconds
    %    of each trial's CF, the ONLY window where a fused/Hybrid classifier
    %    can actually differ from pure MI (bayesian_fuse.m forces alpha=0 --
    %    fused==MI exactly -- from cvsa_influence onward, so the ALL-frames
    %    AUC in 13a/13b dilutes any early advantage with that later, by-
    %    construction-identical portion of every trial). Skipped if no
    %    subject's roc_summary.mat has scores_win yet. ─────────────────────
    if any(has_roc_win)
        fig14a = figure('Name', 'Group Analysis — Windowed ROC Across Subjects', 'Color', 'w', ...
                       'NumberTitle', 'off', 'Visible', fig_vis);
        set(fig14a, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

        tl14a = tiledlayout(fig14a, 1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
        last_roc_pi_w = find(any(~isnan(roc_auc_win_subj), 1), 1, 'last');
        subj_legend_h_w = gobjects(0);
        subj_legend_labels_w = {};
        for pi = 1:3
            if all(isnan(roc_auc_win_subj(:,pi))), continue; end
            axp = nexttile(tl14a, pi); hold(axp, 'on');
            plot(axp, [0 1], [0 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
            show_subj_legend = (pi == last_roc_pi_w);
            for si = 1:n_subj
                curve = reshape(roc_tpr_win_subj_grid(si,pi,:), 1, numel(ROC_FPR_GRID));
                if all(isnan(curve)), continue; end
                hplot = plot(axp, ROC_FPR_GRID, curve, 'Color', subj_roc_colors(si,:), ...
                             'LineWidth', 1.2, 'HandleVisibility', 'off');
                if show_subj_legend
                    subj_legend_h_w(end+1)      = hplot;      %#ok<AGROW>
                    subj_legend_labels_w{end+1} = subjects{si}; %#ok<AGROW>
                end
            end
            m  = roc_mean_tpr_win(pi,:);
            se = roc_sem_tpr_win(pi,:);
            fill(axp, [ROC_FPR_GRID, fliplr(ROC_FPR_GRID)], [m+se, fliplr(m-se)], cols3{pi}, ...
                 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(axp, ROC_FPR_GRID, m, '-', 'Color', cols3{pi}, 'LineWidth', 3.5, ...
                 'DisplayName', 'cross-subject mean \pm SEM');
            xlabel(axp, 'False Positive Rate'); ylabel(axp, 'True Positive Rate');
            xlim(axp, [0 1]); ylim(axp, [0 1]); axis(axp, 'square'); grid(axp, 'on');
            n_v = sum(~isnan(roc_auc_win_subj(:,pi)));
            title(axp, sprintf('%s WINDOW ONLY (n=%d subj)\nmean AUC=%.3f', ...
                  par_labels{pi}, n_v, mean(roc_auc_win_subj(:,pi),'omitnan')), ...
                  'FontWeight', 'bold', 'FontSize', 9);
            legend(axp, 'Location', 'southeast', 'FontSize', 8);   % just the mean +/- SEM entry
        end

        ax14a = nexttile(tl14a, 4); hold(ax14a, 'on');
        for pi = 1:3
            v = roc_auc_win_subj(:,pi);
            m = mean(v, 'omitnan'); se = std(v, 'omitnan') / sqrt(max(1, sum(~isnan(v))));
            bar(ax14a, pi, m, 0.6, 'FaceColor', cols3{pi}, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar(ax14a, pi, m, se, 'k', 'LineWidth', 1.2, 'CapSize', 6);
            text(ax14a, pi, m + se + 0.015, sprintf('%.3f', m), 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
        ya_top = min(1, max(roc_auc_win_subj(:), [], 'omitnan') + 0.22);
        draw_sig_bracket_local(ax14a, 1, 3, ya_top,        sprintf('Hybrid vs MI: %s (d=%+.2f)', stars(p_aucw_hyb_mi), d_aucw_hyb_mi_cohen));
        draw_sig_bracket_local(ax14a, 2, 3, ya_top - 0.07, sprintf('Hybrid vs CVSA: %s (d=%+.2f)', stars(p_aucw_hyb_cvs), d_aucw_hyb_cvs_cohen));
        draw_sig_bracket_local(ax14a, 1, 2, ya_top - 0.14, sprintf('MI vs CVSA: %s (d=%+.2f)', stars(p_aucw_mi_cvs), d_aucw_mi_cvs_cohen));
        set(ax14a, 'XTick', 1:3, 'XTickLabel', par_labels, 'YLim', [0.4, ya_top+0.10]);
        ylabel(ax14a, 'AUC');
        title(ax14a, sprintf('Pairwise AUC, WINDOW ONLY (n=%d subj)', n_subj), 'FontWeight', 'bold', 'FontSize', 9);
        grid(ax14a, 'on');

        if ~isempty(subj_legend_h_w)
            lgd14a = legend(subj_legend_h_w, subj_legend_labels_w, 'FontSize', 7, 'NumColumns', 2);
            lgd14a.Layout.Tile = 5;
            lgd14a.Title.String = 'Subjects';
        end

        sgtitle(fig14a, sprintf(['Group Analysis — Windowed Classifier ROC/AUC: first %.1fs of CF ONLY (n=%d subjects)\n' ...
                'Isolates the ONLY window where P_fused can differ from P_MI\n' ...
                '(bayesian_fuse.m: alpha=0, fused==MI exactly, from cvsa_influence onward)\n' ...
                'ALL-CF-frames version: 13a_%s_roc_analysis.svg'], ...
                mean_window_s, n_subj, group_tag), 'Interpreter', 'none');

        saveas(fig14a, fullfile(out_dir, sprintf('14a_%s_roc_window.svg', group_tag)), 'svg');
        if ~SHOW_FIGURES, close(fig14a); end

        % ── Figure 14b: WINDOWED calibration & sensitivity -- exact mirror
        %    of Fig 13b, same window as 14a. ─────────────────────────────
        fig14b = figure('Name', 'Group Analysis — Windowed Calibration & Sensitivity', 'Color', 'w', ...
                       'NumberTitle', 'off', 'Visible', fig_vis);
        set(fig14b, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

        ax14b1 = subplot(1, 3, 1); hold(ax14b1, 'on');
        for pi = 1:3
            v = brier_win_subj(:,pi);
            m = mean(v, 'omitnan'); se = std(v, 'omitnan') / sqrt(max(1, sum(~isnan(v))));
            bar(ax14b1, pi, m, 0.6, 'FaceColor', cols3{pi}, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar(ax14b1, pi, m, se, 'k', 'LineWidth', 1.2, 'CapSize', 6);
            text(ax14b1, pi, m + se + 0.008, sprintf('%.3f', m), 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
        yb_top_w = min(0.5, max(brier_win_subj(:), [], 'omitnan') + 0.10);
        draw_sig_bracket_local(ax14b1, 1, 3, yb_top_w,        sprintf('Hybrid vs MI: %s (d=%+.2f)', stars(p_brierw_hyb_mi), d_brierw_hyb_mi_cohen));
        draw_sig_bracket_local(ax14b1, 2, 3, yb_top_w - 0.035, sprintf('Hybrid vs CVSA: %s (d=%+.2f)', stars(p_brierw_hyb_cvs), d_brierw_hyb_cvs_cohen));
        draw_sig_bracket_local(ax14b1, 1, 2, yb_top_w - 0.07,  sprintf('MI vs CVSA: %s (d=%+.2f)', stars(p_brierw_mi_cvs), d_brierw_mi_cvs_cohen));
        set(ax14b1, 'XTick', 1:3, 'XTickLabel', par_labels, 'YLim', [0, yb_top_w+0.05]);
        ylabel(ax14b1, 'Brier score (LOWER = better calibrated)');
        title(ax14b1, sprintf('Calibration, WINDOW ONLY (n=%d subj)', n_subj), 'FontWeight', 'bold', 'FontSize', 8);
        grid(ax14b1, 'on');

        ax14b2 = subplot(1, 3, 2); hold(ax14b2, 'on');
        for pi = 1:3
            v = dprime_win_subj(:,pi);
            m = mean(v, 'omitnan'); se = std(v, 'omitnan') / sqrt(max(1, sum(~isnan(v))));
            bar(ax14b2, pi, m, 0.6, 'FaceColor', cols3{pi}, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar(ax14b2, pi, m, se, 'k', 'LineWidth', 1.2, 'CapSize', 6);
            text(ax14b2, pi, m + se + 0.05, sprintf('%.2f', m), 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
        yd_top_w = max(dprime_win_subj(:), [], 'omitnan') + 0.9;
        draw_sig_bracket_local(ax14b2, 1, 3, yd_top_w,       sprintf('Hybrid vs MI: %s (d=%+.2f)', stars(p_dprimew_hyb_mi), d_dprimew_hyb_mi_cohen));
        draw_sig_bracket_local(ax14b2, 2, 3, yd_top_w - 0.3, sprintf('Hybrid vs CVSA: %s (d=%+.2f)', stars(p_dprimew_hyb_cvs), d_dprimew_hyb_cvs_cohen));
        draw_sig_bracket_local(ax14b2, 1, 2, yd_top_w - 0.6, sprintf('MI vs CVSA: %s (d=%+.2f)', stars(p_dprimew_mi_cvs), d_dprimew_mi_cvs_cohen));
        set(ax14b2, 'XTick', 1:3, 'XTickLabel', par_labels, 'YLim', [min(0,min(dprime_win_subj(:),[],'omitnan')-0.2), yd_top_w+0.4]);
        yline(ax14b2, 0, 'k:', 'HandleVisibility', 'off');
        ylabel(ax14b2, 'd-prime (HIGHER = better separated)');
        title(ax14b2, sprintf('Sensitivity, WINDOW ONLY (n=%d subj)', n_subj), 'FontWeight', 'bold', 'FontSize', 8);
        grid(ax14b2, 'on');

        ax14b3 = subplot(1, 3, 3); hold(ax14b3, 'on');
        plot(ax14b3, [0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', 'perfect calibration');
        for pi = 1:3
            conf_mean_w = squeeze(mean(calib_conf_win_subj(:,pi,:), 1, 'omitnan'))';
            obs_mean_w  = squeeze(mean(calib_obs_win_subj(:,pi,:),  1, 'omitnan'))';
            n_contrib_w = squeeze(sum(~isnan(calib_obs_win_subj(:,pi,:)), 1))';
            ok_w = n_contrib_w > 0;
            plot(ax14b3, conf_mean_w(ok_w), obs_mean_w(ok_w), '-o', 'Color', cols3{pi}, 'MarkerFaceColor', cols3{pi}, ...
                 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', par_labels{pi});
        end
        xlabel(ax14b3, 'mean predicted P(target) in bin ("confidence")');
        ylabel(ax14b3, 'observed fraction target=1 ("accuracy")');
        xlim(ax14b3, [0 1]); ylim(ax14b3, [0 1]); axis(ax14b3, 'square');
        legend(ax14b3, 'Location', 'best', 'FontSize', 8);
        title(ax14b3, sprintf('Reliability, WINDOW ONLY (%d bins)', CALIB_NBINS), 'FontWeight', 'bold', 'FontSize', 8);
        grid(ax14b3, 'on');

        sgtitle(fig14b, sprintf(['Group Analysis — Windowed Calibration & Sensitivity: first %.1fs of CF ONLY (n=%d subjects)\n' ...
                'Same as 13b_%s_calibration_dprime.svg, restricted to frames before cvsa_influence elapses\n' ...
                'ALL-CF-frames version: 13b_%s_calibration_dprime.svg'], ...
                mean_window_s, n_subj, group_tag, group_tag), 'Interpreter', 'none');

        saveas(fig14b, fullfile(out_dir, sprintf('14b_%s_calibration_dprime_window.svg', group_tag)), 'svg');
        if ~SHOW_FIGURES, close(fig14b); end
    end
else
    fprintf('\n[main_group_analysis] No ROC data found (run main_roc_analysis first) — skipping Fig 5.\n');
end

% ── Figure 6: cross-subject correlates of real Hybrid performance (bonus) ──
% Five scatter panels: does a subject's counterfactual advantage / fusion
% mechanism / neurophysiological grounding / CVSA-only "quality" / raw
% signal quality predict how well they actually do in real Hybrid sessions?
% The 4th panel directly supports the "CVSA helps even when imperfect"
% thesis: if the Hybrid-MI advantage stays positive even for subjects with
% weak CVSA-only accuracy, adding CVSA is not merely helping when it
% happens to be good. The 5th panel is a confound check: is a subject's
% real accuracy actually explained by how clean their raw EEG was (dead/
% noisy channel rate during CF, main_session_overview.m Fig 11), rather
% than by the BCI method itself?
fig6 = figure('Name', 'Group Analysis — Performance Correlates', 'Color', 'w', ...
              'NumberTitle', 'off', 'Visible', fig_vis);
set(fig6, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

ax_c1 = subplot(2, 3, 1);
scatter_with_fit(ax_c1, 100*d_mi_subj, 100*d_real_hyb_mi, subjects, r_cf_vs_real, p_cf_vs_real, ...
    '\Delta acc, counterfactual (Hybrid-MI, %)', '\Delta acc, real (Hybrid-MI, %)', ...
    'Counterfactual advantage predicts real advantage?', COL_hybrid);

ax_c2 = subplot(2, 3, 2);
scatter_with_fit(ax_c2, fus_adv_subj, 100*real_acc(:,3), subjects, r_fus_vs_real, p_fus_vs_real, ...
    'fusion advantage  mean(P_{fused}-P_{MI})', 'real Hybrid HIT rate (%)', ...
    'Fusion mechanism predicts real accuracy?', COL_hybrid);

ax_c3 = subplot(2, 3, 3);
scatter_with_fit(ax_c3, erd_r_subj(:,3), 100*real_acc(:,3), subjects, r_erd_vs_real, p_erd_vs_real, ...
    'Hybrid ERD/ERS vs CSP-weight (Pearson r)', 'real Hybrid HIT rate (%)', ...
    'Neurophysiological grounding predicts real accuracy?', COL_hybrid);

ax_c4 = subplot(2, 3, 4);
scatter_with_fit(ax_c4, 100*cf_acc_dec(:,3), 100*d_mi_subj_dec, subjects, r_cvsa_quality, p_cvsa_quality, ...
    'CVSA-only decided-trial accuracy (%)', '\Delta acc, counterfactual (Hybrid-MI, %)', ...
    'Does the Hybrid advantage survive weak CVSA?', COL_cvsa_only);
yline(ax_c4, 0, 'k-', 'HandleVisibility', 'off');

ax_c5 = subplot(2, 3, 5);
scatter_with_fit(ax_c5, 100*sigqc_rate_overall, 100*real_acc(:,3), subjects, r_sigqc_vs_real, p_sigqc_vs_real, ...
    '% trials with dead/noisy channel during CF', 'real Hybrid HIT rate (%)', ...
    'Confound check: is accuracy just raw signal quality?', [0.55 0.35 0.15]);

sgtitle(fig6, sprintf(['Group Analysis — Cross-Subject Correlates of Real Hybrid Performance (n=%d subjects)\n' ...
        'each point = one subject; permutation test on Pearson r'], n_subj), 'Interpreter', 'none');

saveas(fig6, fullfile(out_dir, sprintf('15_%s_performance_correlates.svg', group_tag)), 'svg');
if ~SHOW_FIGURES, close(fig6); end

% ── Figure 17: EXECUTIVE VERDICT -- is Hybrid actually worth it? ───────────
% Consolidates every Hybrid-vs-MI-only / Hybrid-vs-CVSA-only paired test
% already computed above (real + counterfactual accuracy in both
% conventions, speed/TTH, ITR, classifier AUC) into ONE scorecard, instead
% of needing to piece the answer together across ~10 separate figures.
% Green = Hybrid significantly better (p<0.05, correct direction -- lower
% is better for TTH, higher for everything else); red = Hybrid
% significantly WORSE; grey = not significant. Every cell already used the
% SAME paired sign-flip test + Cohen's d as its source figure -- nothing
% new is computed here, this is purely a synthesis/aggregation.
if exist('d_auc_hyb_mi', 'var')
    auc_mi_d = d_auc_hyb_mi;   auc_mi_p = p_auc_hyb_mi;   auc_mi_dcohen = d_auc_hyb_mi_cohen;
    auc_cv_d = d_auc_hyb_cvs;  auc_cv_p = p_auc_hyb_cvs;  auc_cv_dcohen = d_auc_hyb_cvs_cohen;
else
    auc_mi_d = NaN; auc_mi_p = NaN; auc_mi_dcohen = NaN;
    auc_cv_d = NaN; auc_cv_p = NaN; auc_cv_dcohen = NaN;
end

verdict_rows = {
    'Real accuracy (TIMEOUT=fail)',  100*d_real_hyb_mi,     p_real_hyb_mi,     d_real_hyb_mi_cohen,     100*d_real_hyb_cvs,     p_real_hyb_cvs,     d_real_hyb_cvs_cohen,     true,  '%+.1f%%';
    'Real accuracy (decided)',       100*d_real_hyb_mi_dec, p_real_hyb_mi_dec, d_real_hyb_mi_dec_cohen, 100*d_real_hyb_cvs_dec, p_real_hyb_cvs_dec, d_real_hyb_cvs_dec_cohen, true,  '%+.1f%%';
    'Counterfactual acc (TIMEOUT=fail)', 100*d_mi_subj,     p_mi_grp,          d_mi_grp_cohen,          100*d_cvs_subj,         p_cvs_grp,          d_cvs_grp_cohen,          true,  '%+.1f%%';
    'Counterfactual acc (decided)',  100*d_mi_subj_dec,     p_mi_grp_dec,      d_mi_grp_dec_cohen,      100*d_cvs_subj_dec,     p_cvs_grp_dec,      d_cvs_grp_dec_cohen,      true,  '%+.1f%%';
    'Speed / TTH, real (s)',         d_tth_hyb_mi,          p_tth_hyb_mi,      d_tth_hyb_mi_cohen,      d_tth_hyb_cvs,          p_tth_hyb_cvs,      d_tth_hyb_cvs_cohen,      false, '%+.2fs';
    'Speed / TTH, counterfactual (s)', d_cftth_hyb_mi,      p_cftth_hyb_mi,    d_cftth_hyb_mi_cohen,    d_cftth_hyb_cvs,        p_cftth_hyb_cvs,    d_cftth_hyb_cvs_cohen,    false, '%+.2fs';
    'ITR, real (bits/min)',          d_itr_hyb_mi,          p_itr_hyb_mi,      d_itr_hyb_mi_cohen,      d_itr_hyb_cvs,          p_itr_hyb_cvs,      d_itr_hyb_cvs_cohen,      true,  '%+.2f';
    'ITR, counterfactual (bits/min)', d_cfitr_hyb_mi,       p_cfitr_hyb_mi,    d_cfitr_hyb_mi_cohen,    d_cfitr_hyb_cvs,        p_cfitr_hyb_cvs,    d_cfitr_hyb_cvs_cohen,    true,  '%+.2f';
    'Classifier ROC AUC',            auc_mi_d,              auc_mi_p,          auc_mi_dcohen,           auc_cv_d,               auc_cv_p,           auc_cv_dcohen,            true,  '%+.3f';
};
n_vrows = size(verdict_rows, 1);

verdict_code = nan(n_vrows, 2);   % col 1 = vs MI-only, col 2 = vs CVSA-only
verdict_txt  = cell(n_vrows, 2);
for vr = 1:n_vrows
    [verdict_code(vr,1), verdict_txt{vr,1}] = verdict_cell_local(verdict_rows{vr,2}, verdict_rows{vr,3}, verdict_rows{vr,4}, verdict_rows{vr,8}, verdict_rows{vr,9});
    [verdict_code(vr,2), verdict_txt{vr,2}] = verdict_cell_local(verdict_rows{vr,5}, verdict_rows{vr,6}, verdict_rows{vr,7}, verdict_rows{vr,8}, verdict_rows{vr,9});
end

fig17 = figure('Name', 'Group Analysis — Executive Verdict', 'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig17, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
ax17 = axes(fig17);
verdict_cmap = [0.85 0.30 0.25; 0.85 0.85 0.85; 0.30 0.65 0.30];   % red / grey / green, for codes -1/0/+1
imagesc(ax17, verdict_code, [-1 1]);
colormap(ax17, verdict_cmap);
set(ax17, 'XTick', 1:2, 'XTickLabel', {'Hybrid vs MI-only', 'Hybrid vs CVSA-only'}, 'FontSize', 10, ...
          'YTick', 1:n_vrows, 'YTickLabel', verdict_rows(:,1), 'TickLabelInterpreter', 'none');
for vr = 1:n_vrows
    for vc = 1:2
        txt_col = 'k';
        text(ax17, vc, vr, verdict_txt{vr,vc}, 'HorizontalAlignment', 'center', ...
             'FontSize', 9, 'Color', txt_col, 'FontWeight', 'bold');
    end
end

n_win_mi  = sum(verdict_code(:,1) == 1); n_lose_mi  = sum(verdict_code(:,1) == -1); n_ns_mi  = sum(verdict_code(:,1) == 0);
n_win_cv  = sum(verdict_code(:,2) == 1); n_lose_cv  = sum(verdict_code(:,2) == -1); n_ns_cv  = sum(verdict_code(:,2) == 0);
fprintf('\n══════ EXECUTIVE VERDICT: is Hybrid worth it? (n=%d subjects) ══════\n', n_subj);
fprintf('  vs MI-only  : significantly BETTER on %d/%d metrics, WORSE on %d/%d, n.s. on %d/%d\n', ...
        n_win_mi, n_vrows, n_lose_mi, n_vrows, n_ns_mi, n_vrows);
fprintf('  vs CVSA-only: significantly BETTER on %d/%d metrics, WORSE on %d/%d, n.s. on %d/%d\n', ...
        n_win_cv, n_vrows, n_lose_cv, n_vrows, n_ns_cv, n_vrows);
fprintf('  (green = Hybrid significantly better, p<0.05 correct-direction paired sign-flip test; red = significantly worse; grey = n.s.)\n');
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

title(ax17, sprintf(['Group Analysis — Executive Verdict: is Hybrid worth it? (n=%d subjects)\n' ...
        'vs MI-only: %d/%d metrics significantly better, %d/%d worse, %d/%d n.s.  |  ' ...
        'vs CVSA-only: %d/%d better, %d/%d worse, %d/%d n.s.\n' ...
        'every cell reuses the paired sign-flip test + Cohen''s d already computed for its own source figure above -- nothing new computed here'], ...
        n_subj, n_win_mi, n_vrows, n_lose_mi, n_vrows, n_ns_mi, n_vrows, n_win_cv, n_vrows, n_lose_cv, n_vrows, n_ns_cv, n_vrows), ...
        'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'none');

saveas(fig17, fullfile(out_dir, sprintf('17_%s_verdict_scorecard.svg', group_tag)), 'svg');
if ~SHOW_FIGURES, close(fig17); end

% ── Figure 18: does CALIBRATION-time signal quality predict how much a
%    subject benefits from Hybrid? "Is the group Hybrid helps a special
%    subgroup -- specifically the subjects who already had weak/noisy
%    signal at calibration time (BEFORE any evaluation/Hybrid session even
%    happened)?" Uses topo_erders.m's own ERD/ERS class-discrimination
%    metric (|ERD(c1)-ERD(c2)|, averaged over that subject's CSP-selected
%    channels/bands), computed from the CALIBRATION folder specifically --
%    entirely independent of evaluation/Hybrid data, so no circularity.
%    A NEGATIVE correlation (weaker calibration signal -> BIGGER advantage)
%    would support a compensatory effect; near-zero/positive would not.
%
%    Two advantage definitions shown side by side (columns), since it is
%    genuinely ambiguous which one answers "how much did Hybrid help":
%      ALL TRIALS (top row)   -- the practical, real-world bottom line:
%        counterfactual Hybrid-vs-unimodal decided-accuracy advantage
%        (d_mi_subj/d_cvs_subj), computed over the WHOLE session regardless
%        of when each trial resolved. This is what a subject actually
%        experiences using the system.
%      cvsa_influence WINDOW ONLY (bottom row) -- the mechanistic question:
%        integrator buffer advantage (buf_adv_mi_win_subj/buf_adv_cvs_win_subj)
%        restricted to the ONLY window where the fused signal can actually
%        differ from pure MI/CVSA (bayesian_fuse.m forces alpha=0, fused==MI
%        exactly, from cvsa_influence onward) -- isolates the fusion
%        algorithm's own contribution from the later, by-construction-
%        identical part of every trial.
%    If both rows agree (same sign/significance), that is more convincing
%    than either alone; if they disagree, the difference is itself
%    informative (e.g. an all-trials effect not present in the window-only
%    one could reflect something other than the fusion mechanism itself).
flat_erd_calib = flat_erd(strcmp({flat_erd.session}, 'calibration'));
calib_discrim_mi   = nan(n_subj, 1);
calib_discrim_cvsa = nan(n_subj, 1);
for si = 1:n_subj
    mask_mi = strcmp({flat_erd_calib.subject}, subjects{si}) & strcmpi({flat_erd_calib.band_origin}, 'MI');
    if any(mask_mi)
        rows = flat_erd_calib(mask_mi);
        calib_discrim_mi(si) = mean([rows.discrimination], 'omitnan');
    end
    mask_cv = strcmp({flat_erd_calib.subject}, subjects{si}) & strcmpi({flat_erd_calib.band_origin}, 'CVSA');
    if any(mask_cv)
        rows = flat_erd_calib(mask_cv);
        calib_discrim_cvsa(si) = mean([rows.discrimination], 'omitnan');
    end
end
n_calib_subj = sum(~isnan(calib_discrim_mi) | ~isnan(calib_discrim_cvsa));
r_calib_mi_all = NaN; p_calib_mi_all = NaN; r_calib_cv_all = NaN; p_calib_cv_all = NaN;
r_calib_mi_win = NaN; p_calib_mi_win = NaN; r_calib_cv_win = NaN; p_calib_cv_win = NaN;

if n_calib_subj > 0
    [r_calib_mi_all,  p_calib_mi_all]  = pearson_perm_test_local(calib_discrim_mi,   d_mi_subj,           n_perm);
    [r_calib_cv_all,  p_calib_cv_all]  = pearson_perm_test_local(calib_discrim_cvsa, d_cvs_subj,          n_perm);
    [r_calib_mi_win,  p_calib_mi_win]  = pearson_perm_test_local(calib_discrim_mi,   buf_adv_mi_win_subj, n_perm);
    [r_calib_cv_win,  p_calib_cv_win]  = pearson_perm_test_local(calib_discrim_cvsa, buf_adv_cvs_win_subj,n_perm);

    fig18 = figure('Name', 'Group Analysis — Calibration Baseline vs Hybrid Advantage', 'Color', 'w', ...
                   'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig18, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    ax18a = subplot(2,2,1);
    scatter_with_fit(ax18a, calib_discrim_mi, 100*d_mi_subj, subjects, r_calib_mi_all, p_calib_mi_all, ...
        'Calibration MI ERD discrimination |ERD(c1)-ERD(c2)|', '\Delta acc, counterfactual (Hybrid-MI, %, ALL trials)', ...
        'ALL TRIALS: does weak calibration MI signal predict MI advantage?', COL_mi_only);

    ax18b = subplot(2,2,2);
    scatter_with_fit(ax18b, calib_discrim_cvsa, 100*d_cvs_subj, subjects, r_calib_cv_all, p_calib_cv_all, ...
        'Calibration CVSA lateralization |ERD(c1)-ERD(c2)|', '\Delta acc, counterfactual (Hybrid-CVSA, %, ALL trials)', ...
        'ALL TRIALS: does weak calibration CVSA signal predict CVSA advantage?', COL_cvsa_only);

    ax18c = subplot(2,2,3);
    scatter_with_fit(ax18c, calib_discrim_mi, buf_adv_mi_win_subj, subjects, r_calib_mi_win, p_calib_mi_win, ...
        'Calibration MI ERD discrimination |ERD(c1)-ERD(c2)|', '\Delta buffer, cvsa\_influence WINDOW ONLY (Hybrid-MI)', ...
        'WINDOW ONLY: same question, isolated to where fusion can differ from MI', COL_mi_only);

    ax18d = subplot(2,2,4);
    scatter_with_fit(ax18d, calib_discrim_cvsa, buf_adv_cvs_win_subj, subjects, r_calib_cv_win, p_calib_cv_win, ...
        'Calibration CVSA lateralization |ERD(c1)-ERD(c2)|', '\Delta buffer, cvsa\_influence WINDOW ONLY (Hybrid-CVSA)', ...
        'WINDOW ONLY: same question, isolated to where fusion can differ from CVSA', COL_cvsa_only);

    sgtitle(fig18, sprintf(['Group Analysis — Does Calibration-Time Signal Quality Predict the Hybrid Advantage? (n=%d subjects with calibration data)\n' ...
            'Calibration ERD/ERS discrimination is measured BEFORE any evaluation/Hybrid session -- no circularity with the outcome being predicted\n' ...
            'NEGATIVE r = subjects with WEAKER calibration signal get a BIGGER advantage from Hybrid (compensatory); near-zero/positive = no such pattern'], ...
            n_calib_subj), 'Interpreter', 'none');

    fprintf('\n══════ Does calibration-time signal quality predict the Hybrid advantage? (n=%d subjects with calibration data) ══════\n', n_calib_subj);
    fprintf('  ALL TRIALS   : calib MI vs Hyb-MI advantage   r=%+.2f  p=%.4f  %s\n', r_calib_mi_all, p_calib_mi_all, stars(p_calib_mi_all));
    fprintf('  ALL TRIALS   : calib CVSA vs Hyb-CVSA advantage r=%+.2f  p=%.4f  %s\n', r_calib_cv_all, p_calib_cv_all, stars(p_calib_cv_all));
    fprintf('  WINDOW ONLY  : calib MI vs Hyb-MI buffer adv.   r=%+.2f  p=%.4f  %s\n', r_calib_mi_win, p_calib_mi_win, stars(p_calib_mi_win));
    fprintf('  WINDOW ONLY  : calib CVSA vs Hyb-CVSA buffer adv. r=%+.2f  p=%.4f  %s\n', r_calib_cv_win, p_calib_cv_win, stars(p_calib_cv_win));
    fprintf('  (negative r = weaker calibration signal -> bigger Hybrid advantage = compensatory pattern)\n');
    fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

    saveas(fig18, fullfile(out_dir, sprintf('18_%s_calibration_baseline_vs_advantage.svg', group_tag)), 'svg');
    if ~SHOW_FIGURES, close(fig18); end
else
    fprintf('\n[main_group_analysis] Fig 18 skipped: no subject has calibration-folder topo_erders data (run topo_erders.m on calibration/ first).\n');
end

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
group_summary.rescue_frac_inf  = rescue_frac_inf;  group_summary.grand_rescue_frac_inf  = grand_rescue_frac_inf;   % [n_subj x 4] / [1x4], cols: rescued/hurt/both-ok/both-bad
group_summary.rescue_frac_all  = rescue_frac_all;  group_summary.grand_rescue_frac_all  = grand_rescue_frac_all;
group_summary.rescue_frac_post = rescue_frac_post; group_summary.grand_rescue_frac_post = grand_rescue_frac_post;
group_summary.n_trials_post_subj = n_trials_post_subj; group_summary.n_trials_probs_subj = n_trials_probs_subj;
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
group_summary.p_bufadvwin_naive_group    = p_bufadvwin_naive_grp;       % plain sign-flip test on buf_adv_mi_win_subj (floor-limited, companion to the meta p above)
group_summary.d_bufadvwin_naive_cohen    = d_bufadvwin_naive_cohen;
group_summary.buf_adv_cvs_win_subj       = buf_adv_cvs_win_subj;        % per-subject mean Hybrid-CVSA buffer advantage, cvsa_influence window
% Fig 18: calibration-time (independent, pre-evaluation) ERD discrimination vs Hybrid advantage
group_summary.calib_discrim_mi   = calib_discrim_mi;    % [n_subj x 1], NaN = no calibration topo_erders data for that subject
group_summary.calib_discrim_cvsa = calib_discrim_cvsa;
group_summary.r_calib_mi_all  = r_calib_mi_all;  group_summary.p_calib_mi_all  = p_calib_mi_all;   % ALL-trials counterfactual advantage
group_summary.r_calib_cv_all  = r_calib_cv_all;  group_summary.p_calib_cv_all  = p_calib_cv_all;
group_summary.r_calib_mi_win  = r_calib_mi_win;  group_summary.p_calib_mi_win  = p_calib_mi_win;   % cvsa_influence-window-only buffer advantage
group_summary.r_calib_cv_win  = r_calib_cv_win;  group_summary.p_calib_cv_win  = p_calib_cv_win;
group_summary.fraccorr_raw_subj = fraccorr_raw_subj; group_summary.fraccorr_int_subj = fraccorr_int_subj;   % [n_subj x 3] Hybrid/MI-only/CVSA-only
group_summary.p_fcr_mi = p_fcr_mi; group_summary.p_fcr_cvs = p_fcr_cvs; group_summary.d_fcr_mi_cohen = d_fcr_mi_cohen; group_summary.d_fcr_cvs_cohen = d_fcr_cvs_cohen;
group_summary.p_fci_mi = p_fci_mi; group_summary.p_fci_cvs = p_fci_cvs; group_summary.d_fci_mi_cohen = d_fci_mi_cohen; group_summary.d_fci_cvs_cohen = d_fci_cvs_cohen;
group_summary.real_tth = real_tth; group_summary.real_T_all = real_T_all;   % [n_subj x 3] MI/CVSA/Hybrid
group_summary.grand_tth = grand_tth; group_summary.grand_T_all = grand_T_all;
group_summary.cf_tth_subj = cf_tth_subj; group_summary.cf_T_all = cf_T_all; % [n_subj x 3] Hybrid/MI-only/CVSA-only
group_summary.grand_cf_tth = grand_cf_tth; group_summary.grand_cf_T_all = grand_cf_T_all;
group_summary.real_itr_bpt = real_itr_bpt; group_summary.real_itr_bpm = real_itr_bpm;
group_summary.grand_real_itr_bpt = grand_real_itr_bpt; group_summary.grand_real_itr_bpm = grand_real_itr_bpm;
group_summary.cf_itr_bpt = cf_itr_bpt; group_summary.cf_itr_bpm = cf_itr_bpm;
group_summary.grand_cf_itr_bpt = grand_cf_itr_bpt; group_summary.grand_cf_itr_bpm = grand_cf_itr_bpm;
group_summary.p_tth_hyb_mi = p_tth_hyb_mi; group_summary.p_tth_hyb_cvs = p_tth_hyb_cvs;
group_summary.d_tth_hyb_mi_cohen = d_tth_hyb_mi_cohen; group_summary.d_tth_hyb_cvs_cohen = d_tth_hyb_cvs_cohen;
group_summary.p_cftth_hyb_mi = p_cftth_hyb_mi; group_summary.p_cftth_hyb_cvs = p_cftth_hyb_cvs;
group_summary.d_cftth_hyb_mi_cohen = d_cftth_hyb_mi_cohen; group_summary.d_cftth_hyb_cvs_cohen = d_cftth_hyb_cvs_cohen;
group_summary.p_itr_hyb_mi = p_itr_hyb_mi; group_summary.p_itr_hyb_cvs = p_itr_hyb_cvs;
group_summary.d_itr_hyb_mi_cohen = d_itr_hyb_mi_cohen; group_summary.d_itr_hyb_cvs_cohen = d_itr_hyb_cvs_cohen;
group_summary.p_cfitr_hyb_mi = p_cfitr_hyb_mi; group_summary.p_cfitr_hyb_cvs = p_cfitr_hyb_cvs;
group_summary.d_cfitr_hyb_mi_cohen = d_cfitr_hyb_mi_cohen; group_summary.d_cfitr_hyb_cvs_cohen = d_cfitr_hyb_cvs_cohen;
group_summary.keepup_mi_early = keepup_mi_early; group_summary.keepup_mi_late = keepup_mi_late;   % [n_subj x 1]
group_summary.keepup_cv_early = keepup_cv_early; group_summary.keepup_cv_late = keepup_cv_late;
group_summary.n_early_subj = n_early_subj; group_summary.n_late_subj = n_late_subj;
group_summary.n_hit_subj_oc = n_hit_subj;       % Hybrid's own HIT trials = n_early_subj + n_late_subj (the denominator below)
group_summary.n_total_subj_oc = n_total_subj;   % ALL Hybrid trials incl. MISS/TIMEOUT -- context only
group_summary.pct_early_of_trials = pct_early_of_trials; group_summary.pct_late_of_trials = pct_late_of_trials;   % % of n_hit_subj_oc (sums to 100%)
group_summary.grand_pct_early = grand_pct_early; group_summary.grand_pct_late = grand_pct_late;
group_summary.grand_keepup_mi_early = grand_keepup_mi_early; group_summary.grand_keepup_mi_late = grand_keepup_mi_late;
group_summary.grand_keepup_cv_early = grand_keepup_cv_early; group_summary.grand_keepup_cv_late = grand_keepup_cv_late;
group_summary.p_keepup_mi = p_keepup_mi; group_summary.p_keepup_cv = p_keepup_cv;
group_summary.d_keepup_mi_cohen = d_keepup_mi_cohen; group_summary.d_keepup_cv_cohen = d_keepup_cv_cohen;
% "Cost of fusion" (Fig 07b): on Hybrid's own MISS/TIMEOUT trials, would MI-only/CVSA-only alone have HIT instead?
group_summary.rescue_mi_all = rescue_mi_all; group_summary.rescue_cv_all = rescue_cv_all;     % [n_subj x 1]
group_summary.rescue_mi_miss = rescue_mi_miss; group_summary.rescue_cv_miss = rescue_cv_miss;
group_summary.rescue_mi_to = rescue_mi_to; group_summary.rescue_cv_to = rescue_cv_to;
group_summary.n_fail_subj = n_fail_subj; group_summary.n_miss_subj = n_miss_subj; group_summary.n_to_subj = n_to_subj;
group_summary.grand_rescue_mi_all = grand_rescue_mi_all; group_summary.grand_rescue_cv_all = grand_rescue_cv_all;
group_summary.grand_rescue_mi_miss = grand_rescue_mi_miss; group_summary.grand_rescue_cv_miss = grand_rescue_cv_miss;
group_summary.grand_rescue_mi_to = grand_rescue_mi_to; group_summary.grand_rescue_cv_to = grand_rescue_cv_to;
group_summary.p_rescue_all = p_rescue_all; group_summary.p_rescue_miss = p_rescue_miss; group_summary.p_rescue_to = p_rescue_to;
group_summary.d_rescue_all_cohen = d_rescue_all_cohen; group_summary.d_rescue_miss_cohen = d_rescue_miss_cohen; group_summary.d_rescue_to_cohen = d_rescue_to_cohen;
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
group_summary.p_auc_hyb_mi = p_auc_hyb_mi; group_summary.p_auc_hyb_cvs = p_auc_hyb_cvs; group_summary.p_auc_mi_cvs = p_auc_mi_cvs;   % pairwise AUC comparisons
group_summary.d_auc_hyb_mi_cohen = d_auc_hyb_mi_cohen; group_summary.d_auc_hyb_cvs_cohen = d_auc_hyb_cvs_cohen; group_summary.d_auc_mi_cvs_cohen = d_auc_mi_cvs_cohen;
group_summary.brier_subj = brier_subj; group_summary.dprime_subj = dprime_subj;   % [n_subj x 3] MI/CVSA/Hybrid
group_summary.p_brier_hyb_mi = p_brier_hyb_mi; group_summary.p_brier_hyb_cvs = p_brier_hyb_cvs; group_summary.p_brier_mi_cvs = p_brier_mi_cvs;
group_summary.d_brier_hyb_mi_cohen = d_brier_hyb_mi_cohen; group_summary.d_brier_hyb_cvs_cohen = d_brier_hyb_cvs_cohen; group_summary.d_brier_mi_cvs_cohen = d_brier_mi_cvs_cohen;
group_summary.p_dprime_hyb_mi = p_dprime_hyb_mi; group_summary.p_dprime_hyb_cvs = p_dprime_hyb_cvs; group_summary.p_dprime_mi_cvs = p_dprime_mi_cvs;
group_summary.d_dprime_hyb_mi_cohen = d_dprime_hyb_mi_cohen; group_summary.d_dprime_hyb_cvs_cohen = d_dprime_hyb_cvs_cohen; group_summary.d_dprime_mi_cvs_cohen = d_dprime_mi_cvs_cohen;
if any(has_roc_win)
    group_summary.roc_auc_win_subj = roc_auc_win_subj; group_summary.brier_win_subj = brier_win_subj; group_summary.dprime_win_subj = dprime_win_subj;   % [n_subj x 3] MI/CVSA/Hybrid
    group_summary.window_s_subj = window_s_subj;
    group_summary.p_aucw_hyb_mi = p_aucw_hyb_mi; group_summary.p_aucw_hyb_cvs = p_aucw_hyb_cvs; group_summary.p_aucw_mi_cvs = p_aucw_mi_cvs;
    group_summary.d_aucw_hyb_mi_cohen = d_aucw_hyb_mi_cohen; group_summary.d_aucw_hyb_cvs_cohen = d_aucw_hyb_cvs_cohen; group_summary.d_aucw_mi_cvs_cohen = d_aucw_mi_cvs_cohen;
    group_summary.p_brierw_hyb_mi = p_brierw_hyb_mi; group_summary.p_brierw_hyb_cvs = p_brierw_hyb_cvs; group_summary.p_brierw_mi_cvs = p_brierw_mi_cvs;
    group_summary.d_brierw_hyb_mi_cohen = d_brierw_hyb_mi_cohen; group_summary.d_brierw_hyb_cvs_cohen = d_brierw_hyb_cvs_cohen; group_summary.d_brierw_mi_cvs_cohen = d_brierw_mi_cvs_cohen;
    group_summary.p_dprimew_hyb_mi = p_dprimew_hyb_mi; group_summary.p_dprimew_hyb_cvs = p_dprimew_hyb_cvs; group_summary.p_dprimew_mi_cvs = p_dprimew_mi_cvs;
    group_summary.d_dprimew_hyb_mi_cohen = d_dprimew_hyb_mi_cohen; group_summary.d_dprimew_hyb_cvs_cohen = d_dprimew_hyb_cvs_cohen; group_summary.d_dprimew_mi_cvs_cohen = d_dprimew_mi_cvs_cohen;
end
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
group_summary.sigqc_rate           = sigqc_rate;            % [n_subj x 3] MI/CVSA/Hybrid: fraction of trials with a dead/noisy channel during CF
group_summary.sigqc_n_flagged      = sigqc_n_flagged;       % [n_subj x 3] raw counts
group_summary.sigqc_n_total        = sigqc_n_total;         % [n_subj x 3] raw counts (denominator)
group_summary.grand_sigqc_rate     = grand_sigqc_rate;      % [1x3] trial-pooled across all subjects
group_summary.sigqc_rate_overall   = sigqc_rate_overall;    % [n_subj x 1] pooled across paradigms -- used by well_vs_bad_comparison
group_summary.r_sigqc_vs_real  = r_sigqc_vs_real;     % cross-subject Pearson r: signal-quality flagged-trial rate vs real Hybrid accuracy
group_summary.p_sigqc_vs_real  = p_sigqc_vs_real;
% Fig 17 executive verdict: -1/0/+1 (worse/n.s./better) per metric (rows,
% see verdict_rows labels) x [vs MI-only, vs CVSA-only] -- pure synthesis
% of tests already saved individually above, kept here too for convenience.
group_summary.verdict_row_labels = verdict_rows(:,1);
group_summary.verdict_code       = verdict_code;   % [n_metrics x 2]
save(fullfile(out_dir, sprintf('group_summary_%s.mat', group_tag)), 'group_summary');
fprintf('\nSaved group_summary_%s.mat and figures to %s\n', group_tag, out_dir);

end % main_group_analysis

% ── Local helpers ────────────────────────────────────────────────────────────

function plot_subject_grand_bars(ax, val_mat, grand_vec, subjects, labels, cols, bw, is_pct, show_labels, sig_pairs)
% PLOT_SUBJECT_GRAND_BARS  One hand-rolled box (mean +/- SD, whiskers to the
%   TRUE min/max -- no Statistics Toolbox, no 1.5*IQR truncation/outlier
%   convention) per CONDITION (one per column of val_mat), pooling ALL
%   subjects; every subject's own value is overlaid as a jittered dot so you
%   can see exactly where each subject falls, not just the summary. The
%   trial-pooled GRAND value (grand_vec -- NOT the mean of the per-subject
%   dots, can differ when subjects contribute unequal trial counts) is
%   marked with a black diamond on top of its condition's box.
%   subjects/bw/show_labels are kept as parameters for call-site
%   compatibility with every existing caller but are not used by this
%   box-plot rendering (subjects are identified by their dot, not a per-
%   subject x position; show_labels' per-bar numeric text isn't needed once
%   the box+dots already show the distribution and console tables already
%   print exact per-subject numbers).
%   sig_pairs (optional): cell array of {i, j, p} triples -- draws a
%   significance bracket + stars between condition columns i and j, using a
%   p-value ALREADY computed by the caller (e.g. the same paired sign-flip
%   test reported in that figure's console table/title elsewhere) -- no new
%   statistics computed here, purely a visual addition. Stacked above the
%   highest box/whisker in the panel, YLim extended to fit them.
    if nargin < 8, is_pct = true; end
    if nargin < 9, show_labels = false; end %#ok<INUSD>
    if nargin < 10, sig_pairs = {}; end
    n_cond = size(val_mat, 2);
    scale = 1; if is_pct, scale = 100; end
    box_hw = 0.28;
    hold(ax, 'on');
    overall_max = -Inf;
    for k = 1:n_cond
        v = scale * val_mat(:,k);
        v_ok = v(~isnan(v));
        if isempty(v_ok), continue; end
        mu = mean(v_ok); sd = std(v_ok);
        lo = mu - sd; hi = mu + sd;
        mn = min(v_ok); mx = max(v_ok);
        overall_max = max(overall_max, mx);

        plot(ax, [k k], [mn lo], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(ax, [k k], [hi mx], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(ax, [k-0.12, k+0.12], [mn mn], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(ax, [k-0.12, k+0.12], [mx mx], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        patch(ax, [k-box_hw, k+box_hw, k+box_hw, k-box_hw], [lo lo hi hi], cols{k}, ...
              'FaceAlpha', 0.35, 'EdgeColor', 'k', 'LineWidth', 1.2, 'DisplayName', labels{k});
        plot(ax, [k-box_hw, k+box_hw], [mu mu], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');

        jit = (rand(numel(v),1) - 0.5) * 0.36;
        scatter(ax, k + jit, v, 24, cols{k}*0.55, 'filled', 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.5, 'MarkerFaceAlpha', 0.8, 'HandleVisibility', 'off');

        gv = scale * grand_vec(k);
        if ~isnan(gv)
            plot(ax, k, gv, 'd', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', ...
                 'MarkerSize', 9, 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end
    set(ax, 'XTick', 1:n_cond, 'XTickLabel', labels, 'XLim', [0.4, n_cond+0.6]);
    if is_pct, set(ax, 'YLim', [0, 110]); end

    if ~isempty(sig_pairs)
        if isinf(overall_max), overall_max = 0; end
        step = 0.09 * max(overall_max, 1);
        y_top = overall_max + 2.5*step;
        for spi = 1:numel(sig_pairs)
            pr = sig_pairs{spi};
            draw_sig_bracket_local(ax, pr{1}, pr{2}, y_top - (spi-1)*step, stars(pr{3}));
        end
        if is_pct
            set(ax, 'YLim', [0, y_top + step]);
        else
            yl = ylim(ax);
            ylim(ax, [yl(1), y_top + step]);
        end
    end
end

function plot_delta_boxplot(ax, delta_mat, grand_delta_vec, labels, cols, unit_scale, p_vals)
% PLOT_DELTA_BOXPLOT  Same hand-rolled box convention as
%   plot_subject_grand_bars (mean +/- SD box, whiskers to the true min/max,
%   jittered per-subject dots, black diamond = trial-pooled GRAND), but for
%   DELTA/advantage columns (e.g. Hybrid-MI, Hybrid-CVSA) rather than raw
%   per-condition values -- directly shows whether an advantage is UNIFORM
%   across the cohort (tight box, dots close together) or HETEROGENEOUS
%   (wide box, dots spanning positive to negative -- helps some subjects a
%   lot, others not at all or even negatively). delta_mat: [n_subj x
%   n_deltas]. unit_scale: 100 for fraction->%, 1 for already-scaled units
%   (seconds, bits/min, AUC points). grand_delta_vec may be [] to skip the
%   GRAND marker (e.g. when no natural trial-pooled equivalent exists).
%   p_vals (optional): [1 x n_deltas], the SAME paired sign-flip p already
%   reported in that figure's title/console -- draws a star above each box
%   marking whether THAT delta is significantly different from zero (one-
%   sample style, not a between-box bracket -- no new statistic computed).
    if nargin < 6, unit_scale = 100; end
    if nargin < 7, p_vals = []; end
    n_delta = size(delta_mat, 2);
    box_hw = 0.28;
    hold(ax, 'on');
    yline(ax, 0, 'k-', 'HandleVisibility', 'off');
    for k = 1:n_delta
        v = unit_scale * delta_mat(:,k);
        v_ok = v(~isnan(v));
        if isempty(v_ok), continue; end
        mu = mean(v_ok); sd = std(v_ok);
        lo = mu - sd; hi = mu + sd;
        mn = min(v_ok); mx = max(v_ok);

        plot(ax, [k k], [mn lo], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(ax, [k k], [hi mx], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(ax, [k-0.12, k+0.12], [mn mn], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(ax, [k-0.12, k+0.12], [mx mx], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        patch(ax, [k-box_hw, k+box_hw, k+box_hw, k-box_hw], [lo lo hi hi], cols{k}, ...
              'FaceAlpha', 0.35, 'EdgeColor', 'k', 'LineWidth', 1.2, 'DisplayName', labels{k});
        plot(ax, [k-box_hw, k+box_hw], [mu mu], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');

        jit = (rand(numel(v),1) - 0.5) * 0.36;
        scatter(ax, k + jit, v, 24, cols{k}*0.55, 'filled', 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.5, 'MarkerFaceAlpha', 0.8, 'HandleVisibility', 'off');

        if ~isempty(grand_delta_vec) && numel(grand_delta_vec) >= k && ~isnan(grand_delta_vec(k))
            gv = unit_scale * grand_delta_vec(k);
            plot(ax, k, gv, 'd', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', ...
                 'MarkerSize', 9, 'LineWidth', 1, 'HandleVisibility', 'off');
        end

        if ~isempty(p_vals) && numel(p_vals) >= k && ~isnan(p_vals(k))
            span = max(mx - mn, eps);
            text(ax, k, mx + 0.08*span, stars(p_vals(k)), 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold', 'FontSize', 10);
        end
    end
    set(ax, 'XTick', 1:n_delta, 'XTickLabel', labels, 'XLim', [0.4, n_delta+0.6]);
end

function subj = subject_of(fpath, root_dir)
    rel = erase(fpath, [char(root_dir) filesep]);
    parts = strsplit(rel, filesep);
    subj = parts{1};
end

function sess = session_of(fpath)
% SESSION_OF  'evaluation' | 'calibration' | 'unknown', inferred from the
%   folder path itself (the .mat carries no session-type field of its own)
%   -- run_subject_analysis.m calls topo_erders() once on the evaluation/
%   folder and once on the sibling calibration/ folder, so the folder name
%   is the only place this distinction lives. Same helper as
%   group_topo_erders.m, duplicated here rather than shared (this file's
%   own convention for small helpers).
    parts = strsplit(fpath, filesep);
    if any(strcmpi(parts, 'evaluation'))
        sess = 'evaluation';
    elseif any(strcmpi(parts, 'calibration'))
        sess = 'calibration';
    else
        sess = 'unknown';
    end
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

function draw_sig_bracket_local(ax, x1, x2, y, label)
% DRAW_SIG_BRACKET_LOCAL  Horizontal significance bracket between bar
%   positions x1 and x2 at height y, with a stars/n.s. label centred above.
    line(ax, [x1 x1 x2 x2], [y-0.01 y y y-0.01], 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(ax, (x1+x2)/2, y+0.005, label, 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
end

function plot_group_marker(ax, n_subj, deltas, col, pval, x_offset)
% PLOT_GROUP_MARKER  Group mean +- SEM errorbar with significance stars,
%   placed just right of a per-subject delta scatter at x = n_subj+x_offset.
%   Shared by Fig 2c's four panels (TTH/ITR x real/counterfactual) to avoid
%   repeating the same errorbar+text block eight times. Text offset is
%   relative to that series' own SEM (not a fixed pixel/data offset), so it
%   reads correctly whether the panel's y-axis is seconds or bits/min.
    m  = mean(deltas, 'omitnan');
    se = std(deltas, 'omitnan') / sqrt(max(1, sum(~isnan(deltas))));
    errorbar(ax, n_subj+x_offset, m, se, 'o', 'Color', col, 'MarkerFaceColor', col, 'LineWidth', 1.5, 'CapSize', 6);
    if isnan(se) || se == 0, se_txt = abs(m)*0.1 + eps; else, se_txt = se; end
    text(ax, n_subj+x_offset, m + se_txt*sign(m+eps)*1.4, stars(pval), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

function [brier, dprime] = brier_dprime_local(scores, labels)
% BRIER_DPRIME_LOCAL  Brier score (mean((P(target)-label)^2), LOWER=better
%   calibrated) and d-prime (z(hit rate)-z(false alarm rate) at the P=0.5
%   operating point, HIGHER=better separated) on a pooled (score,label)
%   frame set. Shared by the ALL-frames and windowed (cvsa_influence-
%   restricted) computations so both use the identical formula.
    lbl = logical(labels);
    brier = mean((scores - double(lbl)).^2, 'omitnan');
    hr  = mean(scores(lbl)  > 0.5, 'omitnan');
    far = mean(scores(~lbl) > 0.5, 'omitnan');
    dprime = NaN;
    if ~isnan(hr) && ~isnan(far)
        dprime = norminv_local(hr) - norminv_local(far);   % clamped away from 0/1 inside norminv_local
    end
end

function plot_acc_vs_timeout_scatter(ax, acc_dec, to_rate, cols, labels)
% PLOT_ACC_VS_TIMEOUT_SCATTER  One point per subject per group (column of
%   acc_dec/to_rate): x = timeout rate (%), y = decided-trial accuracy (%).
%   Decided-trial accuracy alone can look great on very few decisions (high
%   timeout rate) -- this panel puts both numbers on one plot instead of two
%   separate bar charts, so "good decisions but rarely reached" (top-right)
%   reads differently from "good and frequent decisions" (top-left). A big
%   star per group marks that group's own subject-pooled mean position.
    hold(ax, 'on');
    n_grp = size(acc_dec, 2);
    for g = 1:n_grp
        x = 100*to_rate(:,g); y = 100*acc_dec(:,g);
        scatter(ax, x, y, 55, cols{g}, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.75, 'DisplayName', labels{g});
        mx = mean(x, 'omitnan'); my = mean(y, 'omitnan');
        if ~isnan(mx) && ~isnan(my)
            plot(ax, mx, my, 'p', 'MarkerSize', 16, 'MarkerFaceColor', cols{g}, 'MarkerEdgeColor', 'k', ...
                 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    yline(ax, 50, 'k:', 'HandleVisibility', 'off');
    xlabel(ax, 'timeout rate (%)');
    ylabel(ax, 'accuracy on decided trials (%)');
    legend(ax, 'Location', 'best', 'FontSize', 8);
    title(ax, {'Decision quality vs. how often a decision is reached at all', ...
          '(large star = group mean; top-left = frequent AND good decisions)'}, 'FontWeight', 'bold', 'FontSize', 8);
    grid(ax, 'on');
end

function pct = pooled_frac_4(c1, c2, c3, c4)
% POOLED_FRAC_4  Pool 4 category counts (cell arrays of per-file scalars,
%   possibly containing empty entries for files without that field -- MATLAB
%   silently skips [] when concatenated with [c1{:}]) into fractions of the
%   grand total. Returns NaN(1,4) if every category is empty/zero.
    counts = [sum([c1{:}], 'omitnan'), sum([c2{:}], 'omitnan'), sum([c3{:}], 'omitnan'), sum([c4{:}], 'omitnan')];
    tot = sum(counts);
    if tot == 0, pct = nan(1,4); else, pct = counts / tot; end
end

function m = weighted_mean_cols(vals, weights)
% WEIGHTED_MEAN_COLS  Column-wise weighted mean; vals and weights are the
%   SAME [n x k] shape (each column weighted by its own matching count --
%   e.g. TTH weighted by n_hit, T-miss by n_miss). NaN entries (or entries
%   with zero weight) contribute nothing to either the sum or the
%   denominator, so an empty category doesn't drag the mean toward zero.
    w = double(weights);
    w(isnan(vals) | isnan(w)) = 0;
    v = vals; v(isnan(v)) = 0;
    denom = sum(w, 1);
    m = sum(v .* w, 1) ./ max(denom, eps);
    m(denom == 0) = NaN;
end

function T = pool3_weighted(v1, n1, v2, n2, v3, n3)
% POOL3_WEIGHTED  Combine three category means (e.g. TTH/T-miss/T-timeout)
%   into one overall mean-time-to-ANY-outcome, weighted by each category's
%   own trial count. NaN/zero-weight categories contribute nothing. Works
%   element-wise on any matching-size inputs (scalar row or full matrix).
    v1(isnan(v1) | n1==0) = 0;
    v2(isnan(v2) | n2==0) = 0;
    v3(isnan(v3) | n3==0) = 0;
    denom = n1 + n2 + n3;
    T = (v1.*n1 + v2.*n2 + v3.*n3) ./ max(denom, eps);
    T(denom==0) = NaN;
end

function B = decided_rate_itr_local(acc_dec, to_rate, N_CLS)
% DECIDED_RATE_ITR_LOCAL  bits/trial = (decided-trial rate) x (bits/trial
%   from DECIDED-trials-only accuracy) -- see the long comment above this
%   function's call site for why this, and not itr_bits_per_trial_local(acc)
%   directly on the TIMEOUT=fail accuracy, is the correct formula whenever
%   a stream has a non-negligible timeout rate. A 100%-timeout stream
%   (to_rate=1, acc_dec undefined/NaN) correctly returns 0, not NaN.
    if isnan(to_rate), B = NaN; return; end
    decided_rate = 1 - to_rate;
    if decided_rate <= 0, B = 0; return; end
    if isnan(acc_dec), B = NaN; return; end
    B = decided_rate * itr_bits_per_trial_local(acc_dec, N_CLS);
end

function B = itr_bits_per_trial_local(P, N)
% ITR_BITS_PER_TRIAL_LOCAL  Wolpaw et al. bits/trial for an N-class forced
%   choice at accuracy P: B = log2(N) + P*log2(P) + (1-P)*log2((1-P)/(N-1)).
%   P=0 and P=1 edge cases handled explicitly (0*log2(0) := 0). Same formula
%   as main_hybrid_advantage_integ.m's itr_bits_per_trial (duplicated locally
%   rather than shared, matching this codebase's convention for small
%   hand-rolled stats helpers, e.g. sign_flip_test_local in that file).
    if isnan(P) || isnan(N) || N < 2, B = NaN; return; end
    B = log2(N);
    if P > 0, B = B + P * log2(P); end
    if P < 1, B = B + (1-P) * log2((1-P) / (N-1)); end
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

function [verdict_code, txt] = verdict_cell_local(delta_vec, p, dcohen, higher_better, delta_fmt)
% VERDICT_CELL_LOCAL  Turns an already-computed (delta vector, p, Cohen's d)
%   triple into a -1/0/+1 verdict code (Hybrid worse/n.s./better, correcting
%   for metrics where LOWER is actually better, e.g. TTH) plus a display
%   string, for Fig 17's executive-verdict scorecard.
    md = mean(delta_vec, 'omitnan');
    if isnan(p) || isnan(md)
        verdict_code = NaN;
        txt = 'n/a';
        return;
    end
    eff = md; if ~higher_better, eff = -md; end
    if p < 0.05 && eff > 0
        verdict_code = 1;
    elseif p < 0.05 && eff < 0
        verdict_code = -1;
    else
        verdict_code = 0;
    end
    txt = sprintf([delta_fmt '\n%s (d=%+.2f)'], md, stars(p), dcohen);
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

