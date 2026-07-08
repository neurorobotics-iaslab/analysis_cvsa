%% MAIN_ROC_ANALYSIS  Classifier-level ROC/AUC from the raw sLDA output (pre-integrator).
%
%   Complements main_threshold_sweep.m: that script sweeps the INTEGRATOR's
%   buffer thresholds and is deliberately NOT a classic ROC (two independent
%   per-class thresholds + a TIMEOUT outcome do not reduce to a single binary
%   label). This script instead evaluates the CLASSIFIER itself: the raw
%   per-frame sLDA probability P(class 1), pooled across CF frames, swept
%   over every possible scalar threshold in [0,1] -- a genuine ROC/AUC,
%   entirely independent of the integrator/buffer/threshold machinery.
%
%   Reuses trials(t).raw from integrate_signal.m, which already carries the
%   paradigm-correct per-frame signal:
%     MI paradigm     -> raw MI sLDA P(c)
%     CVSA paradigm   -> raw CVSA sLDA P(c)
%     Hybrid paradigm -> the cosine-annealed Bayesian-fused P(c)
%                        (bayesian_fuse.m), treated here as the output of a
%                        single hypothetical classifier.
%   Label per frame = the trial's target class (frame-level binary label,
%   constant within a trial). Artifact-flagged frames and the N_PRE reset
%   frame are excluded from pooling: artifact frames reflect EEG corrupted
%   upstream of the classifier, not a genuine class-conditional sample, and
%   are exactly the frames the real system ignores (artifact gate freezes
%   the integrator on them).
%
%   Mirrors main_threshold_sweep.m's file-selection/grouping structure:
%   works across a mix of MI/CVSA/Hybrid GDFs in one run, one own ROC image
%   per paradigm present, plus one per individual GDF file. UNLIKE
%   threshold_sweep's pooled "ALL" image, a ROC pooling raw scores from
%   different classifiers (MI-only vs CVSA-only vs Hybrid-fused) would not
%   be a meaningful single AUC -- each paradigm is its own classifier on its
%   own probability scale, so blindly pooling frames would be invalid.
%   Instead "general" here means a comparison figure overlaying each
%   paradigm's OWN ROC curve (each with its own AUC) on one plot -- the
%   useful cross-paradigm view for a paper, without conflating classifiers.
%
%   Called by run_subject_analysis.m as part of the standard per-session
%   pipeline (all GDFs in an evaluation/ folder). roc_summary.mat's
%   roc_by_group entries additionally carry the raw pooled (score, label)
%   pairs behind each paradigm's curve (not just the derived fpr/tpr), so
%   main_group_analysis.m can validly RE-POOL a subject's sessions together
%   (same deployed classifier/calibration across a subject's days) into a
%   per-subject ROC, then macro-average those per-subject curves across the
%   cohort (interpolated onto a common FPR grid) -- pooling raw scores
%   ACROSS subjects/classifiers directly would not be valid, exactly as
%   discussed for the "ALL paradigms" case above.
%
%   Output, under <gdf_dir>/analysis_results/roc_analysis/:
%     roc_<paradigm>.svg          one per paradigm present
%     roc_ALL_paradigms.svg       overlay of each paradigm's own ROC curve
%     per_file/roc_<basename>.svg one per selected GDF
%     roc_summary.mat             fpr/tpr/thr/auc + n per group and per file
%                                  (+ raw pooled scores/labels per group)

function main_roc_analysis(gdf_dir, gdf_names, show_figures)
%   Callable as a function: main_roc_analysis(gdf_dir, gdf_names, show_figures)
%   With no arguments, shows the GUI file picker (interactive mode).

if nargin < 3, show_figures = false; end
SHOW_FIGURES = show_figures;

% --- Make every subfolder visible to the simulator --------------------------
this_dir = '/home/paolo/bci_vr_ws/src/analysis_bci/matlab_simulation';
addpath(this_dir);
addpath(fullfile(this_dir, 'io'));
addpath(fullfile(this_dir, 'processing'));
addpath(fullfile(this_dir, 'artifacts'));
addpath(fullfile(this_dir, 'classifier'));
addpath(fullfile(this_dir, 'integrator'));
addpath(fullfile(this_dir, 'plotting'));
addpath(fullfile(this_dir, 'utils'));

% --- GUI or batch file selection -------------------------------------------
if nargin < 1 || isempty(gdf_dir)
    default_dir = '/home/paolo/bci_vr_ws/recordings';
    [gdf_names, gdf_dir] = uigetfile({'*.gdf', 'GDF recordings (*.gdf)'}, ...
                                    'Select GDF recording(s)', default_dir, 'MultiSelect', 'on');
    if isequal(gdf_names, 0)
        error('main_roc_analysis:cancel', 'No GDF selected.');
    end
elseif nargin < 2 || isempty(gdf_names)
    f = dir(fullfile(gdf_dir, '*.gdf'));
    gdf_names = {f.name};
    if isempty(gdf_names), error('main_roc_analysis:nofiles', 'No GDF files in %s', gdf_dir); end
end
if ischar(gdf_names), gdf_names = {gdf_names}; end
n_files = numel(gdf_names);

all_trials = struct([]);   % pooled across files: .raw, .artifact, .target_class, .n_pre, .n_cf, .paradigm, .basename

for file_idx = 1:n_files
gdf_path = fullfile(gdf_dir, gdf_names{file_idx});
fprintf('\n[%d/%d] %s\n', file_idx, n_files, gdf_names{file_idx});

%% --- Load GDF + YAML ---------------------------------------------------
[signal, header, basename] = load_gdf(gdf_path);
[params, ~] = load_params_yaml(gdf_path);

paradigm   = params.integrator.paradigm;
fs         = double(params.acquisition.samplerate);
framerate  = double(params.acquisition.framerate);
chunk_size = round(fs / framerate);
if abs(fs - header.SampleRate) > 1e-3
    fs = header.SampleRate;
    chunk_size = round(fs / framerate);
end

bufsize_proc = double(params.RingBufferCfg.params.size);
bufsize_art  = double(params.RingBufferCfgArtifact.params.size);
eog_names    = to_strcell(params.CarCfg.params.EOG_ch_names);

do_car_mi = true;
if isfield(params, 'processing_fbcsp_mi')
    do_car_mi = logical(params.processing_fbcsp_mi.do_car);
    nchannels = double(params.processing_fbcsp_mi.nchannels);
end
do_car_cvsa = true;
if isfield(params, 'processing_fbcsp_cvsa')
    do_car_cvsa = logical(params.processing_fbcsp_cvsa.do_car);
    nchannels = double(params.processing_fbcsp_cvsa.nchannels);
end

signal = signal(:, 1:nchannels);   % strip ACC/non-EEG channels

%% --- Load CSP + sLDA per used paradigm ---------------------------------
use_mi   = ismember(paradigm, {'mi', 'hybrid'});
use_cvsa = ismember(paradigm, {'cvsa', 'hybrid'});
csp_mi   = []; slda_mi   = [];
csp_cvsa = []; slda_cvsa = [];
if use_mi
    csp_mi   = load_csp(params, 'mi');
    slda_mi  = load_slda(params, 'mi');
end
if use_cvsa
    csp_cvsa  = load_csp(params, 'cvsa');
    slda_cvsa = load_slda(params, 'cvsa');
end

%% --- Apply processing (one stream per paradigm) ------------------------
features_mi = []; header_mi = [];
features_cv = []; header_cv = [];
info_proc   = [];
if use_mi
    proc_cfg_mi = struct('samplerate', fs, 'chunk_size', chunk_size, ...
                         'bufsize', bufsize_proc, 'filter_order', 4, ...
                         'do_car', do_car_mi, 'eog_names', {eog_names});
    [features_mi, header_mi, info_proc] = apply_processing(signal, header, csp_mi, proc_cfg_mi);
end
if use_cvsa
    proc_cfg_cvsa = struct('samplerate', fs, 'chunk_size', chunk_size, ...
                           'bufsize', bufsize_proc, 'filter_order', 4, ...
                           'do_car', do_car_cvsa, 'eog_names', {eog_names});
    [features_cv, header_cv, info2] = apply_processing(signal, header, csp_cvsa, proc_cfg_cvsa);
    if isempty(info_proc), info_proc = info2; end
end

%% --- Apply artifact detection on raw signal ----------------------------
art_cfg = params.ArtifactCfg.params;
art_cfg.EOG_ch_names = to_strcell(art_cfg.EOG_ch_names);
cfg_art = struct('samplerate', fs, 'chunk_size', chunk_size, 'bufsize_artifact', bufsize_art);
[art_flags, ~] = detect_artifacts(signal, header, art_cfg, cfg_art);

%% --- Apply sLDA over the whole feature stream --------------------------
p_mi_aligned   = []; if use_mi,   p_mi_aligned   = apply_slda(features_mi, slda_mi,   csp_mi.bands  ); end
p_cvsa_aligned = []; if use_cvsa, p_cvsa_aligned = apply_slda(features_cv, slda_cvsa, csp_cvsa.bands); end

%% --- Integrate per trial (only trials(t).raw / .artifact are needed here;
%%     buffer/threshold fields are ignored -- this script does not evaluate
%%     the integrator at all, see main_threshold_sweep.m for that.) --------
int_cfg = params.integrator;
if ~isfield(int_cfg, 'increment'),            int_cfg.increment = 1; end
if ~isfield(int_cfg, 'thresholds_rejection'), int_cfg.thresholds_rejection = []; end
if ~isfield(int_cfg, 'cvsa_influence'),       int_cfg.cvsa_influence = 2.5; end
if ~isfield(int_cfg, 'thresholds'),           int_cfg.thresholds = params.training_node.thresholds; end

if use_mi,   header_chunks = header_mi;
else,        header_chunks = header_cv;
end
header_chunks.framerate = framerate;

trials = integrate_signal(p_mi_aligned, p_cvsa_aligned, art_flags, ...
                          header_chunks, int_cfg, paradigm);

%% --- Pool usable trials (need a target class and at least one CF chunk)
for t = 1:numel(trials)
    tr = trials(t);
    if isnan(tr.target_class) || tr.n_cf == 0, continue; end
    idx_acc = numel(all_trials) + 1;
    all_trials(idx_acc).raw          = tr.raw; %#ok<AGROW>
    all_trials(idx_acc).artifact     = tr.artifact;
    all_trials(idx_acc).target_class = tr.target_class;
    all_trials(idx_acc).n_pre        = tr.n_pre;
    all_trials(idx_acc).n_cf         = tr.n_cf;
    all_trials(idx_acc).paradigm     = paradigm;
    all_trials(idx_acc).basename     = basename;
end

end % file_idx loop

if isempty(all_trials)
    error('main_roc_analysis:no_trials', 'No usable trials found across the selected file(s).');
end

n_trials = numel(all_trials);
fprintf('\nPooled %d trials across %d file(s).\n', n_trials, n_files);

out_dir      = fullfile(gdf_dir, 'analysis_results', 'roc_analysis');
out_dir_file = fullfile(out_dir, 'per_file');
if ~exist(out_dir, 'dir'),      mkdir(out_dir);      end
if ~exist(out_dir_file, 'dir'), mkdir(out_dir_file); end

fig_vis = 'off'; if SHOW_FIGURES, fig_vis = 'on'; end

COL_mi     = [0.85 0.30 0.10];
COL_cvsa   = [0.10 0.60 0.30];
COL_hybrid = [0.18 0.45 0.75];
par_colors = containers.Map({'mi', 'cvsa', 'hybrid'}, {COL_mi, COL_cvsa, COL_hybrid});

%% --- One ROC image per paradigm present ------------------------------------
paradigms_present = unique({all_trials.paradigm}, 'stable');

roc_by_group = struct('name', {}, 'n_trials', {}, 'n_frames', {}, ...
                       'fpr', {}, 'tpr', {}, 'thr', {}, 'auc', {}, ...
                       'scores', {}, 'labels', {});
for g = 1:numel(paradigms_present)
    name = paradigms_present{g};
    trials_g = all_trials(strcmp({all_trials.paradigm}, name));
    [scores, labels] = pool_group_frames(trials_g);
    [fpr, tpr, thr, auc_val] = compute_roc_curve(scores, labels);

    roc_by_group(g).name     = name; %#ok<AGROW>
    roc_by_group(g).n_trials = numel(trials_g);
    roc_by_group(g).n_frames = numel(scores);
    roc_by_group(g).fpr      = fpr;
    roc_by_group(g).tpr      = tpr;
    roc_by_group(g).thr      = thr;
    roc_by_group(g).auc      = auc_val;
    % Raw pooled (score, label) pairs, kept (not just the derived curve) so
    % main_group_analysis can validly RE-POOL a subject's sessions together
    % (same deployed classifier/calibration across a subject's days) before
    % computing that subject's own ROC -- pooling raw scores ACROSS subjects
    % would not be valid (different classifiers/calibrations), which is why
    % only THIS per-subject-usable field is kept, not one for the cross-
    % subject step.
    roc_by_group(g).scores   = single(scores);
    roc_by_group(g).labels   = logical(labels);

    print_roc(name, numel(trials_g), numel(scores), scores, labels, auc_val);

    col = COL_mi;
    if isKey(par_colors, name), col = par_colors(name); end
    out_file = fullfile(out_dir, sprintf('roc_%s.svg', name));
    save_roc_figure(out_file, upper(name), numel(trials_g), numel(scores), ...
        fpr, tpr, auc_val, scores, labels, col, fig_vis, SHOW_FIGURES);
end

%% --- Comparison figure: each paradigm's OWN ROC overlaid (NOT a pooled curve)
fig_cmp = figure('Name', 'ROC comparison across paradigms', 'Color', 'w', ...
                 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig_cmp, 'Units', 'normalized', 'Position', [0 0 1 1]);
ax = axes(fig_cmp); hold(ax, 'on');
plot(ax, [0 1], [0 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
for g = 1:numel(roc_by_group)
    name = roc_by_group(g).name;
    col = COL_mi;
    if isKey(par_colors, name), col = par_colors(name); end
    plot(ax, roc_by_group(g).fpr, roc_by_group(g).tpr, '-', 'Color', col, 'LineWidth', 2, ...
         'DisplayName', sprintf('%s  (AUC=%.3f, n=%d trials)', upper(name), roc_by_group(g).auc, roc_by_group(g).n_trials));
end
xlabel(ax, 'False Positive Rate'); ylabel(ax, 'True Positive Rate');
xlim(ax, [0 1]); ylim(ax, [0 1]); axis(ax, 'square'); grid(ax, 'on');
legend(ax, 'Location', 'southeast', 'FontSize', 10);
title(ax, sprintf(['ROC comparison across paradigms (n=%d files)\n' ...
        'each curve is its OWN classifier (raw sLDA, or fused P for Hybrid) -- NOT a pooled curve'], n_files), ...
        'Interpreter', 'none', 'FontWeight', 'bold');

out_file_cmp = fullfile(out_dir, 'roc_ALL_paradigms.svg');
saveas(fig_cmp, out_file_cmp, 'svg');
fprintf('Saved %s\n', out_file_cmp);
if ~SHOW_FIGURES, close(fig_cmp); end

%% --- One image per individual GDF file -----------------------------------
basenames_present = unique({all_trials.basename}, 'stable');
roc_by_file = struct('basename', {}, 'paradigm', {}, 'n_trials', {}, 'n_frames', {}, ...
                      'fpr', {}, 'tpr', {}, 'thr', {}, 'auc', {});
for f = 1:numel(basenames_present)
    bn = basenames_present{f};
    trials_f = all_trials(strcmp({all_trials.basename}, bn));
    par_f    = trials_f(1).paradigm;
    [scores, labels] = pool_group_frames(trials_f);
    [fpr, tpr, thr, auc_val] = compute_roc_curve(scores, labels);

    roc_by_file(f).basename = bn; %#ok<AGROW>
    roc_by_file(f).paradigm = par_f;
    roc_by_file(f).n_trials = numel(trials_f);
    roc_by_file(f).n_frames = numel(scores);
    roc_by_file(f).fpr      = fpr;
    roc_by_file(f).tpr      = tpr;
    roc_by_file(f).thr      = thr;
    roc_by_file(f).auc      = auc_val;

    print_roc(bn, numel(trials_f), numel(scores), scores, labels, auc_val);

    col = COL_mi;
    if isKey(par_colors, par_f), col = par_colors(par_f); end
    out_file = fullfile(out_dir_file, sprintf('roc_%s.svg', bn));
    save_roc_figure(out_file, sprintf('%s (%s)', bn, upper(par_f)), numel(trials_f), numel(scores), ...
        fpr, tpr, auc_val, scores, labels, col, fig_vis, SHOW_FIGURES);
end

%% --- Save .mat for reuse (e.g. future cross-subject aggregation) ---------
roc_summary = struct('roc_by_group', roc_by_group, 'roc_by_file', roc_by_file, ...
    'n_trials', n_trials, 'n_files', n_files);
save(fullfile(out_dir, 'roc_summary.mat'), 'roc_summary');
fprintf('Saved %s\n', fullfile(out_dir, 'roc_summary.mat'));

end % main_roc_analysis

% ── Local helpers (must be after all script statements) ────────────────────

function [scores, labels] = pool_group_frames(trials_g)
% POOL_GROUP_FRAMES  Flatten CF-window frames (excluding the N_PRE reset
%   frame, artifact-flagged frames, and NaN raw values -- e.g. hybrid frames
%   before both MI and CVSA streams are available) across a set of trials
%   into a single (score, label) list. score = P(class 1) for that frame;
%   label = true if the trial's target class is class 1.
    scores = [];
    labels = [];
    for it = 1:numel(trials_g)
        tr = trials_g(it);
        cf_range = (tr.n_pre + 1):(tr.n_pre + tr.n_cf);
        p1  = tr.raw(cf_range, 1);
        art = tr.artifact(cf_range);
        valid = ~art & ~isnan(p1);
        n_valid = sum(valid);
        if n_valid == 0, continue; end
        scores = [scores; p1(valid)]; %#ok<AGROW>
        labels = [labels; repmat(tr.target_class == 1, n_valid, 1)]; %#ok<AGROW>
    end
end

function print_roc(label, n_trials, n_frames, scores, labels, auc_val)
% PRINT_ROC  Console summary: AUC + sensitivity/specificity at the natural
%   p=0.5 operating point (matches main_session_overview's frame-accuracy
%   convention: fraction of CF frames where argmax(P) == target).
    if n_frames == 0
        fprintf('[%s] n_trials=%d  n_frames=0 -- skipped\n', label, n_trials);
        return;
    end
    labels = logical(labels);
    tpr_half = sum(scores >= 0.5 & labels)  / max(1, sum(labels));
    fpr_half = sum(scores >= 0.5 & ~labels) / max(1, sum(~labels));
    fprintf('[%s] n_trials=%d  n_frames=%d  AUC=%.3f  |  @p=0.5: sensitivity=%.0f%%  specificity=%.0f%%\n', ...
            label, n_trials, n_frames, auc_val, 100*tpr_half, 100*(1-fpr_half));
end

function save_roc_figure(out_file, label, n_trials, n_frames, fpr, tpr, auc_val, ...
                          scores, labels, col, fig_vis, show_figures)
% SAVE_ROC_FIGURE  One self-contained image: ROC curve with chance diagonal
%   and the p=0.5 natural operating point marked.
    fig = figure('Name', sprintf('ROC -- %s', label), 'Color', 'w', ...
                 'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig, 'Units', 'normalized', 'Position', [0 0 1 1]);
    ax = axes(fig); hold(ax, 'on');
    plot(ax, [0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', 'chance');
    plot(ax, fpr, tpr, '-', 'Color', col, 'LineWidth', 2.5, 'DisplayName', 'ROC');

    labels = logical(labels);
    tpr_half = sum(scores >= 0.5 & labels)  / max(1, sum(labels));
    fpr_half = sum(scores >= 0.5 & ~labels) / max(1, sum(~labels));
    plot(ax, fpr_half, tpr_half, 'o', 'MarkerSize', 10, 'MarkerFaceColor', col, ...
         'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'operating point p=0.5');

    xlabel(ax, 'False Positive Rate'); ylabel(ax, 'True Positive Rate');
    xlim(ax, [0 1]); ylim(ax, [0 1]); axis(ax, 'square'); grid(ax, 'on');
    legend(ax, 'Location', 'southeast', 'FontSize', 10);
    title(ax, sprintf('%s  --  AUC=%.3f  (n=%d trials, %d frames)', label, auc_val, n_trials, n_frames), ...
          'Interpreter', 'none', 'FontWeight', 'bold');

    saveas(fig, out_file, 'svg');
    fprintf('Saved %s\n', out_file);
    if ~show_figures, close(fig); end
end
