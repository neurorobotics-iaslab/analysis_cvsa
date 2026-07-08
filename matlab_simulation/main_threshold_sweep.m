%% MAIN_THRESHOLD_SWEEP  HIT / TIMEOUT rate as a function of both classes'
%   integrator hit-thresholds, swept independently over a fixed grid.
%
%   Re-evaluates each trial's ALREADY-SIMULATED integrator buffer
%   (trials(t).integrated, produced once by integrate_signal.m) against
%   every (threshold_class1, threshold_class2) combination on a grid from
%   TH_START to TH_END (fixed extremes, so every saved image spans the same
%   axis range regardless of the recording's real threshold) in TH_STEP
%   steps -- no re-running of FBCSP/sLDA/integrator per threshold, since
%   the buffer dynamics (Buffer.cpp leaky integrator) do not depend on the
%   hit-detection thresholds at all, only the pass/fail check does.
%
%   Outcome rule per (th1, th2), mirrors Training.cpp is_target_hit /
%   main_hybrid_advantage_integ.m's classify_trial_outcome: within the CF
%   window, the first class i whose integrated[i] >= th(i) - 5e-3 is
%   reached "wins" that trial -- HIT if that class is the cue's target,
%   MISS if it is the other class, TIMEOUT if neither threshold is ever
%   reached. This is deliberately NOT a classic ROC/AUC (there is no single
%   binary positive/negative label once two independent per-class
%   thresholds and a TIMEOUT outcome are in play) -- instead it is a 2D
%   heatmap over the (th1, th2) plane, which is the natural generalisation
%   here.
%
%   Works across a mix of MI/CVSA/Hybrid GDFs selected in one run. Each
%   paradigm actually present gets its OWN saved image (HIT rate + TIMEOUT
%   rate side by side), plus a pooled "ALL" image -- the buffer/threshold
%   logic is identical across paradigms, but MI/CVSA/Hybrid can behave very
%   differently (e.g. Hybrid's CVSA-fusion window), so they are never
%   silently merged into one figure. Additionally, every individual GDF
%   file gets its own image, saved in a "per_file/" subfolder.
%
%   Output, under <gdf_dir>/analysis_results/threshold_sweep/:
%     threshold_sweep_<paradigm>.svg   one per paradigm present
%     threshold_sweep_ALL.svg          pooled across all paradigms
%     per_file/threshold_sweep_<basename>.svg   one per selected GDF
%     threshold_sweep_summary.mat      all grids + real thresholds

function main_threshold_sweep(gdf_dir, gdf_names, show_figures)
%   Callable as a function: main_threshold_sweep(gdf_dir, gdf_names, show_figures)
%   With no arguments, shows the GUI file picker (interactive mode).

% --- Sweep configuration (fixed grid extremes, independent of paradigm/file)
TH_START     = 0.55;
TH_END       = 1.00;
TH_STEP      = 0.05;
if nargin < 3, show_figures = false; end
SHOW_FIGURES = show_figures;

TH_VALS = TH_START:TH_STEP:TH_END;   % same grid used for both classes' axes

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
        error('main_threshold_sweep:cancel', 'No GDF selected.');
    end
elseif nargin < 2 || isempty(gdf_names)
    f = dir(fullfile(gdf_dir, '*.gdf'));
    gdf_names = {f.name};
    if isempty(gdf_names), error('main_threshold_sweep:nofiles', 'No GDF files in %s', gdf_dir); end
end
if ischar(gdf_names), gdf_names = {gdf_names}; end
n_files = numel(gdf_names);

all_trials = struct([]);   % pooled across files: .integrated, .target_class, .n_pre, .n_cf, .paradigm, .basename
real_thresholds_per_file = struct([]);   % .basename, .thresholds -- so per-file plots mark their own real threshold

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

%% --- Integrate per trial (buffer dynamics do not depend on int_cfg.thresholds)
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

%% --- Record this file's real recording-time thresholds -----------------
thr_file = to_vec(int_cfg.thresholds);
idx_thr = numel(real_thresholds_per_file) + 1;
real_thresholds_per_file(idx_thr).basename   = basename; %#ok<AGROW>
real_thresholds_per_file(idx_thr).paradigm   = paradigm;
real_thresholds_per_file(idx_thr).thresholds = thr_file;

%% --- Pool usable trials (need a target class and at least one CF chunk)
for t = 1:numel(trials)
    tr = trials(t);
    if isnan(tr.target_class) || tr.n_cf == 0, continue; end
    idx_acc = numel(all_trials) + 1;
    all_trials(idx_acc).integrated   = tr.integrated; %#ok<AGROW>
    all_trials(idx_acc).target_class = tr.target_class;
    all_trials(idx_acc).n_pre        = tr.n_pre;
    all_trials(idx_acc).n_cf         = tr.n_cf;
    all_trials(idx_acc).paradigm     = paradigm;
    all_trials(idx_acc).basename     = basename;
end

end % file_idx loop

if isempty(all_trials)
    error('main_threshold_sweep:no_trials', 'No usable trials found across the selected file(s).');
end

n_trials = numel(all_trials);
fprintf('\nPooled %d trials across %d file(s).\n', n_trials, n_files);

out_dir      = fullfile(gdf_dir, 'analysis_results', 'threshold_sweep');
out_dir_file = fullfile(out_dir, 'per_file');
if ~exist(out_dir, 'dir'),      mkdir(out_dir);      end
if ~exist(out_dir_file, 'dir'), mkdir(out_dir_file); end

fig_vis = 'off'; if SHOW_FIGURES, fig_vis = 'on'; end

%% --- One image per paradigm + a pooled "ALL" image -----------------------
paradigms_present = unique({all_trials.paradigm}, 'stable');
group_names = [paradigms_present, {'ALL'}];

sweep_by_group = struct('name', {}, 'n', {}, 'acc_grid', {}, 'to_rate_grid', {});
for g = 1:numel(group_names)
    name = group_names{g};
    if strcmp(name, 'ALL')
        trials_g = all_trials;
        thr_g = mean_thresholds(real_thresholds_per_file, {real_thresholds_per_file.basename});
    else
        trials_g = all_trials(strcmp({all_trials.paradigm}, name));
        files_g  = {real_thresholds_per_file(strcmp({real_thresholds_per_file.paradigm}, name)).basename};
        thr_g = mean_thresholds(real_thresholds_per_file, files_g);
    end
    [acc_grid, to_rate_grid] = compute_sweep_grids(trials_g, TH_VALS, TH_VALS);

    sweep_by_group(g).name         = name; %#ok<AGROW>
    sweep_by_group(g).n            = numel(trials_g);
    sweep_by_group(g).acc_grid     = acc_grid;
    sweep_by_group(g).to_rate_grid = to_rate_grid;

    print_best(name, numel(trials_g), acc_grid, to_rate_grid, TH_VALS, thr_g);

    out_file = fullfile(out_dir, sprintf('threshold_sweep_%s.svg', name));
    save_sweep_figure(out_file, name, numel(trials_g), acc_grid, to_rate_grid, ...
        TH_VALS, TH_VALS, thr_g, fig_vis, SHOW_FIGURES);
end

%% --- One image per individual GDF file -----------------------------------
basenames_present = unique({all_trials.basename}, 'stable');
sweep_by_file = struct('basename', {}, 'paradigm', {}, 'n', {}, 'acc_grid', {}, 'to_rate_grid', {});
for f = 1:numel(basenames_present)
    bn = basenames_present{f};
    trials_f = all_trials(strcmp({all_trials.basename}, bn));
    par_f    = trials_f(1).paradigm;
    thr_f    = real_thresholds_per_file(strcmp({real_thresholds_per_file.basename}, bn)).thresholds;

    [acc_grid, to_rate_grid] = compute_sweep_grids(trials_f, TH_VALS, TH_VALS);

    sweep_by_file(f).basename     = bn; %#ok<AGROW>
    sweep_by_file(f).paradigm     = par_f;
    sweep_by_file(f).n            = numel(trials_f);
    sweep_by_file(f).acc_grid     = acc_grid;
    sweep_by_file(f).to_rate_grid = to_rate_grid;

    print_best(bn, numel(trials_f), acc_grid, to_rate_grid, TH_VALS, thr_f);

    out_file = fullfile(out_dir_file, sprintf('threshold_sweep_%s.svg', bn));
    save_sweep_figure(out_file, sprintf('%s (%s)', bn, par_f), numel(trials_f), ...
        acc_grid, to_rate_grid, TH_VALS, TH_VALS, thr_f, fig_vis, SHOW_FIGURES);
end

%% --- Save .mat for reuse --------------------------------------------------
threshold_sweep_summary = struct('th_vals', TH_VALS, ...
    'sweep_by_group', sweep_by_group, 'sweep_by_file', sweep_by_file, ...
    'real_thresholds_per_file', real_thresholds_per_file, ...
    'n_trials', n_trials, 'n_files', n_files);
save(fullfile(out_dir, 'threshold_sweep_summary.mat'), 'threshold_sweep_summary');
fprintf('Saved %s\n', fullfile(out_dir, 'threshold_sweep_summary.mat'));


end % main_threshold_sweep

% ── Local helpers (must be after all script statements) ────────────────────

function thr = mean_thresholds(real_thresholds_per_file, basenames)
% MEAN_THRESHOLDS  Average real threshold pair across a set of files (by
%   basename) -- used as the marker position for pooled/paradigm images
%   that can span files with slightly different real thresholds.
    mask = ismember({real_thresholds_per_file.basename}, basenames);
    thr_mat = cell2mat({real_thresholds_per_file(mask).thresholds}');
    thr = mean(thr_mat, 1);
end

function print_best(label, n, acc_grid, to_rate_grid, th_vals, real_thr)
% PRINT_BEST  Console summary: best HIT-rate cell + outcome at the (nearest
%   grid point to the) real recording-time threshold.
    if n == 0
        fprintf('[%s] n=0 -- skipped\n', label);
        return;
    end
    [best_acc, best_idx] = max(acc_grid(:));
    [bi1, bi2] = ind2sub(size(acc_grid), best_idx);
    [~, ri1] = min(abs(th_vals - real_thr(1)));
    [~, ri2] = min(abs(th_vals - real_thr(2)));
    fprintf(['[%s] n=%d  best HIT rate %.0f%% at threshold=[%.2f, %.2f]  ' ...
              '(real~=[%.2f, %.2f] -> HIT rate %.0f%%, TIMEOUT rate %.0f%%)\n'], ...
            label, n, 100*best_acc, th_vals(bi1), th_vals(bi2), ...
            real_thr(1), real_thr(2), 100*acc_grid(ri1, ri2), 100*to_rate_grid(ri1, ri2));
end

function [acc_grid, to_rate_grid] = compute_sweep_grids(trials_g, th1_vals, th2_vals)
% COMPUTE_SWEEP_GRIDS  HIT-rate and TIMEOUT-rate grids over (th1, th2) for
%   one group of trials, reusing each trial's already-simulated buffer
%   trace (no pipeline re-run per threshold).
    n1 = numel(th1_vals); n2 = numel(th2_vals);
    n  = numel(trials_g);
    hit_count = zeros(n1, n2);
    to_count  = zeros(n1, n2);

    for it = 1:n
        tr = trials_g(it);
        cf_range = (tr.n_pre + 1):(tr.n_pre + tr.n_cf);
        cmax1 = cummax(tr.integrated(cf_range, 1));
        cmax2 = cummax(tr.integrated(cf_range, 2));

        % First CF-frame index (1-based, within cf_range) where each class'
        % running max reaches each candidate threshold; Inf if never
        % reached. Equivalent to the first frame where the raw value itself
        % crosses the threshold, since cummax is non-decreasing.
        idx1 = first_cross_idx(cmax1, th1_vals);   % [1 x n1]
        idx2 = first_cross_idx(cmax2, th2_vals);   % [1 x n2]

        idx1_grid = repmat(idx1(:),  1, n2);
        idx2_grid = repmat(idx2(:)', n1, 1);

        timeout_grid = isinf(idx1_grid) & isinf(idx2_grid);
        % Tie-break favours class 1, matching Training.cpp's if/elseif
        % order (checks class 1 before class 2 at every frame).
        winner1_grid = (~timeout_grid) & (idx1_grid <= idx2_grid);
        winner2_grid = (~timeout_grid) & ~winner1_grid;

        if tr.target_class == 1
            hit_count = hit_count + winner1_grid;
        else
            hit_count = hit_count + winner2_grid;
        end
        to_count = to_count + timeout_grid;
    end

    if n == 0
        acc_grid = nan(n1, n2); to_rate_grid = nan(n1, n2);
    else
        acc_grid     = hit_count / n;
        to_rate_grid = to_count  / n;
    end
end

function idx = first_cross_idx(cmax_sig, th_vals)
% FIRST_CROSS_IDX  For each candidate threshold, the first frame index
%   (1-based, within the CF window) where the running max of the integrated
%   signal reaches threshold-5e-3 (same tolerance as Training.cpp), or Inf
%   if the threshold is never reached.
    idx = inf(1, numel(th_vals));
    for v = 1:numel(th_vals)
        k = find(cmax_sig >= th_vals(v) - 5e-3, 1, 'first');
        if ~isempty(k), idx(v) = k; end
    end
end

function save_sweep_figure(out_file, label, n, acc_grid, to_rate_grid, ...
                            th1_vals, th2_vals, real_thr, fig_vis, show_figures)
% SAVE_SWEEP_FIGURE  One self-contained image: HIT rate + TIMEOUT rate
%   heatmaps side by side for a single group (paradigm/ALL/file).
    fig = figure('Name', sprintf('Threshold sweep -- %s', label), 'Color', 'w', ...
                 'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig, 'Units', 'normalized', 'Position', [0 0 1 1]);

    ax1 = subplot(1, 2, 1);
    plot_sweep_heatmap(ax1, th1_vals, th2_vals, 100*acc_grid, real_thr, parula(256), 'HIT rate (%)');

    ax2 = subplot(1, 2, 2);
    plot_sweep_heatmap(ax2, th1_vals, th2_vals, 100*to_rate_grid, real_thr, flipud(autumn(256)), 'TIMEOUT rate (%)');

    sgtitle(fig, sprintf(['%s  (n=%d trials)\n' ...
            'red square = real recording-time thresholds [%.2f, %.2f]'], ...
            label, n, real_thr(1), real_thr(2)), 'Interpreter', 'none');

    saveas(fig, out_file, 'svg');
    fprintf('Saved %s\n', out_file);
    if ~show_figures, close(fig); end
end

function plot_sweep_heatmap(ax, th1_vals, th2_vals, grid_pct, real_thr, cmap, label)
% PLOT_SWEEP_HEATMAP  imagesc heatmap over the (th1, th2) grid with every
%   swept value as an axis tick, per-cell percentage annotated (grid is at
%   most ~10x10, so text stays readable), and the real-threshold cell
%   marked with a red square.
    imagesc(ax, 1:numel(th2_vals), 1:numel(th1_vals), grid_pct);
    set(ax, 'YDir', 'normal');
    colormap(ax, cmap);
    cb = colorbar(ax); ylabel(cb, label);
    set(ax, 'XTick', 1:numel(th2_vals), ...
            'XTickLabel', arrayfun(@(v) sprintf('%.2f', v), th2_vals, 'UniformOutput', false));
    set(ax, 'YTick', 1:numel(th1_vals), ...
            'YTickLabel', arrayfun(@(v) sprintf('%.2f', v), th1_vals, 'UniformOutput', false));
    xlabel(ax, 'threshold class 2'); ylabel(ax, 'threshold class 1');
    title(ax, label, 'FontWeight', 'bold');
    for i = 1:numel(th1_vals)
        for j = 1:numel(th2_vals)
            text(ax, j, i, sprintf('%.0f', grid_pct(i, j)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 7, ...
                 'Color', pick_text_color(grid_pct(i, j)));
        end
    end
    [~, i1] = min(abs(th1_vals - real_thr(1)));
    [~, i2] = min(abs(th2_vals - real_thr(2)));
    rectangle(ax, 'Position', [i2-0.5, i1-0.5, 1, 1], 'EdgeColor', 'r', 'LineWidth', 2.5);
    axis(ax, 'tight');
end

function c = pick_text_color(v)
% PICK_TEXT_COLOR  Black text on light cells, white on dark cells.
    if v > 50, c = [1 1 1]; else, c = [0 0 0]; end
end
