%% MAIN_HYBRID_ADVANTAGE  Offline analysis of HYBRID BCI evaluation recordings.
%
%   Pipeline identical to main_simulate.m up to the trials variable.
%   Hybrid-only: aborts if paradigm ~= 'hybrid'.
%
%   For each trial computes (over valid, non-NaN CF frames):
%     mean P_MI(target class)    — raw MI sLDA output
%     mean P_CVSA(target class)  — raw CVSA sLDA output
%     mean P_fused(target class) — Bayesian LOP output (integrator input)
%     mean buffer(target class)  — leaky-WTA integrator state
%
%   Prints a per-trial diagnostic table and a per-outcome summary.
%   Shows one figure: per-trial scatter coloured by outcome (HIT/MISS/TO).

clear; clc; close all;

% --- Display options -------------------------------------------------------
SHOW_FIGURES = false;   % true: figures pop up on screen; false: created hidden (export only)
if SHOW_FIGURES, fig_vis = 'on'; else, fig_vis = 'off'; end

% --- Make every subfolder visible to the simulator -----------------------
this_dir = '/home/paolo/bci_vr_ws/src/analysis_bci/matlab_simulation';
addpath(this_dir);
addpath(fullfile(this_dir, 'io'));
addpath(fullfile(this_dir, 'processing'));
addpath(fullfile(this_dir, 'artifacts'));
addpath(fullfile(this_dir, 'classifier'));
addpath(fullfile(this_dir, 'integrator'));
addpath(fullfile(this_dir, 'plotting'));
addpath(fullfile(this_dir, 'utils'));

% --- GUI: pick GDF(s); everything else flows from the sibling YAML(s) ---
default_dir = '/home/paolo/bci_vr_ws/recordings';

[gdf_names, gdf_dir] = uigetfile({'*.gdf', 'GDF recordings (*.gdf)'}, ...
                                'Select GDF recording(s)', default_dir, 'MultiSelect', 'on');
if isequal(gdf_names, 0)
    error('main_simulate:cancel', 'No GDF selected.');
end
if ischar(gdf_names), gdf_names = {gdf_names}; end
n_files = numel(gdf_names);

for file_idx = 1:n_files
gdf_path = fullfile(gdf_dir, gdf_names{file_idx});
fprintf('\n[%d/%d] %s\n', file_idx, n_files, gdf_names{file_idx});

%% --- Load GDF --------------------------------------------------------
[signal, header, basename] = load_gdf(gdf_path);

%% --- Load parameters YAML -----------------------------
[params, ~] = load_params_yaml(gdf_path);

paradigm  = params.integrator.paradigm;
if ~strcmp(paradigm, 'hybrid')
    fprintf('  skipping %s: paradigm=%s (hybrid-only)\n', basename, paradigm);
    continue;
end
fs        = double(params.acquisition.samplerate);
framerate = double(params.acquisition.framerate);
chunk_size = round(fs / framerate);
if abs(fs - header.SampleRate) > 1e-3
    log_step('main_simulate: YAML samplerate=%.1f != GDF samplerate=%.1f -> using GDF', ...
             fs, header.SampleRate);
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

signal = signal(:,1:nchannels); % remove the last 3 values: associated with ACC values

log_step('main_simulate: paradigm=%s, fs=%g, framerate=%g, chunk=%d, bufproc=%d, bufart=%d, nchannels=%d', ...
         paradigm, fs, framerate, chunk_size, bufsize_proc, bufsize_art, nchannels);

%% --- Load CSP + sLDA per used paradigm --------------------------------
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

%% --- Apply processing (one stream per paradigm) ----------------------
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

%% --- Apply artifact detection on raw signal --------------------------
art_cfg = params.ArtifactCfg.params;
art_cfg.EOG_ch_names = to_strcell(art_cfg.EOG_ch_names);
cfg_art = struct('samplerate', fs, 'chunk_size', chunk_size, 'bufsize_artifact', bufsize_art);
[art_flags, info_art] = detect_artifacts(signal, header, art_cfg, cfg_art);

%% --- Apply sLDA over the whole feature stream ------------------------
%   Note on alignment: features (apply_processing) and art_flags
%   (detect_artifacts) both live on the same chunk-index axis (k = 1..n_chunks).
%   The artifact ringbuf fills first (bufsize_art / chunk_size chunks earlier
%   than the processing ringbuf), so art_flags are already "valid" by the
%   time features become valid -- the +bufsize_proc/chunk offset is implicit
%   in starting integration only at the first 781, which is well past both
%   buffers' fill points.
n_chunks    = info_proc.n_chunks;
first_valid = info_proc.first_valid_chunk;
art_lead    = max(0, round((bufsize_proc - bufsize_art) / chunk_size));

p_mi_aligned   = []; if use_mi,   p_mi_aligned   = apply_slda(features_mi, slda_mi,   csp_mi.bands  ); end
p_cvsa_aligned = []; if use_cvsa, p_cvsa_aligned = apply_slda(features_cv, slda_cvsa, csp_cvsa.bands); end

log_step('main_simulate: streams aligned (n_chunks=%d, first_valid_proc=%d, art_lead=%d chunks)', ...
         n_chunks, first_valid, art_lead);

%% --- Integrate per trial ---------------------------------------------
int_cfg = params.integrator;
% Pull the dynamic_reconfigure-able fields with safe defaults
if ~isfield(int_cfg, 'increment'),               int_cfg.increment = 1; end
if ~isfield(int_cfg, 'thresholds_rejection'),    int_cfg.thresholds_rejection = []; end
if ~isfield(int_cfg, 'cvsa_influence'),          int_cfg.cvsa_influence = 2.5; end
if ~isfield(int_cfg, 'thresholds'),              int_cfg.thresholds = params.training_node.thresholds; end

% header_chunks for trial extents (use whichever paradigm produced a header)
if use_mi,   header_chunks = header_mi;
else,        header_chunks = header_cv;
end
header_chunks.framerate = framerate;

trials = integrate_signal(p_mi_aligned, p_cvsa_aligned, art_flags, ...
                          header_chunks, int_cfg, paradigm);

%% --- Read REAL outcomes from GDF events (897=HIT, 898=MISS, 899=TIMEOUT)
HIT_CODE_SIM = 897;  MISS_CODE_SIM = 898;  TIMEOUT_CODE_SIM = 899;  CF_CODE_SIM = 781;
POS_ev = header.EVENT.POS;  TYP_ev = header.EVENT.TYP;
cf_samp = POS_ev(TYP_ev == CF_CODE_SIM);
n_trials = numel(trials);
trial_outcome_real = zeros(1, n_trials);   % 0 = not found
for t = 1:min(n_trials, numel(cf_samp))
    after = POS_ev > cf_samp(t);
    idx   = find((TYP_ev==HIT_CODE_SIM | TYP_ev==MISS_CODE_SIM | TYP_ev==TIMEOUT_CODE_SIM) & after, 1);
    if ~isempty(idx), trial_outcome_real(t) = TYP_ev(idx); end
end
n_hit_real  = sum(trial_outcome_real == HIT_CODE_SIM);
n_miss_real = sum(trial_outcome_real == MISS_CODE_SIM);
n_to_real   = sum(trial_outcome_real == TIMEOUT_CODE_SIM);
log_step('main_simulate: GDF outcomes -> HIT=%d  MISS=%d  TIMEOUT=%d  (sim PASS=%d)', ...
         n_hit_real, n_miss_real, n_to_real, sum([trials.pass]));

%% ── Analysis of trials ───────────────────────────────────────────────────
% For each trial (split by target class):
%   mean P(target class)  — from MI, CVSA, fused streams — over valid CF frames
%   frame accuracy        — fraction of valid CF frames where P(target) > 0.5
%   mean buffer           — leaky-WTA integrator state (target class)

classes    = to_vec(int_cfg.classes);
n_cls      = numel(classes);
thresholds = to_vec(int_cfg.thresholds);
p_rest     = to_vec(int_cfg.init_val); p_rest = p_rest(1);
cls_names  = arrayfun(@(x) num2str(x), classes, 'UniformOutput', false);
cvsa_inf = int_cfg.cvsa_influence;
trial_outcome_real = trial_outcome_real(:);   % ensure column vector (same orientation as trial_class)

% Per-trial scalars (target class, valid non-NaN CF frames only)
mean_P_mi       = nan(n_trials, 1);
mean_P_cvsa     = nan(n_trials, 1);
mean_P_cvsa_inf = nan(n_trials, 1);
mean_P_fus      = nan(n_trials, 1);
mean_buf        = nan(n_trials, 1);
acc_mi          = nan(n_trials, 1);
acc_cvsa        = nan(n_trials, 1);
acc_cvsa_inf    = nan(n_trials, 1);
acc_fus         = nan(n_trials, 1);
acc_buf         = nan(n_trials, 1);   % fraction CF frames where buf(target) > p_rest
trial_class     = nan(n_trials, 1);

fprintf('\n  Trial diagnostics  (m=mean P(target), a=frame acc):\n');
fprintf('  %-4s  %-5s  %-4s  %-6s %-5s  %-6s %-5s  %-6s %-5s  %-6s %-5s  %-6s %-5s\n', ...
        '#', 'cue', 'res', 'mMI', 'aMI', 'mCV', 'aCV', 'mCVi', 'aCVi', 'mFus', 'aFus', 'mBuf', 'aBuf');

for t = 1:n_trials
    c  = trials(t).target_class;
    if isnan(c) || c < 1 || c > n_cls, continue; end
    trial_class(t) = c;
    np = trials(t).n_pre;
    nc = trials(t).n_cf;

    pm  = trials(t).p_mi      (np+1 : np+nc, :);   % [nc x 2]
    pc  = trials(t).p_cvsa    (np+1 : np+nc, :);   % [nc x 2]
    pf  = trials(t).raw       (np+1 : np+nc, :);   % [nc x 2]  fused LOP
    buf = trials(t).integrated(np+1 : np+nc, c);   % [nc x 1]  target class

    vld = ~any(isnan(pm) | isnan(pc) | isnan(pf), 2);
    vld_inf = vld;
    if ~any(vld), continue; end
    if size(vld, 1) > cvsa_inf*framerate
        vld_inf(cvsa_inf*framerate:end) = false;
    end

    mean_P_mi(t)       = mean(pm(vld, c));
    mean_P_cvsa(t)     = mean(pc(vld, c));
    mean_P_cvsa_inf(t) = mean(pc(vld_inf, c));
    mean_P_fus(t)      = mean(pf(vld, c));
    acc_mi(t)          = mean(pm(vld, c) > 0.5);
    acc_cvsa(t)        = mean(pc(vld, c) > 0.5);
    acc_cvsa_inf(t)    = mean(pc(vld_inf, c) > 0.5);
    acc_fus(t)         = mean(pf(vld, c) > 0.5);
    vld_buf            = ~isnan(buf);
    mean_buf(t)        = mean(buf(vld_buf));
    acc_buf(t)         = mean(buf(vld_buf) > p_rest);

    switch trial_outcome_real(t)
        case HIT_CODE_SIM,     oc_str = 'HIT';
        case MISS_CODE_SIM,    oc_str = 'MISS';
        case TIMEOUT_CODE_SIM, oc_str = 'TO';
        otherwise,             oc_str = '?';
    end
    fprintf('  %-4d  %-5d  %-4s  %-6.3f %-5.2f  %-6.3f %-5.2f  %-6.3f %-5.2f  %-6.3f %-5.2f  %-6.3f %-5.2f\n', ...
            t, trials(t).onset_code, oc_str, ...
            mean_P_mi(t), acc_mi(t), ...
            mean_P_cvsa(t), acc_cvsa(t), ...
            mean_P_cvsa_inf(t), acc_cvsa_inf(t), ...
            mean_P_fus(t), acc_fus(t), ...
            mean_buf(t), acc_buf(t));
end

% ── Per-class summary ─────────────────────────────────────────────────────
fprintf('\n  Per-class summary  (mean over trials of that class):\n');
fprintf('  %-8s  %-5s  %-9s  %-6s  %-9s  %-6s  %-12s  %-8s  %-9s  %-6s  %-9s  %-7s\n', ...
        'class', 'n', 'meanP_MI', 'accMI', 'meanP_CV', 'accCV', 'meanCVinf', 'accCVinf', 'meanP_fus', 'accFus', 'meanBuf', 'accBuf');
for c = 1:n_cls
    mask = trial_class == c;
    if ~any(mask), continue; end
    fprintf('  %-8s  %-5d  %-9.3f  %-6.2f  %-9.3f  %-6.2f  %-12.3f  %-8.2f  %-9.3f  %-6.2f  %-9.3f  %-7.2f\n', ...
            cls_names{c}, sum(mask), ...
            mean(mean_P_mi(mask),       'omitnan'), mean(acc_mi(mask),       'omitnan'), ...
            mean(mean_P_cvsa(mask),     'omitnan'), mean(acc_cvsa(mask),     'omitnan'), ...
            mean(mean_P_cvsa_inf(mask), 'omitnan'), mean(acc_cvsa_inf(mask), 'omitnan'), ...
            mean(mean_P_fus(mask),      'omitnan'), mean(acc_fus(mask),      'omitnan'), ...
            mean(mean_buf(mask),        'omitnan'), mean(acc_buf(mask),       'omitnan'));
end

% ── Per-outcome summary ───────────────────────────────────────────────────
fprintf('\n  Per-outcome summary  (mean over trials with that outcome):\n');
fprintf('  %-8s  %-5s  %-9s  %-6s  %-9s  %-6s  %-12s  %-8s  %-9s  %-6s  %-9s  %-7s\n', ...
        'outcome', 'n', 'meanP_MI', 'accMI', 'meanP_CV', 'accCV', 'meanCVinf', 'accCVinf', 'meanP_fus', 'accFus', 'meanBuf', 'accBuf');
for oc_info = {{HIT_CODE_SIM,'HIT'},{MISS_CODE_SIM,'MISS'},{TIMEOUT_CODE_SIM,'TIMEOUT'}}
    oc = oc_info{1}{1};  lbl = oc_info{1}{2};
    mask = trial_outcome_real == oc;
    if ~any(mask), continue; end
    fprintf('  %-8s  %-5d  %-9.3f  %-6.2f  %-9.3f  %-6.2f  %-12.3f  %-8.2f  %-9.3f  %-6.2f  %-9.3f  %-7.2f\n', ...
            lbl, sum(mask), ...
            mean(mean_P_mi(mask),       'omitnan'), mean(acc_mi(mask),       'omitnan'), ...
            mean(mean_P_cvsa(mask),     'omitnan'), mean(acc_cvsa(mask),     'omitnan'), ...
            mean(mean_P_cvsa_inf(mask), 'omitnan'), mean(acc_cvsa_inf(mask), 'omitnan'), ...
            mean(mean_P_fus(mask),      'omitnan'), mean(acc_fus(mask),      'omitnan'), ...
            mean(mean_buf(mask),        'omitnan'), mean(acc_buf(mask),       'omitnan'));
end
fprintf('  Thresholds: cls1=%.2f  cls2=%.2f  |  cvsa_influence=%.1f s\n', thresholds(1), thresholds(2), cvsa_inf);

%% ── Shared colours / markers for Figures 1–2 ────────────────────────────
COL_HIT  = [0.15 0.65 0.15];
COL_MISS = [0.90 0.50 0.10];
COL_TO   = [0.80 0.15 0.15];
COL_MI   = [0.85 0.30 0.10];
COL_CV   = [0.15 0.35 0.80];
COL_FUS  = [0.20 0.50 0.80];
COL_BUF  = [0.35 0.20 0.65];
mk_cls   = {'o', 's'};   % circle = class 1, square = class 2
mk_sz    = [65, 60];

%% ── Figure 1: per-trial MEAN of P(target class) ──────────────────────────
mean_data   = {mean_P_mi, mean_P_cvsa, mean_P_cvsa_inf, mean_P_fus, mean_buf};
mean_labels = {'mean P_{MI}', 'mean P_{CVSA}', ...
               sprintf('mean P_{CVSA}(%.1fs)', cvsa_inf), ...
               'mean P_{fused}', 'mean Buffer'};
mean_cols   = {COL_MI, COL_CV, COL_CV, COL_FUS, COL_BUF};
mean_ylims  = {[0.2,1.0],[0.2,1.0],[0.2,1.0],[0.2,1.0],[p_rest-0.05,1.0]};
mean_yref   = {0.5, 0.5, 0.5, 0.5, p_rest};

fig1 = figure('Name', sprintf('Mean P per trial — %s', basename), ...
              'Color', 'w', 'NumberTitle', 'off', 'Position', [40 40 1600 420], 'Visible', fig_vis);
set(fig1, 'Units','normalized', 'OuterPosition',[0 0 1 1]);
ax_m = gobjects(1,5);
for sp = 1:5
    ax_m(sp) = subplot(1,5,sp);
    draw_scatter_panel(ax_m(sp), mean_data{sp}, trial_class, trial_outcome_real, n_trials, n_cls, ...
                       mk_cls, mk_sz, COL_HIT, COL_MISS, COL_TO, mean_yref{sp}, mean_ylims{sp}, ...
                       'trial #', mean_labels{sp}, mean_cols{sp});
end
add_scatter_legend(ax_m(1), COL_HIT, COL_MISS, COL_TO, n_hit_real, n_miss_real, n_to_real, cls_names);
sgtitle(fig1, sprintf('%s  |  HIT=%d  MISS=%d  TO=%d  — mean over CF frames', ...
        basename, n_hit_real, n_miss_real, n_to_real), 'FontSize', 10, 'Interpreter', 'none');

%% ── Figure 2: per-trial ACCURACY (fraction frames > chance) ─────────────
acc_data   = {acc_mi, acc_cvsa, acc_cvsa_inf, acc_fus, acc_buf};
acc_labels = {'acc P_{MI} > 0.5', 'acc P_{CVSA} > 0.5', ...
              sprintf('acc P_{CVSA} > 0.5  (%.1fs)', cvsa_inf), ...
              'acc P_{fused} > 0.5', 'acc Buffer > p_{rest}'};
acc_cols   = {COL_MI, COL_CV, COL_CV, COL_FUS, COL_BUF};
acc_ylims  = repmat({[0, 1]}, 1, 5);
acc_yref   = {0.5, 0.5, 0.5, 0.5, 0.5};

fig2 = figure('Name', sprintf('Frame accuracy per trial — %s', basename), ...
              'Color', 'w', 'NumberTitle', 'off', 'Position', [60 60 1600 420], 'Visible', fig_vis);
set(fig2, 'Units','normalized', 'OuterPosition',[0 0 1 1]);
ax_a = gobjects(1,5);
for sp = 1:5
    ax_a(sp) = subplot(1,5,sp);
    draw_scatter_panel(ax_a(sp), acc_data{sp}, trial_class, trial_outcome_real, n_trials, n_cls, ...
                       mk_cls, mk_sz, COL_HIT, COL_MISS, COL_TO, acc_yref{sp}, acc_ylims{sp}, ...
                       'trial #', acc_labels{sp}, acc_cols{sp});
end
add_scatter_legend(ax_a(1), COL_HIT, COL_MISS, COL_TO, n_hit_real, n_miss_real, n_to_real, cls_names);
sgtitle(fig2, sprintf('%s  |  HIT=%d  MISS=%d  TO=%d  — frame accuracy (P > chance)', ...
        basename, n_hit_real, n_miss_real, n_to_real), 'FontSize', 10, 'Interpreter', 'none');

%% ── Figure 3: does CVSA fusion help? ────────────────────────────────────────

% ── Per-trial: CF duration (plotted separately in Figure 4) ────────────────
dur_s = nan(n_trials, 1);
for t = 1:n_trials
    c = trial_class(t);
    if isnan(c) || c < 1 || c > n_cls, continue; end
    dur_s(t) = trials(t).n_cf * chunk_size / fs;
end
max_dur = max(dur_s, [], 'omitnan');

% ── Per-trial time course of P_MI(target) vs P_fused(target), per class ────
max_nc = 0;
for t = 1:n_trials
    if ~isnan(trial_class(t)), max_nc = max(max_nc, trials(t).n_cf); end
end
time_axis = (0:max_nc-1) * chunk_size / fs;

mi_tc  = cell(1, n_cls);
fus_tc = cell(1, n_cls);
nc_cls = zeros(1, n_cls);   % longest CF (in chunks) among trials of each class
for c = 1:n_cls
    mi_tc{c}  = nan(n_trials, max_nc);
    fus_tc{c} = nan(n_trials, max_nc);
end

% ── Per-trial CVSA-fusion advantage + frame-level rescue/hurt counts,
%    evaluated over the CVSA-influence window where alpha > 0 ──────────────
fus_adv = nan(n_trials, 1);   % mean(P_fused - P_MI) on target class
n_resc     = 0;   % frames: P_MI(target)<=0.5 -> P_fused(target)>0.5  (CVSA rescued)
n_hurt_fr  = 0;   % frames: P_MI(target)>0.5  -> P_fused(target)<=0.5 (CVSA hurt)
n_both_ok  = 0;   % frames: both correct
n_both_bad = 0;   % frames: both wrong

for t = 1:n_trials
    c = trial_class(t);
    if isnan(c) || c < 1 || c > n_cls, continue; end
    np = trials(t).n_pre;
    nc = trials(t).n_cf;
    pm = trials(t).p_mi(np+1:np+nc, c);
    pf = trials(t).raw (np+1:np+nc, c);
    vld = ~isnan(pm) & ~isnan(pf);
    if ~any(vld), continue; end

    mi_tc{c}(t, 1:nc)  = pm;
    fus_tc{c}(t, 1:nc) = pf;
    nc_cls(c) = max(nc_cls(c), nc);

    vld_inf = vld;
    if numel(vld) > cvsa_inf * framerate
        vld_inf(round(cvsa_inf*framerate)+1:end) = false;
    end
    if any(vld_inf)
        fus_adv(t) = mean(pf(vld_inf) - pm(vld_inf));

        mi_ok  = pm(vld_inf) > 0.5;
        fus_ok = pf(vld_inf) > 0.5;
        n_resc     = n_resc     + sum(~mi_ok &  fus_ok);
        n_hurt_fr  = n_hurt_fr  + sum( mi_ok & ~fus_ok);
        n_both_ok  = n_both_ok  + sum( mi_ok &  fus_ok);
        n_both_bad = n_both_bad + sum(~mi_ok & ~fus_ok);
    end
end

n_def  = sum(~isnan(fus_adv));
n_help = sum(fus_adv > 0);
n_hurt = sum(fus_adv < 0);
fprintf('\n  CVSA-fusion advantage  (mean[P_fused(target) - P_MI(target)] over first %.1fs of CF):\n', cvsa_inf);
fprintf('  helped (>0): %d/%d   hurt (<0): %d/%d   mean delta: %+.3f\n', ...
        n_help, n_def, n_hurt, n_def, mean(fus_adv, 'omitnan'));

net_fr = n_resc - n_hurt_fr;
if net_fr > 0,     verdict_fr = 'CVSA helped';
elseif net_fr < 0, verdict_fr = 'CVSA hurt';
else,              verdict_fr = 'neutral';
end
fprintf('\n  CVSA-fusion frame-level effect  (within first %.1fs of CF):\n', cvsa_inf);
fprintf('  rescued (MI wrong -> fused correct): %d\n', n_resc);
fprintf('  hurt    (MI correct -> fused wrong): %d\n', n_hurt_fr);
fprintf('  both correct: %d   both wrong: %d\n', n_both_ok, n_both_bad);
fprintf('  net effect: %+d frames  ->  %s\n', net_fr, verdict_fr);

fig3 = figure('Name', sprintf('Does CVSA fusion help? — %s', basename), ...
              'Color', 'w', 'NumberTitle', 'off', 'Position', [80 80 1300 700], 'Visible', fig_vis);
set(fig3, 'Units','normalized', 'OuterPosition',[0 0 1 1]);

% ── Panel (1,1): frame-level rescue/hurt bar chart ──────────────────────────
ax1 = subplot(2,2,1);  hold(ax1,'on');
bar_cats = categorical({'rescued','hurt','both correct','both wrong'});
bar_cats = reordercats(bar_cats, {'rescued','hurt','both correct','both wrong'});
bar_vals = [n_resc, n_hurt_fr, n_both_ok, n_both_bad];
bar_cols = [0.15 0.65 0.15; 0.80 0.15 0.15; 0.55 0.55 0.55; 0.85 0.85 0.85];
b1 = bar(ax1, bar_cats, bar_vals, 'FaceColor', 'flat');
b1.CData = bar_cols;
for i = 1:numel(bar_vals)
    text(ax1, i, bar_vals(i), sprintf(' %d', bar_vals(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
end
ylabel(ax1, 'CF frames (within CVSA-infl. window)', 'FontSize', 9);
title(ax1, sprintf('CVSA effect on frame correctness (first %.1fs)\nnet = %+d frames  ->  %s', ...
      cvsa_inf, net_fr, verdict_fr), 'FontSize', 10, 'FontWeight', 'bold');
grid(ax1, 'on');

% ── Panel (1,2): per-trial CVSA-fusion advantage metric ─────────────────────
ax2 = subplot(2,2,2);  hold(ax2,'on');
for t = 1:n_trials
    if isnan(fus_adv(t)), continue; end
    c = trial_class(t);  if isnan(c)||c<1||c>n_cls, c=1; end
    switch trial_outcome_real(t)
        case HIT_CODE_SIM,     col = COL_HIT;
        case MISS_CODE_SIM,    col = COL_MISS;
        case TIMEOUT_CODE_SIM, col = COL_TO;
        otherwise,             col = [0.5 0.5 0.5];
    end
    scatter(ax2, t, fus_adv(t), mk_sz(c), col, 'filled', ...
            'Marker', mk_cls{c}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
end
yline(ax2, 0, 'k-', 'LineWidth', 1.2, 'HandleVisibility','off');
set(ax2, 'XLim',[0,n_trials+1], 'XTick',1:n_trials);
xlabel(ax2,'trial #','FontSize',9);
ylabel(ax2, 'mean(P_{fused} - P_{MI})  on target class', 'FontSize', 9);
title(ax2, sprintf('CVSA-fusion advantage (first %.1fs)  —  helped %d/%d, hurt %d/%d', ...
      cvsa_inf, n_help, n_def, n_hurt, n_def), 'FontSize', 10, 'FontWeight', 'bold');
grid(ax2,'on');

% ── Panels (2,1)/(2,2): P_MI(target) vs P_fused(target) over CF time, per class ──
for c = 1:n_cls
    ax = subplot(2,2,2+c);  hold(ax,'on');
    L = nc_cls(c);
    if L == 0
        text(ax, 0.5, 0.5, 'no trials', 'Units','normalized', 'HorizontalAlignment','center');
        title(ax, sprintf('Class %s  (n=0)', cls_names{c}), 'FontSize', 10, 'FontWeight', 'bold');
        continue;
    end
    ta    = time_axis(1:L);
    m_mi  = mean(mi_tc{c}(:,1:L),  1, 'omitnan');
    m_fus = mean(fus_tc{c}(:,1:L), 1, 'omitnan');
    n_eff = sum(~isnan(mi_tc{c}(:,1:L)), 1);
    s_mi  = std(mi_tc{c}(:,1:L),  0, 1, 'omitnan') ./ sqrt(max(n_eff,1));
    s_fus = std(fus_tc{c}(:,1:L), 0, 1, 'omitnan') ./ sqrt(max(n_eff,1));

    fill(ax, [ta, fliplr(ta)], [m_mi-s_mi, fliplr(m_mi+s_mi)], ...
         COL_MI, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill(ax, [ta, fliplr(ta)], [m_fus-s_fus, fliplr(m_fus+s_fus)], ...
         COL_FUS, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, ta, m_mi,  '-', 'Color', COL_MI,  'LineWidth', 2, 'DisplayName', 'P_{MI}(target)');
    plot(ax, ta, m_fus, '-', 'Color', COL_FUS, 'LineWidth', 2, 'DisplayName', 'P_{fused}(target)');
    yline(ax, 0.5, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    xline(ax, cvsa_inf, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
          'DisplayName', sprintf('CVSA infl. ends (%.1fs)', cvsa_inf));
    set(ax, 'XLim', [0, max(ta)], 'YLim', [0, 1]);
    xlabel(ax, 'time from CF onset (s)', 'FontSize', 9);
    ylabel(ax, 'P(target class)', 'FontSize', 9);
    title(ax, sprintf('Class %s  (n=%d)', cls_names{c}, sum(trial_class==c & ~isnan(fus_adv))), ...
          'FontSize', 10, 'FontWeight', 'bold');
    legend(ax, 'FontSize', 8, 'Location', 'best');
    grid(ax, 'on');
end

sgtitle(fig3, sprintf('%s  |  HIT=%d  MISS=%d  TO=%d  — does CVSA fusion help?', ...
        basename, n_hit_real, n_miss_real, n_to_real), 'FontSize', 10, 'Interpreter', 'none');

%% ── Figure 4: CF trial duration + CVSA influence over time, by trial length ─
fig4 = figure('Name', sprintf('CF duration & CVSA time-influence — %s', basename), ...
              'Color', 'w', 'NumberTitle', 'off', 'Position', [120 120 950 420], 'Visible', fig_vis);
set(fig4, 'Units','normalized', 'OuterPosition',[0 0 1 1]);

% ── Panel 1: CF trial duration per trial (colored by outcome) ──────────────
ax4a = subplot(1,2,1);  hold(ax4a,'on');
for t = 1:n_trials
    if isnan(dur_s(t)), continue; end
    c = trial_class(t);  if isnan(c)||c<1||c>n_cls, c=1; end
    switch trial_outcome_real(t)
        case HIT_CODE_SIM,     col = COL_HIT;
        case MISS_CODE_SIM,    col = COL_MISS;
        case TIMEOUT_CODE_SIM, col = COL_TO;
        otherwise,             col = [0.5 0.5 0.5];
    end
    scatter(ax4a, t, dur_s(t), mk_sz(c), col, 'filled', ...
            'Marker', mk_cls{c}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
end
for oc = {HIT_CODE_SIM, MISS_CODE_SIM, TIMEOUT_CODE_SIM}
    mask = trial_outcome_real == oc{1};
    v = dur_s(mask);  v = v(~isnan(v));
    if isempty(v), continue; end
    switch oc{1}, case HIT_CODE_SIM, lc=COL_HIT; case MISS_CODE_SIM, lc=COL_MISS; case TIMEOUT_CODE_SIM, lc=COL_TO; end
    yline(ax4a, mean(v), '--', 'Color', lc, 'LineWidth', 1.5, 'HandleVisibility','off');
end
yline(ax4a, cvsa_inf, 'k:', 'LineWidth', 1.2, 'DisplayName', sprintf('cvsa_{influence} (%.1fs)', cvsa_inf));
h4 = gobjects(1,3);
h4(1) = scatter(ax4a,NaN,NaN,50,COL_HIT, 'filled','o','MarkerEdgeColor','k','DisplayName',sprintf('HIT (%d)', n_hit_real));
h4(2) = scatter(ax4a,NaN,NaN,50,COL_MISS,'filled','o','MarkerEdgeColor','k','DisplayName',sprintf('MISS (%d)',n_miss_real));
h4(3) = scatter(ax4a,NaN,NaN,50,COL_TO,  'filled','o','MarkerEdgeColor','k','DisplayName',sprintf('TO (%d)',  n_to_real));
legend(ax4a, h4, 'FontSize', 7, 'Location', 'best');
set(ax4a, 'XLim',[0,n_trials+1], 'XTick',1:n_trials, 'YLim',[0, max_dur*1.1]);
xlabel(ax4a,'trial #','FontSize',9);  ylabel(ax4a,'duration (s)','FontSize',9);
title(ax4a, sprintf('%s  —  CF trial duration', basename), 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter','none');
grid(ax4a,'on');

% ── Panel 2: mean P_MI -> P_fused for HIT trials, by hit-timing category ───
% MISS/TIMEOUT trials are excluded. For HIT trials, dur_s(t) is the time
% from CF onset to hit, so it splits trials into those that hit within the
% CVSA-influence window vs those that hit after it. For each group, frames
% within the first cvsa_inf seconds of CF are pooled (same window as Fig 5).
hit_cat_labels = {sprintf('HIT within %.1fs', cvsa_inf), sprintf('HIT after %.1fs', cvsa_inf)};
hit_cat_short  = {sprintf('hit < %.1fs', cvsa_inf), sprintf('hit \\geq %.1fs', cvsa_inf)};
hit_cat_cols   = [COL_CV; COL_BUF];

hit_pmi  = {[], []};
hit_pfus = {[], []};

for t = 1:n_trials
    if trial_outcome_real(t) ~= HIT_CODE_SIM, continue; end
    c = trial_class(t);
    if isnan(c) || c < 1 || c > n_cls, continue; end
    if isnan(dur_s(t)), continue; end
    np = trials(t).n_pre;  nc = trials(t).n_cf;
    pm = trials(t).p_mi(np+1:np+nc, c);
    pf = trials(t).raw (np+1:np+nc, c);
    vld = ~isnan(pm) & ~isnan(pf);

    n_inf = min(nc, round(cvsa_inf * framerate));
    idx = find(vld(1:n_inf));
    if isempty(idx), continue; end

    if dur_s(t) < cvsa_inf, k = 1; else, k = 2; end
    hit_pmi{k}  = [hit_pmi{k};  pm(idx)];
    hit_pfus{k} = [hit_pfus{k}; pf(idx)];
end

n_hit_fr     = [numel(hit_pmi{1}), numel(hit_pmi{2})];
mean_hit_mi  = nan(1,2);  mean_hit_fus = nan(1,2);
sem_hit_mi   = nan(1,2);  sem_hit_fus  = nan(1,2);
mean_hit_d   = nan(1,2);
for k = 1:2
    if n_hit_fr(k) == 0, continue; end
    mean_hit_mi(k)  = mean(hit_pmi{k});
    mean_hit_fus(k) = mean(hit_pfus{k});
    sem_hit_mi(k)   = std(hit_pmi{k})  / sqrt(n_hit_fr(k));
    sem_hit_fus(k)  = std(hit_pfus{k}) / sqrt(n_hit_fr(k));
    mean_hit_d(k)   = mean(hit_pfus{k} - hit_pmi{k});
end

n_hit_within = sum(trial_outcome_real == HIT_CODE_SIM & ~isnan(trial_class(:)) & dur_s <  cvsa_inf);
n_hit_after  = sum(trial_outcome_real == HIT_CODE_SIM & ~isnan(trial_class(:)) & dur_s >= cvsa_inf);
n_hit_trials = [n_hit_within, n_hit_after];

fprintf('\n  HIT trials, P_MI -> P_fused within first %.1fs of CF, by hit timing:\n', cvsa_inf);
for k = 1:2
    if n_hit_fr(k) == 0
        fprintf('  %-18s trials=%3d   (no valid frames)\n', hit_cat_labels{k}, n_hit_trials(k));
        continue;
    end
    fprintf('  %-18s trials=%3d   frames=%5d   P_MI=%.3f -> P_fused=%.3f   (delta=%+.3f)\n', ...
            hit_cat_labels{k}, n_hit_trials(k), n_hit_fr(k), mean_hit_mi(k), mean_hit_fus(k), mean_hit_d(k));
end

ax4b = subplot(1,2,2); hold(ax4b,'on');
for k = 1:2
    if n_hit_fr(k) == 0, continue; end
    plot(ax4b, [k-0.08, k+0.08], [mean_hit_mi(k), mean_hit_fus(k)], '-', ...
         'Color', hit_cat_cols(k,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    errorbar(ax4b, k-0.08, mean_hit_mi(k),  sem_hit_mi(k),  'o', 'Color', hit_cat_cols(k,:), ...
             'MarkerFaceColor', hit_cat_cols(k,:), 'CapSize', 3, 'HandleVisibility', 'off');
    errorbar(ax4b, k+0.08, mean_hit_fus(k), sem_hit_fus(k), 's', 'Color', hit_cat_cols(k,:), ...
             'MarkerFaceColor', 'w', 'LineWidth', 1.2, 'CapSize', 3, 'HandleVisibility', 'off');
    text(ax4b, k, max(mean_hit_mi(k),mean_hit_fus(k)) + 0.04, sprintf('%+.3f', mean_hit_d(k)), ...
         'HorizontalAlignment', 'center', 'FontSize', 8);
end
yline(ax4b, 0.5, 'k:', 'HandleVisibility', 'off');
h4b = gobjects(1,2);
h4b(1) = scatter(ax4b, NaN, NaN, 40, [0.3 0.3 0.3], 'filled', 'o', 'DisplayName', 'P_{MI}');
h4b(2) = scatter(ax4b, NaN, NaN, 40, [0.3 0.3 0.3], 'Marker', 's', 'DisplayName', 'P_{fused}');
legend(ax4b, h4b, 'FontSize', 8, 'Location', 'best');
set(ax4b, 'XLim', [0.5,2.5], 'YLim', [0,1.1], 'XTick', 1:2, 'XTickLabel', ...
    {sprintf('%s (n=%d)', hit_cat_short{1}, n_hit_within), sprintf('%s (n=%d)', hit_cat_short{2}, n_hit_after)});
ylabel(ax4b, 'P(target class)', 'FontSize', 9);
title(ax4b, sprintf('HIT trials: mean P_{MI} \\rightarrow P_{fused}\n(within first %.1fs of CF)', cvsa_inf), ...
      'FontSize', 10, 'FontWeight', 'bold');
grid(ax4b, 'on');

%% ── Figure 5: magnitude of the CVSA-fusion effect, frame-by-frame ──────────
% Per-frame data within the CVSA-influence window (alpha > 0), all trials.

fr_pmi   = [];
fr_pcvsa = [];
fr_pfus  = [];
fr_alpha = [];
fr_class = [];

for t = 1:n_trials
    c = trial_class(t);
    if isnan(c) || c < 1 || c > n_cls, continue; end
    np = trials(t).n_pre;
    nc = trials(t).n_cf;
    pm = trials(t).p_mi  (np+1:np+nc, c);
    pc = trials(t).p_cvsa(np+1:np+nc, c);
    pf = trials(t).raw   (np+1:np+nc, c);
    vld = ~isnan(pm) & ~isnan(pc) & ~isnan(pf);

    n_inf = min(nc, round(cvsa_inf * framerate));
    idx = find(vld(1:n_inf));
    if isempty(idx), continue; end

    t_sec = (idx-1) * chunk_size / fs;
    alpha = 0.5 * (1 + cos(pi * min(t_sec, cvsa_inf) / cvsa_inf));

    fr_pmi   = [fr_pmi;   pm(idx)];
    fr_pcvsa = [fr_pcvsa; pc(idx)];
    fr_pfus  = [fr_pfus;  pf(idx)];
    fr_alpha = [fr_alpha; alpha];
    fr_class = [fr_class; repmat(c, numel(idx), 1)];
end

fr_delta = fr_pfus - fr_pmi;     % signed effect of fusion on P(target class)

% ── Classify every frame by MI/CVSA agreement on the target class ──────────
mi_ok   = fr_pmi   > 0.5;
cvsa_ok = fr_pcvsa > 0.5;

fr_cat = zeros(size(fr_pmi));
fr_cat(mi_ok  & cvsa_ok)  = 1;   % agree, both correct
fr_cat(mi_ok  & ~cvsa_ok) = 2;   % MI correct, CVSA wrong  -> "cost?"
fr_cat(~mi_ok & cvsa_ok)  = 3;   % MI wrong, CVSA correct  -> "rescue?"
fr_cat(~mi_ok & ~cvsa_ok) = 4;   % agree, both wrong

cat_labels  = {'agree: both correct', 'MI correct, CVSA wrong', 'MI wrong, CVSA correct', 'agree: both wrong'};
cat_short   = {'agree-ok', 'MI ok / CVSA no', 'MI no / CVSA ok', 'agree-wrong'};
cat_cols    = [0.20 0.60 0.20; 0.85 0.55 0.10; 0.20 0.40 0.85; 0.70 0.20 0.20];

n_fr     = numel(fr_delta);
n_cat    = zeros(1,4);
mean_mi  = nan(1,4);
mean_fus = nan(1,4);
sem_mi   = nan(1,4);
sem_fus  = nan(1,4);
mean_d   = nan(1,4);

fprintf('\n  CVSA vs MI agreement breakdown  (within first %.1fs of CF, %d frames):\n', cvsa_inf, n_fr);
for k = 1:4
    m = fr_cat == k;
    n_cat(k) = sum(m);
    if n_cat(k) == 0, continue; end
    mean_mi(k)  = mean(fr_pmi(m));
    mean_fus(k) = mean(fr_pfus(m));
    sem_mi(k)   = std(fr_pmi(m))  / sqrt(n_cat(k));
    sem_fus(k)  = std(fr_pfus(m)) / sqrt(n_cat(k));
    mean_d(k)   = mean(fr_delta(m));
    fprintf('  %-24s n=%5d (%4.1f%%)   P_MI=%.3f -> P_fused=%.3f   (delta=%+.3f)\n', ...
            cat_labels{k}, n_cat(k), 100*n_cat(k)/n_fr, mean_mi(k), mean_fus(k), mean_d(k));
end
fprintf('  Rescue effect (MI wrong, CVSA correct): delta=%+.3f (n=%d)\n', mean_d(3), n_cat(3));
fprintf('  Cost effect   (MI correct, CVSA wrong): delta=%+.3f (n=%d)\n', mean_d(2), n_cat(2));

fig5 = figure('Name', sprintf('CVSA/MI agreement & fusion effect — %s', basename), ...
              'Color', 'w', 'NumberTitle', 'off', 'Position', [80 80 1300 420], 'Visible', fig_vis);
set(fig5, 'Units','normalized', 'OuterPosition',[0 0 1 1]);

% ── Panel 1: P_MI(target) vs P_fused(target), per CF frame, by agreement ────
ax1 = subplot(1,3,1); hold(ax1,'on');
for k = 1:4
    m = fr_cat == k;
    if ~any(m), continue; end
    scatter(ax1, fr_pmi(m), fr_pfus(m), 10, cat_cols(k,:), 'filled', ...
            'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none', ...
            'DisplayName', sprintf('%s (n=%d)', cat_labels{k}, n_cat(k)));
end
plot(ax1, [0 1], [0 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(ax1, 0.5, 'k:', 'HandleVisibility', 'off');
yline(ax1, 0.5, 'k:', 'HandleVisibility', 'off');
set(ax1, 'XLim', [0,1], 'YLim', [0,1]);
xlabel(ax1, 'P_{MI}(target)', 'FontSize', 9);
ylabel(ax1, 'P_{fused}(target)', 'FontSize', 9);
title(ax1, 'Per-frame effect, colored by MI/CVSA agreement', 'FontSize', 10, 'FontWeight', 'bold');
legend(ax1, 'FontSize', 7, 'Location', 'best');
grid(ax1, 'on'); axis(ax1, 'square');

% ── Panel 2: mean P_MI -> P_fused per agreement category (slopegraph) ───────
ax2 = subplot(1,3,2); hold(ax2,'on');
for k = 1:4
    if n_cat(k) == 0, continue; end
    plot(ax2, [k-0.08, k+0.08], [mean_mi(k), mean_fus(k)], '-', ...
         'Color', cat_cols(k,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    errorbar(ax2, k-0.08, mean_mi(k),  sem_mi(k),  'o', 'Color', cat_cols(k,:), ...
             'MarkerFaceColor', cat_cols(k,:), 'CapSize', 3, 'HandleVisibility', 'off');
    errorbar(ax2, k+0.08, mean_fus(k), sem_fus(k), 's', 'Color', cat_cols(k,:), ...
             'MarkerFaceColor', 'w', 'LineWidth', 1.2, 'CapSize', 3, 'HandleVisibility', 'off');
    text(ax2, k, max(mean_mi(k),mean_fus(k)) + 0.04, sprintf('%+.3f', mean_d(k)), ...
         'HorizontalAlignment', 'center', 'FontSize', 8);
end
yline(ax2, 0.5, 'k:', 'HandleVisibility', 'off');
h2 = gobjects(1,2);
h2(1) = scatter(ax2, NaN, NaN, 40, [0.3 0.3 0.3], 'filled', 'o', 'DisplayName', 'P_{MI}');
h2(2) = scatter(ax2, NaN, NaN, 40, [0.3 0.3 0.3], 'Marker', 's', 'DisplayName', 'P_{fused}');
legend(ax2, h2, 'FontSize', 8, 'Location', 'best');
set(ax2, 'XLim', [0.5,4.5], 'YLim', [0,1.1], 'XTick', 1:4, 'XTickLabel', cat_short);
xtickangle(ax2, 12);
ylabel(ax2, 'P(target class)', 'FontSize', 9);
title(ax2, 'Mean P_{MI} \rightarrow P_{fused}  per agreement category', 'FontSize', 10, 'FontWeight', 'bold');
grid(ax2, 'on');

% ── Panel 3: rescue vs cost effect, scaling with CVSA weight alpha(t) ───────
ax3 = subplot(1,3,3); hold(ax3,'on');
m3 = fr_cat == 3;   % rescue candidates: MI wrong, CVSA correct
m2 = fr_cat == 2;   % cost candidates:   MI correct, CVSA wrong
scatter(ax3, fr_alpha(m3), fr_delta(m3), 12, cat_cols(3,:), 'filled', 'MarkerFaceAlpha', 0.35, ...
        'DisplayName', sprintf('rescue: MI wrong, CVSA right (n=%d)', n_cat(3)));
scatter(ax3, fr_alpha(m2), fr_delta(m2), 12, cat_cols(2,:), 'filled', 'MarkerFaceAlpha', 0.35, ...
        'DisplayName', sprintf('cost: MI right, CVSA wrong (n=%d)', n_cat(2)));
yline(ax3, 0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
set(ax3, 'XLim', [0,1]);
xlabel(ax3, '\alpha(t)  (CVSA weight in fusion)', 'FontSize', 9);
ylabel(ax3, 'P_{fused}-P_{MI}  on target class', 'FontSize', 9);
title(ax3, sprintf('Rescue vs cost effect vs \\alpha\nrescue mean=%+.3f | cost mean=%+.3f', mean_d(3), mean_d(2)), ...
      'FontSize', 10, 'FontWeight', 'bold');
legend(ax3, 'FontSize', 8, 'Location', 'best');
grid(ax3, 'on');

sgtitle(fig5, sprintf('%s  |  HIT=%d  MISS=%d  TO=%d  — does CVSA help, broken down by MI/CVSA agreement (first %.1fs of CF)', ...
        basename, n_hit_real, n_miss_real, n_to_real, cvsa_inf), 'FontSize', 10, 'Interpreter', 'none');

%% --- Save figures -------------------------------------------------------
out_dir = fullfile(gdf_dir, 'analysis_results', 'advantage_hybrid');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
saveas(fig1, fullfile(out_dir, sprintf('advantage_%s_meanP.svg',          basename)), 'svg');
saveas(fig2, fullfile(out_dir, sprintf('advantage_%s_frame_accuracy.svg', basename)), 'svg');
saveas(fig3, fullfile(out_dir, sprintf('advantage_%s_cvsa_fusion.svg',    basename)), 'svg');
saveas(fig4, fullfile(out_dir, sprintf('advantage_%s_cf_duration.svg',    basename)), 'svg');
saveas(fig5, fullfile(out_dir, sprintf('advantage_%s_agreement.svg',      basename)), 'svg');
fprintf('Saved figures to %s\n', out_dir);
close([fig1, fig2, fig3, fig4, fig5]);

end % file_idx loop


% ── Local helpers (must be after all script statements) ───────────────────

function draw_scatter_panel(ax, dat, trial_class, trial_outcome_real, n_trials, n_cls, ...
                             mk_cls, mk_sz, COL_HIT, COL_MISS, COL_TO, yref, ylim_sp, ...
                             xlabel_str, title_str, title_col)
    hold(ax, 'on');
    for t = 1:n_trials
        val = dat(t);  if isnan(val), continue; end
        c = trial_class(t);  if isnan(c)||c<1||c>n_cls, c=1; end
        switch trial_outcome_real(t)
            case 897, col = COL_HIT;
            case 898, col = COL_MISS;
            case 899, col = COL_TO;
            otherwise, col = [0.5 0.5 0.5];
        end
        scatter(ax, t, val, mk_sz(c), col, 'filled', ...
                'Marker', mk_cls{c}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end
    for oc = {897, 898, 899}
        mask = trial_outcome_real == oc{1};
        v = dat(mask);  v = v(~isnan(v));
        if isempty(v), continue; end
        switch oc{1}
            case 897, lc = COL_HIT;
            case 898, lc = COL_MISS;
            case 899, lc = COL_TO;
        end
        yline(ax, mean(v), '--', 'Color', lc, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    yline(ax, yref, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    set(ax, 'XLim', [0,n_trials+1], 'YLim', ylim_sp, 'XTick', 1:n_trials);
    xlabel(ax, xlabel_str, 'FontSize', 9);
    title(ax, title_str, 'FontSize', 10, 'FontWeight', 'bold', 'Color', title_col);
    grid(ax, 'on');
end

function add_scatter_legend(ax_leg, COL_HIT, COL_MISS, COL_TO, n_hit, n_miss, n_to, cls_names)
    h = gobjects(1, 5);
    h(1) = scatter(ax_leg, NaN, NaN, 55, COL_HIT,  'filled', 'o', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('HIT (%d)',  n_hit));
    h(2) = scatter(ax_leg, NaN, NaN, 55, COL_MISS, 'filled', 'o', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('MISS (%d)', n_miss));
    h(3) = scatter(ax_leg, NaN, NaN, 55, COL_TO,   'filled', 'o', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('TO (%d)',   n_to));
    h(4) = scatter(ax_leg, NaN, NaN, 55, [0.5 0.5 0.5], 'filled', 'o', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('cls1 (%s)', cls_names{1}));
    h(5) = scatter(ax_leg, NaN, NaN, 55, [0.5 0.5 0.5], 'filled', 's', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('cls2 (%s)', cls_names{2}));
    legend(ax_leg, h, 'FontSize', 8, 'Location', 'southoutside', 'NumColumns', 5);
end
