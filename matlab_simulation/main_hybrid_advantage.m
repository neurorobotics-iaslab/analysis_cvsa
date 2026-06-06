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

% --- GUI: pick the GDF; everything else flows from the sibling YAML -----
default_dir = '/home/paolo/bci_vr_ws/recordings';

[gdf_name, gdf_dir] = uigetfile({'*.gdf', 'GDF recordings (*.gdf)'}, ...
                                'Select a GDF recording', default_dir);
if isequal(gdf_name, 0)
    error('main_simulate:cancel', 'No GDF selected.');
end
gdf_path = fullfile(gdf_dir, gdf_name);

%% --- Load GDF --------------------------------------------------------
[signal, header, basename] = load_gdf(gdf_path);

%% --- Load parameters YAML -----------------------------
[params, ~] = load_params_yaml(gdf_path);

paradigm  = params.integrator.paradigm;
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
p_rest     = to_vec(int_cfg.init_val)(1);
cls_names  = arrayfun(@(x) num2str(x), classes, 'UniformOutput', false);

% Per-trial scalars (target class, valid non-NaN CF frames only)
mean_P_mi   = nan(n_trials, 1);
mean_P_cvsa = nan(n_trials, 1);
mean_P_fus  = nan(n_trials, 1);
mean_buf    = nan(n_trials, 1);
acc_mi      = nan(n_trials, 1);   % fraction frames P_MI(target)   > 0.5
acc_cvsa    = nan(n_trials, 1);   % fraction frames P_CVSA(target) > 0.5
acc_fus     = nan(n_trials, 1);   % fraction frames P_fused(target) > 0.5
trial_class = nan(n_trials, 1);

fprintf('\n  Trial diagnostics:\n');
fprintf('  %-5s  %-8s  %-6s  %-9s  %-6s  %-9s  %-6s  %-9s  %-6s  %-8s\n', ...
        '#', 'cue', 'result', 'meanP_MI', 'accMI', 'meanP_CV', 'accCV', ...
        'meanP_fus', 'accFus', 'meanBuf');

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
    if ~any(vld), continue; end

    mean_P_mi(t)   = mean(pm(vld, c));
    mean_P_cvsa(t) = mean(pc(vld, c));
    mean_P_fus(t)  = mean(pf(vld, c));
    acc_mi(t)      = mean(pm(vld, c) > 0.5);
    acc_cvsa(t)    = mean(pc(vld, c) > 0.5);
    acc_fus(t)     = mean(pf(vld, c) > 0.5);
    vld_buf        = ~isnan(buf);
    mean_buf(t)    = mean(buf(vld_buf));

    switch trial_outcome_real(t)
        case HIT_CODE_SIM,     oc_str = 'HIT';
        case MISS_CODE_SIM,    oc_str = 'MISS';
        case TIMEOUT_CODE_SIM, oc_str = 'TO';
        otherwise,             oc_str = '?';
    end
    fprintf('  %-5d  %-8d  %-6s  %-9.3f  %-6.2f  %-9.3f  %-6.2f  %-9.3f  %-6.2f  %-8.3f\n', ...
            t, trials(t).onset_code, oc_str, ...
            mean_P_mi(t), acc_mi(t), ...
            mean_P_cvsa(t), acc_cvsa(t), ...
            mean_P_fus(t), acc_fus(t), ...
            mean_buf(t));
end

% ── Per-class summary ─────────────────────────────────────────────────────
fprintf('\n  Per-class summary  (mean over trials of that class):\n');
fprintf('  %-8s  %-5s  %-9s  %-6s  %-9s  %-6s  %-9s  %-6s  %-8s\n', ...
        'class', 'n', 'meanP_MI', 'accMI', 'meanP_CV', 'accCV', 'meanP_fus', 'accFus', 'meanBuf');
for c = 1:n_cls
    mask = trial_class == c;
    if ~any(mask), continue; end
    fprintf('  %-8s  %-5d  %-9.3f  %-6.2f  %-9.3f  %-6.2f  %-9.3f  %-6.2f  %-8.3f\n', ...
            cls_names{c}, sum(mask), ...
            mean(mean_P_mi(mask),   'omitnan'), mean(acc_mi(mask),   'omitnan'), ...
            mean(mean_P_cvsa(mask), 'omitnan'), mean(acc_cvsa(mask), 'omitnan'), ...
            mean(mean_P_fus(mask),  'omitnan'), mean(acc_fus(mask),  'omitnan'), ...
            mean(mean_buf(mask),    'omitnan'));
end

% ── Per-outcome summary ───────────────────────────────────────────────────
fprintf('\n  Per-outcome summary  (mean over trials with that outcome):\n');
fprintf('  %-8s  %-5s  %-9s  %-6s  %-9s  %-6s  %-9s  %-6s  %-8s\n', ...
        'outcome', 'n', 'meanP_MI', 'accMI', 'meanP_CV', 'accCV', 'meanP_fus', 'accFus', 'meanBuf');
for oc_info = {{HIT_CODE_SIM,'HIT'},{MISS_CODE_SIM,'MISS'},{TIMEOUT_CODE_SIM,'TIMEOUT'}}
    oc = oc_info{1}{1};  lbl = oc_info{1}{2};
    mask = trial_outcome_real == oc;
    if ~any(mask), continue; end
    fprintf('  %-8s  %-5d  %-9.3f  %-6.2f  %-9.3f  %-6.2f  %-9.3f  %-6.2f  %-8.3f\n', ...
            lbl, sum(mask), ...
            mean(mean_P_mi(mask),   'omitnan'), mean(acc_mi(mask),   'omitnan'), ...
            mean(mean_P_cvsa(mask), 'omitnan'), mean(acc_cvsa(mask), 'omitnan'), ...
            mean(mean_P_fus(mask),  'omitnan'), mean(acc_fus(mask),  'omitnan'), ...
            mean(mean_buf(mask),    'omitnan'));
end
fprintf('  Thresholds: cls1=%.2f  cls2=%.2f\n', thresholds(1), thresholds(2));

%% ── Figure 1: per-trial mean P(target class) — scatter coloured by outcome
COL_HIT  = [0.15 0.65 0.15];
COL_MISS = [0.90 0.50 0.10];
COL_TO   = [0.80 0.15 0.15];
COL_MI   = [0.85 0.30 0.10];
COL_CV   = [0.15 0.35 0.80];
COL_FUS  = [0.20 0.50 0.80];
COL_BUF  = [0.35 0.20 0.65];
mk_cls   = {'o', 's'};   % circle = class 1, square = class 2
mk_sz    = [65, 60];

stream_data   = {mean_P_mi, mean_P_cvsa, mean_P_fus, mean_buf};
stream_labels = {'P_{MI}(class asked)', 'P_{CVSA}(class asked)', ...
                 'P_{fused}(class asked)', 'Buffer(class asked)'};
stream_cols   = {COL_MI, COL_CV, COL_FUS, COL_BUF};
ylims_sp      = {[0.2, 1.0], [0.2, 1.0], [0.2, 1.0], [p_rest-0.05, 1.0]};
yref_sp       = {0.5, 0.5, 0.5, p_rest};

fig1 = figure('Name', sprintf('Mean P per trial — %s', basename), ...
              'Color', 'w', 'NumberTitle', 'off', 'Position', [40 40 1400 420]);
ax_sp = gobjects(1, 4);

for sp = 1:4
    ax = subplot(1, 4, sp);  ax_sp(sp) = ax;  hold(ax, 'on');
    dat = stream_data{sp};

    for t = 1:n_trials
        val = dat(t);  if isnan(val), continue; end
        c = trial_class(t);  if isnan(c)||c<1||c>n_cls, c=1; end
        switch trial_outcome_real(t)
            case HIT_CODE_SIM,     col = COL_HIT;
            case MISS_CODE_SIM,    col = COL_MISS;
            case TIMEOUT_CODE_SIM, col = COL_TO;
            otherwise,             col = [0.5 0.5 0.5];
        end
        scatter(ax, t, val, mk_sz(c), col, 'filled', ...
                'Marker', mk_cls{c}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end

    % Dashed mean per outcome
    for oc = {HIT_CODE_SIM, MISS_CODE_SIM, TIMEOUT_CODE_SIM}
        mask = trial_outcome_real == oc{1};
        v = dat(mask);  v = v(~isnan(v));
        if isempty(v), continue; end
        switch oc{1}
            case HIT_CODE_SIM,     lc = COL_HIT;
            case MISS_CODE_SIM,    lc = COL_MISS;
            case TIMEOUT_CODE_SIM, lc = COL_TO;
        end
        yline(ax, mean(v), '--', 'Color', lc, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end

    yline(ax, yref_sp{sp}, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    if sp == 4
        for ci = 1:n_cls
            yline(ax, thresholds(ci), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, ...
                  'Label', sprintf('thr%s=%.2f', cls_names{ci}, thresholds(ci)), ...
                  'FontSize', 7, 'HandleVisibility', 'off');
        end
    end

    set(ax, 'XLim', [0,n_trials+1], 'YLim', ylims_sp{sp}, 'XTick', 1:n_trials);
    xlabel(ax, 'trial #', 'FontSize', 9);
    title(ax, stream_labels{sp}, 'FontSize', 10, 'FontWeight', 'bold', 'Color', stream_cols{sp});
    grid(ax, 'on');
end

% Legend on first panel
ax_leg = ax_sp(1);
h = gobjects(1, 5);
h(1) = scatter(ax_leg, NaN, NaN, 55, COL_HIT,  'filled', 'o', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('HIT (%d)',  n_hit_real));
h(2) = scatter(ax_leg, NaN, NaN, 55, COL_MISS, 'filled', 'o', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('MISS (%d)', n_miss_real));
h(3) = scatter(ax_leg, NaN, NaN, 55, COL_TO,   'filled', 'o', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('TO (%d)',   n_to_real));
h(4) = scatter(ax_leg, NaN, NaN, 55, [0.5 0.5 0.5], 'filled', 'o', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('cls1 (%s)', cls_names{1}));
h(5) = scatter(ax_leg, NaN, NaN, 55, [0.5 0.5 0.5], 'filled', 's', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('cls2 (%s)', cls_names{2}));
legend(ax_leg, h, 'FontSize', 8, 'Location', 'southoutside', 'NumColumns', 5);
sgtitle(fig1, sprintf('%s  |  HIT=%d  MISS=%d  TO=%d', basename, n_hit_real, n_miss_real, n_to_real), ...
        'FontSize', 10, 'Interpreter', 'none');

%% ── Figure 2: frame accuracy per class — grouped bars (HIT vs MISS+TO)
fig2 = figure('Name', sprintf('Frame accuracy per class — %s', basename), ...
              'Color', 'w', 'NumberTitle', 'off', 'Position', [80 80 400*n_cls 400]);

stream_acc    = {acc_mi, acc_cvsa, acc_fus};
stream_xlbls  = {'MI', 'CVSA', 'fused'};
stream_fcols  = {COL_MI, COL_CV, COL_FUS};
n_streams     = numel(stream_acc);

for c = 1:n_cls
    ax = subplot(1, n_cls, c);  hold(ax, 'on');

    mask_cls  = trial_class == c;
    mask_hit  = mask_cls & (trial_outcome_real == HIT_CODE_SIM);
    mask_miss = mask_cls & (trial_outcome_real == MISS_CODE_SIM | trial_outcome_real == TIMEOUT_CODE_SIM);

    x = 1:n_streams;
    bw = 0.25;
    for s = 1:n_streams
        v_hit  = stream_acc{s}(mask_hit);
        v_miss = stream_acc{s}(mask_miss);
        v_all  = stream_acc{s}(mask_cls);
        mh  = mean(v_hit,  'omitnan');
        mm  = mean(v_miss, 'omitnan');
        mall= mean(v_all,  'omitnan');

        if ~isnan(mh)
            bar(ax, x(s)-bw, mh,   bw*0.9, 'FaceColor', COL_HIT,  'EdgeColor', 'k', 'HandleVisibility', 'off');
        end
        if ~isnan(mm)
            bar(ax, x(s),    mm,   bw*0.9, 'FaceColor', COL_MISS, 'EdgeColor', 'k', 'HandleVisibility', 'off');
        end
        if ~isnan(mall)
            bar(ax, x(s)+bw, mall, bw*0.9, 'FaceColor', stream_fcols{s}, 'EdgeColor', 'k', ...
                'FaceAlpha', 0.55, 'HandleVisibility', 'off');
        end
    end

    yline(ax, 0.5, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    set(ax, 'XTick', x, 'XTickLabel', stream_xlbls, 'YLim', [0, 1]);
    ylabel(ax, 'frame accuracy  (P > 0.5)', 'FontSize', 9);
    title(ax, sprintf('Class %s  (n=%d)', cls_names{c}, sum(mask_cls)), 'FontSize', 11);
    grid(ax, 'on');

    if c == 1
        bar(ax, NaN, NaN, 'FaceColor', COL_HIT,  'EdgeColor', 'k', 'DisplayName', sprintf('HIT (n=%d)',       sum(mask_hit)));
        bar(ax, NaN, NaN, 'FaceColor', COL_MISS, 'EdgeColor', 'k', 'DisplayName', sprintf('MISS+TO (n=%d)',   sum(mask_miss)));
        bar(ax, NaN, NaN, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'k', 'FaceAlpha', 0.55, 'DisplayName', 'overall');
        legend(ax, 'FontSize', 8, 'Location', 'south');
    end
end

sgtitle(fig2, sprintf('%s  |  frame accuracy per stream per class', basename), ...
        'FontSize', 10, 'Interpreter', 'none');
