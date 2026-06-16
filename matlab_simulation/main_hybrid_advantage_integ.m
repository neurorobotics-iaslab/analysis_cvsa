%% MAIN_HYBRID_VS_MI  Counterfactual: hybrid (fused) vs MI-only vs CVSA-only.
%
%   Loads one or more HYBRID evaluation GDFs (skips non-hybrid files) and
%   replays the pipeline exactly like main_simulate.m up to p_mi_aligned /
%   p_cvsa_aligned / art_flags / int_cfg. From that SAME pair of raw sLDA
%   streams, the leaky-WTA integrator (integrate_signal) is then run THREE
%   times with the SAME buffer/threshold parameters (int_cfg), changing only
%   which signal feeds the buffer:
%
%     - 'hybrid' : Bayesian-fused MI+CVSA (as actually run online)
%     - 'mi'     : raw MI sLDA output only
%     - 'cvsa'   : raw CVSA sLDA output only
%
%   This answers: "given the exact same classifier outputs, would the trial
%   have ended differently (HIT/MISS/TIMEOUT, and at what time) if the
%   integrator had been driven by MI alone or CVSA alone instead of the
%   fused signal?"
%
%   Per-trial simulated outcome is derived from the integrator buffer itself
%   (not from GDF events, which only reflect the actual hybrid run):
%     HIT     : buffer[target_class] >= thresholds[target_class] first
%     MISS    : buffer[other_class]  >= thresholds[other_class]  first
%     TIMEOUT : neither threshold reached within the CF window
%
%   Figures (per file, saved under per_file/):
%     1. plot_trials_streams — one panel per trial: P(c1) control signal for
%        all three streams + thresholds; background = REAL hybrid outcome.
%     2. Per-file advantage analysis (2x2): mean P(target) per trial,
%        time-to-outcome per trial, simulated-outcome heatmap (trial x
%        stream), outcome counts per stream.
%
%   Figure (once, saved at the top level): aggregate summary across ALL
%   trials/files (2x3): outcome counts, accuracy, mean time-to-HIT/MISS/
%   TIMEOUT per stream, and a Hybrid-vs-MI / Hybrid-vs-CVSA "who hit?"
%   breakdown.

clear; clc; close all;

% --- Display options -------------------------------------------------------
SHOW_FIGURES = false;   % true: figures pop up on screen; false: created hidden (export only)
if SHOW_FIGURES, fig_vis = 'on'; else, fig_vis = 'off'; end

OC_STR = {'?', 'HIT', 'MISS', 'TO'};   % outcome_code (0..3) -> label

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

% --- GUI: pick hybrid GDF(s); everything else flows from the sibling YAML(s)
default_dir = '/home/paolo/bci_vr_ws/recordings';

[gdf_names, gdf_dir] = uigetfile({'*.gdf', 'GDF recordings (*.gdf)'}, ...
                                'Select HYBRID GDF recording(s)', default_dir, 'MultiSelect', 'on');
if isequal(gdf_names, 0)
    error('main_hybrid_vs_mi:cancel', 'No GDF selected.');
end
if ischar(gdf_names), gdf_names = {gdf_names}; end
n_files = numel(gdf_names);

out_dir       = fullfile(gdf_dir, 'analysis_results', 'hybrid_vs_unimodal');
out_dir_files = fullfile(out_dir, 'per_file');
if ~exist(out_dir, 'dir'),       mkdir(out_dir);       end
if ~exist(out_dir_files, 'dir'), mkdir(out_dir_files); end

% --- Aggregate accumulators (across all files/trials) ----------------------
all_oc_hyb  = zeros(0, 2);   % [code, t_event_s]
all_oc_mi   = zeros(0, 2);
all_oc_cvsa = zeros(0, 2);
all_meanP_hyb  = [];
all_meanP_mi   = [];
all_meanP_cvsa = [];
all_meanBuf_hyb  = [];   % mean integrator output for the TARGET class over CF (per trial)
all_meanBuf_mi   = [];
all_meanBuf_cvsa = [];
all_buf_hyb  = {};       % full integrator-output trajectory for the TARGET class over CF (per trial)
all_buf_mi   = {};
all_buf_cvsa = {};
all_trial_class = [];
all_file_label   = {};

for file_idx = 1:n_files
gdf_path = fullfile(gdf_dir, gdf_names{file_idx});
fprintf('\n[%d/%d] %s\n', file_idx, n_files, gdf_names{file_idx});

%% --- Load GDF --------------------------------------------------------
[signal, header, basename] = load_gdf(gdf_path);

%% --- Load parameters YAML -----------------------------
[params, ~] = load_params_yaml(gdf_path);

paradigm = params.integrator.paradigm;
if ~strcmp(paradigm, 'hybrid')
    fprintf('  skipping %s: paradigm=%s (hybrid-only)\n', basename, paradigm);
    continue;
end

fs        = double(params.acquisition.samplerate);
framerate = double(params.acquisition.framerate);
chunk_size = round(fs / framerate);
if abs(fs - header.SampleRate) > 1e-3
    log_step('main_hybrid_vs_mi: YAML samplerate=%.1f != GDF samplerate=%.1f -> using GDF', ...
             fs, header.SampleRate);
    fs = header.SampleRate;
    chunk_size = round(fs / framerate);
end

bufsize_proc = double(params.RingBufferCfg.params.size);
bufsize_art  = double(params.RingBufferCfgArtifact.params.size);
eog_names    = to_strcell(params.CarCfg.params.EOG_ch_names);

do_car_mi   = logical(params.processing_fbcsp_mi.do_car);
do_car_cvsa = logical(params.processing_fbcsp_cvsa.do_car);
nchannels   = double(params.processing_fbcsp_mi.nchannels);

signal = signal(:,1:nchannels); % remove the last 3 values: associated with ACC values

log_step('main_hybrid_vs_mi: paradigm=%s, fs=%g, framerate=%g, chunk=%d, bufproc=%d, bufart=%d, nchannels=%d', ...
         paradigm, fs, framerate, chunk_size, bufsize_proc, bufsize_art, nchannels);

%% --- Load CSP + sLDA for both MI and CVSA -----------------------------
csp_mi   = load_csp(params, 'mi');
slda_mi  = load_slda(params, 'mi');
csp_cvsa = load_csp(params, 'cvsa');
slda_cvsa = load_slda(params, 'cvsa');

%% --- Apply processing (one stream per modality) -----------------------
proc_cfg_mi = struct('samplerate', fs, 'chunk_size', chunk_size, ...
                     'bufsize', bufsize_proc, 'filter_order', 4, ...
                     'do_car', do_car_mi, 'eog_names', {eog_names});
[features_mi, header_mi, info_proc] = apply_processing(signal, header, csp_mi, proc_cfg_mi);

proc_cfg_cvsa = struct('samplerate', fs, 'chunk_size', chunk_size, ...
                       'bufsize', bufsize_proc, 'filter_order', 4, ...
                       'do_car', do_car_cvsa, 'eog_names', {eog_names});
[features_cv, header_cv, ~] = apply_processing(signal, header, csp_cvsa, proc_cfg_cvsa);

%% --- Apply artifact detection on raw signal ----------------------------
art_cfg = params.ArtifactCfg.params;
art_cfg.EOG_ch_names = to_strcell(art_cfg.EOG_ch_names);
cfg_art = struct('samplerate', fs, 'chunk_size', chunk_size, 'bufsize_artifact', bufsize_art);
[art_flags, ~] = detect_artifacts(signal, header, art_cfg, cfg_art);

%% --- Apply sLDA over the whole feature stream --------------------------
p_mi_aligned   = apply_slda(features_mi, slda_mi,   csp_mi.bands);
p_cvsa_aligned = apply_slda(features_cv, slda_cvsa, csp_cvsa.bands);

log_step('main_hybrid_vs_mi: streams aligned (n_chunks=%d, first_valid_proc=%d)', ...
         info_proc.n_chunks, info_proc.first_valid_chunk);

%% --- Integrator config (shared by all three counterfactual runs) ------
int_cfg = params.integrator;
if ~isfield(int_cfg, 'increment'),               int_cfg.increment = 1; end
if ~isfield(int_cfg, 'thresholds_rejection'),    int_cfg.thresholds_rejection = []; end
if ~isfield(int_cfg, 'cvsa_influence'),          int_cfg.cvsa_influence = 2.5; end
if ~isfield(int_cfg, 'thresholds'),              int_cfg.thresholds = params.training_node.thresholds; end

header_chunks = header_mi;
header_chunks.framerate = framerate;

classes    = to_vec(int_cfg.classes);
n_cls      = numel(classes);
thresholds = to_vec(int_cfg.thresholds);

%% --- Run the integrator three times: hybrid / MI-only / CVSA-only ------
trials_hyb  = integrate_signal(p_mi_aligned, p_cvsa_aligned, art_flags, header_chunks, int_cfg, 'hybrid');
trials_mi   = integrate_signal(p_mi_aligned, p_cvsa_aligned, art_flags, header_chunks, int_cfg, 'mi');
trials_cvsa = integrate_signal(p_mi_aligned, p_cvsa_aligned, art_flags, header_chunks, int_cfg, 'cvsa');
n_trials = numel(trials_hyb);

%% --- Read REAL outcomes from GDF events (897=HIT, 898=MISS, 899=TIMEOUT)
HIT_CODE_SIM = 897;  MISS_CODE_SIM = 898;  TIMEOUT_CODE_SIM = 899;  CF_CODE_SIM = 781;
POS_ev = header.EVENT.POS;  TYP_ev = header.EVENT.TYP;
cf_samp = POS_ev(TYP_ev == CF_CODE_SIM);
trial_outcome_real = zeros(1, n_trials);
for t = 1:min(n_trials, numel(cf_samp))
    after = POS_ev > cf_samp(t);
    idx   = find((TYP_ev==HIT_CODE_SIM | TYP_ev==MISS_CODE_SIM | TYP_ev==TIMEOUT_CODE_SIM) & after, 1);
    if ~isempty(idx), trial_outcome_real(t) = TYP_ev(idx); end
end

%% --- Per-trial: simulated outcome (code + time) and mean P(target) -----
oc = struct('hyb', zeros(n_trials,2), 'mi', zeros(n_trials,2), 'cvsa', zeros(n_trials,2));
meanP   = struct('hyb', nan(n_trials,1), 'mi', nan(n_trials,1), 'cvsa', nan(n_trials,1));
meanBuf = struct('hyb', nan(n_trials,1), 'mi', nan(n_trials,1), 'cvsa', nan(n_trials,1));
buf_hyb_cell  = cell(n_trials, 1);   % integrator output for the target class over CF, per trial
buf_mi_cell   = cell(n_trials, 1);
buf_cvsa_cell = cell(n_trials, 1);
trial_class = nan(n_trials, 1);

fprintf('\n  Trial  class  REAL    HYB            MI-only        CVSA-only\n');
for t = 1:n_trials
    c = trials_hyb(t).target_class;
    if isnan(c) || c < 1 || c > n_cls
        continue;
    end
    trial_class(t) = c;

    [oc.hyb(t,1),  oc.hyb(t,2)]  = classify_trial_outcome(trials_hyb(t),  thresholds, framerate);
    [oc.mi(t,1),   oc.mi(t,2)]   = classify_trial_outcome(trials_mi(t),   thresholds, framerate);
    [oc.cvsa(t,1), oc.cvsa(t,2)] = classify_trial_outcome(trials_cvsa(t), thresholds, framerate);

    np = trials_hyb(t).n_pre;  nc = trials_hyb(t).n_cf;
    meanP.hyb(t)  = mean(trials_hyb(t).raw (np+1:np+nc, c), 'omitnan');
    meanP.mi(t)   = mean(trials_mi(t).raw  (np+1:np+nc, c), 'omitnan');
    meanP.cvsa(t) = mean(trials_cvsa(t).raw(np+1:np+nc, c), 'omitnan');

    % Control-signal (integrator output) for the CUED class -- this is what
    % actually drives the VR feedback and what crosses thresholds(c).
    buf_hyb_cell{t}  = trials_hyb(t).integrated (np+1:np+nc, c);
    buf_mi_cell{t}   = trials_mi(t).integrated  (np+1:np+nc, c);
    buf_cvsa_cell{t} = trials_cvsa(t).integrated(np+1:np+nc, c);
    meanBuf.hyb(t)  = mean(buf_hyb_cell{t},  'omitnan');
    meanBuf.mi(t)   = mean(buf_mi_cell{t},   'omitnan');
    meanBuf.cvsa(t) = mean(buf_cvsa_cell{t}, 'omitnan');

    switch trial_outcome_real(t)
        case HIT_CODE_SIM,     real_str = 'HIT';
        case MISS_CODE_SIM,    real_str = 'MISS';
        case TIMEOUT_CODE_SIM, real_str = 'TO';
        otherwise,             real_str = '?';
    end
    fprintf('  %-5d  %-5d  %-6s  %-4s (%4.1fs)  %-4s (%4.1fs)  %-4s (%4.1fs)\n', ...
            t, c, real_str, ...
            OC_STR{oc.hyb(t,1)+1},  oc.hyb(t,2), ...
            OC_STR{oc.mi(t,1)+1},   oc.mi(t,2), ...
            OC_STR{oc.cvsa(t,1)+1}, oc.cvsa(t,2));
end

%% --- Figure 1: per-trial control-signal comparison ---------------------
fig1 = plot_trials_streams(trials_hyb, trials_mi, trials_cvsa, int_cfg, framerate, ...
                            basename, trial_outcome_real, oc, SHOW_FIGURES);
if ~isempty(fig1)
    saveas(fig1, fullfile(out_dir_files, sprintf('signals_%s.svg', basename)), 'svg');
    if ~SHOW_FIGURES, close(fig1); end
end

%% --- Figure 2: per-file advantage analysis (2x3) ------------------------
COL_HYB  = [0.05 0.35 0.85];
COL_MI   = [0.85 0.40 0.10];
COL_CV   = [0.49 0.18 0.56];
COL_HIT  = [0.15 0.65 0.15];
COL_MISS = [0.80 0.15 0.15];   % MISS (898): wrong class threshold reached -> red
COL_TO   = [0.95 0.65 0.10];   % TIMEOUT (899): no threshold reached -> orange/yellow
COL_UNK  = [0.65 0.65 0.65];
OC_COLS  = [COL_UNK; COL_HIT; COL_MISS; COL_TO];   % index = code + 1

valid_t = find(~isnan(trial_class));
nv      = numel(valid_t);
M = [oc.hyb(valid_t,1), oc.mi(valid_t,1), oc.cvsa(valid_t,1)];
x_ticks = 1:nv;

fig2 = figure('Name', sprintf('Hybrid vs MI/CVSA-only — per-trial analysis — %s', basename), ...
              'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Panel (1,1): mean raw sLDA P(target class) per trial — classifier signal quality
ax1 = subplot(2,3,1); hold(ax1, 'on');
bar(ax1, [meanP.hyb(valid_t), meanP.mi(valid_t), meanP.cvsa(valid_t)]);
yline(ax1, 0.5, 'k:', 'HandleVisibility', 'off');
set(ax1, 'XTick', x_ticks, 'XTickLabel', valid_t, 'YLim', [0,1]);
xlabel(ax1, 'trial #'); ylabel(ax1, 'mean P(target) over CF');
legend(ax1, {'Hybrid','MI-only','CVSA-only'}, 'Location', 'best', 'FontSize', 8);
title(ax1, 'Mean sLDA P(target class) over CF', 'FontWeight', 'bold');
grid(ax1, 'on');

% Panel (1,2): mean integrator (control signal) for the TARGET class per trial
%   This is the leaky-buffer value for the cued class — what drives the VR
%   object toward the target threshold.  Higher = buffer moved more toward
%   the win threshold for the class the subject was trying to hit.
ax2 = subplot(2,3,2); hold(ax2, 'on');
bar(ax2, [meanBuf.hyb(valid_t), meanBuf.mi(valid_t), meanBuf.cvsa(valid_t)]);
yline(ax2, thresholds(1) - 5e-3, 'k--', 'thr', 'HandleVisibility', 'off');
yline(ax2, to_vec(int_cfg.init_val(1)), 'k:', 'HandleVisibility', 'off');
set(ax2, 'XTick', x_ticks, 'XTickLabel', valid_t, 'YLim', [0,1]);
xlabel(ax2, 'trial #'); ylabel(ax2, 'mean buffer(target class)');
legend(ax2, {'Hybrid','MI-only','CVSA-only'}, 'Location', 'best', 'FontSize', 8);
title(ax2, 'Mean integrator output for TARGET class (control signal level)', 'FontWeight', 'bold');
grid(ax2, 'on');

% Panel (1,3): time to outcome per trial, colored by outcome
ax3 = subplot(2,3,3); hold(ax3, 'on');
bar(ax3, [oc.hyb(valid_t,2), oc.mi(valid_t,2), oc.cvsa(valid_t,2)]);
set(ax3, 'XTick', x_ticks, 'XTickLabel', valid_t);
xlabel(ax3, 'trial #'); ylabel(ax3, 'time to outcome (s)');
legend(ax3, {'Hybrid','MI-only','CVSA-only'}, 'Location', 'best', 'FontSize', 8);
title(ax3, 'Time to outcome (HIT / MISS / TIMEOUT)', 'FontWeight', 'bold');
grid(ax3, 'on');

% Panel (2,1): simulated-outcome heatmap (trial x stream)
ax4 = subplot(2,3,4);
imagesc(ax4, M, [-0.5, 3.5]);
colormap(ax4, OC_COLS);
set(ax4, 'XTick', 1:3, 'XTickLabel', {'Hybrid','MI-only','CVSA-only'}, ...
         'YTick', 1:nv, 'YTickLabel', valid_t, 'YDir', 'normal');
xlabel(ax4, 'stream'); ylabel(ax4, 'trial #');
title(ax4, 'Simulated outcome per stream', 'FontWeight', 'bold');
for ti = 1:nv
    for ss = 1:3
        text(ax4, ss, ti, OC_STR{M(ti,ss)+1}, 'HorizontalAlignment', 'center', 'FontSize', 7);
    end
end

% Panel (2,2): outcome counts stacked bar per stream
ax5 = subplot(2,3,5); hold(ax5, 'on');
counts = zeros(3,3);
for ss = 1:3
    counts(ss,1) = sum(M(:,ss)==1);
    counts(ss,2) = sum(M(:,ss)==2);
    counts(ss,3) = sum(M(:,ss)==3);
end
b2 = bar(ax5, counts, 'stacked');
b2(1).FaceColor = COL_HIT; b2(2).FaceColor = COL_MISS; b2(3).FaceColor = COL_TO;
set(ax5, 'XTick', 1:3, 'XTickLabel', {'Hybrid','MI-only','CVSA-only'});
ylabel(ax5, '# trials');
legend(ax5, {'HIT','MISS','TIMEOUT'}, 'Location', 'best', 'FontSize', 8);
title(ax5, 'Outcome counts per stream', 'FontWeight', 'bold');
grid(ax5, 'on');

% Panel (2,3): per-trial control-signal advantage (Hybrid - MI-only) and
%   (Hybrid - CVSA-only) for the cued class.  Positive = hybrid integrator
%   stayed closer to the target threshold for that trial's cued class.
ax6 = subplot(2,3,6); hold(ax6, 'on');
d_mi  = meanBuf.hyb(valid_t) - meanBuf.mi(valid_t);
d_cvs = meanBuf.hyb(valid_t) - meanBuf.cvsa(valid_t);
bar_x = (1:nv)' + [-0.2, 0.2];
bar(ax6, bar_x(:,1), d_mi,  0.35, 'FaceColor', COL_MI, 'EdgeColor', 'none');
bar(ax6, bar_x(:,2), d_cvs, 0.35, 'FaceColor', COL_CV, 'EdgeColor', 'none');
yline(ax6, 0, 'k-', 'HandleVisibility', 'off');
set(ax6, 'XTick', 1:nv, 'XTickLabel', valid_t);
xlabel(ax6, 'trial #');
ylabel(ax6, '\Delta mean buffer(target): Hybrid - other');
legend(ax6, {'Hybrid - MI-only','Hybrid - CVSA-only'}, 'Location', 'best', 'FontSize', 8);
title(ax6, 'Per-trial control-signal advantage (positive = fusion helps target class)', 'FontWeight', 'bold');
grid(ax6, 'on');

sgtitle(fig2, sprintf('%s — Hybrid vs MI-only vs CVSA-only (counterfactual, same buffer/thresholds)', ...
        basename), 'Interpreter', 'none');

saveas(fig2, fullfile(out_dir_files, sprintf('advantage_%s.svg', basename)), 'svg');
if ~SHOW_FIGURES, close(fig2); end

%% --- Accumulate for the aggregate figure --------------------------------
all_oc_hyb      = [all_oc_hyb;  oc.hyb(valid_t,:)];
all_oc_mi       = [all_oc_mi;   oc.mi(valid_t,:)];
all_oc_cvsa     = [all_oc_cvsa; oc.cvsa(valid_t,:)];
all_meanP_hyb   = [all_meanP_hyb;  meanP.hyb(valid_t)];
all_meanP_mi    = [all_meanP_mi;   meanP.mi(valid_t)];
all_meanP_cvsa  = [all_meanP_cvsa; meanP.cvsa(valid_t)];
all_meanBuf_hyb  = [all_meanBuf_hyb;  meanBuf.hyb(valid_t)];
all_meanBuf_mi   = [all_meanBuf_mi;   meanBuf.mi(valid_t)];
all_meanBuf_cvsa = [all_meanBuf_cvsa; meanBuf.cvsa(valid_t)];
all_buf_hyb  = [all_buf_hyb;  buf_hyb_cell(valid_t)];
all_buf_mi   = [all_buf_mi;   buf_mi_cell(valid_t)];
all_buf_cvsa = [all_buf_cvsa; buf_cvsa_cell(valid_t)];
all_trial_class = [all_trial_class; trial_class(valid_t)];
all_file_label  = [all_file_label; repmat({basename}, numel(valid_t), 1)];

end % file_idx loop

%% ═══════════════════════════════════════════════════════════════════════
%  AGGREGATE SUMMARY ACROSS ALL FILES/TRIALS
%  ═══════════════════════════════════════════════════════════════════════
n_tot = size(all_oc_hyb, 1);
if n_tot == 0
    fprintf('\nNo hybrid trials found across the selected file(s) — nothing to summarize.\n');
    return;
end

COL_HYB  = [0.05 0.35 0.85];
COL_MI   = [0.85 0.40 0.10];
COL_CV   = [0.49 0.18 0.56];
COL_HIT  = [0.15 0.65 0.15];
COL_MISS = [0.80 0.15 0.15];   % MISS (898): wrong class threshold reached -> red
COL_TO   = [0.95 0.65 0.10];   % TIMEOUT (899): no threshold reached -> orange/yellow
streams_oc   = {all_oc_hyb, all_oc_mi, all_oc_cvsa};
stream_names = {'Hybrid','MI-only','CVSA-only'};
stream_cols  = {COL_HYB, COL_MI, COL_CV};

% --- counts + accuracy per stream -----------------------------------------
n_hit  = zeros(1,3); n_miss = zeros(1,3); n_to = zeros(1,3); acc = zeros(1,3); acc_se = zeros(1,3);
t_hit_mean  = nan(1,3); t_hit_sem  = nan(1,3);
t_miss_mean = nan(1,3); t_miss_sem = nan(1,3);
t_to_mean   = nan(1,3); t_to_sem   = nan(1,3);

fprintf('\n══════════════════ Aggregate summary (%d trials, %d file(s)) ══════════════════\n', n_tot, n_files);
fprintf('  %-10s %5s %5s %5s  %6s   %8s   %8s   %8s\n', ...
        'stream', 'HIT', 'MISS', 'TO', 'acc%', 'TTH(s)', 'Tmiss(s)', 'Tto(s)');
for s = 1:3
    codes = streams_oc{s}(:,1);
    times = streams_oc{s}(:,2);
    n_hit(s)  = sum(codes==1);
    n_miss(s) = sum(codes==2);
    n_to(s)   = sum(codes==3);
    acc(s)    = n_hit(s) / n_tot;
    acc_se(s) = sqrt(acc(s) * (1-acc(s)) / n_tot);

    v = times(codes==1); if ~isempty(v), t_hit_mean(s)  = mean(v); t_hit_sem(s)  = std(v)/sqrt(numel(v)); end
    v = times(codes==2); if ~isempty(v), t_miss_mean(s) = mean(v); t_miss_sem(s) = std(v)/sqrt(numel(v)); end
    v = times(codes==3); if ~isempty(v), t_to_mean(s)   = mean(v); t_to_sem(s)   = std(v)/sqrt(numel(v)); end

    fprintf('  %-10s %5d %5d %5d  %5.1f%%   %8s   %8s   %8s\n', ...
            stream_names{s}, n_hit(s), n_miss(s), n_to(s), 100*acc(s), ...
            fmt_mean_sem(t_hit_mean(s),  t_hit_sem(s)), ...
            fmt_mean_sem(t_miss_mean(s), t_miss_sem(s)), ...
            fmt_mean_sem(t_to_mean(s),   t_to_sem(s)));
end

% --- "who hits?" breakdown: Hybrid vs MI-only, Hybrid vs CVSA-only ---------
cmp_labels = {'both HIT','only Hybrid HIT','only other HIT','neither HIT'};
hyb_hit = all_oc_hyb(:,1) == 1;
mi_hit  = all_oc_mi(:,1)  == 1;
cv_hit  = all_oc_cvsa(:,1) == 1;

cmp_mi  = [sum(hyb_hit & mi_hit), sum(hyb_hit & ~mi_hit), sum(~hyb_hit & mi_hit), sum(~hyb_hit & ~mi_hit)];
cmp_cv  = [sum(hyb_hit & cv_hit), sum(hyb_hit & ~cv_hit), sum(~hyb_hit & cv_hit), sum(~hyb_hit & ~cv_hit)];

fprintf('\n  Hybrid vs MI-only    : %s\n', fmt_cmp(cmp_labels, cmp_mi));
fprintf('  Hybrid vs CVSA-only  : %s\n', fmt_cmp(cmp_labels, cmp_cv));
fprintf('  Hybrid net advantage vs MI-only   : %+d trials (rescued - cost)\n', cmp_mi(2) - cmp_mi(3));
fprintf('  Hybrid net advantage vs CVSA-only : %+d trials (rescued - cost)\n', cmp_cv(2) - cmp_cv(3));

% --- per-class breakdown ----------------------------------------------------
classes_all = unique(all_trial_class(~isnan(all_trial_class)))';
fprintf('\n  Per-class accuracy (HIT rate):\n');
fprintf('  %-8s %5s %8s %8s %8s\n', 'class', 'n', 'Hybrid%', 'MI-only%', 'CVSA-only%');
for c = classes_all
    m = all_trial_class == c;
    fprintf('  %-8d %5d %8.1f %8.1f %8.1f\n', c, sum(m), ...
            100*mean(all_oc_hyb(m,1)==1), 100*mean(all_oc_mi(m,1)==1), 100*mean(all_oc_cvsa(m,1)==1));
end
fprintf('═══════════════════════════════════════════════════════════════════════════════\n');

% --- Paired control-signal advantage: sign-flip permutation test -----------
N_PERM = 2000;
d_mi_all  = all_meanBuf_hyb - all_meanBuf_mi;
d_cvs_all = all_meanBuf_hyb - all_meanBuf_cvsa;
p_buf_mi  = sign_flip_test(d_mi_all,  N_PERM);
p_buf_cvs = sign_flip_test(d_cvs_all, N_PERM);
fprintf('\n  Mean integrator output (target class) — paired sign-flip permutation test:\n');
fprintf('    Hybrid vs MI-only   : delta=%.4f  p=%.4f  %s\n', mean(d_mi_all,  'omitnan'), p_buf_mi,  p_to_stars(p_buf_mi));
fprintf('    Hybrid vs CVSA-only : delta=%.4f  p=%.4f  %s\n', mean(d_cvs_all, 'omitnan'), p_buf_cvs, p_to_stars(p_buf_cvs));
fprintf('    (two-sided, %d permutations, single comparison)\n', N_PERM);

%% --- Aggregate figure (3x3) ----------------------------------------------
fig3 = figure('Name', 'Hybrid vs MI-only vs CVSA-only — aggregate summary', ...
              'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig3, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

acc_delta_mi  = 100*(acc(1)-acc(2));
acc_delta_cvs = 100*(acc(1)-acc(3));

% Panel (1,1): outcome counts stacked bar
ax1 = subplot(3,3,1); hold(ax1, 'on');
counts_agg = [n_hit; n_miss; n_to]';
b_agg = bar(ax1, counts_agg, 'stacked');
b_agg(1).FaceColor = COL_HIT; b_agg(2).FaceColor = COL_MISS; b_agg(3).FaceColor = COL_TO;
for s = 1:3
    text(ax1, s, n_tot*1.05, sprintf('H%d M%d T%d', n_hit(s), n_miss(s), n_to(s)), ...
         'HorizontalAlignment', 'center', 'FontSize', 8);
end
set(ax1, 'XTick', 1:3, 'XTickLabel', stream_names, 'YLim', [0, n_tot*1.18]);
ylabel(ax1, '# trials'); legend(ax1, {'HIT','MISS','TIMEOUT'}, 'Location', 'best', 'FontSize', 8);
title(ax1, sprintf('Outcome counts (n=%d)', n_tot), 'FontWeight', 'bold');
grid(ax1, 'on');

% Panel (1,2): accuracy — title states the finding directly
ax2 = subplot(3,3,2); hold(ax2, 'on');
for s = 1:3
    bar(ax2, s, 100*acc(s), 0.6, 'FaceColor', stream_cols{s}, 'EdgeColor', 'k');
    errorbar(ax2, s, 100*acc(s), 100*acc_se(s), 'k', 'LineWidth', 1.2, 'CapSize', 6);
    text(ax2, s, 100*acc(s)+4, sprintf('%.0f%%', 100*acc(s)), 'HorizontalAlignment', 'center', 'FontSize', 9);
end
yline(ax2, 50, 'k--', 'HandleVisibility', 'off');
set(ax2, 'XTick', 1:3, 'XTickLabel', stream_names, 'YLim', [0,108]);
ylabel(ax2, 'HIT rate (%)');
title(ax2, sprintf('Accuracy — Hyb %.0f%%,  MI %.0f%% (%+.0f%%),  CVSA %.0f%% (%+.0f%%)', ...
     100*acc(1), 100*acc(2), acc_delta_mi, 100*acc(3), acc_delta_cvs), 'FontWeight', 'bold');
grid(ax2, 'on');

% Panel (1,3): mean time-to-HIT — title states speed advantage
ax3 = subplot(3,3,3); hold(ax3, 'on');
plot_mean_sem_bars(ax3, t_hit_mean, t_hit_sem, n_hit, stream_names, stream_cols);
if ~isnan(t_hit_mean(1)) && ~isnan(t_hit_mean(2))
    title(ax3, sprintf('Mean TTH — Hybrid faster by %.2fs vs MI', t_hit_mean(2)-t_hit_mean(1)), ...
         'FontWeight', 'bold');
else
    title(ax3, 'Mean time to HIT', 'FontWeight', 'bold');
end
ylabel(ax3, 'time (s)'); grid(ax3, 'on');

% Panel (2,1): mean time-to-MISS
ax4 = subplot(3,3,4); hold(ax4, 'on');
plot_mean_sem_bars(ax4, t_miss_mean, t_miss_sem, n_miss, stream_names, stream_cols);
title(ax4, 'Mean time to MISS (898 — wrong class threshold)', 'FontWeight', 'bold');
ylabel(ax4, 'time (s)'); grid(ax4, 'on');

% Panel (2,2): mean time-to-TIMEOUT (restored)
ax_to = subplot(3,3,5); hold(ax_to, 'on');
plot_mean_sem_bars(ax_to, t_to_mean, t_to_sem, n_to, stream_names, stream_cols);
title(ax_to, 'Mean time to TIMEOUT (899 — no threshold reached)', 'FontWeight', 'bold');
ylabel(ax_to, 'time (s)'); grid(ax_to, 'on');

% Panel (2,3): "who hits?" — Hybrid vs MI-only
%   "Hyb only" = fusion RESCUED the trial (MI-only would have missed/timed out)
%   "MI only"  = fusion COST the trial (MI-only would have hit, fusion didn't)
cmp_labels_short = {'Both HIT', 'Hyb only', 'MI only', 'Neither'};
ax5 = subplot(3,3,6); hold(ax5, 'on');
b_mi2 = bar(ax5, cmp_mi(:), 0.6);
b_mi2.FaceColor = 'flat';
b_mi2.CData = [0.15 0.65 0.15; COL_HYB; COL_MI; 0.65 0.65 0.65];
for i = 1:4
    text(ax5, i, cmp_mi(i)+0.3, num2str(cmp_mi(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end
set(ax5, 'XTick', 1:4, 'XTickLabel', cmp_labels_short);
ylabel(ax5, '# trials');
net_mi = cmp_mi(2) - cmp_mi(3);
title(ax5, sprintf('Hybrid vs MI-only — rescued=%d, cost=%d, net=%+d trials', ...
      cmp_mi(2), cmp_mi(3), net_mi), 'FontWeight', 'bold');
grid(ax5, 'on');

% Panel (3,1): "who hits?" — Hybrid vs CVSA-only
cmp_labels_cvs = {'Both HIT', 'Hyb only', 'CVSA only', 'Neither'};
ax6 = subplot(3,3,7); hold(ax6, 'on');
b_cv2 = bar(ax6, cmp_cv(:), 0.6);
b_cv2.FaceColor = 'flat';
b_cv2.CData = [0.15 0.65 0.15; COL_HYB; COL_CV; 0.65 0.65 0.65];
for i = 1:4
    text(ax6, i, cmp_cv(i)+0.3, num2str(cmp_cv(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end
set(ax6, 'XTick', 1:4, 'XTickLabel', cmp_labels_cvs);
ylabel(ax6, '# trials');
net_cvs = cmp_cv(2) - cmp_cv(3);
title(ax6, sprintf('Hybrid vs CVSA-only — rescued=%d, cost=%d, net=%+d trials', ...
      cmp_cv(2), cmp_cv(3), net_cvs), 'FontWeight', 'bold');
grid(ax6, 'on');

% Panel (3,2): mean integrator output for TARGET class per stream + sign test
ax7 = subplot(3,3,8); hold(ax7, 'on');
mb_data = {all_meanBuf_hyb, all_meanBuf_mi, all_meanBuf_cvsa};
mb_mean = cellfun(@(x) mean(x, 'omitnan'), mb_data);
mb_sem  = cellfun(@(x) std(x,  'omitnan') / sqrt(sum(~isnan(x))), mb_data);
for s = 1:3
    bar(ax7, s, mb_mean(s), 0.6, 'FaceColor', stream_cols{s}, 'EdgeColor', 'k');
    errorbar(ax7, s, mb_mean(s), mb_sem(s), 'k', 'LineWidth', 1.2, 'CapSize', 6);
end
yline(ax7, thresholds(1) - 5e-3, 'k--', 'win thr', 'HandleVisibility', 'off');
yline(ax7, 0.5, 'k:', 'p\_rest', 'HandleVisibility', 'off');
y_brk = max(mb_mean + mb_sem) + 0.04;
draw_sig_bracket(ax7, 1, 2, y_brk,        p_to_stars(p_buf_mi));
draw_sig_bracket(ax7, 1, 3, y_brk + 0.07, p_to_stars(p_buf_cvs));
set(ax7, 'XTick', 1:3, 'XTickLabel', stream_names, 'YLim', [0.45, y_brk + 0.16]);
ylabel(ax7, 'mean buffer(target class)');
title(ax7, sprintf('Mean control signal (cued class)\nHyb-MI: %s (p=%.3f), Hyb-CVSA: %s (p=%.3f)', ...
     p_to_stars(p_buf_mi), p_buf_mi, p_to_stars(p_buf_cvs), p_buf_cvs), 'FontWeight', 'bold');
grid(ax7, 'on');

% Panel (3,3): per-trial advantage distribution (histogram of Hybrid - unimodal)
ax8 = subplot(3,3,9); hold(ax8, 'on');
histogram(ax8, d_mi_all,  'BinWidth', 0.02, ...
          'FaceColor', COL_MI, 'FaceAlpha', 0.55, 'EdgeColor', 'none');
histogram(ax8, d_cvs_all, 'BinWidth', 0.02, ...
          'FaceColor', COL_CV, 'FaceAlpha', 0.55, 'EdgeColor', 'none');
xline(ax8, 0, 'k--', 'HandleVisibility', 'off');
mn_d_mi_tot  = mean(d_mi_all,  'omitnan');
mn_d_cvs_tot = mean(d_cvs_all, 'omitnan');
xline(ax8, mn_d_mi_tot,  '-', 'Color', COL_MI,  'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(ax8, mn_d_cvs_tot, '-', 'Color', COL_CV,  'LineWidth', 1.5, 'HandleVisibility', 'off');
n_above_mi  = sum(d_mi_all  > 0);
n_above_cvs = sum(d_cvs_all > 0);
legend(ax8, {sprintf('Hyb-MI  mean=%+.3f  %s  (%d/%d>0)', mn_d_mi_tot,  p_to_stars(p_buf_mi),  n_above_mi,  n_tot), ...
             sprintf('Hyb-CVSA mean=%+.3f %s  (%d/%d>0)', mn_d_cvs_tot, p_to_stars(p_buf_cvs), n_above_cvs, n_tot)}, ...
       'Location', 'best', 'FontSize', 8);
xlabel(ax8, '\Delta mean buffer(target)  [Hybrid - unimodal]');
ylabel(ax8, '# trials');
title(ax8, sprintf('Per-trial control-signal advantage\n(positive = Hybrid pushed buffer higher for cued class)'), ...
     'FontWeight', 'bold');
grid(ax8, 'on');

sgtitle(fig3, sprintf('Hybrid vs MI-only vs CVSA-only — %d trials, %d file(s) | Accuracy: Hyb %+.0f%% vs MI, Hyb %+.0f%% vs CVSA', ...
        n_tot, n_files, acc_delta_mi, acc_delta_cvs), 'Interpreter', 'none');

saveas(fig3, fullfile(out_dir, 'summary_hybrid_vs_unimodal.svg'), 'svg');
if ~SHOW_FIGURES, close(fig3); end

%% --- Figure 4: time-resolved significance (buffer trajectory) ------------
%
%   For each time step k from CF onset, aligns all trials and computes
%   mean +- SEM of buffer(target_class) per stream, plus a paired sign-flip
%   permutation test (Hybrid vs MI-only, Hybrid vs CVSA-only) using only
%   the trials still running at that step.
%
%   IMPORTANT: these are POINTWISE tests with NO correction for multiple
%   comparisons.  With n~40 trials and ~100 time bins, isolated p<0.05
%   frames are likely false positives.  Look for SUSTAINED runs of low
%   p-values as cluster-level evidence.  This is an exploratory analysis.

max_nc   = max(cellfun(@numel, all_buf_hyb));
B_hyb    = nan(n_tot, max_nc);
B_mi_mat = nan(n_tot, max_nc);
B_cv_mat = nan(n_tot, max_nc);
for i = 1:n_tot
    nc = numel(all_buf_hyb{i});
    B_hyb(i,   1:nc) = all_buf_hyb{i};
    B_mi_mat(i,1:nc) = all_buf_mi{i};
    B_cv_mat(i,1:nc) = all_buf_cvsa{i};
end
t_ax    = (0:max_nc-1) / framerate;
n_valid = sum(~isnan(B_hyb), 1);
MIN_N   = max(6, round(n_tot * 0.15));

mn_hyb = mean(B_hyb,    1, 'omitnan');
mn_mi  = mean(B_mi_mat, 1, 'omitnan');
mn_cv  = mean(B_cv_mat, 1, 'omitnan');
se_hyb = std(B_hyb,    [], 1, 'omitnan') ./ sqrt(max(n_valid, 1));
se_mi  = std(B_mi_mat, [], 1, 'omitnan') ./ sqrt(max(n_valid, 1));
se_cv  = std(B_cv_mat, [], 1, 'omitnan') ./ sqrt(max(n_valid, 1));

N_PERM_T = 500;
p_t_mi = nan(1, max_nc);
p_t_cv = nan(1, max_nc);
fprintf('\n  Computing temporal significance (%d time bins x %d perms)...\n', max_nc, N_PERM_T);
for k = 1:max_nc
    d_mi_k = B_hyb(:,k) - B_mi_mat(:,k);
    d_cv_k = B_hyb(:,k) - B_cv_mat(:,k);
    ok_mi  = ~isnan(d_mi_k);  ok_cv = ~isnan(d_cv_k);
    if sum(ok_mi) >= MIN_N, p_t_mi(k) = sign_flip_test(d_mi_k(ok_mi), N_PERM_T); end
    if sum(ok_cv) >= MIN_N, p_t_cv(k) = sign_flip_test(d_cv_k(ok_cv), N_PERM_T); end
end
fprintf('  p<0.05 frames (uncorrected): Hybrid vs MI=%d/%d, Hybrid vs CVSA=%d/%d\n', ...
        sum(p_t_mi < 0.05, 'omitnan'), sum(~isnan(p_t_mi)), ...
        sum(p_t_cv < 0.05, 'omitnan'), sum(~isnan(p_t_cv)));
fprintf('  NOTE: pointwise only, no cluster correction (n=%d trials, exploratory)\n', n_tot);

fig4 = figure('Name', 'Temporal significance — buffer(target class)', ...
              'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig4, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Subplot 1: mean +- SEM buffer(target class) per stream vs time
ax_t1 = subplot(3,1,1); hold(ax_t1, 'on');
fill(ax_t1, [t_ax, fliplr(t_ax)], [mn_hyb+se_hyb, fliplr(mn_hyb-se_hyb)], COL_HYB, ...
     'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(ax_t1, [t_ax, fliplr(t_ax)], [mn_mi+se_mi,   fliplr(mn_mi-se_mi)],   COL_MI, ...
     'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(ax_t1, [t_ax, fliplr(t_ax)], [mn_cv+se_cv,   fliplr(mn_cv-se_cv)],   COL_CV, ...
     'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax_t1, t_ax, mn_hyb, '-',  'Color', COL_HYB, 'LineWidth', 2.5);
plot(ax_t1, t_ax, mn_mi,  '--', 'Color', COL_MI,  'LineWidth', 1.8);
plot(ax_t1, t_ax, mn_cv,  ':',  'Color', COL_CV,  'LineWidth', 1.8);
yline(ax_t1, thresholds(1) - 5e-3, 'k--', 'win thr', 'HandleVisibility', 'off');
yline(ax_t1, 0.5, 'k:', 'HandleVisibility', 'off');
if isfield(int_cfg, 'cvsa_influence')
    xline(ax_t1, int_cfg.cvsa_influence, '--', 'Color', [0.55 0.55 0.55], ...
          'HandleVisibility', 'off');
    text(ax_t1, int_cfg.cvsa_influence + 0.05, 0.52, '\alpha\rightarrow0', ...
         'FontSize', 8, 'Color', [0.55 0.55 0.55]);
end
% grey patches where fewer than MIN_N trials are still running
lo_idx = find(n_valid < MIN_N & n_valid > 0);
for k = lo_idx
    patch(ax_t1, t_ax(k) + [-1 1 1 -1]/(2*framerate), [0 0 1.05 1.05], ...
          [0.82 0.82 0.82], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
          'HandleVisibility', 'off');
end
legend(ax_t1, {'Hybrid (fused)','MI-only','CVSA-only'}, 'Location', 'best', 'FontSize', 9);
xlabel(ax_t1, 't from CF onset (s)');
ylabel(ax_t1, 'buffer(target class)  mean +- SEM');
title(ax_t1, sprintf('Integrator output for TARGET class — aligned to CF onset (n=%d trials total)', n_tot), ...
     'FontWeight', 'bold');
grid(ax_t1, 'on');

% Compute per-trial delta trajectories for the two comparisons.
% D_mi(i,k)  = hybrid buffer(target) - MI-only buffer(target) at frame k, trial i.
% Positive = Hybrid pushed the buffer higher toward the win threshold.
D_mi = B_hyb - B_mi_mat;
D_cv = B_hyb - B_cv_mat;
mn_d_mi_t = mean(D_mi, 1, 'omitnan');
se_d_mi_t = std(D_mi, [], 1, 'omitnan') ./ sqrt(max(n_valid, 1));
mn_d_cv_t = mean(D_cv, 1, 'omitnan');
se_d_cv_t = std(D_cv, [], 1, 'omitnan') ./ sqrt(max(n_valid, 1));

% Helper: shade significant runs on a delta-trajectory axes.
% UniformOutput=false prevents cellfun from trying to concatenate patch handles.
shade_sig = @(ax, p_vec, col) cellfun(@(seg) patch(ax, ...
    [t_ax(seg(1)) t_ax(seg(end)) t_ax(seg(end)) t_ax(seg(1))], ...
    [-1 -1 1 1], col, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off'), ...
    find_runs(p_vec < 0.05), 'UniformOutput', false);

% Subplot 2: DELTA = Hybrid - MI-only over time
%   Positive = Hybrid pushed the control signal higher for the cued class.
%   Coloured shading = frames where sign-flip test p<0.05 (exploratory).
ax_t2 = subplot(3,1,2); hold(ax_t2, 'on');
shade_sig(ax_t2, p_t_mi, COL_MI);
fill(ax_t2, [t_ax, fliplr(t_ax)], [mn_d_mi_t+se_d_mi_t, fliplr(mn_d_mi_t-se_d_mi_t)], ...
     COL_MI, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax_t2, t_ax, mn_d_mi_t, '-', 'Color', COL_MI, 'LineWidth', 2.2);
yline(ax_t2, 0, 'k-', 'HandleVisibility', 'off');
if isfield(int_cfg, 'cvsa_influence')
    xline(ax_t2, int_cfg.cvsa_influence, '--', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
    text(ax_t2, int_cfg.cvsa_influence + 0.05, 0, '\alpha\rightarrow0', ...
         'FontSize', 8, 'Color', [0.55 0.55 0.55], 'VerticalAlignment', 'bottom');
end
[pk_mi, pk_mi_k] = max(mn_d_mi_t);
if ~isnan(pk_mi) && pk_mi > 0
    text(ax_t2, t_ax(pk_mi_k), pk_mi + se_d_mi_t(pk_mi_k) + 0.005, ...
         sprintf('peak +%.3f\n@ %.1fs', pk_mi, t_ax(pk_mi_k)), ...
         'FontSize', 8, 'Color', COL_MI, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
lo_idx_d = find(n_valid < MIN_N & n_valid > 0);
for k = lo_idx_d
    patch(ax_t2, t_ax(k) + [-1 1 1 -1]/(2*framerate), [-1 -1 1 1], ...
          [0.82 0.82 0.82], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
end
yl_mi_v = abs([mn_d_mi_t + se_d_mi_t, mn_d_mi_t - se_d_mi_t]);
yl_mi = max(yl_mi_v(~isnan(yl_mi_v))) * 1.5;
if isempty(yl_mi) || ~isfinite(yl_mi) || yl_mi < 0.02, yl_mi = 0.1; end
set(ax_t2, 'YLim', [-yl_mi, yl_mi]);
xlabel(ax_t2, 't from CF onset (s)');
ylabel(ax_t2, '\Delta buffer(target)  [Hybrid - MI-only]');
title(ax_t2, sprintf('Hybrid vs MI-only advantage over time — mean ± SEM  (shading = p<0.05, %d perms, exploratory)', N_PERM_T), ...
     'FontWeight', 'bold');
grid(ax_t2, 'on');

% Subplot 3: DELTA = Hybrid - CVSA-only over time
ax_t3 = subplot(3,1,3); hold(ax_t3, 'on');
shade_sig(ax_t3, p_t_cv, COL_CV);
fill(ax_t3, [t_ax, fliplr(t_ax)], [mn_d_cv_t+se_d_cv_t, fliplr(mn_d_cv_t-se_d_cv_t)], ...
     COL_CV, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax_t3, t_ax, mn_d_cv_t, '-', 'Color', COL_CV, 'LineWidth', 2.2);
yline(ax_t3, 0, 'k-', 'HandleVisibility', 'off');
if isfield(int_cfg, 'cvsa_influence')
    xline(ax_t3, int_cfg.cvsa_influence, '--', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
end
[pk_cv, pk_cv_k] = max(mn_d_cv_t);
if ~isnan(pk_cv) && pk_cv > 0
    text(ax_t3, t_ax(pk_cv_k), pk_cv + se_d_cv_t(pk_cv_k) + 0.005, ...
         sprintf('peak +%.3f\n@ %.1fs', pk_cv, t_ax(pk_cv_k)), ...
         'FontSize', 8, 'Color', COL_CV, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
for k = lo_idx_d
    patch(ax_t3, t_ax(k) + [-1 1 1 -1]/(2*framerate), [-1 -1 1 1], ...
          [0.82 0.82 0.82], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
end
yl_cv_v = abs([mn_d_cv_t + se_d_cv_t, mn_d_cv_t - se_d_cv_t]);
yl_cv = max(yl_cv_v(~isnan(yl_cv_v))) * 1.5;
if isempty(yl_cv) || ~isfinite(yl_cv) || yl_cv < 0.02, yl_cv = 0.1; end
set(ax_t3, 'YLim', [-yl_cv, yl_cv]);
xlabel(ax_t3, 't from CF onset (s)');
ylabel(ax_t3, '\Delta buffer(target)  [Hybrid - CVSA-only]');
title(ax_t3, sprintf('Hybrid vs CVSA-only advantage over time — mean ± SEM  (shading = p<0.05, %d perms, exploratory)', N_PERM_T), ...
     'FontWeight', 'bold');
grid(ax_t3, 'on');

sgtitle(fig4, sprintf(['Hybrid advantage in control signal over time | n=%d trials | ' ...
        'zero line = no advantage | grey = n<%d trials (unreliable) | shading = p<0.05 uncorrected'], ...
        n_tot, MIN_N), 'Interpreter', 'none');

saveas(fig4, fullfile(out_dir, 'temporal_significance.svg'), 'svg');
if ~SHOW_FIGURES, close(fig4); end

%% --- Effect size + bootstrap CI (no Stats Toolbox needed) ----------------
%   Cohen's d (paired): d = mean(delta) / std(delta) — for paired data this
%   is equivalent to the one-sample t-test effect size on the differences.
d_mi_cohen  = mean(d_mi_all,  'omitnan') / std(d_mi_all,  'omitnan');
d_cvs_cohen = mean(d_cvs_all, 'omitnan') / std(d_cvs_all, 'omitnan');
fprintf('\n  Cohen''s d (paired buffer advantage, target class):\n');
fprintf('    Hybrid vs MI-only   : d=%.3f\n', d_mi_cohen);
fprintf('    Hybrid vs CVSA-only : d=%.3f\n', d_cvs_cohen);
fprintf('    (|d|>0.2 small, >0.5 medium, >0.8 large)\n');

%   Bootstrap 95%% CI on accuracy per stream and on accuracy deltas.
n_boot = 2000;
boot_acc = zeros(n_boot, 3);
for b = 1:n_boot
    idx = randi(n_tot, n_tot, 1);
    boot_acc(b,:) = [mean(all_oc_hyb(idx,1)==1), ...
                     mean(all_oc_mi(idx,1)==1), ...
                     mean(all_oc_cvsa(idx,1)==1)];
end
boot_sorted  = sort(boot_acc, 1);
lo_b = max(1,      round(0.025 * n_boot));
hi_b = min(n_boot, round(0.975 * n_boot));
ci_lo_acc = boot_sorted(lo_b, :);
ci_hi_acc = boot_sorted(hi_b, :);
boot_d_mi  = sort(boot_acc(:,1) - boot_acc(:,2));
boot_d_cvs = sort(boot_acc(:,1) - boot_acc(:,3));
fprintf('\n  Accuracy bootstrap 95%%%% CI (n_boot=%d, percentile method):\n', n_boot);
for s = 1:3
    fprintf('    %-10s: %.0f%%  [%.0f%%, %.0f%%]\n', stream_names{s}, ...
            100*acc(s), 100*ci_lo_acc(s), 100*ci_hi_acc(s));
end
fprintf('  Accuracy delta bootstrap 95%%%% CI:\n');
fprintf('    Hybrid - MI-only   : %+.0f%%  [%+.0f%%, %+.0f%%]\n', ...
        100*(acc(1)-acc(2)), 100*boot_d_mi(lo_b),  100*boot_d_mi(hi_b));
fprintf('    Hybrid - CVSA-only : %+.0f%%  [%+.0f%%, %+.0f%%]\n', ...
        100*(acc(1)-acc(3)), 100*boot_d_cvs(lo_b), 100*boot_d_cvs(hi_b));

%% --- Figure 5: deep advantage analysis -----------------------------------
%
%   Four panels designed to make the strongest possible case that hybrid
%   outperforms MI-only and CVSA-only on this data:
%
%   (1,1) Cumulative HIT fraction over time — shows hybrid hits MORE trials
%         AND hits them FASTER (curve shifts left and up).
%   (1,2) Per-trial scatter: meanBuf_mi vs meanBuf_hyb (target class).
%         Points above the diagonal = fusion pushed the buffer higher for
%         the cued class.  Colour = simulated hybrid outcome.
%   (2,1) Same scatter for Hybrid vs CVSA-only.
%   (2,2) Accuracy bar with bootstrap 95%% CI (more honest for small n than
%         binomial SE).  CIs that don't overlap with Hybrid suggest a
%         reliable advantage.

max_t_cf = max([all_oc_hyb(:,2); all_oc_mi(:,2); all_oc_cvsa(:,2)]);
t_fine   = linspace(0, max_t_cf, 300);
all_ocs  = {all_oc_hyb, all_oc_mi, all_oc_cvsa};
cum_hit  = zeros(3, numel(t_fine));
for s = 1:3
    for ki = 1:numel(t_fine)
        cum_hit(s, ki) = mean(all_ocs{s}(:,1) == 1 & all_ocs{s}(:,2) <= t_fine(ki));
    end
end
% per-class cumulative (used for dashed thin lines)
classes_all_u = unique(all_trial_class(~isnan(all_trial_class)))';
cum_hit_cls = zeros(3, numel(classes_all_u), numel(t_fine));
for s = 1:3
    for ci2 = 1:numel(classes_all_u)
        mc = all_trial_class == classes_all_u(ci2);
        if ~any(mc), continue; end
        for ki = 1:numel(t_fine)
            cum_hit_cls(s, ci2, ki) = ...
                sum(all_ocs{s}(mc,1)==1 & all_ocs{s}(mc,2)<=t_fine(ki)) / n_tot;
        end
    end
end

fig5 = figure('Name', 'Hybrid vs MI/CVSA — deep advantage analysis', ...
              'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig5, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Panel (1,1): cumulative HIT fraction over time
ax_c = subplot(2,2,1); hold(ax_c, 'on');
ls_cls = {'--', ':'};
cls_labels = arrayfun(@(c) sprintf('class %d', c), classes_all_u, 'UniformOutput', false);
% per-class thin lines first (background)
for ci2 = 1:numel(classes_all_u)
    for s = 1:3
        plot(ax_c, t_fine, squeeze(cum_hit_cls(s,ci2,:)), ...
             ls_cls{min(ci2,2)}, 'Color', stream_cols{s}, ...
             'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end
% overall thick lines
h_lines = gobjects(1,3);
for s = 1:3
    h_lines(s) = plot(ax_c, t_fine, cum_hit(s,:), '-', ...
                      'Color', stream_cols{s}, 'LineWidth', 2.5);
end
if isfield(int_cfg, 'cvsa_influence')
    xline(ax_c, int_cfg.cvsa_influence, '--', 'Color', [0.55 0.55 0.55], ...
          'HandleVisibility', 'off');
end
% Final accuracy labels and advantage annotation
for s = 1:3
    text(ax_c, max_t_cf, cum_hit(s,end), sprintf(' %.0f%%', 100*cum_hit(s,end)), ...
         'Color', stream_cols{s}, 'FontSize', 9, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
end
gap_mi_cum  = cum_hit(1,end) - cum_hit(2,end);
gap_cvs_cum = cum_hit(1,end) - cum_hit(3,end);
text(ax_c, max_t_cf*0.05, 0.97, ...
     sprintf('Hybrid advantage:\nvs MI-only:   %+.0f%%\nvs CVSA-only: %+.0f%%', ...
             100*gap_mi_cum, 100*gap_cvs_cum), ...
     'FontSize', 9, 'FontWeight', 'bold', 'Color', COL_HYB, ...
     'VerticalAlignment', 'top', 'BackgroundColor', [1 1 1 0.75], 'EdgeColor', [0.7 0.7 0.7]);
xlabel(ax_c, 't from CF onset (s)'); ylabel(ax_c, 'fraction of all trials with HIT');
legend(ax_c, h_lines, stream_names, 'Location', 'best', 'FontSize', 9);
title(ax_c, sprintf(['Cumulative HIT fraction over time (n=%d trials)\n' ...
      'curve shifts LEFT = hits faster; ends HIGHER = more accurate'], n_tot), ...
      'FontWeight', 'bold');
set(ax_c, 'YLim', [0, 1]); grid(ax_c, 'on');

% Pre-define outcome color matrix to avoid invalid bracket-index syntax.
OUTCOME_COLS = [COL_HIT; COL_MISS; COL_TO];   % rows: HIT, MISS, TO

% Panel (1,2): per-trial scatter Hybrid vs MI-only (target-class buffer)
ax_s1 = subplot(2,2,2); hold(ax_s1, 'on');
hyb_oc_codes = all_oc_hyb(:,1);
scatter_markers = {'o', 's', 'd'};   % HIT, MISS, TO
scatter_labels  = {'HIT (hyb)', 'MISS (hyb)', 'TO (hyb)'};
for g = 1:3
    m_g = hyb_oc_codes == g;
    if ~any(m_g), continue; end
    scatter(ax_s1, all_meanBuf_mi(m_g), all_meanBuf_hyb(m_g), 55, ...
            scatter_markers{g}, 'MarkerFaceColor', OUTCOME_COLS(g,:), ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.75);
end
ax_lim1 = [min([all_meanBuf_mi; all_meanBuf_hyb])-0.01, ...
            max([all_meanBuf_mi; all_meanBuf_hyb])+0.01];
plot(ax_s1, ax_lim1, ax_lim1, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
text(ax_s1, ax_lim1(2)-0.02, ax_lim1(2)-0.04, 'y=x', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
fill(ax_s1, [ax_lim1(1) ax_lim1(2) ax_lim1(2) ax_lim1(1)], ...
           [ax_lim1(1) ax_lim1(1) ax_lim1(2) ax_lim1(2)], ...
     COL_HYB, 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
n_above_s1 = sum(all_meanBuf_hyb > all_meanBuf_mi);
text(ax_s1, ax_lim1(1)+0.01, ax_lim1(2)-0.02, ...
     sprintf('%d/%d above diag\n= Hybrid higher', n_above_s1, n_tot), ...
     'FontSize', 9, 'Color', COL_HYB, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
legend(ax_s1, scatter_labels, 'Location', 'southeast', 'FontSize', 8);
xlabel(ax_s1, 'mean buffer(target)  MI-only');
ylabel(ax_s1, 'mean buffer(target)  Hybrid');
title(ax_s1, sprintf('Hybrid vs MI-only — per-trial control signal\n(%d/%d trials: Hybrid pushed buffer higher for cued class)', ...
     n_above_s1, n_tot), 'FontWeight', 'bold');
set(ax_s1, 'XLim', ax_lim1, 'YLim', ax_lim1); axis(ax_s1, 'equal'); grid(ax_s1, 'on');

% Panel (2,1): per-trial scatter Hybrid vs CVSA-only
ax_s2 = subplot(2,2,3); hold(ax_s2, 'on');
for g = 1:3
    m_g = hyb_oc_codes == g;
    if ~any(m_g), continue; end
    scatter(ax_s2, all_meanBuf_cvsa(m_g), all_meanBuf_hyb(m_g), 55, ...
            scatter_markers{g}, 'MarkerFaceColor', OUTCOME_COLS(g,:), ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.75);
end
ax_lim2 = [min([all_meanBuf_cvsa; all_meanBuf_hyb])-0.01, ...
            max([all_meanBuf_cvsa; all_meanBuf_hyb])+0.01];
plot(ax_s2, ax_lim2, ax_lim2, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
text(ax_s2, ax_lim2(2)-0.02, ax_lim2(2)-0.04, 'y=x', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
fill(ax_s2, [ax_lim2(1) ax_lim2(2) ax_lim2(2) ax_lim2(1)], ...
           [ax_lim2(1) ax_lim2(1) ax_lim2(2) ax_lim2(2)], ...
     COL_HYB, 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
n_above_s2 = sum(all_meanBuf_hyb > all_meanBuf_cvsa);
text(ax_s2, ax_lim2(1)+0.01, ax_lim2(2)-0.02, ...
     sprintf('%d/%d above diag\n= Hybrid higher', n_above_s2, n_tot), ...
     'FontSize', 9, 'Color', COL_HYB, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
legend(ax_s2, scatter_labels, 'Location', 'southeast', 'FontSize', 8);
xlabel(ax_s2, 'mean buffer(target)  CVSA-only');
ylabel(ax_s2, 'mean buffer(target)  Hybrid');
title(ax_s2, sprintf('Hybrid vs CVSA-only — per-trial control signal\n(%d/%d trials: Hybrid pushed buffer higher for cued class)', ...
     n_above_s2, n_tot), 'FontWeight', 'bold');
set(ax_s2, 'XLim', ax_lim2, 'YLim', ax_lim2); axis(ax_s2, 'equal'); grid(ax_s2, 'on');

% Panel (2,2): accuracy bars with bootstrap 95%% CI
ax_b = subplot(2,2,4); hold(ax_b, 'on');
for s = 1:3
    bar(ax_b, s, 100*acc(s), 0.6, 'FaceColor', stream_cols{s}, 'EdgeColor', 'k');
    errorbar(ax_b, s, 100*acc(s), ...
             100*(acc(s)-ci_lo_acc(s)), 100*(ci_hi_acc(s)-acc(s)), ...
             'k', 'LineWidth', 1.5, 'CapSize', 8);
    text(ax_b, s, 100*ci_hi_acc(s)+2, ...
         sprintf('%.0f%%\n[%.0f, %.0f]', 100*acc(s), 100*ci_lo_acc(s), 100*ci_hi_acc(s)), ...
         'HorizontalAlignment', 'center', 'FontSize', 8);
end
yline(ax_b, 50, 'k--', 'HandleVisibility', 'off');
set(ax_b, 'XTick', 1:3, 'XTickLabel', stream_names, 'YLim', [0, 108]);
ylabel(ax_b, 'HIT rate (%)');
title(ax_b, sprintf('Accuracy with bootstrap 95%%%% CI\n(n_boot=%d — asymmetric CI honest for small n)', n_boot), ...
     'FontWeight', 'bold');
grid(ax_b, 'on');

sgtitle(fig5, sprintf(['Hybrid advantage — deeper analysis (%d trials, %d file(s))\n' ...
        'Cohen''s d: Hyb-MI=%.2f, Hyb-CVSA=%.2f | ' ...
        'Accuracy delta [95%%%%CI]: Hyb-MI=%+.0f%% [%+.0f,+%.0f], ' ...
        'Hyb-CVSA=%+.0f%% [%+.0f,+%.0f]'], ...
        n_tot, n_files, d_mi_cohen, d_cvs_cohen, ...
        100*(acc(1)-acc(2)), 100*boot_d_mi(lo_b),  100*boot_d_mi(hi_b), ...
        100*(acc(1)-acc(3)), 100*boot_d_cvs(lo_b), 100*boot_d_cvs(hi_b)), ...
        'Interpreter', 'none');

saveas(fig5, fullfile(out_dir, 'advantage_deep.svg'), 'svg');
fprintf('\nSaved all figures to %s\n', out_dir);
if ~SHOW_FIGURES, close(fig5); end


% ── Local helpers (must be after all script statements) ────────────────────

function [outcome_code, t_event_s] = classify_trial_outcome(tr, thresholds, framerate)
% CLASSIFY_TRIAL_OUTCOME  Offline HIT/MISS/TIMEOUT classification from the
%   integrator buffer itself (matches Training.cpp is_target_hit: first
%   class i for which integrated[i] >= thresholds[i] - 5e-3 wins).
%     outcome_code: 0=unknown (no target class), 1=HIT, 2=MISS, 3=TIMEOUT
%     t_event_s   : time from CF onset to the event (s); for TIMEOUT this
%                   is the CF window duration.
    HIT = 1; MISS = 2; TIMEOUT = 3;
    cf_range = (tr.n_pre+1):(tr.n_pre+tr.n_cf);
    t_event_s = (tr.n_cf - 1) / framerate;
    if isnan(tr.target_class)
        outcome_code = 0;
        return;
    end
    outcome_code = TIMEOUT;
    for idx = 1:numel(cf_range)
        k = cf_range(idx);
        if tr.integrated(k,1) >= thresholds(1) - 5e-3
            winner = 1;
        elseif tr.integrated(k,2) >= thresholds(2) - 5e-3
            winner = 2;
        else
            continue;
        end
        if winner == tr.target_class, outcome_code = HIT; else, outcome_code = MISS; end
        t_event_s = (idx - 1) / framerate;
        return;
    end
end

function s = fmt_mean_sem(m, sem)
% FMT_MEAN_SEM  "mean+-sem" string, or "--" if no data.
    if isnan(m), s = '--'; else, s = sprintf('%.2f+-%.2f', m, sem); end
end

function s = fmt_cmp(labels, vals)
% FMT_CMP  "label1=n1, label2=n2, ..." string for console summaries.
    parts = cell(1, numel(labels));
    for i = 1:numel(labels)
        parts{i} = sprintf('%s=%d', labels{i}, vals(i));
    end
    s = strjoin(parts, ', ');
end

function plot_mean_sem_bars(ax, means, sems, ns, stream_names, stream_cols)
% PLOT_MEAN_SEM_BARS  One bar per stream with mean+-SEM error bar and an
%   "n=" annotation; streams with no data (NaN mean) are skipped.
    for s = 1:numel(means)
        if isnan(means(s)), continue; end
        bar(ax, s, means(s), 0.6, 'FaceColor', stream_cols{s}, 'EdgeColor', 'k');
        errorbar(ax, s, means(s), sems(s), 'k', 'LineWidth', 1.2, 'CapSize', 6);
        text(ax, s, means(s) + sems(s) + 0.05, sprintf('n=%d', ns(s)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    set(ax, 'XTick', 1:numel(stream_names), 'XTickLabel', stream_names);
end

function p = sign_flip_test(diffs, n_perm)
% SIGN_FLIP_TEST  Two-sided paired permutation test (sign-flip null) on
%   H0: mean(diffs) = 0.  No Statistics Toolbox required.
%   Returns NaN if fewer than 2 valid (non-NaN) observations.
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

function s = p_to_stars(p)
% P_TO_STARS  Convert a p-value to a significance-star string.
    if isnan(p),      s = '';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,              s = 'n.s.';
    end
end

function draw_sig_bracket(ax, x1, x2, y, label)
% DRAW_SIG_BRACKET  Horizontal significance bracket between two bar positions.
    if isempty(label) || strcmp(label, ''), return; end
    tick = 0.005;
    plot(ax, [x1 x1 x2 x2], [y-tick y y y-tick], 'k-', 'HandleVisibility', 'off');
    text(ax, (x1+x2)/2, y + tick, label, ...
         'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

function segs = find_runs(mask)
% FIND_RUNS  Return cell array of contiguous-true index vectors from a logical mask.
    mask = mask(:)';
    segs = {};
    k = 1;
    while k <= numel(mask)
        if mask(k)
            k2 = k;
            while k2 <= numel(mask) && mask(k2), k2 = k2 + 1; end
            segs{end+1} = k:(k2-1); %#ok<AGROW>
            k = k2;
        else
            k = k + 1;
        end
    end
end
