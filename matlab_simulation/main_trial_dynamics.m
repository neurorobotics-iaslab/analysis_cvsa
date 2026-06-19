%% MAIN_TRIAL_DYNAMICS  Within-session learning/fatigue trends across trial order.
%
%   Loads one or more evaluation GDFs of ANY paradigm (MI/CVSA/Hybrid mixed,
%   paradigm inferred from filename like main_session_overview). Uses ONLY
%   real GDF events (897/898/899 relative to 781) — no CSP/sLDA/artifact
%   simulation, no companion YAML required — so it is fast and works even
%   when the calibration model files are missing.
%
%   For each file, trial index within that session is normalised to
%   fractional progress (0 = first trial, 1 = last trial), so sessions of
%   different length can be pooled within the same paradigm. This tests
%   whether performance changes systematically over the COURSE of a single
%   session — learning (improves) vs fatigue (degrades) — separately for
%   each paradigm.
%
%   Three figures (one per paradigm found, MI/CVSA/Hybrid):
%     Fig 1 — accuracy in the first half vs second half of each session
%             (paired per-file dots + per-paradigm mean), sign-flip test on
%             the paired delta (second half - first half)
%     Fig 2 — time-to-HIT vs normalised trial position (all HIT trials,
%             pooled across files of that paradigm), linear fit, Pearson r
%             with a permutation-based p-value (shuffle trial position)
%     Fig 3 — accuracy across quartiles of session progress, per paradigm
%
%   Console: per-file half-split table; per-paradigm pooled statistics.
%
%   Saves trial_dynamics_summary.mat (flat per-trial table: file, paradigm,
%   trial_idx, frac, outcome_code, dt) under
%   <gdf_dir>/analysis_results/trial_dynamics/.

clear; clc; close all;

SHOW_FIGURES = false;
if SHOW_FIGURES, fig_vis = 'on'; else, fig_vis = 'off'; end

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir,'io'), fullfile(this_dir,'utils'));

HIT_EV = 897; MISS_EV = 898; TO_EV = 899; CF_EV = 781;

default_dir = '/home/paolo/bci_vr_ws/recordings';
[gdf_names, gdf_dir] = uigetfile({'*.gdf','GDF files (*.gdf)'}, ...
    'Select evaluation GDF file(s)', default_dir, 'MultiSelect', 'on');
if isequal(gdf_names, 0), error('main_trial_dynamics:cancel', 'No file selected.'); end
if ischar(gdf_names), gdf_names = {gdf_names}; end
n_files = numel(gdf_names);

COL = struct('mi', [0.85 0.30 0.10], 'cvsa', [0.10 0.60 0.30], 'hybrid', [0.18 0.45 0.75]);

% ── Flat per-trial table across all files ───────────────────────────────────
flat = struct('file_label', {}, 'paradigm', {}, 'trial_idx', {}, 'n_trials', {}, ...
              'frac', {}, 'outcome_code', {}, 'dt', {});

fprintf('\n══════════════════ Per-file half-split accuracy ══════════════════\n');
fprintf('  %-30s  %-7s  %8s  %8s  %8s\n', 'file', 'paradigm', 'n_trials', '1st-half', '2nd-half');

for fi = 1:n_files
    gdf_path = fullfile(gdf_dir, gdf_names{fi});
    [~, basename] = fileparts(gdf_path);
    paradigm = detect_paradigm(basename);

    [~, header, ~] = load_gdf(gdf_path);
    fs  = header.SampleRate;
    POS = header.EVENT.POS;
    TYP = header.EVENT.TYP;

    cf_pos = POS(TYP == CF_EV);
    n_cf   = numel(cf_pos);
    if n_cf < 2, continue; end

    outcomes = zeros(1, n_cf);
    dt_vals  = nan(1, n_cf);
    for t = 1:n_cf
        after = POS > cf_pos(t);
        idx = find((TYP==HIT_EV | TYP==MISS_EV | TYP==TO_EV) & after, 1);
        if isempty(idx), continue; end
        outcomes(t) = TYP(idx);
        dt_vals(t)  = (POS(idx) - cf_pos(t)) / fs;
    end

    for t = 1:n_cf
        idx = numel(flat) + 1;
        flat(idx).file_label   = basename;
        flat(idx).paradigm     = paradigm;
        flat(idx).trial_idx    = t;
        flat(idx).n_trials     = n_cf;
        flat(idx).frac         = (t-1) / max(1, n_cf-1);
        flat(idx).outcome_code = outcomes(t);
        flat(idx).dt           = NaN;
        if outcomes(t) == HIT_EV, flat(idx).dt = dt_vals(t); end
    end

    half1 = 1:floor(n_cf/2);
    half2 = (floor(n_cf/2)+1):n_cf;
    acc1 = mean(outcomes(half1) == HIT_EV);
    acc2 = mean(outcomes(half2) == HIT_EV);
    fprintf('  %-30s  %-7s  %8d  %7.0f%%  %7.0f%%\n', basename, paradigm, n_cf, 100*acc1, 100*acc2);
end
fprintf('═══════════════════════════════════════════════════════════════════\n');

pars = {'mi', 'cvsa', 'hybrid'};
par_labels = {'MI', 'CVSA', 'Hybrid'};
par_present = ismember(pars, unique({flat.paradigm}));

out_dir = fullfile(gdf_dir, 'analysis_results', 'trial_dynamics');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

n_perm = 2000;

% ═══════════════════════════════════════════════════════════════════════════
%  FIG 1 — first-half vs second-half accuracy, per file, per paradigm
% ═══════════════════════════════════════════════════════════════════════════
fig1 = figure('Name', 'Trial Dynamics — Within-session accuracy (1st vs 2nd half)', ...
              'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig1, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

fprintf('\n══════════════════ Within-session accuracy trend (per paradigm) ══════════════════\n');
n_par_present = sum(par_present);
sp = 0;
for pi = 1:3
    if ~par_present(pi), continue; end
    sp = sp + 1;
    files_p = unique({flat(strcmp({flat.paradigm}, pars{pi})).file_label});
    n_fp = numel(files_p);
    d1 = nan(1, n_fp); d2 = nan(1, n_fp);
    for fi = 1:n_fp
        mask = strcmp({flat.file_label}, files_p{fi});
        rows = flat(mask);
        n_t  = rows(1).n_trials;
        half1 = [rows.trial_idx] <= floor(n_t/2);
        half2 = ~half1;
        d1(fi) = mean([rows(half1).outcome_code] == HIT_EV);
        d2(fi) = mean([rows(half2).outcome_code] == HIT_EV);
    end
    delta = d2 - d1;
    p_val = sign_flip_test_local(delta, n_perm);

    ax = subplot(1, n_par_present, sp); hold(ax, 'on');
    for fi = 1:n_fp
        plot(ax, [1 2], 100*[d1(fi), d2(fi)], '-o', 'Color', COL.(pars{pi}), ...
             'MarkerFaceColor', COL.(pars{pi}), 'LineWidth', 1, 'MarkerSize', 6);
    end
    m1 = mean(d1, 'omitnan'); m2 = mean(d2, 'omitnan');
    plot(ax, [1 2], 100*[m1 m2], '-', 'Color', 'k', 'LineWidth', 2.5);
    scatter(ax, [1 2], 100*[m1 m2], 80, 'k', 'filled', 'Marker', 's');
    set(ax, 'XTick', [1 2], 'XTickLabel', {'1st half', '2nd half'}, 'XLim', [0.7, 2.3], 'YLim', [0, 105]);
    ylabel(ax, 'HIT rate (%)');
    title(ax, sprintf('%s (n=%d files)\n\\Delta=%+.1f%%  p=%.3f  %s', ...
          par_labels{pi}, n_fp, 100*mean(delta,'omitnan'), p_val, stars_local(p_val)), 'FontWeight', 'bold');
    grid(ax, 'on');

    fprintf('  %-8s  n_files=%2d  mean 1st-half=%.0f%%  mean 2nd-half=%.0f%%  delta=%+.1f%%  p=%.4f  %s\n', ...
            par_labels{pi}, n_fp, 100*m1, 100*m2, 100*mean(delta,'omitnan'), p_val, stars_local(p_val));
end
fprintf('  (paired sign-flip permutation across files, %d perms; positive delta = 2nd half better = learning)\n', n_perm);
fprintf('═══════════════════════════════════════════════════════════════════════════════════\n');
sgtitle(fig1, 'Within-session accuracy: 1st half vs 2nd half (learning if 2nd > 1st, fatigue if 2nd < 1st)', ...
        'Interpreter', 'none');
saveas(fig1, fullfile(out_dir, 'trial_dynamics_half_split.svg'), 'svg');
if ~SHOW_FIGURES, close(fig1); end

% ═══════════════════════════════════════════════════════════════════════════
%  FIG 2 — TTH vs normalised trial position, per paradigm
% ═══════════════════════════════════════════════════════════════════════════
fig2 = figure('Name', 'Trial Dynamics — TTH vs session progress', ...
              'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

fprintf('\n══════════════════ TTH vs session-progress trend (per paradigm) ══════════════════\n');
sp = 0;
for pi = 1:3
    if ~par_present(pi), continue; end
    sp = sp + 1;
    mask = strcmp({flat.paradigm}, pars{pi}) & ~isnan([flat.dt]);
    frac_v = [flat(mask).frac];
    dt_v   = [flat(mask).dt];

    ax = subplot(1, n_par_present, sp); hold(ax, 'on');
    scatter(ax, frac_v, dt_v, 35, COL.(pars{pi}), 'filled', 'MarkerFaceAlpha', 0.55);
    r_val = NaN; p_perm = NaN;
    if numel(frac_v) >= 3
        p_fit = polyfit(frac_v, dt_v, 1);
        xfit = linspace(0, 1, 50);
        plot(ax, xfit, polyval(p_fit, xfit), 'k--', 'LineWidth', 1.5);
        C = corrcoef(frac_v, dt_v);
        r_val = C(1,2);
        p_perm = corr_perm_test_local(frac_v, dt_v, r_val, n_perm);
    end
    set(ax, 'XLim', [0,1]);
    xlabel(ax, 'normalised trial position (0=first, 1=last)');
    ylabel(ax, 'time to HIT (s)');
    title(ax, sprintf('%s  (n=%d HIT trials)\nr=%.2f  p=%.3f  %s', ...
          par_labels{pi}, numel(frac_v), r_val, p_perm, stars_local(p_perm)), 'FontWeight', 'bold');
    grid(ax, 'on');

    fprintf('  %-8s  n_hit_trials=%4d  r=%+.3f  p=%.4f  %s  (negative r = gets faster = learning)\n', ...
            par_labels{pi}, numel(frac_v), r_val, p_perm, stars_local(p_perm));
end
fprintf('  (Pearson r, permutation test shuffling trial position, %d perms)\n', n_perm);
fprintf('═══════════════════════════════════════════════════════════════════════════════════\n');
sgtitle(fig2, 'Time-to-HIT vs normalised session progress (negative trend = speeding up = learning)', ...
        'Interpreter', 'none');
saveas(fig2, fullfile(out_dir, 'trial_dynamics_tth_trend.svg'), 'svg');
if ~SHOW_FIGURES, close(fig2); end

% ═══════════════════════════════════════════════════════════════════════════
%  FIG 3 — accuracy by quartile of session progress, per paradigm
% ═══════════════════════════════════════════════════════════════════════════
fig3 = figure('Name', 'Trial Dynamics — Accuracy by session quartile', ...
              'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig3, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
ax = axes(fig3); hold(ax, 'on'); %#ok<LAXES>

q_edges = [0, 0.25, 0.5, 0.75, 1.0001];
q_labels = {'Q1 (0-25%)', 'Q2 (25-50%)', 'Q3 (50-75%)', 'Q4 (75-100%)'};
bw = 0.2;
present_idx = find(par_present);
for k = 1:numel(present_idx)
    pi = present_idx(k);
    mask_p = strcmp({flat.paradigm}, pars{pi});
    frac_p = [flat(mask_p).frac];
    oc_p   = [flat(mask_p).outcome_code];
    q_acc = nan(1, 4);
    for q = 1:4
        m = frac_p >= q_edges(q) & frac_p < q_edges(q+1);
        if any(m), q_acc(q) = mean(oc_p(m) == HIT_EV); end
    end
    x = (1:4) + (k - (numel(present_idx)+1)/2) * bw;
    bar(ax, x, 100*q_acc, bw*0.9, 'FaceColor', COL.(pars{pi}), 'EdgeColor', 'k', ...
        'DisplayName', par_labels{pi});
end
set(ax, 'XTick', 1:4, 'XTickLabel', q_labels, 'YLim', [0, 105]);
ylabel(ax, 'HIT rate (%)');
legend(ax, 'Location', 'best', 'FontSize', 9);
title(ax, 'Accuracy across session quartiles, per paradigm', 'FontWeight', 'bold');
grid(ax, 'on');
sgtitle(fig3, 'Accuracy by quartile of session progress', 'Interpreter', 'none');
saveas(fig3, fullfile(out_dir, 'trial_dynamics_quartiles.svg'), 'svg');
if ~SHOW_FIGURES, close(fig3); end

% ── Save flat per-trial table ───────────────────────────────────────────────
trial_dynamics_summary = flat;
save(fullfile(out_dir, 'trial_dynamics_summary.mat'), 'trial_dynamics_summary');
fprintf('\nSaved figures and trial_dynamics_summary.mat to %s\n', out_dir);


% ── Local helpers ────────────────────────────────────────────────────────────

function par = detect_paradigm(fname)
    s = lower(fname);
    if contains(s,'hybrid'),   par = 'hybrid';
    elseif contains(s,'cvsa'), par = 'cvsa';
    elseif contains(s,'mi'),   par = 'mi';
    else,                      par = 'unknown';
    end
end

function s = stars_local(p)
    if isnan(p),      s = '';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,              s = 'n.s.';
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

function p = corr_perm_test_local(x, y, r_obs, n_perm)
    n = numel(x);
    if n < 3, p = NaN; return; end
    cnt = 0;
    for i = 1:n_perm
        y_shuf = y(randperm(n));
        C = corrcoef(x, y_shuf);
        if abs(C(1,2)) >= abs(r_obs) - 1e-12
            cnt = cnt + 1;
        end
    end
    p = (cnt + 1) / (n_perm + 1);
end
