%% MAIN_TRIAL_DYNAMICS  Within-session learning/fatigue trends across trial order.
%
%   QUESTION THIS ANSWERS: within a single evaluation session, do LATER
%   trials go better or worse than EARLIER trials? E.g. a 20-trial MI
%   session: is trial 18 more likely to be a HIT than trial 2? Faster?
%   This is a trial-ORDER effect (learning/fatigue across the session),
%   NOT the within-trial ERD/ERS time-course (that is topo_erders.m).
%
%   WHY THIS MATTERS BEYOND ITSELF: it is a validity/confound check for
%   every other cross-paradigm comparison in this package. If MI/CVSA/Hybrid
%   sessions are recorded in a fixed order within a day and performance
%   drifts across that order, an apparent "paradigm X is worse" finding
%   elsewhere (main_session_overview, main_validate_counterfactual,
%   main_group_analysis) could really be a time-of-day/fatigue effect
%   confounded with recording order — Fig 4 below checks exactly this.
%
%   Loads one or more evaluation GDFs of ANY paradigm (MI/CVSA/Hybrid mixed,
%   paradigm inferred from filename like main_session_overview). Uses ONLY
%   real GDF events (897/898/899 relative to 781) — no CSP/sLDA/artifact
%   simulation, no companion YAML required — so it is fast and works even
%   when the calibration model files are missing.
%
%   Each file's trial index is normalised to fractional session progress
%   (0 = first trial, 1 = last trial), so sessions of different length pool
%   correctly within the same paradigm.
%
%   Four figures:
%     Fig 1 — accuracy TREND across normalised session progress, one line
%             per paradigm (5 bins, mean +- binomial SE), overlaid in one
%             plot. Significance = point-biserial correlation between trial
%             position and HIT/not-HIT (every trial, not the binned means),
%             with a permutation p-value. This is the main, most readable
%             answer to "do trials get better or worse across the session?"
%     Fig 2 — first half vs second half of each session (paired per-file
%             dots + per-paradigm mean), sign-flip test on the paired delta.
%             A simple, file-level cross-check of Fig 1's trend.
%     Fig 3 — time-to-HIT vs normalised trial position (HIT trials only),
%             linear fit, Pearson r with a permutation p-value — same
%             question but on SPEED instead of accuracy.
%     Fig 4 — accuracy vs CHRONOLOGICAL recording order across ALL selected
%             files regardless of paradigm (timestamp parsed from the
%             bag_bci filename convention). Also prints the mean order
%             index per paradigm — if these are far apart, a paradigm
%             comparison done elsewhere may be confounded with recording
%             order/fatigue rather than reflecting a true paradigm effect.
%             Skipped if a timestamp can't be parsed from every filename.
%
%   Console: per-file half-split table; per-paradigm pooled statistics for
%   Figs 1-3; chronological order table + per-paradigm mean order index for
%   Fig 4.
%
%   Saves trial_dynamics_summary.mat (flat per-trial table: file, paradigm,
%   trial_idx, frac, outcome_code, dt; plus file_info: file-level chronological
%   table with timestamp and accuracy, used for Fig 4) under
%   <gdf_dir>/analysis_results/trial_dynamics/.

function main_trial_dynamics(gdf_dir, gdf_names, show_figures)
%   Callable as a function: main_trial_dynamics(gdf_dir, gdf_names, show_figures)
%   With no arguments, shows the GUI file picker (interactive mode).

if nargin < 3, show_figures = false; end
SHOW_FIGURES = show_figures;
if SHOW_FIGURES, fig_vis = 'on'; else, fig_vis = 'off'; end

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir,'io'), fullfile(this_dir,'utils'));

HIT_EV = 897; MISS_EV = 898; TO_EV = 899; CF_EV = 781;

if nargin < 1 || isempty(gdf_dir)
    default_dir = '/home/paolo/bci_vr_ws/recordings';
    [gdf_names, gdf_dir] = uigetfile({'*.gdf','GDF files (*.gdf)'}, ...
        'Select evaluation GDF file(s)', default_dir, 'MultiSelect', 'on');
    if isequal(gdf_names, 0), error('main_trial_dynamics:cancel', 'No file selected.'); end
elseif nargin < 2 || isempty(gdf_names)
    f = dir(fullfile(gdf_dir, '*.gdf'));
    gdf_names = {f.name};
    if isempty(gdf_names), error('main_trial_dynamics:nofiles', 'No GDF files in %s', gdf_dir); end
end
if ischar(gdf_names), gdf_names = {gdf_names}; end
n_files = numel(gdf_names);

COL = struct('mi', [0.85 0.30 0.10], 'cvsa', [0.10 0.60 0.30], 'hybrid', [0.18 0.45 0.75]);

% ── Flat per-trial table across all files ───────────────────────────────────
flat = struct('file_label', {}, 'paradigm', {}, 'trial_idx', {}, 'n_trials', {}, ...
              'frac', {}, 'outcome_code', {}, 'dt', {});

% ── File-level table: chronological order across ALL files, any paradigm ───
%   (catches a confound the within-file analysis below cannot: if one
%   paradigm is always recorded earlier/later in the day than the others,
%   any apparent paradigm effect could really be a time-of-day/fatigue effect)
file_info = struct('file_label', {}, 'paradigm', {}, 'timestamp', {}, 'acc', {});

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

    idx_fi = numel(file_info) + 1;
    file_info(idx_fi).file_label = basename;
    file_info(idx_fi).paradigm   = paradigm;
    file_info(idx_fi).timestamp  = parse_timestamp(basename);
    file_info(idx_fi).acc        = mean(outcomes == HIT_EV);
end
fprintf('═══════════════════════════════════════════════════════════════════\n');

pars = {'mi', 'cvsa', 'hybrid'};
par_labels = {'MI', 'CVSA', 'Hybrid'};
par_present = ismember(pars, unique({flat.paradigm}));
present_idx = find(par_present);

out_dir = fullfile(gdf_dir, 'analysis_results', 'trial_dynamics');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

n_perm = 2000;
N_BINS = 5;
bin_edges  = linspace(0, 1.0001, N_BINS+1);
bin_centers = (bin_edges(1:end-1) + min(bin_edges(2:end),1)) / 2;

% ═══════════════════════════════════════════════════════════════════════════
%  FIG 1 — accuracy TREND across normalised session progress (main figure)
% ═══════════════════════════════════════════════════════════════════════════
fig1 = figure('Name', 'Trial Dynamics — Accuracy trend across session progress', ...
              'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig1, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
ax1 = axes(fig1); hold(ax1, 'on'); %#ok<LAXES>

fprintf('\n══════════════════ Accuracy trend vs session progress (per paradigm) ══════════════════\n');
legend_h = gobjects(1, numel(present_idx));
for k = 1:numel(present_idx)
    pi = present_idx(k);
    mask_p = strcmp({flat.paradigm}, pars{pi});
    frac_p = [flat(mask_p).frac];
    is_hit = double([flat(mask_p).outcome_code] == HIT_EV);

    r_val = NaN; p_val = NaN;
    if numel(frac_p) >= 3
        C = corrcoef(frac_p, is_hit);
        r_val = C(1,2);
        p_val = corr_perm_test_local(frac_p, is_hit, r_val, n_perm);
    end

    bin_acc = nan(1, N_BINS); bin_se = nan(1, N_BINS); bin_n = zeros(1, N_BINS);
    for b = 1:N_BINS
        m = frac_p >= bin_edges(b) & frac_p < bin_edges(b+1);
        bin_n(b) = sum(m);
        if bin_n(b) > 0
            bin_acc(b) = mean(is_hit(m));
            bin_se(b)  = sqrt(bin_acc(b)*(1-bin_acc(b)) / bin_n(b));
        end
    end

    errorbar(ax1, bin_centers, 100*bin_acc, 100*bin_se, '-o', 'Color', COL.(pars{pi}), ...
             'MarkerFaceColor', COL.(pars{pi}), 'LineWidth', 2, 'MarkerSize', 7, 'CapSize', 6);
    legend_h(k) = plot(ax1, NaN, NaN, '-o', 'Color', COL.(pars{pi}), 'MarkerFaceColor', COL.(pars{pi}), ...
         'LineWidth', 2, 'DisplayName', sprintf('%s  (r=%+.2f, p=%.3f %s)', ...
         par_labels{pi}, r_val, p_val, stars_local(p_val)));

    fprintf('  %-8s  n_trials=%4d  r=%+.3f  p=%.4f  %s  (negative r = degrades, positive r = improves)\n', ...
            par_labels{pi}, numel(frac_p), r_val, p_val, stars_local(p_val));
end
yline(ax1, 50, 'k:', 'chance', 'HandleVisibility', 'off');
set(ax1, 'XLim', [0,1], 'YLim', [0,105]);
xlabel(ax1, 'normalised session progress (0=first trial, 1=last trial)');
ylabel(ax1, 'HIT rate (%)  —  mean \pm binomial SE per bin');
legend(ax1, legend_h, 'Location', 'best', 'FontSize', 9);
title(ax1, 'Accuracy trend across the session, per paradigm (5 bins)', 'FontWeight', 'bold');
grid(ax1, 'on');
fprintf('  (point-biserial r between trial position and HIT/not-HIT, permutation test, %d perms)\n', n_perm);
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

sgtitle(fig1, 'Within-session accuracy trend — positive r = learning, negative r = fatigue', 'Interpreter', 'none');
saveas(fig1, fullfile(out_dir, 'trial_dynamics_trend.svg'), 'svg');
if ~SHOW_FIGURES, close(fig1); end

% ═══════════════════════════════════════════════════════════════════════════
%  FIG 2 — first-half vs second-half accuracy, per file, per paradigm
%  (simple file-level cross-check of the Fig 1 trend)
% ═══════════════════════════════════════════════════════════════════════════
fig2 = figure('Name', 'Trial Dynamics — Within-session accuracy (1st vs 2nd half)', ...
              'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

fprintf('\n══════════════════ Half-split accuracy cross-check (per paradigm) ══════════════════\n');
n_par_present = numel(present_idx);
sp = 0;
for k = 1:n_par_present
    pi = present_idx(k);
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
sgtitle(fig2, 'Within-session accuracy: 1st half vs 2nd half (learning if 2nd > 1st, fatigue if 2nd < 1st)', ...
        'Interpreter', 'none');
saveas(fig2, fullfile(out_dir, 'trial_dynamics_half_split.svg'), 'svg');
if ~SHOW_FIGURES, close(fig2); end

% ═══════════════════════════════════════════════════════════════════════════
%  FIG 3 — TTH vs normalised trial position, per paradigm (speed dimension)
% ═══════════════════════════════════════════════════════════════════════════
fig3 = figure('Name', 'Trial Dynamics — TTH vs session progress', ...
              'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
set(fig3, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

fprintf('\n══════════════════ TTH vs session-progress trend (per paradigm) ══════════════════\n');
sp = 0;
for k = 1:n_par_present
    pi = present_idx(k);
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
sgtitle(fig3, 'Time-to-HIT vs normalised session progress (negative trend = speeding up = learning)', ...
        'Interpreter', 'none');
saveas(fig3, fullfile(out_dir, 'trial_dynamics_tth_trend.svg'), 'svg');
if ~SHOW_FIGURES, close(fig3); end

% ═══════════════════════════════════════════════════════════════════════════
%  FIG 4 — accuracy vs CHRONOLOGICAL recording order, across ALL files
%  (confound check: is one paradigm always recorded earlier/later in the day?)
% ═══════════════════════════════════════════════════════════════════════════
has_ts = ~any(isnan([file_info.timestamp]));
fprintf('\n══════════════════ Chronological order vs accuracy (confound check) ══════════════════\n');
if ~has_ts || numel(file_info) < 3
    fprintf('  Skipped: could not parse a YYYYMMDD.HHMMSS timestamp from every filename, or too few files.\n');
    fprintf('  Expected pattern: <subject>.<YYYYMMDD>.<HHMMSS>.<...>.gdf (bag_bci naming convention).\n');
else
    [~, ord] = sort([file_info.timestamp]);
    file_info = file_info(ord);
    n_fi = numel(file_info);
    order_idx = 1:n_fi;
    acc_v = [file_info.acc];

    C = corrcoef(order_idx, acc_v);
    r_ord = C(1,2);
    p_ord = corr_perm_test_local(order_idx, acc_v, r_ord, n_perm);

    fprintf('  %-4s  %-30s  %-7s  %8s\n', '#', 'file', 'paradigm', 'acc');
    for i = 1:n_fi
        fprintf('  %-4d  %-30s  %-7s  %7.0f%%\n', i, file_info(i).file_label, file_info(i).paradigm, 100*file_info(i).acc);
    end
    fprintf('  Recording-order trend: r=%+.3f  p=%.4f  %s  (negative r = later recordings score worse = fatigue across the whole block)\n', ...
            r_ord, p_ord, stars_local(p_ord));

    fprintf('\n  Mean recording-order index per paradigm (flags confounding if very different):\n');
    for pi = 1:3
        m = strcmp({file_info.paradigm}, pars{pi});
        if ~any(m), continue; end
        fprintf('    %-8s  mean order=%.1f / %d  (n=%d files)\n', par_labels{pi}, mean(order_idx(m)), n_fi, sum(m));
    end
    fprintf('  If these means are far apart, any MI/CVSA/Hybrid comparison elsewhere may be confounded with recording order.\n');

    fig4 = figure('Name', 'Trial Dynamics — Accuracy vs recording order', ...
                  'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
    set(fig4, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    ax4 = axes(fig4); hold(ax4, 'on'); %#ok<LAXES>
    mk = {'o','s','^'};
    for pi = 1:3
        m = strcmp({file_info.paradigm}, pars{pi});
        if ~any(m), continue; end
        scatter(ax4, order_idx(m), 100*acc_v(m), 80, COL.(pars{pi}), 'filled', ...
                'Marker', mk{pi}, 'MarkerEdgeColor', 'k', 'DisplayName', par_labels{pi});
    end
    p_fit = polyfit(order_idx, 100*acc_v, 1);
    plot(ax4, [1 n_fi], polyval(p_fit, [1 n_fi]), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    set(ax4, 'XTick', order_idx, 'YLim', [0,105]);
    xlabel(ax4, 'recording order (chronological, across all paradigms)');
    ylabel(ax4, 'file accuracy (%)');
    legend(ax4, 'Location', 'best', 'FontSize', 9);
    title(ax4, sprintf('Accuracy vs chronological recording order  (r=%+.2f, p=%.3f %s)', ...
          r_ord, p_ord, stars_local(p_ord)), 'FontWeight', 'bold');
    grid(ax4, 'on');
    sgtitle(fig4, 'Confound check: does performance drift across the recording session, regardless of paradigm?', ...
            'Interpreter', 'none');
    saveas(fig4, fullfile(out_dir, 'trial_dynamics_recording_order.svg'), 'svg');
    if ~SHOW_FIGURES, close(fig4); end
end
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

% ── Save flat per-trial table + file-level chronological table ─────────────
trial_dynamics_summary = flat;
save(fullfile(out_dir, 'trial_dynamics_summary.mat'), 'trial_dynamics_summary', 'file_info');
fprintf('\nSaved figures and trial_dynamics_summary.mat to %s\n', out_dir);


end % main_trial_dynamics

% ── Local helpers ────────────────────────────────────────────────────────────

function ts = parse_timestamp(fname)
% PARSE_TIMESTAMP  Numeric YYYYMMDDHHMMSS from a "<subject>.<YYYYMMDD>.<HHMMSS>.<...>"
%   filename (bag_bci naming convention). NaN if the pattern isn't found.
    parts = strsplit(fname, '.');
    ts = NaN;
    for i = 1:numel(parts)-1
        if numel(parts{i}) == 8 && all(isstrprop(parts{i}, 'digit')) && ...
           numel(parts{i+1}) == 6 && all(isstrprop(parts{i+1}, 'digit'))
            ts = str2double([parts{i}, parts{i+1}]);
            return;
        end
    end
end

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
