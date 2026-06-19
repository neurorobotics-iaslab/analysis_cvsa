%% MAIN_VALIDATE_COUNTERFACTUAL  Real between-session accuracy vs offline simulation.
%
%   Loads two .mat files:
%     1. session_summary.mat      — from main_session_overview (3 paradigm sessions)
%     2. counterfactual_summary.mat — from main_hybrid_advantage_integ (hybrid data)
%
%   Compares, per paradigm:
%     real MI session acc  vs  simulated MI-only (from hybrid data)
%     real CVSA session acc vs simulated CVSA-only (from hybrid data)
%     real Hybrid session acc vs simulated Hybrid  [sanity check]
%
%   Agreement near the identity line validates the counterfactual: the offline
%   estimate "what would have happened in MI-only mode" is representative of an
%   actual dedicated MI session. Large divergence → dual-task cognitive cost.
%
%   TTH comparison follows the same logic.
%
%   Stream index in counterfactual_summary: 1=Hybrid, 2=MI-only, 3=CVSA-only.
%
%   Saves:
%     counterfactual_validation.svg  — bar chart + scatter + TTH comparison

clear; clc; close all;

SHOW_FIGURES = true;
if SHOW_FIGURES, fig_vis = 'on'; else, fig_vis = 'off'; end

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir,'io'), fullfile(this_dir,'utils'));

% ── Load session summary ────────────────────────────────────────────────────
[f1, d1] = uigetfile({'*.mat','MAT files (*.mat)'}, ...
    'Select session_summary.mat  (from main_session_overview)', ...
    '/home/paolo/bci_vr_ws/recordings');
if isequal(f1, 0), error('No file selected.'); end
s1 = load(fullfile(d1, f1));
summary_file = s1.summary_file;

% ── Load counterfactual summary ─────────────────────────────────────────────
[f2, d2] = uigetfile({'*.mat','MAT files (*.mat)'}, ...
    'Select counterfactual_summary.mat  (from main_hybrid_advantage_integ)', ...
    '/home/paolo/bci_vr_ws/recordings');
if isequal(f2, 0), error('No file selected.'); end
s2 = load(fullfile(d2, f2));
cf = s2.counterfactual;

% ── Aggregate real performance per paradigm ─────────────────────────────────
pars       = {'mi',  'cvsa', 'hybrid'};
par_labels = {'MI',  'CVSA', 'Hybrid'};
real_acc   = nan(1, 3);
real_tth   = nan(1, 3);
real_n     = zeros(1, 3);

for pi = 1:3
    mask = strcmp({summary_file.paradigm}, pars{pi});
    if ~any(mask), continue; end
    nh = sum([summary_file(mask).n_hit]);
    nt = sum([summary_file(mask).n_hit] + [summary_file(mask).n_miss] + [summary_file(mask).n_to]);
    real_acc(pi) = nh / max(1, nt);
    real_n(pi)   = nt;
    tth_all = [];
    for fi = find(mask)
        tth_all = [tth_all, summary_file(fi).tth_vals]; %#ok<AGROW>
    end
    if ~isempty(tth_all), real_tth(pi) = mean(tth_all); end
end

% ── Map counterfactual streams to paradigm order [MI, CVSA, Hybrid] ─────────
%   cf.stream_names = {'Hybrid','MI-only','CVSA-only'} → indices 1,2,3
sim_acc     = [cf.acc(2),       cf.acc(3),       cf.acc(1)];
sim_ci_lo   = [cf.ci_lo_acc(2), cf.ci_lo_acc(3), cf.ci_lo_acc(1)];
sim_ci_hi   = [cf.ci_hi_acc(2), cf.ci_hi_acc(3), cf.ci_hi_acc(1)];
sim_tth     = [cf.t_hit_mean(2), cf.t_hit_mean(3), cf.t_hit_mean(1)];
sim_tth_sem = [cf.t_hit_sem(2),  cf.t_hit_sem(3),  cf.t_hit_sem(1)];

% ── Console summary ─────────────────────────────────────────────────────────
fprintf('\n══════════════════ Counterfactual validation ══════════════════\n');
fprintf('  %-8s  %9s  %9s  %8s  %8s  %8s\n', ...
        'paradigm', 'real acc', 'sim acc', 'delta', 'n_real', 'n_sim');
for pi = 1:3
    r_str = '--';  s_str = '--';  d_str = '--';
    if ~isnan(real_acc(pi)), r_str = sprintf('%.1f%%', 100*real_acc(pi)); end
    if ~isnan(sim_acc(pi)),  s_str = sprintf('%.1f%%', 100*sim_acc(pi)); end
    if ~isnan(real_acc(pi)) && ~isnan(sim_acc(pi))
        d_str = sprintf('%+.1f%%', 100*(real_acc(pi)-sim_acc(pi)));
    end
    fprintf('  %-8s  %9s  %9s  %8s  %8d  %8d\n', ...
            par_labels{pi}, r_str, s_str, d_str, real_n(pi), cf.n_trials);
end
fprintf('\n  TTH comparison (HIT trials):\n');
fprintf('  %-8s  %9s  %9s\n', 'paradigm', 'real TTH', 'sim TTH');
for pi = 1:3
    r_str = '--';  s_str = '--';
    if ~isnan(real_tth(pi)), r_str = sprintf('%.2fs', real_tth(pi)); end
    if ~isnan(sim_tth(pi)),  s_str = sprintf('%.2fs', sim_tth(pi)); end
    fprintf('  %-8s  %9s  %9s\n', par_labels{pi}, r_str, s_str);
end
fprintf('\n  Interpretation:\n');
fprintf('  |delta| < 10%% → simulation is a valid proxy for a real dedicated session.\n');
fprintf('  |delta| > 10%% → dual-task condition during hybrid session degraded the signal.\n');
fprintf('═══════════════════════════════════════════════════════════════\n\n');

% ── Colors ──────────────────────────────────────────────────────────────────
COL_REAL = {[0.85 0.30 0.10], [0.10 0.60 0.30], [0.18 0.45 0.75]};
COL_SIM  = {[0.93 0.70 0.60], [0.65 0.88 0.75], [0.65 0.78 0.93]};

analysis_results_dir = fileparts(fileparts(d2));   % .../analysis_results/hybrid_advantage_integ/ -> .../analysis_results/
out_dir = fullfile(analysis_results_dir, 'validate_counterfactual');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% ── Figure: 3-panel layout ──────────────────────────────────────────────────
%   Panel 1 — grouped bars: real (solid) vs simulated (lighter) per paradigm
%   Panel 2 — scatter real vs simulated (3 points) + identity line
%   Panel 3 — TTH comparison (same layout as panel 1)

fig = figure('Name', 'Counterfactual Validation', 'Color', 'w', ...
             'NumberTitle', 'off', 'Visible', fig_vis);
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

bw = 0.35;
x_pos = [1 2 3];

% ── Panel 1: accuracy ───────────────────────────────────────────────────────
ax1 = subplot(1, 3, 1); hold(ax1, 'on');
for pi = 1:3
    if ~isnan(real_acc(pi))
        bar(ax1, x_pos(pi)-bw/2, 100*real_acc(pi), bw, ...
            'FaceColor', COL_REAL{pi}, 'EdgeColor', 'k', 'LineWidth', 1.2);
        text(ax1, x_pos(pi)-bw/2, 100*real_acc(pi)+2, ...
             sprintf('%.0f%%', 100*real_acc(pi)), ...
             'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    end
    if ~isnan(sim_acc(pi))
        bar(ax1, x_pos(pi)+bw/2, 100*sim_acc(pi), bw, ...
            'FaceColor', COL_SIM{pi}, 'EdgeColor', 'k', 'LineWidth', 1.2, ...
            'LineStyle', '--');
        errorbar(ax1, x_pos(pi)+bw/2, 100*sim_acc(pi), ...
                 100*(sim_acc(pi)-sim_ci_lo(pi)), 100*(sim_ci_hi(pi)-sim_acc(pi)), ...
                 'k', 'LineWidth', 1.2, 'CapSize', 6);
        text(ax1, x_pos(pi)+bw/2, 100*sim_ci_hi(pi)+2.5, ...
             sprintf('%.0f%%', 100*sim_acc(pi)), ...
             'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    if ~isnan(real_acc(pi)) && ~isnan(sim_acc(pi))
        delta = 100*(real_acc(pi) - sim_acc(pi));
        text(ax1, x_pos(pi), max(100*real_acc(pi), 100*sim_ci_hi(pi))+8, ...
             sprintf('\\Delta%+.0f%%', delta), ...
             'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.35 0.35 0.35]);
    end
end
b_r = bar(ax1, NaN, NaN, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k', ...
          'DisplayName', 'Real session');
b_s = bar(ax1, NaN, NaN, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'k', ...
          'LineStyle', '--', 'DisplayName', 'Simulated (from hybrid data)');
legend(ax1, [b_r, b_s], 'Location', 'south', 'FontSize', 8);
yline(ax1, 50, '--k', 'FontSize', 7, 'HandleVisibility', 'off');
set(ax1, 'XTick', x_pos, 'XTickLabel', par_labels, 'YLim', [0, 118]);
ylabel(ax1, 'HIT rate (%)');
title(ax1, {'Accuracy', 'solid = real session;  dashed = simulated'}, ...
      'FontWeight', 'bold');
grid(ax1, 'on');

% ── Panel 2: scatter real vs simulated (identity line = perfect agreement) ──
ax2 = subplot(1, 3, 2); hold(ax2, 'on');
valid = ~isnan(real_acc) & ~isnan(sim_acc);
if any(valid)
    all_v = [100*real_acc(valid), 100*sim_acc(valid)];
    lims  = [max(0, min(all_v)-10), min(108, max(all_v)+10)];

    fill(ax2, [lims(1) lims(2) lims(2) lims(1)], ...
              [lims(1)-10 lims(2)-10 lims(2)+10 lims(1)+10], ...
         [0.88 0.88 0.88], 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
    plot(ax2, lims, lims, 'k--', 'LineWidth', 1.5, 'DisplayName', 'y = x  (perfect agreement)');
    text(ax2, lims(2)-3, lims(2)+2, '±10%', 'FontSize', 7, 'Color', [0.5 0.5 0.5]);

    for pi = 1:3
        if ~valid(pi), continue; end
        scatter(ax2, 100*sim_acc(pi), 100*real_acc(pi), 140, COL_REAL{pi}, ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        text(ax2, 100*sim_acc(pi)+1.5, 100*real_acc(pi), par_labels{pi}, ...
             'FontSize', 10, 'FontWeight', 'bold', 'Color', COL_REAL{pi});
    end

    legend(ax2, 'Location', 'southeast', 'FontSize', 8);
    set(ax2, 'XLim', lims, 'YLim', lims); axis(ax2, 'equal');
end
xlabel(ax2, 'Simulated accuracy (%)  [from hybrid data]');
ylabel(ax2, 'Real session accuracy (%)');
title(ax2, {'Simulated vs Real — accuracy', 'points near diagonal: simulation valid'}, ...
      'FontWeight', 'bold');
grid(ax2, 'on'); box(ax2, 'on');

% ── Panel 3: TTH comparison ─────────────────────────────────────────────────
ax3 = subplot(1, 3, 3); hold(ax3, 'on');
any_tth = false;
for pi = 1:3
    if ~isnan(real_tth(pi))
        bar(ax3, x_pos(pi)-bw/2, real_tth(pi), bw, ...
            'FaceColor', COL_REAL{pi}, 'EdgeColor', 'k', 'LineWidth', 1.2);
        text(ax3, x_pos(pi)-bw/2, real_tth(pi)+0.05, sprintf('%.1fs', real_tth(pi)), ...
             'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
        any_tth = true;
    end
    if ~isnan(sim_tth(pi))
        bar(ax3, x_pos(pi)+bw/2, sim_tth(pi), bw, ...
            'FaceColor', COL_SIM{pi}, 'EdgeColor', 'k', 'LineWidth', 1.2, ...
            'LineStyle', '--');
        if ~isnan(sim_tth_sem(pi)) && sim_tth_sem(pi) > 0
            errorbar(ax3, x_pos(pi)+bw/2, sim_tth(pi), sim_tth_sem(pi), ...
                     'k', 'LineWidth', 1.2, 'CapSize', 6);
        end
        text(ax3, x_pos(pi)+bw/2, sim_tth(pi)+0.05, sprintf('%.1fs', sim_tth(pi)), ...
             'HorizontalAlignment', 'center', 'FontSize', 9);
        any_tth = true;
    end
end
if ~any_tth
    text(ax3, 2, 0.5, 'No HIT trials available', 'HorizontalAlignment', 'center');
end
set(ax3, 'XTick', x_pos, 'XTickLabel', par_labels);
ylabel(ax3, 'Mean TTH (s)');
title(ax3, {'Time to HIT', 'solid = real;  dashed = simulated'}, 'FontWeight', 'bold');
grid(ax3, 'on');

% ── Sgtitle with per-paradigm delta summary ──────────────────────────────────
delta_str = '';
for pi = 1:3
    if ~isnan(real_acc(pi)) && ~isnan(sim_acc(pi))
        delta_str = [delta_str, sprintf('%s: %+.1f%%  ', par_labels{pi}, ...
                     100*(real_acc(pi)-sim_acc(pi)))]; %#ok<AGROW>
    end
end
sgtitle(fig, sprintf(['Counterfactual Validation  |  n_sim=%d trials  |  ' ...
        'Accuracy delta (real − simulated):  %s\n' ...
        '|Δ|<10%%: simulation valid   |Δ|>10%%: dual-task cost'], ...
        cf.n_trials, strtrim(delta_str)), 'Interpreter', 'none', 'FontSize', 11);

saveas(fig, fullfile(out_dir, 'counterfactual_validation.svg'), 'svg');
fprintf('Saved counterfactual_validation.svg to %s\n', out_dir);
if ~SHOW_FIGURES, close(fig); end
