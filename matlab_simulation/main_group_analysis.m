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
%
%   Subject ID is inferred from the folder structure: the first path
%   component directly under the selected root folder.
%
%   Multiple .mat files for the same subject (e.g. several evaluation
%   sessions) are pooled: hit/miss/to counts are summed, accuracy is
%   recomputed from the pooled counts; counterfactual accuracy is pooled
%   weighted by n_trials; fusion-mechanism and ERD/ERS metrics are averaged
%   (frame counts summed).
%
%   Four figures:
%     Fig 1 — per-subject real-session accuracy per paradigm (MI/CVSA/Hybrid),
%             grouped bars + group mean +- SEM dashed lines.
%     Fig 2 — per-subject counterfactual accuracy (Hybrid/MI-only/CVSA-only)
%             and Hybrid-MI / Hybrid-CVSA accuracy deltas, with a group-level
%             sign-flip permutation test across SUBJECTS (not trials) — the
%             statistically meaningful unit once more than one subject exists.
%     Fig 3 — per-subject CVSA-fusion advantage (mean P_fused-P_MI) and
%             rescue-vs-cost frame-level effect, with the same group-level
%             sign-flip test across subjects — tests whether the fusion
%             mechanism described per-file by main_hybrid_advantage_probs
%             is consistent across the whole cohort, not just one subject.
%     Fig 4 — per-subject ERD/ERS class discrimination and CSP-weight
%             correlation (Pearson r), per task, with a group-level
%             sign-flip test on r across subjects — tests whether the
%             neurophysiological grounding shown per-file by topo_erders.m
%             holds across the whole cohort.
%
%   Console: per-subject table + group-level paired statistics.
%
%   Saves group_summary.mat (per-subject aggregated table) and
%   group_analysis.svg under <root>/group_analysis/.

clear; clc; close all;

SHOW_FIGURES = true;
if SHOW_FIGURES, fig_vis = 'on'; else, fig_vis = 'off'; end

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir,'utils'));

root_dir = uigetdir('/home/paolo/bci_vr_ws/recordings', 'Select root recordings folder');
if isequal(root_dir, 0), error('No folder selected.'); end

session_files = dir(fullfile(root_dir, '**', 'session_overview', 'session_summary.mat'));
cf_files      = dir(fullfile(root_dir, '**', 'hybrid_advantage_integ', 'counterfactual_summary.mat'));
probs_files   = dir(fullfile(root_dir, '**', 'hybrid_advantage_probs', 'hybrid_advantage_probs_summary.mat'));
erd_files     = dir(fullfile(root_dir, '**', 'topo_erders', 'topo_erders_summary.mat'));

fprintf('Found %d session_summary.mat, %d counterfactual_summary.mat, %d hybrid_advantage_probs_summary.mat, %d topo_erders_summary.mat under %s\n', ...
        numel(session_files), numel(cf_files), numel(probs_files), numel(erd_files), root_dir);
if isempty(session_files) && isempty(cf_files) && isempty(probs_files) && isempty(erd_files)
    error('main_group_analysis:nodata', 'No .mat summaries found. Run main_session_overview / main_hybrid_advantage_integ / main_hybrid_advantage_probs / topo_erders first.');
end

% ── Flatten real-session per-file records, tagged with subject ─────────────
flat_real = struct('subject', {}, 'paradigm', {}, 'n_hit', {}, 'n_miss', {}, 'n_to', {}, 'tth_vals', {});
for i = 1:numel(session_files)
    fpath = fullfile(session_files(i).folder, session_files(i).name);
    subj  = subject_of(fpath, root_dir);
    s = load(fpath); sf = s.summary_file;
    for k = 1:numel(sf)
        idx = numel(flat_real) + 1;
        flat_real(idx).subject  = subj;
        flat_real(idx).paradigm = sf(k).paradigm;
        flat_real(idx).n_hit    = sf(k).n_hit;
        flat_real(idx).n_miss   = sf(k).n_miss;
        flat_real(idx).n_to     = sf(k).n_to;
        flat_real(idx).tth_vals = sf(k).tth_vals;
    end
end

% ── Flatten counterfactual per-file records, tagged with subject ───────────
flat_cf = struct('subject', {}, 'acc', {}, 'n_trials', {}, 't_hit_mean', {});
for i = 1:numel(cf_files)
    fpath = fullfile(cf_files(i).folder, cf_files(i).name);
    subj  = subject_of(fpath, root_dir);
    s = load(fpath); cf = s.counterfactual;
    idx = numel(flat_cf) + 1;
    flat_cf(idx).subject    = subj;
    flat_cf(idx).acc        = cf.acc;          % [1x3] Hybrid, MI-only, CVSA-only
    flat_cf(idx).n_trials   = cf.n_trials;
    flat_cf(idx).t_hit_mean = cf.t_hit_mean;   % [1x3]
end

% ── Flatten hybrid_advantage_probs per-file records, tagged with subject ───
flat_probs = struct('subject', {}, 'fus_adv_mean', {}, 'n_rescued_fr', {}, 'n_hurt_fr', {}, ...
                     'rescue_delta', {}, 'cost_delta', {});
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

subjects = unique([{flat_real.subject}, {flat_cf.subject}, {flat_probs.subject}, {flat_erd.subject}]);
n_subj   = numel(subjects);
pars     = {'mi', 'cvsa', 'hybrid'};
par_labels = {'MI', 'CVSA', 'Hybrid'};

% ── Per-subject real accuracy/TTH per paradigm (pooled across files) ───────
real_acc = nan(n_subj, 3);
real_n   = zeros(n_subj, 3);
real_tth = nan(n_subj, 3);
for si = 1:n_subj
    for pi = 1:3
        mask = strcmp({flat_real.subject}, subjects{si}) & strcmp({flat_real.paradigm}, pars{pi});
        if ~any(mask), continue; end
        nh = sum([flat_real(mask).n_hit]);
        nt = sum([flat_real(mask).n_hit] + [flat_real(mask).n_miss] + [flat_real(mask).n_to]);
        real_acc(si,pi) = nh / max(1, nt);
        real_n(si,pi)   = nt;
        tth_all = [flat_real(mask).tth_vals];
        if ~isempty(tth_all), real_tth(si,pi) = mean(tth_all); end
    end
end

% ── Per-subject counterfactual accuracy (pooled, weighted by n_trials) ─────
cf_acc = nan(n_subj, 3);   % Hybrid, MI-only, CVSA-only
cf_n   = zeros(n_subj, 1);
for si = 1:n_subj
    mask = strcmp({flat_cf.subject}, subjects{si});
    if ~any(mask), continue; end
    rows = flat_cf(mask);
    ntot = sum([rows.n_trials]);
    acc_w = zeros(1,3);
    for r = 1:numel(rows)
        acc_w = acc_w + rows(r).acc * rows(r).n_trials;
    end
    cf_acc(si,:) = acc_w / max(1, ntot);
    cf_n(si) = ntot;
end

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

% ── Console: per-subject table ──────────────────────────────────────────────
fprintf('\n══════════════════ Per-subject real-session accuracy ══════════════════\n');
fprintf('  %-10s  %8s  %8s  %8s   %8s  %8s  %8s\n', 'subject', 'MI%', 'CVSA%', 'Hybrid%', 'n_MI', 'n_CVSA', 'n_Hyb');
for si = 1:n_subj
    fprintf('  %-10s  %8s  %8s  %8s   %8d  %8d  %8d\n', subjects{si}, ...
            fmt_pct(real_acc(si,1)), fmt_pct(real_acc(si,2)), fmt_pct(real_acc(si,3)), ...
            real_n(si,1), real_n(si,2), real_n(si,3));
end

fprintf('\n══════════════════ Per-subject counterfactual accuracy ══════════════════\n');
fprintf('  %-10s  %8s  %8s  %8s   %8s  %8s\n', 'subject', 'Hybrid%', 'MI-only%', 'CVSA-only%', 'n_trials', 'Hyb-MI');
for si = 1:n_subj
    if isnan(cf_acc(si,1)), continue; end
    fprintf('  %-10s  %8s  %8s  %8s   %8d  %+7.1f%%\n', subjects{si}, ...
            fmt_pct(cf_acc(si,1)), fmt_pct(cf_acc(si,2)), fmt_pct(cf_acc(si,3)), ...
            cf_n(si), 100*(cf_acc(si,1)-cf_acc(si,2)));
end

% ── Group-level paired statistics across SUBJECTS ───────────────────────────
d_mi_subj  = cf_acc(:,1) - cf_acc(:,2);
d_cvs_subj = cf_acc(:,1) - cf_acc(:,3);
n_perm = 2000;
p_mi_grp  = sign_flip_test_local(d_mi_subj,  n_perm);
p_cvs_grp = sign_flip_test_local(d_cvs_subj, n_perm);
fprintf('\n══════════════════ Group-level counterfactual advantage (n=%d subjects) ══════════════════\n', sum(~isnan(d_mi_subj)));
fprintf('  Hybrid - MI-only   : mean delta=%+.1f%%  p=%.4f  %s\n', ...
        100*mean(d_mi_subj,'omitnan'),  p_mi_grp,  stars(p_mi_grp));
fprintf('  Hybrid - CVSA-only : mean delta=%+.1f%%  p=%.4f  %s\n', ...
        100*mean(d_cvs_subj,'omitnan'), p_cvs_grp, stars(p_cvs_grp));
fprintf('  (two-sided sign-flip permutation, %d perms; subject is the statistical unit)\n', n_perm);
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

% ── Group-level CVSA-fusion mechanism across SUBJECTS ───────────────────────
p_fusadv_grp = sign_flip_test_local(fus_adv_subj, n_perm);
d_resc_cost  = rescue_subj - cost_subj;
p_resc_grp   = sign_flip_test_local(d_resc_cost, n_perm);
fprintf('\n══════════════════ Group-level CVSA-fusion mechanism (n=%d subjects) ══════════════════\n', sum(~isnan(fus_adv_subj)));
fprintf('  %-10s  %10s  %10s  %10s  %12s\n', 'subject', 'fus_adv', 'rescue_d', 'cost_d', 'net_resc(fr)');
for si = 1:n_subj
    if isnan(fus_adv_subj(si)), continue; end
    fprintf('  %-10s  %+9.3f  %+9.3f  %+9.3f  %+11d\n', subjects{si}, ...
            fus_adv_subj(si), rescue_subj(si), cost_subj(si), net_rescue_subj(si));
end
fprintf('  mean fusion advantage (P_fused-P_MI) : %+.3f  p=%.4f  %s\n', ...
        mean(fus_adv_subj, 'omitnan'), p_fusadv_grp, stars(p_fusadv_grp));
fprintf('  mean (rescue - cost) delta           : %+.3f  p=%.4f  %s\n', ...
        mean(d_resc_cost, 'omitnan'), p_resc_grp, stars(p_resc_grp));
fprintf('  (two-sided sign-flip permutation, %d perms; subject is the statistical unit)\n', n_perm);
fprintf('═══════════════════════════════════════════════════════════════════════════════════════\n');

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

% ── Colors (consistent with main_session_overview / main_hybrid_advantage_integ) ─
COL_mi     = [0.85 0.30 0.10];
COL_cvsa   = [0.10 0.60 0.30];
COL_hybrid = [0.18 0.45 0.75];
COL_mi_only   = [0.85 0.40 0.10];
COL_cvsa_only = [0.49 0.18 0.56];

out_dir = fullfile(root_dir, 'group_analysis');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% ── Figure 1: per-subject real accuracy per paradigm ───────────────────────
fig1 = figure('Name', 'Group Analysis — Real Session Accuracy', 'Color', 'w', ...
              'NumberTitle', 'off', 'Visible', fig_vis);
set(fig1, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
ax1 = axes(fig1); hold(ax1, 'on'); %#ok<LAXES>

bw = 0.25;
cols3 = {COL_mi, COL_cvsa, COL_hybrid};
for pi = 1:3
    x = (1:n_subj) + (pi-2)*bw;
    v = 100*real_acc(:,pi);
    bar(ax1, x, v, bw*0.9, 'FaceColor', cols3{pi}, 'EdgeColor', 'k', 'LineWidth', 1, ...
        'DisplayName', par_labels{pi});
    m = mean(v, 'omitnan');
    if ~isnan(m)
        plot(ax1, [0.5, n_subj+0.5], [m m], '--', 'Color', cols3{pi}*0.6, ...
             'LineWidth', 1.5, 'HandleVisibility', 'off');
        text(ax1, n_subj+0.55, m, sprintf('\\mu=%.0f%%', m), 'FontSize', 8, ...
             'Color', cols3{pi}*0.6, 'FontWeight', 'bold');
    end
end
yline(ax1, 50, 'k:', 'HandleVisibility', 'off');
set(ax1, 'XTick', 1:n_subj, 'XTickLabel', subjects, 'YLim', [0, 110], 'XLim', [0.4, n_subj+0.9]);
ylabel(ax1, 'HIT rate (%)');
legend(ax1, 'Location', 'south', 'FontSize', 9);
title(ax1, sprintf('Real-session accuracy per subject (n=%d subjects)', n_subj), 'FontWeight', 'bold');
grid(ax1, 'on');
sgtitle(fig1, 'Group Analysis — Real Session Accuracy', 'Interpreter', 'none');

saveas(fig1, fullfile(out_dir, 'group_real_accuracy.svg'), 'svg');
if ~SHOW_FIGURES, close(fig1); end

% ── Figure 2: counterfactual accuracy + group-level advantage ─────────────
fig2 = figure('Name', 'Group Analysis — Counterfactual Advantage', 'Color', 'w', ...
              'NumberTitle', 'off', 'Visible', fig_vis);
set(fig2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

ax2 = subplot(1, 2, 1); hold(ax2, 'on');
cols_cf = {COL_hybrid, COL_mi_only, COL_cvsa_only};
cf_labels = {'Hybrid', 'MI-only', 'CVSA-only'};
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
title(ax3, sprintf('Per-subject advantage  |  group: Hyb-MI=%+.1f%% (%s), Hyb-CVSA=%+.1f%% (%s)', ...
      100*m_mi, stars(p_mi_grp), 100*m_cvs, stars(p_cvs_grp)), 'FontWeight', 'bold', 'FontSize', 9);
grid(ax3, 'on');

sgtitle(fig2, sprintf('Group Analysis — Counterfactual Advantage (n=%d subjects)', n_subj), 'Interpreter', 'none');

saveas(fig2, fullfile(out_dir, 'group_counterfactual_advantage.svg'), 'svg');
if ~SHOW_FIGURES, close(fig2); end

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
title(ax4, sprintf('Per-subject CVSA-fusion advantage  |  group mean=%+.3f (%s)', m_fa, stars(p_fusadv_grp)), ...
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
title(ax5, sprintf('Rescue vs cost per subject  |  group (rescue-cost)=%+.3f (%s)', ...
      mean(d_resc_cost,'omitnan'), stars(p_resc_grp)), 'FontWeight', 'bold', 'FontSize', 9);
grid(ax5, 'on');

sgtitle(fig3, sprintf('Group Analysis — CVSA-Fusion Mechanism (n=%d subjects)', n_subj), 'Interpreter', 'none');

saveas(fig3, fullfile(out_dir, 'group_fusion_mechanism.svg'), 'svg');
if ~SHOW_FIGURES, close(fig3); end

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

% ── Save group summary .mat ────────────────────────────────────────────────
group_summary = struct();
group_summary.subjects = subjects;
group_summary.real_acc = real_acc;   % [n_subj x 3] MI, CVSA, Hybrid
group_summary.real_n   = real_n;
group_summary.real_tth = real_tth;
group_summary.cf_acc   = cf_acc;     % [n_subj x 3] Hybrid, MI-only, CVSA-only
group_summary.cf_n     = cf_n;
group_summary.p_mi_group  = p_mi_grp;
group_summary.p_cvsa_group = p_cvs_grp;
group_summary.fus_adv_subj    = fus_adv_subj;
group_summary.rescue_subj     = rescue_subj;
group_summary.cost_subj       = cost_subj;
group_summary.net_rescue_subj = net_rescue_subj;
group_summary.p_fusadv_group  = p_fusadv_grp;
group_summary.p_rescue_group  = p_resc_grp;
group_summary.erd_discrim_subj = erd_discrim_subj;   % [n_subj x 3] MI, CVSA, Hybrid
group_summary.erd_r_subj       = erd_r_subj;
group_summary.p_erd_r_group    = p_r_grp;            % [1x3]
save(fullfile(out_dir, 'group_summary.mat'), 'group_summary');
fprintf('\nSaved group_summary.mat and figures to %s\n', out_dir);


% ── Local helpers ────────────────────────────────────────────────────────────

function subj = subject_of(fpath, root_dir)
    rel = erase(fpath, [char(root_dir) filesep]);
    parts = strsplit(rel, filesep);
    subj = parts{1};
end

function s = fmt_pct(v)
    if isnan(v), s = '--'; else, s = sprintf('%.1f%%', 100*v); end
end

function s = stars(p)
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
