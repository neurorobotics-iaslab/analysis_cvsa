function well_vs_bad_comparison(gs_well, gs_bad, out_dir, show_figures)
% WELL_VS_BAD_COMPARISON  Direct, UNPAIRED comparison of two manually-
%   assigned subject-performance groups ("well" vs "bad" -- different
%   subjects in each, NOT the paired sign-flip tests main_group_analysis.m
%   uses everywhere else for the SAME subjects across streams/paradigms).
%
%   gs_well / gs_bad are the group_summary structs returned by
%   main_group_analysis(root_dir, well_subjects, show_figures, 'well') and
%   (..., bad_subjects, show_figures, 'bad') respectively -- see
%   batch_group_analysis.m, which is the intended caller: it defines the
%   WELL_SUBJECTS/BAD_SUBJECTS cellstr lists, runs main_group_analysis once
%   per group, and passes the two returned summaries here.
%
%   Answers: does the fusion advantage / classifier quality genuinely differ
%   between strong and weak performers -- e.g. does Hybrid fusion compensate
%   weak performers (bigger advantage in "bad") or only pay off for subjects
%   who are already good (bigger advantage in "well")? Tests real Hybrid
%   decided-trial accuracy, counterfactual Hybrid-MI/Hybrid-CVSA advantage,
%   fusion advantage (mean P_fused-P_MI), ROC AUC (Hybrid, MI-only), and the
%   raw signal-quality flagged-trial rate (dead/noisy channel during CF,
%   main_session_overview.m Fig 11) -- the last one is a confound check: are
%   "bad" performers actually just recorded with worse electrode contact?
%   via a two-sample permutation test (pooled-label shuffle) + unpaired
%   Cohen's d -- both hand-rolled, no Statistics Toolbox required.
%
%   Saves 16_well_vs_bad_comparison.svg under out_dir (no group tag: it
%   spans both groups, unlike every other figure in this package).

if show_figures, fig_vis = 'on'; else, fig_vis = 'off'; end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Backward-compatible: gs_well/gs_bad may have been computed by an older
% main_group_analysis.m that didn't save sigqc_rate_overall yet -- render
% that one panel with NaN (skipped stats) rather than erroring.
if isfield(gs_well, 'sigqc_rate_overall') && isfield(gs_bad, 'sigqc_rate_overall')
    sq_well = 100 * gs_well.sigqc_rate_overall;
    sq_bad  = 100 * gs_bad.sigqc_rate_overall;
else
    sq_well = nan(numel(gs_well.subjects), 1);
    sq_bad  = nan(numel(gs_bad.subjects), 1);
end

metrics = {
    'Real Hybrid acc (decided)',       gs_well.real_acc_dec(:,3),                        gs_bad.real_acc_dec(:,3);
    'CF Hybrid-MI advantage',          gs_well.cf_acc_dec(:,1) - gs_well.cf_acc_dec(:,2), gs_bad.cf_acc_dec(:,1) - gs_bad.cf_acc_dec(:,2);
    'CF Hybrid-CVSA advantage',        gs_well.cf_acc_dec(:,1) - gs_well.cf_acc_dec(:,3), gs_bad.cf_acc_dec(:,1) - gs_bad.cf_acc_dec(:,3);
    'Fusion advantage (P_fused-P_MI)', gs_well.fus_adv_subj,                              gs_bad.fus_adv_subj;
    'ROC AUC (Hybrid)',                gs_well.roc_auc_subj(:,3),                         gs_bad.roc_auc_subj(:,3);
    'ROC AUC (MI-only)',               gs_well.roc_auc_subj(:,1),                         gs_bad.roc_auc_subj(:,1);
    'Signal-quality flagged-trial rate (%)', sq_well,                                     sq_bad;
};

fig = figure('Name', 'Group Analysis — Well vs Bad Performers', 'Color', 'w', ...
             'NumberTitle', 'off', 'Visible', fig_vis);
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

fprintf('\n── Well (n=%d: %s) vs Bad (n=%d: %s) performer comparison (unpaired) ──\n', ...
    numel(gs_well.subjects), strjoin(gs_well.subjects, ', '), ...
    numel(gs_bad.subjects), strjoin(gs_bad.subjects, ', '));

for mi = 1:size(metrics, 1)
    ax = subplot(3, 3, mi); hold(ax, 'on');
    wv = metrics{mi, 2}; bv = metrics{mi, 3};
    p = two_sample_perm_test_local(wv, bv);
    d = cohen_d_unpaired_local(wv, bv);

    jitter_w = (rand(numel(wv),1) - 0.5) * 0.15;
    jitter_b = (rand(numel(bv),1) - 0.5) * 0.15;
    scatter(ax, ones(numel(wv),1) + jitter_w, wv, 30, [0.2 0.5 0.2], 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(ax, 2*ones(numel(bv),1) + jitter_b, bv, 30, [0.7 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.6);

    mw = mean(wv, 'omitnan'); sew = std(wv, 'omitnan') / sqrt(max(1, sum(~isnan(wv))));
    mb = mean(bv, 'omitnan'); seb = std(bv, 'omitnan') / sqrt(max(1, sum(~isnan(bv))));
    errorbar(ax, 1, mw, sew, 'ko', 'MarkerFaceColor', [0.2 0.5 0.2], 'MarkerSize', 9, 'LineWidth', 1.5, 'CapSize', 8);
    errorbar(ax, 2, mb, seb, 'ko', 'MarkerFaceColor', [0.7 0.2 0.2], 'MarkerSize', 9, 'LineWidth', 1.5, 'CapSize', 8);

    yr = [min([wv(:);bv(:)],[],'omitnan'), max([wv(:);bv(:)],[],'omitnan')];
    if diff(yr) == 0 || any(isnan(yr)), yr = [min(0,yr(1)-0.1), max(1,yr(2)+0.1)]; end
    pad = 0.12 * max(diff(yr), eps);
    draw_sig_bracket_local(ax, 1, 2, yr(2) + pad, sprintf('%s (d=%+.2f)', stars_local(p), d));

    set(ax, 'XTick', [1 2], 'XTickLabel', {'well', 'bad'}, 'XLim', [0.5 2.5]);
    title(ax, metrics{mi, 1}, 'FontWeight', 'bold', 'FontSize', 9);
    grid(ax, 'on');

    fprintf('  %-32s well=%.3f±%.3f (n=%d)  bad=%.3f±%.3f (n=%d)  p=%.3f %s  d=%+.2f\n', ...
        metrics{mi,1}, mw, sew, sum(~isnan(wv)), mb, seb, sum(~isnan(bv)), p, stars_local(p), d);
end

sgtitle(fig, sprintf(['Group Analysis — "Well" vs "Bad" Performers (manually assigned, unpaired)\n' ...
        'well n=%d, bad n=%d -- test: two-sample permutation on group means + unpaired Cohen''s d\n' ...
        'Tests whether fusion advantage / classifier quality genuinely differs between strong and weak performers'], ...
        numel(gs_well.subjects), numel(gs_bad.subjects)), 'Interpreter', 'none');

saveas(fig, fullfile(out_dir, '16_well_vs_bad_comparison.svg'), 'svg');
if ~show_figures, close(fig); end
fprintf('Saved 16_well_vs_bad_comparison.svg to %s\n', out_dir);

end % well_vs_bad_comparison


function p = two_sample_perm_test_local(x, y, n_perm)
% TWO_SAMPLE_PERM_TEST_LOCAL  Two-sided permutation test for a difference in
%   means between two INDEPENDENT samples (different subjects in each group,
%   unlike the paired sign-flip tests used in main_group_analysis.m for the
%   SAME subjects across streams). Shuffles pooled group labels; no
%   Statistics Toolbox required.
    if nargin < 3, n_perm = 5000; end
    x = x(~isnan(x)); y = y(~isnan(y));
    nx = numel(x); ny = numel(y);
    if nx < 1 || ny < 1, p = NaN; return; end
    obs = mean(x) - mean(y);
    pooled = [x(:); y(:)];
    n = nx + ny;
    count = 0;
    for k = 1:n_perm
        idx = randperm(n);
        xs = pooled(idx(1:nx)); ys = pooled(idx(nx+1:end));
        if abs(mean(xs) - mean(ys)) >= abs(obs) - eps
            count = count + 1;
        end
    end
    p = count / n_perm;
end

function d = cohen_d_unpaired_local(x, y)
% COHEN_D_UNPAIRED_LOCAL  Unpaired Cohen's d (pooled SD) for two independent
%   groups -- companion to two_sample_perm_test_local.
    x = x(~isnan(x)); y = y(~isnan(y));
    nx = numel(x); ny = numel(y);
    if nx < 2 || ny < 2, d = NaN; return; end
    sp = sqrt(((nx-1)*var(x) + (ny-1)*var(y)) / (nx+ny-2));
    if sp == 0, d = NaN; return; end
    d = (mean(x) - mean(y)) / sp;
end

function s = stars_local(p)
% STARS_LOCAL  Same convention as main_group_analysis.m's stars(): duplicated
%   here (not shared) since it is a private local function there.
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
%   Duplicated from main_group_analysis.m (private local function there).
    line(ax, [x1 x1 x2 x2], [y-0.01 y y y-0.01], 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(ax, (x1+x2)/2, y+0.005, label, 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
end
