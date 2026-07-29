%% GROUP_CSP_SLDA_IMPORTANCE  Cross-subject (group-average) CSP + sLDA
%   channel-importance topoplots, per band, per paradigm (MI/CVSA).
%
%   Run this AFTER main_session_overview.m has been run for every subject
%   (it reads the csp_slda_channels.mat file each run saves alongside
%   session_summary.mat; it does not re-run any pipeline). No EEGLAB
%   dependency -- topoplots are drawn with topo_map.m, the same
%   self-contained standard-10-20 renderer used by group_topo_erders.m.
%
%   Recursively scans a root folder (e.g. recordings/) for every
%   <subject>/.../analysis_results/session_overview/csp_slda_channels.mat.
%   Subject ID = first path component under the root. If a subject has more
%   than one matching file (e.g. several recording days), that subject's
%   own files are averaged first (simple mean) before the cross-subject
%   average -- subject is the unit, same convention as group_topo_erders.m
%   and main_group_analysis.m.
%
%   Channel labels are assumed IDENTICAL across all subjects (fixed
%   hardware montage) -- errors clearly on a mismatch. CSP-selected
%   channels differ per subject; each subject's CSP weights are
%   re-expressed in the FULL reference channel space (zero for channels
%   that subject's CSP didn't select) by main_session_overview.m before
%   saving, so the cross-subject mean below is a simple average over that
%   full space -- a channel only some subjects selected is diluted, not
%   excluded (same treatment group_topo_erders.m gives ERD/ERS).
%
%   BAND DIVISION IS KEPT throughout (unlike main_session_overview.m's own
%   per-subject Fig 4/5, which pool CSP/sLDA importance across all bands
%   into one number) -- every figure below has one row per band.
%
%   Per paradigm (MI, CVSA), three outputs under <root>/group_csp_slda_importance/:
%     1. <paradigm>_csp_importance.svg   -- rows = bands, columns = [class-1
%        weight, class-2 weight, selectivity (c1-c2)/(c1+c2)], each a
%        cross-subject-averaged topoplot. Weight panels use jet with clim
%        [0, max] (one-sided, non-negative -- a diverging map would put 0 at
%        one extreme instead of a neutral centre, see topo_erders.m's own
%        MI-band fix); selectivity is genuinely bidirectional and keeps the
%        diverging RdBu (clim symmetric about 0).
%     2. <paradigm>_slda_importance.svg  -- rows = bands, one column: channel
%        importance weighted by |sLDA coef| x |CSP filter| (jet, one-sided,
%        same reasoning).
%     3. <paradigm>_slda_feature_selection.svg -- comp x band heatmap: % of
%        contributing subjects that selected that (component, band) feature,
%        with the mean |sLDA coef| among those who selected it annotated in
%        parenthesis. This IS the natural cross-subject "average" of a
%        per-subject binary-ish selection pattern (main_session_overview.m
%        Fig 5's own per-subject feat_mat heatmap).
%
%   Each topoplot cell (per band, per column) gets its OWN colour scale
%   (CLIM_PCTL-th percentile, default 99, of the cross-subject mean value at
%   that cell -- NOT one shared scale for the whole image, same reasoning as
%   group_topo_erders.m) and a black significance ring (one-sample t-test
%   across subjects, H0: mean=0, p<0.05, hand-rolled via betainc -- no
%   Statistics Toolbox, no small-cohort permutation floor).
%
%   Console prints, per paradigm, how many of the cohort's subjects
%   contributed.

function group_csp_slda_importance(root_dir, subjects_filter, show_figures)
%   Callable as a function:
%     group_csp_slda_importance()                                  % default recordings root, ALL subjects
%     group_csp_slda_importance(root_dir)                          % given root, ALL subjects
%     group_csp_slda_importance(root_dir, {'a1','a2','a3'})        % only these subjects
%     group_csp_slda_importance(root_dir, {'a1','a2','a3'}, true)  % also show figures on screen
%     group_csp_slda_importance([], {'a1','a2'})                   % default root, subset of subjects
%   No GUI folder picker: root_dir defaults to DEFAULT_ROOT below (edit it
%   directly, or pass a path) rather than prompting.

DEFAULT_ROOT = '/home/paolo/bci_vr_ws/recordings';

if nargin < 3, show_figures = false; end
if nargin < 2, subjects_filter = {}; end
if nargin < 1 || isempty(root_dir), root_dir = DEFAULT_ROOT; end
SHOW_FIGURES = show_figures;
if SHOW_FIGURES, fig_vis = 'on'; else, fig_vis = 'off'; end

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'utils'));   % topo_map.m

files = dir(fullfile(root_dir, '**', 'csp_slda_channels.mat'));
fprintf('Found %d csp_slda_channels.mat under %s\n', numel(files), root_dir);
if isempty(files)
    error('group_csp_slda_importance:nodata', ...
          'No csp_slda_channels.mat found. Run main_session_overview.m for each subject first (regenerates this file).');
end

% ── Load all, tag with subject ──────────────────────────────────────────
entries = struct('subject', {}, 'data', {});
for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    s = load(fpath);
    idx = numel(entries) + 1;
    entries(idx).subject = subject_of(fpath, root_dir);
    entries(idx).data    = s.csp_slda_data;
end

% ── Optional subject filter: only include the requested subjects ────────
if ~isempty(subjects_filter)
    keep = ismember({entries.subject}, subjects_filter);
    fprintf('Subject filter requested: %s  (%d/%d loaded entries match)\n', ...
            strjoin(subjects_filter, ', '), sum(keep), numel(entries));
    entries = entries(keep);
    missing = setdiff(subjects_filter, unique({entries.subject}));
    if ~isempty(missing)
        fprintf('  WARNING: no data found for requested subject(s): %s\n', strjoin(missing, ', '));
    end
    if isempty(entries)
        error('group_csp_slda_importance:nosubjects', 'None of the requested subjects were found under %s', root_dir);
    end
end

subjects = unique({entries.subject}, 'stable');
n_subj   = numel(subjects);
fprintf('Subjects found (n=%d): %s\n', n_subj, strjoin(subjects, ', '));

% ── Reference channel labels must match across ALL subjects (fixed
%    hardware montage) -- error clearly rather than silently misaligning. ──
ref_labels = entries(1).data.ref_labels;
for i = 2:numel(entries)
    if ~isequal(entries(i).data.ref_labels, ref_labels)
        error('group_csp_slda_importance:chanmismatch', ...
              'Subject %s has a different channel set than %s -- all subjects must share the same montage.', ...
              entries(i).subject, entries(1).subject);
    end
end
n_ch = numel(ref_labels);

out_root = fullfile(root_dir, 'group_csp_slda_importance');
if ~exist(out_root, 'dir'), mkdir(out_root); end

CLIM_PCTL = 99;

for par_c = {'mi', 'cvsa'}
    pname = par_c{1};
    has_par = arrayfun(@(e) isfield(e.data, pname), entries);
    if ~any(has_par), continue; end

    ref_idx      = find(has_par, 1);
    bands        = entries(ref_idx).data.(pname).bands;          % [n_bands x 2] Hz
    n_bands      = size(bands, 1);
    n_components = entries(ref_idx).data.(pname).n_components;

    subj_w1   = nan(n_bands, n_ch, n_subj);
    subj_w2   = nan(n_bands, n_ch, n_subj);
    subj_sel  = nan(n_bands, n_ch, n_subj);
    subj_slda = nan(n_bands, n_ch, n_subj);
    subj_has  = false(1, n_subj);
    feat_mats = nan(n_components, n_bands, n_subj);
    n_skipped_bandmismatch = 0;

    for si = 1:n_subj
        idx_e = find(strcmp({entries.subject}, subjects{si}));
        w1_list = {}; w2_list = {}; sel_list = {}; slda_list = {}; feat_list = {};
        for ie = idx_e
            if ~isfield(entries(ie).data, pname), continue; end
            pe = entries(ie).data.(pname);
            if size(pe.bands, 1) ~= n_bands
                n_skipped_bandmismatch = n_skipped_bandmismatch + 1;
                continue;   % different band configuration than the reference -- skip rather than misalign rows
            end
            w1_list{end+1}   = pe.w_c1_band;         %#ok<AGROW>
            w2_list{end+1}   = pe.w_c2_band;         %#ok<AGROW>
            sel_list{end+1}  = pe.selectivity_band;  %#ok<AGROW>
            slda_list{end+1} = pe.slda_w_ch_band;    %#ok<AGROW>
            feat_list{end+1} = pe.feat_mat;          %#ok<AGROW>
        end
        if isempty(w1_list), continue; end
        subj_w1(:,:,si)   = mean(cat(3, w1_list{:}),   3, 'omitnan');
        subj_w2(:,:,si)   = mean(cat(3, w2_list{:}),   3, 'omitnan');
        subj_sel(:,:,si)  = mean(cat(3, sel_list{:}),  3, 'omitnan');
        subj_slda(:,:,si) = mean(cat(3, slda_list{:}), 3, 'omitnan');
        feat_mats(:,:,si) = mean(cat(3, feat_list{:}), 3, 'omitnan');
        subj_has(si) = true;
    end
    if ~any(subj_has), continue; end
    n_contrib = sum(subj_has);

    if n_skipped_bandmismatch > 0
        fprintf('  [%s] WARNING: %d file(s) skipped (band configuration differs from the reference subject''s)\n', ...
                upper(pname), n_skipped_bandmismatch);
    end
    fprintf('[%s] CSP/sLDA importance: %d/%d subjects contributed\n', upper(pname), n_contrib, n_subj);

    plot_csp_importance_grid(subj_w1, subj_w2, subj_sel, subj_has, ref_labels, bands, CLIM_PCTL, ...
                              out_root, sprintf('%s_csp_importance', pname), n_contrib, n_subj, fig_vis);

    plot_slda_importance_grid(subj_slda, subj_has, ref_labels, bands, CLIM_PCTL, ...
                               out_root, sprintf('%s_slda_importance', pname), n_contrib, n_subj, fig_vis);

    plot_feature_selection_consistency(feat_mats, subj_has, bands, n_components, ...
                                        out_root, sprintf('%s_slda_feature_selection', pname), n_contrib, n_subj, fig_vis);
end

fprintf('\nSaved group-level CSP/sLDA importance figures to %s\n', out_root);

end % group_csp_slda_importance

% ── Local helpers ────────────────────────────────────────────────────────

function subj = subject_of(fpath, root_dir)
    rel = erase(fpath, [char(root_dir) filesep]);
    parts = strsplit(rel, filesep);
    subj = parts{1};
end

function plot_csp_importance_grid(subj_w1, subj_w2, subj_sel, subj_has, ch_names, bands, CLIM_PCTL, ...
                                   out_dir, prefix, n_contrib, n_subj, fig_vis)
% PLOT_CSP_IMPORTANCE_GRID  Rows = bands, columns = [class-1 weight,
%   class-2 weight, selectivity], cross-subject mean, one topoplot per cell.
    n_bands = size(bands, 1);
    col_labels = {'Class-1 weight', 'Class-2 weight', 'Selectivity (c1-c2)/(c1+c2)'};

    h = figure('Name', sprintf('%s (group, n=%d/%d)', prefix, n_contrib, n_subj), ...
               'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
    set(h, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    tl = tiledlayout(n_bands, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    tl.Position = [0.03 0.08 0.94 0.84];

    for b = 1:n_bands
        for c = 1:3
            ax = nexttile;
            switch c
                case 1, src = subj_w1;
                case 2, src = subj_w2;
                case 3, src = subj_sel;
            end
            vals_subj = squeeze(src(b,:,subj_has));   % [n_ch x n_contrib]
            data = mean(vals_subj, 2, 'omitnan');
            sig  = local_ttest_pvals(vals_subj') < 0.05;   % [1 x n_ch]

            if c == 3
                mx = local_prctile(abs(data), CLIM_PCTL);
                if mx == 0 || isnan(mx), mx = 1; end
                clim_cell = [-mx mx];
                cmap_here = [];
            else
                mx = local_prctile(data, CLIM_PCTL);
                if mx == 0 || isnan(mx), mx = 1; end
                clim_cell = [0 mx];
                cmap_here = jet(256);
            end

            topo_map(ch_names, data, clim_cell, ax, '', true, false, sig, cmap_here);
            if b == 1, title(ax, col_labels{c}, 'FontSize', 9); end
            if c == 1
                ylabel(ax, sprintf('%g-%g Hz', bands(b,1), bands(b,2)), 'Visible', 'on', 'FontWeight', 'bold');
            end
        end
    end

    annotation(h, 'textbox', [0 0.95 1 0.05], 'String', ...
        sprintf('%s  (n=%d/%d subjects)', prefix, n_contrib, n_subj), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');
    annotation(h, 'textbox', [0 0 1 0.04], 'String', ...
        'dot = cross-subject mean;  black ring = significant (p<0.05, one-sample t-test across subjects)', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 9);

    save_path = fullfile(out_dir, sprintf('%s.svg', prefix));
    saveas(h, save_path, 'svg');
    if ~strcmp(fig_vis, 'on'), close(h); end
    fprintf('  Saved: %s\n', save_path);
end

function plot_slda_importance_grid(subj_slda, subj_has, ch_names, bands, CLIM_PCTL, ...
                                    out_dir, prefix, n_contrib, n_subj, fig_vis)
% PLOT_SLDA_IMPORTANCE_GRID  One column, rows = bands: cross-subject mean
%   channel importance weighted by |sLDA coef| x |CSP filter|.
    n_bands = size(bands, 1);

    h = figure('Name', sprintf('%s (group, n=%d/%d)', prefix, n_contrib, n_subj), ...
               'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
    set(h, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    tl = tiledlayout(1, n_bands, 'TileSpacing', 'compact', 'Padding', 'compact');
    tl.Position = [0.03 0.08 0.94 0.84];

    for b = 1:n_bands
        ax = nexttile;
        vals_subj = squeeze(subj_slda(b,:,subj_has));   % [n_ch x n_contrib]
        data = mean(vals_subj, 2, 'omitnan');
        sig  = local_ttest_pvals(vals_subj') < 0.05;

        mx = local_prctile(data, CLIM_PCTL);
        if mx == 0 || isnan(mx), mx = 1; end

        topo_map(ch_names, data, [0 mx], ax, sprintf('%g-%g Hz', bands(b,1), bands(b,2)), true, false, sig, jet(256));
    end

    annotation(h, 'textbox', [0 0.95 1 0.05], 'String', ...
        sprintf('%s  (n=%d/%d subjects)', prefix, n_contrib, n_subj), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');
    annotation(h, 'textbox', [0 0 1 0.04], 'String', ...
        'dot = cross-subject mean |sLDA coef| x |CSP filter|;  black ring = significant (p<0.05, one-sample t-test across subjects)', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 9);

    save_path = fullfile(out_dir, sprintf('%s.svg', prefix));
    saveas(h, save_path, 'svg');
    if ~strcmp(fig_vis, 'on'), close(h); end
    fprintf('  Saved: %s\n', save_path);
end

function plot_feature_selection_consistency(feat_mats, subj_has, bands, n_components, ...
                                             out_dir, prefix, n_contrib, n_subj, fig_vis)
% PLOT_FEATURE_SELECTION_CONSISTENCY  comp x band heatmap: % of contributing
%   subjects who selected that feature, mean |sLDA coef| among selectors
%   annotated in parenthesis.
    n_bands = size(bands, 1);
    contrib_mats = feat_mats(:,:,subj_has);
    sel_rate  = 100 * mean(~isnan(contrib_mats), 3);
    mean_coef = mean(contrib_mats, 3, 'omitnan');

    h = figure('Name', sprintf('%s (group, n=%d/%d)', prefix, n_contrib, n_subj), ...
               'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
    set(h, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    ax = axes(h);
    imagesc(ax, sel_rate, [0 100]);
    colormap(ax, parula(256));
    colorbar(ax);

    band_labels = arrayfun(@(b) sprintf('%g-%g', bands(b,1), bands(b,2)), 1:n_bands, 'UniformOutput', false);
    set(ax, 'XTick', 1:n_bands, 'XTickLabel', band_labels, 'YTick', 1:n_components);
    xlabel(ax, 'band (Hz)'); ylabel(ax, 'CSP component');

    for comp = 1:n_components
        for b = 1:n_bands
            txt_col = 'k'; if sel_rate(comp,b) > 60, txt_col = 'w'; end
            if isnan(mean_coef(comp,b))
                txt = sprintf('%.0f%%', sel_rate(comp,b));
            else
                txt = sprintf('%.0f%%\n(%.2f)', sel_rate(comp,b), mean_coef(comp,b));
            end
            text(ax, b, comp, txt, 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', txt_col);
        end
    end

    title(ax, sprintf('%s -- %% of subjects selecting each feature (n=%d/%d)\nparenthesis = mean |sLDA coef| among those who selected it', ...
          prefix, n_contrib, n_subj), 'Interpreter', 'none', 'FontSize', 10);

    save_path = fullfile(out_dir, sprintf('%s.svg', prefix));
    saveas(h, save_path, 'svg');
    if ~strcmp(fig_vis, 'on'), close(h); end
    fprintf('  Saved: %s\n', save_path);
end

function p = local_ttest_pvals(X)
% LOCAL_TTEST_PVALS  Two-sided one-sample t-test (H0: mean=0) per column, no
%   Statistics Toolbox required (betainc is base MATLAB). X: [n_subj x n_ch],
%   NaN allowed for a missing subject/channel.
    [~, nch] = size(X);
    p = nan(1, nch);
    for c = 1:nch
        x = X(~isnan(X(:,c)), c);
        nx = numel(x);
        if nx < 2, continue; end
        sx = std(x);
        if sx == 0, continue; end
        t  = mean(x) / (sx / sqrt(nx));
        df = nx - 1;
        p(c) = betainc(df / (df + t^2), df/2, 0.5);
    end
end

function v = local_prctile(x, pct)
% LOCAL_PRCTILE  Linear-interpolation percentile, no Statistics Toolbox
%   required (equivalent to MATLAB's default prctile method on a vector).
    x = x(isfinite(x));
    x = sort(x(:));
    n = numel(x);
    if n == 0, v = NaN; return; end
    if n == 1, v = x(1); return; end
    idx = (pct/100) * (n-1) + 1;
    lo = floor(idx); hi = ceil(idx);
    if lo == hi, v = x(lo); else, v = x(lo) + (x(hi)-x(lo)) * (idx-lo); end
end
