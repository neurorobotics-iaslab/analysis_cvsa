%% GROUP_TOPO_ERDERS  Cross-subject (group-average) ERD/ERS topoplots.
%
%   Run this AFTER topo_erders.m has been run for every subject (it reads
%   the topo_erders_channels.mat file each run saves; it does not re-run any
%   pipeline or need EEGLAB itself -- topoplots are drawn with topo_map.m,
%   a self-contained standard-10-20 renderer already used elsewhere in this
%   package, e.g. main_session_overview.m's CSP/sLDA importance figures).
%
%   Recursively scans a root folder (e.g. recordings/) for every
%   <subject>/.../analysis_results/topo_erders/topo_erders_channels.mat.
%   Subject ID = first path component under the root. If a subject has more
%   than one matching file (e.g. several recording days), that subject's
%   own files are averaged first (simple mean) before the cross-subject
%   average -- subject is the unit, matching main_group_analysis.m's
%   convention for every other metric in this package.
%
%   Channel labels are assumed IDENTICAL across all subjects (the LiveAmp
%   montage is fixed hardware, not subject-specific) -- errors clearly if a
%   mismatch is found rather than silently misaligning channel positions.
%
%   For every (paradigm, band-origin, frequency band) found across the
%   cohort, reconstructs the SAME 7-column time-interval topoplot grid as
%   topo_erders.m (Cue+CF, CF only, Cue, then 1s-wide CF bins), averaging
%   each channel's per-time band power across subjects with a simple
%   (unweighted) mean -- every subject counts equally, regardless of trial
%   count, matching how main_group_analysis.m already treats most per-
%   subject metrics (e.g. fusion advantage, ERD-CSP r).
%
%   Two versions per (paradigm, band), same convention as topo_erders.m:
%     - "all_<paradigm>_<band>"  : every channel shown
%     - "sel_<paradigm>_<band>"  : only channels selected by ANY subject's
%                                  CSP for that origin (union across the
%                                  cohort -- continues the same union logic
%                                  topo_erders.m already applies when a
%                                  subject has multiple files), all others
%                                  forced to 0. Skipped if no subject
%                                  contributed a CSP mask for that band.
%   Display mode per band-origin, same as topo_erders.m:
%     - "MI"   : ERD only (values clipped <= 0), two rows (one per cue)
%     - "CVSA" : lateralization index (cue1 - cue2), one row
%   Hybrid paradigm bands carry both origins (MI-origin and CVSA-origin
%   bands are never mixed), so MI and CVSA topoplots for hybrid are always
%   produced and saved separately, same as for the "mi"/"cvsa" paradigms.
%
%   Output: <root>/group_topo_erders/<paradigm>/all_<paradigm>_<origin>_<band>.svg
%   and .../sel_<paradigm>_<origin>_<band>.svg. Console prints, per group,
%   how many of the cohort's subjects contributed and the CSP-mask footprint.

function group_topo_erders(root_dir, subjects_filter, show_figures)
%   Callable as a function:
%     group_topo_erders()                                  % default recordings root, ALL subjects
%     group_topo_erders(root_dir)                           % given root, ALL subjects
%     group_topo_erders(root_dir, {'a1','a2','a3'})         % only these subjects
%     group_topo_erders(root_dir, {'a1','a2','a3'}, true)   % also show figures on screen
%     group_topo_erders([], {'a1','a2'})                    % default root, subset of subjects
%   No GUI folder picker: root_dir defaults to DEFAULT_ROOT below (edit it
%   directly, or pass a path) rather than prompting.

DEFAULT_ROOT = '/home/paolo/bci_vr_ws/recordings';

if nargin < 3, show_figures = false; end
if nargin < 2, subjects_filter = {}; end
if nargin < 1 || isempty(root_dir), root_dir = DEFAULT_ROOT; end
SHOW_FIGURES = show_figures;
if SHOW_FIGURES, fig_vis = 'on'; else, fig_vis = 'off'; end

this_dir = fileparts(mfilename('fullpath'));
ms_dir   = fullfile(this_dir, '..', 'matlab_simulation');
addpath(fullfile(ms_dir, 'utils'));   % topo_map.m

files = dir(fullfile(root_dir, '**', 'topo_erders', 'topo_erders_channels.mat'));
fprintf('Found %d topo_erders_channels.mat under %s\n', numel(files), root_dir);
if isempty(files)
    error('group_topo_erders:nodata', ...
          'No topo_erders_channels.mat found. Run topo_erders.m for each subject first (regenerates this file).');
end

% ── Load all, tag with subject ──────────────────────────────────────────
entries = struct('subject', {}, 'data', {});
for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    s = load(fpath);
    idx = numel(entries) + 1;
    entries(idx).subject = subject_of(fpath, root_dir);
    entries(idx).data    = s.topo_channel_data;
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
        error('group_topo_erders:nosubjects', 'None of the requested subjects were found under %s', root_dir);
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
        error('group_topo_erders:chanmismatch', ...
              'Subject %s has a different channel set than %s -- all subjects must share the same montage.', ...
              entries(i).subject, entries(1).subject);
    end
end
n_ch = numel(ref_labels);

% ── Time axis / interval definitions -- reconstructed exactly as in
%    topo_erders.m, assumed consistent across subjects (fixed epoch design,
%    same resample rate). ────────────────────────────────────────────────
heat_times     = entries(1).data.heat_times;
epoch_limits   = entries(1).data.epoch_limits;
CUE_DURATION_S = entries(1).data.cue_duration_s;
T = epoch_limits(2);
intervals   = [0 T; CUE_DURATION_S T];
titles_cols = {sprintf('Cue+CF (0-%ds)', T), sprintf('CF only (%g-%ds)', CUE_DURATION_S, T)};
intervals(end+1, :) = [0, CUE_DURATION_S];
titles_cols{end+1}  = sprintf('Cue (0-%gs)', CUE_DURATION_S);
s_t = CUE_DURATION_S;
while s_t < T
    e_t = min(s_t + 1, T);
    intervals(end+1, :) = [s_t, e_t];           %#ok<AGROW>
    titles_cols{end+1}  = sprintf('CF (%g-%gs)', s_t, e_t); %#ok<AGROW>
    s_t = e_t;
end
cf_window_idx = heat_times >= CUE_DURATION_S*1000 & heat_times <= epoch_limits(2)*1000;

% ── Group all subjects' band entries into (paradigm, origin, freq) keys ──
keys = {};
for i = 1:numel(entries)
    for b = 1:numel(entries(i).data.bands)
        bd = entries(i).data.bands(b);
        keys{end+1} = sprintf('%s|%s|%.1f|%.1f', bd.paradigm, bd.band_origin, bd.freq_lo, bd.freq_hi); %#ok<AGROW>
    end
end
ukeys = unique(keys, 'stable');
fprintf('\nFound %d (paradigm, band-origin, frequency) groups across the cohort\n', numel(ukeys));

out_root = fullfile(root_dir, 'group_topo_erders');
if ~exist(out_root, 'dir'), mkdir(out_root); end

for kk = 1:numel(ukeys)
    parts    = strsplit(ukeys{kk}, '|');
    paradigm = parts{1};
    origin   = parts{2};
    freq_lo  = str2double(parts{3});
    freq_hi  = str2double(parts{4});

    % ── Per subject: average that subject's own matching entries first
    %    (e.g. several recording days), simple mean ──────────────────────
    subj_c1   = nan(n_ch, numel(heat_times), n_subj);
    subj_c2   = nan(n_ch, numel(heat_times), n_subj);
    subj_mask = false(n_ch, n_subj);
    subj_has  = false(1, n_subj);
    for si = 1:n_subj
        idx_e = find(strcmp({entries.subject}, subjects{si}));
        c1_list = {}; c2_list = {}; mask_list = {};
        for ie = idx_e
            for b = 1:numel(entries(ie).data.bands)
                bd = entries(ie).data.bands(b);
                if strcmp(bd.paradigm, paradigm) && strcmpi(bd.band_origin, origin) && ...
                        abs(bd.freq_lo - freq_lo) < 0.05 && abs(bd.freq_hi - freq_hi) < 0.05
                    c1_list{end+1} = bd.c1; %#ok<AGROW>
                    c2_list{end+1} = bd.c2; %#ok<AGROW>
                    if strcmpi(origin, 'MI'), m = bd.mi_mask; else, m = bd.cvsa_mask; end
                    if ~isempty(m), mask_list{end+1} = m(:); end %#ok<AGROW>
                end
            end
        end
        if isempty(c1_list), continue; end
        subj_c1(:,:,si) = mean(cat(3, c1_list{:}), 3, 'omitnan');
        subj_c2(:,:,si) = mean(cat(3, c2_list{:}), 3, 'omitnan');
        if ~isempty(mask_list)
            subj_mask(:,si) = any(cat(2, mask_list{:}), 2);
        end
        subj_has(si) = true;
    end
    if ~any(subj_has), continue; end

    % ── Cross-subject average: simple (unweighted) mean over subjects,
    %    CSP mask = UNION across subjects (continues the same union logic
    %    topo_erders.m already applies within one subject's own files). ──
    grp_c1   = mean(subj_c1(:,:,subj_has), 3, 'omitnan');
    grp_c2   = mean(subj_c2(:,:,subj_has), 3, 'omitnan');
    grp_mask = any(subj_mask(:, subj_has), 2);
    n_contrib = sum(subj_has);

    band_name = sprintf('%s_%g-%gHz', origin, freq_lo, freq_hi);
    par_dir = fullfile(out_root, paradigm);
    if ~exist(par_dir, 'dir'), mkdir(par_dir); end

    fprintf('  [%-6s] %-4s %5.1f-%5.1f Hz : %d/%d subjects  (CSP mask: %s)\n', paradigm, origin, freq_lo, freq_hi, ...
            n_contrib, n_subj, mask_summary(grp_mask));

    % ── Colour limits: PERCENTILE-based (not max/min), same convention as
    %    topo_erders.m otherwise. A single noisy channel/timepoint can
    %    otherwise stretch the whole scale and wash out every real
    %    difference -- CLIM_PCTL clips that tail; values beyond it just
    %    saturate to the strongest colour instead of being lost. ──────────
    CLIM_PCTL = 95;
    is_mi = strcmpi(origin, 'MI');
    if is_mi
        neg = min(cat(1, grp_c1(:,cf_window_idx), grp_c2(:,cf_window_idx)), 0);
        mx = local_prctile(-neg(:), CLIM_PCTL); if mx == 0 || isnan(mx), mx = 1; end
        clim = [-mx 0];
    else
        d = grp_c1(:,cf_window_idx) - grp_c2(:,cf_window_idx);
        mx = local_prctile(abs(d(:)), CLIM_PCTL); if mx == 0 || isnan(mx), mx = 1; end
        clim = [-mx mx];
    end

    plot_group_topo_grid(subj_c1, subj_c2, subj_has, ref_labels, heat_times, intervals, titles_cols, origin, clim, [], ...
                          par_dir, sprintf('all_%s_%s', paradigm, band_name), n_contrib, n_subj, freq_lo, freq_hi, fig_vis);
    if any(grp_mask)
        plot_group_topo_grid(subj_c1, subj_c2, subj_has, ref_labels, heat_times, intervals, titles_cols, origin, clim, grp_mask, ...
                              par_dir, sprintf('sel_%s_%s', paradigm, band_name), n_contrib, n_subj, freq_lo, freq_hi, fig_vis);
    end
end

fprintf('\nSaved group-level ERD/ERS topoplots to %s\n', out_root);

end % group_topo_erders

% ── Local helpers ────────────────────────────────────────────────────────

function subj = subject_of(fpath, root_dir)
    rel = erase(fpath, [char(root_dir) filesep]);
    parts = strsplit(rel, filesep);
    subj = parts{1};
end

function s = mask_summary(mask)
    if isempty(mask) || ~any(mask), s = 'no CSP channels';
    else, s = sprintf('%d/%d channels (union across cohort)', sum(mask), numel(mask));
    end
end

function plot_group_topo_grid(subj_c1, subj_c2, subj_has, ch_names, heat_times, intervals, titles_cols, origin, clim, mask, ...
                               out_dir, prefix, n_contrib, n_subj, freq_lo, freq_hi, fig_vis)
% PLOT_GROUP_TOPO_GRID  One figure: rows = cue(s), cols = time interval,
%   drawn with topo_map.m (no EEGLAB dependency). Mirrors topo_erders.m's
%   plot_topo_grid, but averaged across subjects upstream of this call.
%
%   subj_c1/subj_c2  [n_ch x n_time x n_subj] per-subject band power (NOT
%                     pre-averaged) -- kept per-subject here, rather than
%                     passing in the group mean, so each cell can run its
%                     OWN per-channel significance test (subject is the
%                     unit) on exactly the same time-window slice it plots.
%   subj_has          [1 x n_subj] logical, which subjects contributed.
%
%   Each cell's displayed value is still the plain cross-subject mean (same
%   number as before); additionally, a black ring is drawn over channels
%   where a one-sample t-test across subjects (H0: mean=0, two-sided,
%   hand-rolled via betainc -- no Statistics Toolbox) reaches p<0.05, so a
%   real, consistent-across-subjects effect can be told apart from noise at
%   a glance. With few subjects this test is naturally conservative (a
%   parametric t-test has no permutation-count floor, unlike a sign-flip
%   test, so it stays usable from n=2 subjects and gets more powerful as
%   the cohort grows).
    num_intervals = size(intervals, 1);
    is_mi = strcmpi(origin, 'MI');
    if is_mi
        num_rows = 2;
        row_labels = {'Cue 1 (ERD)', 'Cue 2 (ERD)'};
    else
        num_rows = 1;
        row_labels = {'Cue 1 - Cue 2 (lateralization)'};
    end
    bg_zero = ~isempty(mask);   % "sel": fill non-CSP channels with 0 for a smooth full-head map

    h = figure('Name', sprintf('%s — %s %g-%g Hz (group, n=%d/%d)', prefix, origin, freq_lo, freq_hi, n_contrib, n_subj), ...
               'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
    set(h, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    tl = tiledlayout(num_rows, num_intervals, 'TileSpacing', 'compact', 'Padding', 'compact');

    for r = 1:num_rows
        for c = 1:num_intervals
            ax = nexttile;
            t_idx = heat_times >= intervals(c,1)*1000 & heat_times < intervals(c,2)*1000;

            % Per-subject channel value for this cell (raw, unclipped) --
            % feeds both the displayed group mean and the significance test.
            subj_vals = nan(numel(ch_names), numel(subj_has));
            for si = find(subj_has)
                if is_mi
                    src = subj_c1(:,:,si); if r == 2, src = subj_c2(:,:,si); end
                    subj_vals(:,si) = mean(src(:, t_idx), 2, 'omitnan');
                else
                    d1 = mean(subj_c1(:, t_idx, si), 2, 'omitnan');
                    d2 = mean(subj_c2(:, t_idx, si), 2, 'omitnan');
                    subj_vals(:,si) = d1 - d2;
                end
            end
            data = mean(subj_vals, 2, 'omitnan');
            if is_mi
                data(data > 0) = 0;   % ERD only (display convention)
            end
            data(isnan(data)) = 0;

            sig_chan = local_ttest_pvals(subj_vals') < 0.05;   % [1 x n_ch], subj_vals' is [n_subj x n_ch]

            if ~isempty(mask)
                names_here = ch_names(mask);
                vals_here  = data(mask);
                sig_here   = sig_chan(mask);
            else
                names_here = ch_names;
                vals_here  = data;
                sig_here   = sig_chan;
            end

            topo_map(names_here, vals_here, clim, ax, '', c == num_intervals, bg_zero, sig_here);
            if r == 1, title(ax, titles_cols{c}, 'FontSize', 9); end
            if c == 1, ylabel(ax, row_labels{r}, 'Visible', 'on', 'FontWeight', 'bold'); end
        end
    end

    title(tl, sprintf('%s — %s %g-%g Hz  (group average, n=%d/%d subjects | ring = p<0.05 cross-subject t-test)', ...
          prefix, origin, freq_lo, freq_hi, n_contrib, n_subj), 'Interpreter', 'none', 'FontWeight', 'bold');

    save_path = fullfile(out_dir, sprintf('%s.svg', prefix));
    saveas(h, save_path, 'svg');
    if ~strcmp(fig_vis, 'on'), close(h); end
    fprintf('    Saved: %s\n', save_path);
end

function p = local_ttest_pvals(X)
% LOCAL_TTEST_PVALS  Two-sided one-sample t-test (H0: mean=0) per column, no
%   Statistics Toolbox required (betainc is base MATLAB). X: [n_subj x n_ch],
%   NaN allowed for a missing subject/channel. Unlike a sign-flip permutation
%   test (2^n_subj achievable p-values), this parametric test has no
%   small-cohort resolution floor -- at the cost of assuming the per-subject
%   values are roughly normally distributed.
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
