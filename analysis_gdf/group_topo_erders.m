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
%   than one matching file for the SAME session type (e.g. several
%   recording days), that subject's own files are averaged first (simple
%   mean) before the cross-subject average -- subject is the unit, matching
%   main_group_analysis.m's convention for every other metric in this
%   package.
%
%   EVALUATION vs CALIBRATION are kept as two entirely separate groups,
%   never pooled together: run_subject_analysis.m calls topo_erders.m once
%   on the evaluation/ folder and once on the sibling calibration/ folder
%   (per day), so a subject can have BOTH an evaluation and a calibration
%   topo_erders_channels.mat for the same paradigm/band -- pooling them
%   would silently mix two different recording contexts (autopilot/real
%   feedback aside, ERD/ERS during calibration can differ from evaluation,
%   e.g. more trials, different day, no closed-loop feedback). Session type
%   ('evaluation'/'calibration'/'unknown') is inferred from the folder path
%   itself (topo_erders_channels.mat carries no such field), since that is
%   the only place this distinction exists on disk.
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
%   Output: <root>/group_topo_erders/<session>/<paradigm>/all_<paradigm>_<origin>_<band>.svg
%   and .../sel_<paradigm>_<origin>_<band>.svg, where <session> is
%   'evaluation' or 'calibration' (or 'unknown' as a graceful fallback).
%   Console prints, per group, how many of the cohort's subjects
%   contributed and the CSP-mask footprint.
%
%   Additionally, for every (session, paradigm, origin) found, one extra
%   summary image ".../all_<paradigm>_<origin>_CF_only_summary.svg" (and
%   "sel_..." if a CSP mask exists) puts the "CF only" column (whole-
%   feedback-period average) of EVERY band side by side, so you don't have
%   to open each band's own 7-column figure just to compare column 2
%   across bands.

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

% ── Load all, tag with subject + session type (evaluation/calibration) ──
entries = struct('subject', {}, 'session', {}, 'data', {});
for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    s = load(fpath);
    idx = numel(entries) + 1;
    entries(idx).subject = subject_of(fpath, root_dir);
    entries(idx).session = session_of(fpath);
    entries(idx).data    = s.topo_channel_data;
end
n_eval  = sum(strcmp({entries.session}, 'evaluation'));
n_calib = sum(strcmp({entries.session}, 'calibration'));
n_unk   = sum(strcmp({entries.session}, 'unknown'));
if n_unk > 0
    fprintf('  by session: %d evaluation, %d calibration, %d unknown\n', n_eval, n_calib, n_unk);
else
    fprintf('  by session: %d evaluation, %d calibration\n', n_eval, n_calib);
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
% ── Group all subjects' band entries into (session, paradigm, origin, freq)
%    keys -- session kept separate so evaluation and calibration are NEVER
%    pooled into the same average (see header comment above). ────────────
keys = {};
for i = 1:numel(entries)
    for b = 1:numel(entries(i).data.bands)
        bd = entries(i).data.bands(b);
        keys{end+1} = sprintf('%s|%s|%s|%.1f|%.1f', entries(i).session, bd.paradigm, bd.band_origin, bd.freq_lo, bd.freq_hi); %#ok<AGROW>
    end
end
ukeys = unique(keys, 'stable');
fprintf('\nFound %d (session, paradigm, band-origin, frequency) groups across the cohort\n', numel(ukeys));

out_root = fullfile(root_dir, 'group_topo_erders');
if ~exist(out_root, 'dir'), mkdir(out_root); end

% Accumulates, per band, the already-computed subject-then-cross-subject
% arrays below -- reused after the main loop to build the "CF only, all
% bands in one image" summary figure without recomputing any averaging.
cf_summary_entries = struct('session', {}, 'paradigm', {}, 'origin', {}, 'band_name', {}, ...
                             'freq_lo', {}, 'freq_hi', {}, 'subj_c1', {}, 'subj_c2', {}, ...
                             'subj_has', {}, 'grp_mask', {}, 'n_contrib', {});

for kk = 1:numel(ukeys)
    parts    = strsplit(ukeys{kk}, '|');
    session  = parts{1};
    paradigm = parts{2};
    origin   = parts{3};
    freq_lo  = str2double(parts{4});
    freq_hi  = str2double(parts{5});

    % ── Per subject: average that subject's own matching entries first
    %    (e.g. several recording days OF THE SAME SESSION TYPE), simple
    %    mean ─────────────────────────────────────────────────────────────
    subj_c1   = nan(n_ch, numel(heat_times), n_subj);
    subj_c2   = nan(n_ch, numel(heat_times), n_subj);
    subj_mask = false(n_ch, n_subj);
    subj_has  = false(1, n_subj);
    for si = 1:n_subj
        idx_e = find(strcmp({entries.subject}, subjects{si}) & strcmp({entries.session}, session));
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

    % ── CSP mask = UNION across subjects (continues the same union logic
    %    topo_erders.m already applies within one subject's own files).
    %    The cross-subject MEAN itself (used to be grp_c1/grp_c2 here) is
    %    now computed per-cell inside plot_group_topo_grid instead, since
    %    that's also where the per-cell colour scale needs it. ───────────
    grp_mask = any(subj_mask(:, subj_has), 2);
    n_contrib = sum(subj_has);

    band_name = sprintf('%s_%g-%gHz', origin, freq_lo, freq_hi);

    cfe = numel(cf_summary_entries) + 1;
    cf_summary_entries(cfe).session   = session;
    cf_summary_entries(cfe).paradigm  = paradigm;
    cf_summary_entries(cfe).origin    = origin;
    cf_summary_entries(cfe).band_name = band_name;
    cf_summary_entries(cfe).freq_lo   = freq_lo;
    cf_summary_entries(cfe).freq_hi   = freq_hi;
    cf_summary_entries(cfe).subj_c1   = subj_c1;
    cf_summary_entries(cfe).subj_c2   = subj_c2;
    cf_summary_entries(cfe).subj_has  = subj_has;
    cf_summary_entries(cfe).grp_mask  = grp_mask;
    cf_summary_entries(cfe).n_contrib = n_contrib;
    par_dir = fullfile(out_root, session, paradigm);
    if ~exist(par_dir, 'dir'), mkdir(par_dir); end

    fprintf('  [%-11s|%-6s] %-4s %5.1f-%5.1f Hz : %d/%d subjects  (CSP mask: %s)\n', session, paradigm, origin, freq_lo, freq_hi, ...
            n_contrib, n_subj, mask_summary(grp_mask));

    % ── Colour limits: PERCENTILE-based (CLIM_PCTL=99, set to 100 to
    %    disable), computed PER TOPOPLOT (per row x column cell) inside
    %    plot_group_topo_grid below -- NOT one shared scale for the whole
    %    image. A single shared scale meant a cell with genuinely low
    %    variance (e.g. a quiet 1s bin) got rendered on a scale sized for
    %    whatever OTHER cell in the same image had the biggest swing,
    %    crushing its own real (if modest) contrast to a single flat
    %    colour -- indistinguishable from "no ERD at all" even when there
    %    was some. Each topoplot now shows its own full colour range,
    %    computed from its own data -- see local_cell_clim() below. ───────
    CLIM_PCTL = 99;

    plot_group_topo_grid(subj_c1, subj_c2, subj_has, ref_labels, heat_times, intervals, titles_cols, origin, CLIM_PCTL, [], ...
                          par_dir, sprintf('all_%s_%s', paradigm, band_name), n_contrib, n_subj, freq_lo, freq_hi, fig_vis, session);
    if any(grp_mask)
        plot_group_topo_grid(subj_c1, subj_c2, subj_has, ref_labels, heat_times, intervals, titles_cols, origin, CLIM_PCTL, grp_mask, ...
                              par_dir, sprintf('sel_%s_%s', paradigm, band_name), n_contrib, n_subj, freq_lo, freq_hi, fig_vis, session);
    end
end

% ── "CF only" summary: ALL BANDS of a given (session, paradigm, origin)
%    side by side in ONE image, instead of opening N separate per-band
%    7-column figures and picking out column 2 ("CF only") from each.
%    Reuses the subj_c1/subj_c2/subj_has/grp_mask already computed above --
%    same cross-subject averaging convention as everywhere else in this
%    file (each subject's own files averaged first, THEN a simple
%    unweighted mean across subjects -- subject is the unit, not a
%    trial-pooled grand average; "media della media", not "media
%    generale") and the same per-cell percentile colour scale +
%    significance-ring convention as plot_group_topo_grid.
if ~isempty(cf_summary_entries)
    grp_keys = arrayfun(@(e) sprintf('%s|%s|%s', e.session, e.paradigm, e.origin), cf_summary_entries, 'UniformOutput', false);
    u_grp_keys = unique(grp_keys, 'stable');
    for gk = 1:numel(u_grp_keys)
        gparts = strsplit(u_grp_keys{gk}, '|');
        g_session  = gparts{1};
        g_paradigm = gparts{2};
        g_origin   = gparts{3};
        band_idx   = find(strcmp(grp_keys, u_grp_keys{gk}));
        par_dir    = fullfile(out_root, g_session, g_paradigm);
        if ~exist(par_dir, 'dir'), mkdir(par_dir); end

        plot_cf_only_summary(cf_summary_entries(band_idx), ref_labels, heat_times, intervals(2,:), g_origin, CLIM_PCTL, false, ...
                              par_dir, sprintf('all_%s_%s_CF_only_summary', g_paradigm, g_origin), n_subj, fig_vis, g_session);

        any_mask = any(arrayfun(@(e) any(e.grp_mask), cf_summary_entries(band_idx)));
        if any_mask
            plot_cf_only_summary(cf_summary_entries(band_idx), ref_labels, heat_times, intervals(2,:), g_origin, CLIM_PCTL, true, ...
                                  par_dir, sprintf('sel_%s_%s_CF_only_summary', g_paradigm, g_origin), n_subj, fig_vis, g_session);
        end
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

function sess = session_of(fpath)
% SESSION_OF  'evaluation' | 'calibration' | 'unknown', inferred from the
%   folder path itself -- topo_erders_channels.mat carries no session-type
%   field of its own. run_subject_analysis.m calls topo_erders() once on
%   the evaluation/ folder and once on the sibling calibration/ folder (per
%   day), so the folder name is the only place this distinction lives.
    parts = strsplit(fpath, filesep);
    if any(strcmpi(parts, 'evaluation'))
        sess = 'evaluation';
    elseif any(strcmpi(parts, 'calibration'))
        sess = 'calibration';
    else
        sess = 'unknown';
    end
end

function s = mask_summary(mask)
    if isempty(mask) || ~any(mask), s = 'no CSP channels';
    else, s = sprintf('%d/%d channels (union across cohort)', sum(mask), numel(mask));
    end
end

function plot_group_topo_grid(subj_c1, subj_c2, subj_has, ch_names, heat_times, intervals, titles_cols, origin, CLIM_PCTL, mask, ...
                               out_dir, prefix, n_contrib, n_subj, freq_lo, freq_hi, fig_vis, session)
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
%   session           'evaluation'|'calibration'|'unknown' -- cosmetic only
%                      (figure Name/title annotation), so an exported SVG
%                      viewed on its own still says which session type it is.
%
%   Each cell (one row x column tile) gets its OWN colour scale, computed
%   from just that cell's own channel values (CLIM_PCTL-th percentile of
%   |value|, see local_cell_clim below) -- NOT one shared scale for the
%   whole image. A shared scale meant a genuinely low-variance cell (e.g. a
%   quiet 1s bin) got rendered on a range sized for whichever OTHER cell in
%   the same image had the biggest swing, crushing its own modest but real
%   contrast down to a single flat colour that looked identical to "no
%   effect at all". Each tile now shows its own full colour range and its
%   own colourbar, so "this tile is saturated red/blue" always means
%   "near this tile's OWN maximum", not "near some other tile's maximum".
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

    h = figure('Name', sprintf('%s — %s %s %g-%g Hz (group, n=%d/%d)', prefix, session, origin, freq_lo, freq_hi, n_contrib, n_subj), ...
               'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
    set(h, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    % Explicit, hand-picked figure position for the tiledlayout itself
    % (normalized units), leaving a fixed top margin for the title and a
    % fixed bottom margin for the dot/ring legend below -- NOT relying on
    % tiledlayout's own automatic title-margin sizing via title(tl,...) +
    % Padding, which repeatedly still let the title text collide with the
    % top row's per-tile titles (e.g. "Cue+CF (0-5s)") regardless of
    % Padding/FontSize. Both the title and the legend are drawn as
    % figure-level annotation() textboxes in their own reserved bands
    % below, so there is no automatic sizing left to get wrong.
    tl = tiledlayout(num_rows, num_intervals, 'TileSpacing', 'compact', 'Padding', 'compact');
    tl.Position = [0.03 0.09 0.94 0.83];

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
                data(data > 0) = 0;   % ERD only (display convention, matches topo_erders.m)
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

            % THIS CELL's own colour scale -- from vals_here only (i.e.
            % after channel masking, so a "sel_" image's scale reflects
            % only the channels actually shown in it, not the full 32).
            if is_mi
                mx_cell = local_prctile(-min(vals_here, 0), CLIM_PCTL);
                if mx_cell == 0 || isnan(mx_cell), mx_cell = 1; end
                clim_cell = [-mx_cell 0];
            else
                mx_cell = local_prctile(abs(vals_here), CLIM_PCTL);
                if mx_cell == 0 || isnan(mx_cell), mx_cell = 1; end
                clim_cell = [-mx_cell mx_cell];
            end

            % jet for MI (ERD-only, one-sided) to visually match the
            % EEGLAB-rendered per-subject topoplots in topo_erders.m, which
            % use EEGLAB's/MATLAB's default colormap rather than the
            % diverging RdBu used elsewhere in this file (e.g. CVSA's
            % lateralization index, still genuinely bidirectional, keeps
            % the default RdBu -- pass [] for it).
            if is_mi, cmap_here = jet(256); else, cmap_here = []; end
            if is_mi, cbar_lbl = '% ERD vs baseline'; else, cbar_lbl = '% ERD/ERS, cue1-cue2'; end
            % show_cbar = true on EVERY tile now (not just the last column):
            % since every tile has its own independent scale, a single
            % shared colourbar at the row's end would only describe that
            % last tile, not the others.
            topo_map(names_here, vals_here, clim_cell, ax, '', true, bg_zero, sig_here, cmap_here, cbar_lbl);
            if r == 1, title(ax, titles_cols{c}, 'FontSize', 9); end
            if c == 1, ylabel(ax, row_labels{r}, 'Visible', 'on', 'FontWeight', 'bold'); end
        end
    end

    % Title and dot/ring legend as figure-level annotation() textboxes, each
    % pinned to its OWN reserved band (see tl.Position above) -- completely
    % independent of the tiledlayout's internal spacing, so neither can ever
    % collide with the top row's per-tile titles no matter the text length.
    annotation(h, 'textbox', [0 0.94 1 0.06], 'String', ...
        sprintf('%s — [%s] %s %g-%g Hz  (n=%d/%d subjects)', prefix, upper(session), origin, freq_lo, freq_hi, n_contrib, n_subj), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');

    annotation(h, 'textbox', [0 0 1 0.05], 'String', ...
        'dot = value;  black ring = statistically significant (p<0.05, one-sample t-test across subjects)', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 9);

    save_path = fullfile(out_dir, sprintf('%s.svg', prefix));
    saveas(h, save_path, 'svg');
    if ~strcmp(fig_vis, 'on'), close(h); end
    fprintf('    Saved: %s\n', save_path);
end

function plot_cf_only_summary(band_entries, ch_names, heat_times, cf_interval, origin, CLIM_PCTL, use_mask, ...
                               out_dir, prefix, n_subj, fig_vis, session)
% PLOT_CF_ONLY_SUMMARY  One figure: rows = cue(s) (or the lateralization
%   row for CVSA), COLUMNS = BAND (not time interval) -- the cross-subject
%   "CF only" topoplot (the whole-feedback-period average -- same cell
%   that is column 2 of plot_group_topo_grid's 7-column grid) for every
%   band found in this (session, paradigm, origin) group, side by side in
%   one image, so you don't have to open N separate per-band figures and
%   pick out column 2 from each.
%
%   band_entries      struct array (one per band) with fields subj_c1/
%                      subj_c2 [n_ch x n_time x n_subj] (NOT pre-averaged
%                      across subjects -- same per-subject arrays already
%                      computed by the main loop in group_topo_erders.m,
%                      passed in here rather than recomputed), subj_has
%                      [1 x n_subj], grp_mask [n_ch x 1], freq_lo/freq_hi,
%                      n_contrib.
%   use_mask           if true, restrict every panel to the channels in
%                      that band's own grp_mask (union across all bands in
%                      this call); if false, show every channel.
%
%   Same per-cell percentile colour scale (CLIM_PCTL-th percentile of
%   |value|) and the same cross-subject one-sample t-test significance
%   ring (p<0.05) as plot_group_topo_grid -- each band's panel gets its OWN
%   colour scale, not one shared across the whole image, for the same
%   reason documented there.
    is_mi = strcmpi(origin, 'MI');
    if is_mi
        num_rows = 2;
        row_labels = {'Cue 1 (ERD)', 'Cue 2 (ERD)'};
    else
        num_rows = 1;
        row_labels = {'Cue 1 - Cue 2 (lateralization)'};
    end
    n_bands = numel(band_entries);

    combined_mask = false(numel(ch_names), 1);
    if use_mask
        for k = 1:n_bands
            combined_mask = combined_mask | band_entries(k).grp_mask(:);
        end
    end
    bg_zero = use_mask;

    t_idx = heat_times >= cf_interval(1)*1000 & heat_times < cf_interval(2)*1000;

    h = figure('Name', sprintf('%s — %s %s CF-only summary (group)', prefix, session, origin), ...
               'Color', 'w', 'NumberTitle', 'off', 'Visible', fig_vis);
    set(h, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    tl = tiledlayout(num_rows, n_bands, 'TileSpacing', 'compact', 'Padding', 'compact');
    tl.Position = [0.03 0.09 0.94 0.83];

    for r = 1:num_rows
        for bk = 1:n_bands
            ax = nexttile;
            be = band_entries(bk);

            subj_vals = nan(numel(ch_names), numel(be.subj_has));
            for si = find(be.subj_has)
                if is_mi
                    src = be.subj_c1(:,:,si); if r == 2, src = be.subj_c2(:,:,si); end
                    subj_vals(:,si) = mean(src(:, t_idx), 2, 'omitnan');
                else
                    d1 = mean(be.subj_c1(:, t_idx, si), 2, 'omitnan');
                    d2 = mean(be.subj_c2(:, t_idx, si), 2, 'omitnan');
                    subj_vals(:,si) = d1 - d2;
                end
            end
            data = mean(subj_vals, 2, 'omitnan');
            if is_mi, data(data > 0) = 0; end
            data(isnan(data)) = 0;

            sig_chan = local_ttest_pvals(subj_vals') < 0.05;

            if use_mask
                names_here = ch_names(combined_mask);
                vals_here  = data(combined_mask);
                sig_here   = sig_chan(combined_mask);
            else
                names_here = ch_names;
                vals_here  = data;
                sig_here   = sig_chan;
            end

            if is_mi
                mx_cell = local_prctile(-min(vals_here, 0), CLIM_PCTL);
                if mx_cell == 0 || isnan(mx_cell), mx_cell = 1; end
                clim_cell = [-mx_cell 0];
                cmap_here = jet(256);   % one-sided ERD-only, see topo_erders.m/plot_group_topo_grid
            else
                mx_cell = local_prctile(abs(vals_here), CLIM_PCTL);
                if mx_cell == 0 || isnan(mx_cell), mx_cell = 1; end
                clim_cell = [-mx_cell mx_cell];
                cmap_here = [];   % topo_map.m default: diverging RdBu, centred at 0
            end
            if is_mi, cbar_lbl = '% ERD vs baseline'; else, cbar_lbl = '% ERD/ERS, cue1-cue2'; end

            topo_map(names_here, vals_here, clim_cell, ax, '', true, bg_zero, sig_here, cmap_here, cbar_lbl);
            if r == 1
                title(ax, sprintf('%g-%g Hz (n=%d/%d)', be.freq_lo, be.freq_hi, be.n_contrib, n_subj), 'FontSize', 9);
            end
            if bk == 1, ylabel(ax, row_labels{r}, 'Visible', 'on', 'FontWeight', 'bold'); end
        end
    end

    annotation(h, 'textbox', [0 0.94 1 0.06], 'String', ...
        sprintf('%s — [%s] %s — "CF only" average, all bands (n subjects per band shown in each panel title)', ...
                prefix, upper(session), origin), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');

    annotation(h, 'textbox', [0 0 1 0.05], 'String', ...
        'dot = value;  black ring = statistically significant (p<0.05, one-sample t-test across subjects)', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 9);

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
