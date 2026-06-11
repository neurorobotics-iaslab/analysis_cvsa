%% TOPO_ERDERS  ERD/ERS topoplots (EEGLAB newtimef) — total / per-file / 1s-steps
%
%   Loads one or more evaluation GDF files (multi-select; files may come
%   from DIFFERENT paradigms — MI/CVSA/Hybrid — as long as they share the
%   same channel montage), preprocesses each with EEGLAB (resample 128 Hz,
%   1-40 Hz bandpass, epoch around the class-onset events, CAR), then
%   computes per-channel ERD/ERS (% power change vs. pre-stim baseline) via
%   newtimef, exactly as in eeglab_gdf.m.
%
%   Each file's class-onset event codes, paradigm name, and CSP channels +
%   FREQUENCY BANDS are detected from its OWN companion rosparam YAML
%   (paradigm-specific event codes, e.g. MI=769/770, CVSA=730/731,
%   Hybrid=750/751).
%
%   FREQUENCY BANDS + DISPLAY MODE per paradigm, taken from the CSP YAML
%   (processing_fbcsp_{mi,cvsa}.CspCfg.params.bands), NOT a fixed list:
%     - paradigm "mi"     -> MI CSP bands,   "MI"   display (see below)
%     - paradigm "cvsa"   -> CVSA CSP bands, "CVSA" display
%     - paradigm "hybrid" -> BOTH: MI CSP bands with "MI" display AND
%                            CVSA CSP bands with "CVSA" display (both
%                            classifiers run together in hybrid mode)
%   If no CSP info is found for a paradigm, DEFAULT_BANDS is used with
%   "MI" display as a fallback.
%
%   Display modes (one figure per band):
%     - "MI"   : ERD only — two rows, one per cue (class), values clipped
%                to <= 0 (only event-related DESYNCHRONISATION is shown).
%     - "CVSA" : Lateralization-index style — a single row showing the
%                difference (cue1 - cue2) of the per-class ERD/ERS maps.
%   Row labels use the actual cue/event codes for that paradigm.
%
%   Output is organised in one sub-folder PER PARADIGM under
%   <gdf_folder>/results_eeglab/<paradigm>/:
%     - "all_<paradigm>_<MI|CVSA>_<band>"  : grand average over all files
%                                            of that paradigm, all channels
%     - "sel_<paradigm>_<MI|CVSA>_<band>"  : same, but with all non-CSP
%                                            electrodes forced to 0 — MI
%                                            bands use the MI CSP channel
%                                            set, CVSA bands use the CVSA
%                                            CSP channel set (no mixing)
%   and a "per_file" sub-folder inside <paradigm>/ with, for each
%   individual file of that paradigm, "all_<filename>_..." /
%   "sel_<filename>_..." using the same colour scale and bands as that
%   file's paradigm.

clear; clc; close all;

%% --- 1. SETUP ----------------------------------------------------------
eeglab_repo = '/home/paolo/Local/Matlab/eeglab';
if ~exist('eeglab', 'file'), addpath(eeglab_repo); end
[ALLEEG, EEG, CURRENTSET] = eeglab('nogui'); %#ok<ASGLU>

this_dir = fileparts(mfilename('fullpath'));
ms_dir   = fullfile(this_dir, '..', 'matlab_simulation');
addpath(ms_dir, fullfile(ms_dir,'io'), fullfile(ms_dir,'utils'));

%% --- 2. PARAMETERS ------------------------------------------------------
DEFAULT_BANDS = [4 8; 8 14; 14 24; 14 30];   % fallback if a paradigm has no CSP info
epoch_limits  = [-2 5];   % [s] relative to class-onset event; CF assumed in [0, epoch_limits(2)]
CUE_DURATION_S = 1.5;     % [s] cue length (Training.cpp duration/cue=1500ms), immediately before CF onset
baseline_ms = [-CUE_DURATION_S*1000, 0];   % baseline window = cue only (excludes earlier fixation)

T = epoch_limits(2);
intervals   = [0 T; 1 T];
titles_cols = {sprintf('Cue+CF (0-%ds)', T), sprintf('CF only (1-%ds)', T)};
for s = 0:T-1
    intervals(end+1, :)  = [s, s+1];          %#ok<AGROW>
    titles_cols{end+1}   = sprintf('%d-%ds', s, s+1); %#ok<AGROW>
end
num_intervals = size(intervals, 1);

%% --- 3. FILE PICKER ------------------------------------------------------
[files, folder] = uigetfile('*.gdf', 'Seleziona uno o più file GDF', 'MultiSelect', 'on');
if isequal(files,0), disp('Operazione annullata'); return; end
if ischar(files), files = {files}; end
n_files = numel(files);

output_dir = fullfile(folder, 'results_eeglab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% --- 4. PRE-SCAN: paradigm, event codes, CSP channels + bands per file ----
file_paradigm    = cell(1, n_files);   % {f}: 'mi'|'cvsa'|'hybrid'|...
file_event_types = cell(1, n_files);   % {f}: {'<low>','<high>'} cue codes
file_csp_info    = cell(1, n_files);   % {f}: struct array('name','MI'|'CVSA','channels',{cellstr},'bands',[Nx2])

for f = 1:n_files
    [~, basename] = fileparts(files{f});
    full_path = fullfile(folder, files{f});
    file_csp_info{f} = struct('name', {}, 'channels', {}, 'bands', {});

    try
        [params_f, ~] = load_params_yaml(full_path);
        classes_v = sort(to_vec(params_f.integrator.classes));
        file_event_types{f} = {num2str(classes_v(1)), num2str(classes_v(2))};
        file_paradigm{f}    = params_f.integrator.paradigm;
        fprintf('[%d/%d] %s: event types %s/%s (paradigm=%s)\n', f, n_files, basename, ...
                file_event_types{f}{1}, file_event_types{f}{2}, file_paradigm{f});

        if isfield(params_f, 'processing_fbcsp_mi')
            csp_mi = load_csp(params_f, 'mi');
            file_csp_info{f}(end+1) = struct('name', 'MI', 'channels', {csp_mi.selected_channels}, 'bands', csp_mi.bands);
        end
        if isfield(params_f, 'processing_fbcsp_cvsa')
            csp_cvsa = load_csp(params_f, 'cvsa');
            file_csp_info{f}(end+1) = struct('name', 'CVSA', 'channels', {csp_cvsa.selected_channels}, 'bands', csp_cvsa.bands);
        end
    catch ME
        warning('topo_erders:noyaml', ...
            'Companion YAML not usable for %s (%s). Falling back to filename-based event/paradigm detection; CSP info skipped for this file.', basename, ME.message);
    end

    if isempty(file_event_types{f})
        if any(contains(lower(basename), {'cvsa_blbr','cvsa_lbrb'}))
            file_event_types{f} = {'730','731'};
            file_paradigm{f}    = 'cvsa';
        elseif contains(lower(basename), 'hybrid')
            file_event_types{f} = {'750','751'};
            file_paradigm{f}    = 'hybrid';
        elseif contains(lower(basename), 'mi')
            file_event_types{f} = {'769','770'};
            file_paradigm{f}    = 'mi';
        else
            file_event_types{f} = {'771','773'};
            file_paradigm{f}    = 'unknown';
            warning('Pattern non riconosciuto per %s. Uso default: 771, 773.', basename);
        end
        fprintf('[%d/%d] %s: event types %s/%s (paradigm=%s, filename heuristic)\n', f, n_files, basename, ...
                file_event_types{f}{1}, file_event_types{f}{2}, file_paradigm{f});
    end
end

%% --- 4b. PER-PARADIGM BANDS (+ display mode) AND GLOBAL FREQ RANGE ---------
paradigms = unique(file_paradigm, 'stable');   % e.g. {'mi','hybrid','cvsa'}
n_par     = numel(paradigms);

par_bands       = cell(1, n_par);   % {p}: [Nx2] bands used for paradigm p
par_band_names  = cell(1, n_par);   % {p}: {N} band name strings (include display mode)
par_band_origin = cell(1, n_par);   % {p}: {N} 'MI'|'CVSA' display mode per band
par_event_types = cell(1, n_par);   % {p}: {'<low>','<high>'} cue codes for that paradigm

for p = 1:n_par
    idx_p = find(strcmp(file_paradigm, paradigms{p}));
    par_event_types{p} = file_event_types{idx_p(1)};

    if strcmpi(paradigms{p}, 'hybrid')
        bnd_mi   = collect_bands(file_csp_info, idx_p, 'MI');
        bnd_cvsa = collect_bands(file_csp_info, idx_p, 'CVSA');
        if isempty(bnd_mi) && isempty(bnd_cvsa)
            bnd_mi = DEFAULT_BANDS;
            warning('topo_erders:nobands', 'No CSP bands found for paradigm "%s" — using default bands (MI display).', paradigms{p});
        end
        bnd    = [bnd_mi; bnd_cvsa];
        origin = [repmat({'MI'}, size(bnd_mi,1), 1); repmat({'CVSA'}, size(bnd_cvsa,1), 1)];
    else
        variant_name = upper(paradigms{p});   % 'MI' or 'CVSA'
        bnd = collect_bands(file_csp_info, idx_p, variant_name);
        if isempty(bnd)
            bnd = DEFAULT_BANDS;
            warning('topo_erders:nobands', 'No CSP bands found for paradigm "%s" — using default bands.', paradigms{p});
            variant_name = 'MI';   % default display mode
        end
        origin = repmat({variant_name}, size(bnd,1), 1);
    end

    par_bands{p}       = bnd;
    par_band_origin{p} = origin;
    par_band_names{p}  = arrayfun(@(i) sprintf('%s_%g_%g', origin{i}, bnd(i,1), bnd(i,2)), ...
                                   1:size(bnd,1), 'UniformOutput', false);
    fprintf('Paradigm "%s": %d band(s): %s\n', paradigms{p}, size(bnd,1), strjoin(par_band_names{p}, ', '));
end

all_bands_global = vertcat(par_bands{:});
freq_range = [max(1, floor(min(all_bands_global(:,1))) - 1), ceil(max(all_bands_global(:,2))) + 1];
fprintf('Global newtimef frequency range: [%g %g] Hz\n', freq_range(1), freq_range(2));

%% --- 5. PER-FILE LOAD + PREPROCESS + EPOCH + ERSP ------------------------
all_p_db   = cell(1, n_files);   % {f}: [nchan x 2 x nfreq x ntime]
file_label = cell(1, n_files);
file_csp   = cell(1, n_files);   % {f}: struct array('name','MI'|'CVSA','mask',logical), per-file CSP channel masks
ref_labels = {};
times = []; freqs = [];
chanlocs_ref = [];

for f = 1:n_files
    [~, basename] = fileparts(files{f});
    file_label{f} = basename;
    fprintf('\n[%d/%d] %s\n', f, n_files, basename);
    full_path = fullfile(folder, files{f});

    % ── Load + build EEG struct ─────────────────────────────────────────
    [s, h] = sload(full_path);
    total_chans  = size(s, 2);
    eeg_chans_idx = 1:(total_chans - 3);
    s_eeg  = s(:, eeg_chans_idx);
    labels_eeg = h.Label(eeg_chans_idx);

    EEG = pop_importdata('dataformat','array', 'nbchan', length(eeg_chans_idx), ...
        'data', s_eeg', 'srate', h.SampleRate, 'setname', basename);
    EEG.chanlocs = struct('labels', labels_eeg);
    EEG = pop_chanedit(EEG, 'lookup', 'standard_1005.elc');

    if isfield(h, 'EVENT')
        for i = 1:length(h.EVENT.POS)
            EEG.event(i).type     = num2str(h.EVENT.TYP(i));
            EEG.event(i).latency  = h.EVENT.POS(i);
            EEG.event(i).duration = h.EVENT.DUR(i);
        end
    end
    EEG = eeg_checkset(EEG);

    if f == 1
        ref_labels = {EEG.chanlocs.labels};
    elseif numel(EEG.chanlocs) ~= numel(ref_labels) || ~isequal({EEG.chanlocs.labels}, ref_labels)
        error('topo_erders:chanmismatch', ...
              'File %s has different channels than %s — select files from the same montage/session.', ...
              basename, files{1});
    end

    % ── CSP channel masks for this file (from pre-scan) ──────────────────
    file_csp{f} = struct('name', {}, 'mask', {});
    for v = 1:numel(file_csp_info{f})
        file_csp{f}(end+1) = struct('name', file_csp_info{f}(v).name, ...
                                     'mask', make_csp_mask(file_csp_info{f}(v).channels, ref_labels));
    end

    event_types_f = file_event_types{f};

    % ── Resample + bandpass ─────────────────────────────────────────────
    EEG = pop_resample(EEG, 128);
    EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'hicutoff', 40);

    % ── Epoch around class onset + baseline removal ─────────────────────
    EEG_ep = pop_epoch(EEG, event_types_f, epoch_limits, 'newname', [basename '_epoched'], 'epochinfo', 'yes');
    EEG_ep = pop_rmbase(EEG_ep, baseline_ms);

    % ── CAR (matches Car.hpp: average over non-EOG channels, subtracted ──
    %        from ALL channels including the excluded ones) ──────────────
    car_excl_lbl = {'Fp1','Fp2'};
    car_all_lbl  = lower(strtrim({EEG_ep.chanlocs.labels}));
    car_excl_idx = find(ismember(car_all_lbl, lower(car_excl_lbl)));
    car_incl_idx = setdiff(1:EEG_ep.nbchan, car_excl_idx);
    car_ref      = mean(EEG_ep.data(car_incl_idx,:,:), 1);
    EEG_ep.data  = EEG_ep.data - car_ref;

    % ── Class indices ────────────────────────────────────────────────────
    idx_cls = {[], []};
    for i = 1:length(EEG_ep.epoch)
        ev_latencies = [EEG_ep.epoch(i).eventlatency{:}];
        zero_idx = find(ev_latencies == 0, 1);
        if isempty(zero_idx), continue; end
        t = EEG_ep.epoch(i).eventtype;
        if iscell(t), t = t{zero_idx}; end
        if any([isequal(t, event_types_f{1}), isequal(t, str2double(event_types_f{1}))])
            idx_cls{1} = [idx_cls{1}, i];
        elseif any([isequal(t, event_types_f{2}), isequal(t, str2double(event_types_f{2}))])
            idx_cls{2} = [idx_cls{2}, i];
        end
    end
    fprintf('  %d epochs total (cls %s: %d, cls %s: %d)\n', ...
            length(EEG_ep.epoch), event_types_f{1}, numel(idx_cls{1}), event_types_f{2}, numel(idx_cls{2}));

    % ── ERSP per channel/class ───────────────────────────────────────────
    p_db = nan(EEG_ep.nbchan, 2, 0, 0);
    for ch = 1:EEG_ep.nbchan
        for c = 1:2
            [ersp, ~, ~, times_f, freqs_f] = newtimef(EEG_ep.data(ch, :, idx_cls{c}), ...
                EEG_ep.pnts, [EEG_ep.times(1) EEG_ep.times(end)], EEG_ep.srate, 0, ...
                'baseline', baseline_ms, 'freqs', freq_range, 'winsize', 128, ...
                'plotersp', 'off', 'plotitc', 'off', 'verbose', 'off');
            if ch == 1 && c == 1
                times = times_f; freqs = freqs_f;
                p_db = nan(EEG_ep.nbchan, 2, numel(freqs), numel(times));
            end
            p_db(ch, c, :, :) = (10.^(ersp/10) - 1) * 100;   % % ERD(-)/ERS(+)
        end
    end
    all_p_db{f} = p_db;

    if f == 1, chanlocs_ref = EEG_ep.chanlocs; end
    fprintf('  ERSP computed: %d channels x %d freqs x %d times\n', EEG_ep.nbchan, numel(freqs), numel(times));
end

%% --- 6. GROUP FILES BY PARADIGM + PER-PARADIGM AVERAGES --------------------
par_p_db      = cell(1, n_par);   % {p}: [nchan x 2 x nfreq x ntime] grand average for that paradigm
par_mi_mask   = cell(1, n_par);   % {p}: union MI CSP mask across that paradigm's files ([] if none)
par_cvsa_mask = cell(1, n_par);   % {p}: union CVSA CSP mask across that paradigm's files ([] if none)
par_lims      = cell(1, n_par);   % {p}: band_lims cell (one per band of par_bands{p})

for p = 1:n_par
    idx_p = find(strcmp(file_paradigm, paradigms{p}));
    par_p_db{p} = mean(cat(5, all_p_db{idx_p}), 5);

    par_mi_mask{p}   = union_mask_for(file_csp, idx_p, 'MI',   numel(ref_labels));
    par_cvsa_mask{p} = union_mask_for(file_csp, idx_p, 'CVSA', numel(ref_labels));

    bnd    = par_bands{p};
    origin = par_band_origin{p};
    lims = cell(size(bnd,1), 1);
    for idx_b = 1:size(bnd,1)
        freq_idx = freqs >= bnd(idx_b,1) & freqs <= bnd(idx_b,2);
        pw = squeeze(mean(par_p_db{p}(:,:,freq_idx,:), 3));   % [nchan x 2 x ntime]
        if strcmpi(origin{idx_b}, 'MI')
            neg = pw; neg(neg > 0) = 0;
            mx = max(-neg(:));
            if mx == 0 || isnan(mx), mx = 1; end
            lims{idx_b} = [-mx 0];
        else   % CVSA: lateralization index = class1 - class2
            d = pw(:,1,:) - pw(:,2,:);
            mx = max(abs(d(:)));
            if mx == 0 || isnan(mx), mx = 1; end
            lims{idx_b} = [-mx mx];
        end
    end
    par_lims{p} = lims;

    fprintf('Paradigm "%s": %d file(s), MI CSP mask: %s, CVSA CSP mask: %s\n', paradigms{p}, numel(idx_p), ...
            mask_summary(par_mi_mask{p}), mask_summary(par_cvsa_mask{p}));
end

%% --- 7. GENERATE TOPOPLOT GRIDS ---------------------------------------------
fprintf('\n--- GENERATING TOPOPLOTS ---\n');

for p = 1:n_par
    par_dir       = fullfile(output_dir, paradigms{p});
    par_dir_files = fullfile(par_dir, 'per_file');
    if ~exist(par_dir, 'dir'), mkdir(par_dir); end
    if ~exist(par_dir_files, 'dir'), mkdir(par_dir_files); end

    bnd    = par_bands{p};
    bnames = par_band_names{p};
    origin = par_band_origin{p};
    lims   = par_lims{p};
    evt    = par_event_types{p};

    masks_none = cell(size(bnd,1), 1);   % all-channel figures: no masking
    masks_sel  = band_masks_for(origin, par_mi_mask{p}, par_cvsa_mask{p});

    % ── Per-paradigm grand average ───────────────────────────────────────
    plot_topo_grid(par_p_db{p}, chanlocs_ref, freqs, times, bnd, bnames, origin, ...
                   lims, intervals, titles_cols, evt, masks_none, par_dir, sprintf('all_%s', paradigms{p}));
    if any(~cellfun(@isempty, masks_sel))
        plot_topo_grid(par_p_db{p}, chanlocs_ref, freqs, times, bnd, bnames, origin, ...
                       lims, intervals, titles_cols, evt, masks_sel, par_dir, sprintf('sel_%s', paradigms{p}));
    end

    % ── Per-file (same colour scale + bands as this paradigm) ────────────
    idx_p = find(strcmp(file_paradigm, paradigms{p}));
    for f = idx_p
        plot_topo_grid(all_p_db{f}, chanlocs_ref, freqs, times, bnd, bnames, origin, ...
                       lims, intervals, titles_cols, evt, masks_none, par_dir_files, sprintf('all_%s', file_label{f}));

        file_mi_mask   = union_mask_for(file_csp, f, 'MI',   numel(ref_labels));
        file_cvsa_mask = union_mask_for(file_csp, f, 'CVSA', numel(ref_labels));
        file_masks_sel = band_masks_for(origin, file_mi_mask, file_cvsa_mask);
        if any(~cellfun(@isempty, file_masks_sel))
            plot_topo_grid(all_p_db{f}, chanlocs_ref, freqs, times, bnd, bnames, origin, ...
                           lims, intervals, titles_cols, evt, file_masks_sel, par_dir_files, sprintf('sel_%s', file_label{f}));
        end
    end
end

fprintf('\nDone. Figures saved under: %s (one sub-folder per paradigm)\n', output_dir);

%% ── Local helpers (must be after all script statements) ────────────────────
function bnd = collect_bands(file_csp_info, idx_files, name)
% COLLECT_BANDS  Unique [Nx2] CSP bands of the given variant ('MI'|'CVSA')
%                across the given files. Empty if none found.
    bnd = [];
    for f = idx_files
        for v = 1:numel(file_csp_info{f})
            if strcmpi(file_csp_info{f}(v).name, name)
                bnd = [bnd; file_csp_info{f}(v).bands]; %#ok<AGROW>
            end
        end
    end
    if ~isempty(bnd), bnd = unique(bnd, 'rows', 'stable'); end
end

function mask = union_mask_for(file_csp, idx_files, name, n_ch)
% UNION_MASK_FOR  Union of CSP channel masks of the given variant
%                 ('MI'|'CVSA') across the given files. [] if none found.
    mask = false(n_ch, 1);
    found = false;
    for f = idx_files
        for v = 1:numel(file_csp{f})
            if strcmpi(file_csp{f}(v).name, name)
                mask = mask | file_csp{f}(v).mask;
                found = true;
            end
        end
    end
    if ~found, mask = []; end
end

function masks = band_masks_for(band_origin, mi_mask, cvsa_mask)
% BAND_MASKS_FOR  Per-band mask cell: MI-origin bands get mi_mask,
%                 CVSA-origin bands get cvsa_mask (either may be []).
    masks = cell(numel(band_origin), 1);
    for i = 1:numel(band_origin)
        if strcmpi(band_origin{i}, 'MI')
            masks{i} = mi_mask;
        else
            masks{i} = cvsa_mask;
        end
    end
end

function s = mask_summary(mask)
    if isempty(mask), s = 'none'; else, s = sprintf('%d/%d channels', sum(mask), numel(mask)); end
end

function mask = make_csp_mask(channels, ref_labels)
% MAKE_CSP_MASK  Logical column mask over ref_labels marking CSP-selected channels.
    sel     = lower(strtrim(channels));
    all_lbl = lower(strtrim(ref_labels));
    mask    = ismember(all_lbl, sel)';
end

function plot_topo_grid(all_p_db, chanlocs, freqs, times, bands, band_names, band_origin, ...
                         band_lims, intervals, titles_cols, event_types_p, band_masks, out_dir, prefix)
% PLOT_TOPO_GRID  One figure per band: cols = time interval.
%   all_p_db      : [nchan x 2 x nfreq x ntime] percent ERD(-)/ERS(+)
%   bands         : [Nx2] frequency bands (Hz) to plot for this call
%   band_origin   : {N} 'MI' (ERD-only, 2 rows = cues) | 'CVSA' (lateralization
%                   index cue1-cue2, 1 row)
%   band_lims     : {N} colour-axis limits per band
%   event_types_p : {2} cue/event codes for this paradigm, used as row labels
%   band_masks    : {N} of [] or logical [nchan x 1] — channels not in mask
%                   are forced to 0

    num_intervals = size(intervals, 1);

    for idx_b = 1:size(bands,1)
        low_f = bands(idx_b,1); high_f = bands(idx_b,2);
        freq_idx = freqs >= low_f & freqs <= high_f;

        pow_band = squeeze(mean(all_p_db(:,:,freq_idx,:), 3));   % [nchan x 2 x ntime]
        mask     = band_masks{idx_b};
        is_mi    = strcmpi(band_origin{idx_b}, 'MI');

        if is_mi
            num_rows = 2;
            row_labels = {sprintf('Cue %s', event_types_p{1}), sprintf('Cue %s', event_types_p{2})};
        else
            num_rows = 1;
            row_labels = {sprintf('Cue %s - Cue %s', event_types_p{1}, event_types_p{2})};
        end

        h = figure('Name', sprintf('%s — %s %g-%g Hz', prefix, band_origin{idx_b}, low_f, high_f), ...
                    'Color','w', 'NumberTitle','off', 'Visible','off');
        set(h, 'Units','normalized', 'OuterPosition', [0 0 1 1]);
        tiledlayout(num_rows, num_intervals, 'TileSpacing','compact', 'Padding','tight');

        for r = 1:num_rows
            for c = 1:num_intervals
                nexttile;
                t_idx = times >= intervals(c,1)*1000 & times < intervals(c,2)*1000;

                if is_mi
                    data = mean(pow_band(:, r, t_idx), 3);
                    data = data(:);
                    data(data > 0) = 0;   % ERD only
                else
                    d1 = mean(pow_band(:, 1, t_idx), 3);
                    d2 = mean(pow_band(:, 2, t_idx), 3);
                    data = d1(:) - d2(:);
                end
                if ~isempty(mask), data = data .* mask(:); end

                topoplot(data, chanlocs, 'style','both', 'maplimits', band_lims{idx_b}, ...
                         'electrodes','labels', 'whitebk','on');
                if r == 1, title(titles_cols{c}, 'FontSize', 9); end
                if c == 1
                    ylabel(row_labels{r}, 'Visible','on', 'FontWeight','bold');
                end
                if c == num_intervals, colorbar; end
            end
        end

        sgtitle(sprintf('%s — %s %g-%g Hz', prefix, band_origin{idx_b}, low_f, high_f), ...
                'Interpreter','none', 'FontWeight','bold');

        save_path = fullfile(out_dir, sprintf('topo_%s_%s.png', prefix, band_names{idx_b}));
        saveas(h, save_path);
        close(h);
        fprintf('  Saved: %s\n', save_path);
    end
end
