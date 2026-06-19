%% TOPO_ERDERS  ERD/ERS topoplots — total / per-file / 1s-steps
%
%   Loads one or more evaluation GDF files (multi-select; files may come
%   from DIFFERENT paradigms — MI/CVSA/Hybrid — as long as they share the
%   same channel montage), preprocesses each with EEGLAB (resample 128 Hz,
%   1-40 Hz bandpass, epoch around the class-onset events, CAR), then
%   computes per-channel, per-trial ERD/ERS (Pfurtscheller-style band power,
%   % change vs. pre-stim baseline; see "CHANNEL x TIME ERD/ERS HEATMAPS"
%   below). Topoplots are time-averages of this SAME per-trial,
%   NaN-masked band power — see "VARIABLE CF DURATION" below — so topoplots
%   and heatmaps are always consistent with each other.
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
%   <gdf_folder>/analysis_results/topo_erders/<paradigm>/:
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
%
%   CHANNEL x TIME ERD/ERS HEATMAPS ("heat_*"), one per band, same folders
%   as the topoplots ("all_<paradigm>"/"sel_<paradigm>" + per_file):
%     - Y axis = channels, X axis = time from cue onset (t=0) to the end
%       of the longest possible continuous-feedback (CF) window
%       (epoch_limits(2)), one row per cue/class (no lateralization-index
%       collapsing here — both classes are always shown separately).
%     - Colour = % ERD(-, blue)/ERS(+, red), diverging colormap centred at
%       0, with symmetric limits = max(abs(value)) over the data actually
%       displayed in that figure.
%     - "all_*" shows every channel; "sel_*" shows only the CSP-selected
%       channels (MI/CVSA paradigms: that paradigm's own CSP channels;
%       Hybrid: union of MI + CVSA CSP channels, for every band).
%     - VARIABLE CF DURATION (evaluation mode): per trial, ERD/ERS is a
%       Pfurtscheller-style band-power percentage relative to that trial's
%       own pre-cue baseline; samples after that trial's actual CF end
%       (taken from the 897/898/899 outcome event latency) are set to NaN
%       before averaging across trials, so the heatmap at each time point
%       is the mean over only the trials whose CF is still running at that
%       time. NaN (no-data) regions are rendered as the figure background.
%       The topoplots ("topo_*") use this same NaN-masked per-trial band
%       power, averaged over each time interval (omitting NaNs); an
%       interval/channel where NO trial of that class still has CF running
%       (e.g. a 4-5s column when every trial ended earlier) is shown as 0.
%       Topoplot colour limits are computed over the same data, restricted
%       to [CUE_DURATION_S, epoch_limits(2)] (CF onset -> max CF end, i.e.
%       the CF window only).
%
%   ERD/ERS vs. CSP WEIGHT ("corr_<paradigm>"), one figure per paradigm,
%   one panel per band that has CSP info: scatter of each channel's class
%   discrimination |ERD/ERS(class1) - ERD/ERS(class2)| during continuous
%   feedback (cue+CUE_DURATION_S -> max CF end, from the same NaN-masked
%   band power) against that channel's CSP weight (sum of |CSP filter
%   coefficients| over components for that band; 0 for channels not
%   selected by the CSP). |c1-c2| (not the per-class average) is used for
%   BOTH MI and CVSA bands: it captures contralateral ERD/ERS patterns
%   (e.g. C4 ERD for class1 / C3 ERD for class2 both appear as large
%   |c1-c2|, whereas averaging the classes would cancel them out) and
%   matches the lateralization index already used for CVSA topoplots.
%   MI-origin bands use the MI CSP, CVSA-origin bands use the CVSA CSP
%   (hybrid: both, one panel each).
%   Pearson correlation coefficient is shown in each panel's title.
%
%   Also saves topo_erders_summary.mat under analysis_results/topo_erders/
%   (one row per paradigm x band: mean ERD/ERS per class over CSP-selected
%   channels during CF, their discrimination, and the CSP-weight Pearson r)
%   for cross-subject aggregation by main_group_analysis.m.

clear; clc; close all;

% --- Display options -------------------------------------------------------
SHOW_FIGURES = false;   % true: figures pop up on screen; false: created hidden (export only)
if SHOW_FIGURES, fig_vis = 'on'; else, fig_vis = 'off'; end

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
intervals   = [0 T; CUE_DURATION_S T];
titles_cols = {sprintf('Cue+CF (0-%ds)', T), sprintf('CF only (%g-%ds)', CUE_DURATION_S, T)};

% Per-column breakdown, CF-aligned: one Cue bin [0, CUE_DURATION_S], then
% 1s-wide CF bins starting at CUE_DURATION_S (last bin shorter if T-CUE_DURATION_S
% is not a whole number).
intervals(end+1, :) = [0, CUE_DURATION_S];
titles_cols{end+1}  = sprintf('Cue (0-%gs)', CUE_DURATION_S);
s = CUE_DURATION_S;
while s < T
    e = min(s + 1, T);
    intervals(end+1, :) = [s, e];          %#ok<AGROW>
    titles_cols{end+1}  = sprintf('CF (%g-%gs)', s, e); %#ok<AGROW>
    s = e;
end
num_intervals = size(intervals, 1);

% ── Channel x time ERD/ERS heatmaps ──────────────────────────────────────
OUTCOME_EVENT_TYPES = {'897','898','899'};      % HIT/MISS/TIMEOUT -> per-trial CF end
ERD_SMOOTH_MS       = 200;                      % [ms] moving-average window for band power
heat_time_range     = [0, epoch_limits(2)*1000]; % [ms] cue onset -> max possible CF end

%% --- 3. FILE PICKER ------------------------------------------------------
[files, folder] = uigetfile('*.gdf', 'Seleziona uno o più file GDF', 'MultiSelect', 'on');
if isequal(files,0), disp('Operazione annullata'); return; end
if ischar(files), files = {files}; end
n_files = numel(files);

output_dir = fullfile(folder, 'analysis_results', 'topo_erders');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% --- 4. PRE-SCAN: paradigm, event codes, CSP channels + bands per file ----
file_paradigm    = cell(1, n_files);   % {f}: 'mi'|'cvsa'|'hybrid'|...
file_event_types = cell(1, n_files);   % {f}: {'<low>','<high>'} cue codes
file_csp_info    = cell(1, n_files);   % {f}: struct array('name','MI'|'CVSA','channels',{cellstr},'bands',[Nx2])

for f = 1:n_files
    [~, basename] = fileparts(files{f});
    full_path = fullfile(folder, files{f});
    file_csp_info{f} = struct('name', {}, 'channels', {}, 'bands', {}, 'csp_matrices', {});

    try
        [params_f, ~] = load_params_yaml(full_path);
        classes_v = sort(to_vec(params_f.integrator.classes));
        file_event_types{f} = {num2str(classes_v(1)), num2str(classes_v(2))};
        file_paradigm{f}    = params_f.integrator.paradigm;
        fprintf('[%d/%d] %s: event types %s/%s (paradigm=%s)\n', f, n_files, basename, ...
                file_event_types{f}{1}, file_event_types{f}{2}, file_paradigm{f});

        if isfield(params_f, 'processing_fbcsp_mi')
            csp_mi = load_csp(params_f, 'mi');
            file_csp_info{f}(end+1) = struct('name', 'MI', 'channels', {csp_mi.selected_channels}, 'bands', csp_mi.bands, 'csp_matrices', {csp_mi.csp_matrices});
        end
        if isfield(params_f, 'processing_fbcsp_cvsa')
            csp_cvsa = load_csp(params_f, 'cvsa');
            file_csp_info{f}(end+1) = struct('name', 'CVSA', 'channels', {csp_cvsa.selected_channels}, 'bands', csp_cvsa.bands, 'csp_matrices', {csp_cvsa.csp_matrices});
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

%% --- 5. PER-FILE LOAD + PREPROCESS + EPOCH + ERD/ERS BAND POWER ----------
all_heat_pow = cell(1, n_files); % {f}: {nbands} struct with .c1/.c2 [nchan x ntime] (% ERD/ERS)
file_label = cell(1, n_files);
file_csp   = cell(1, n_files);   % {f}: struct array('name','MI'|'CVSA','mask',logical), per-file CSP channel masks
ref_labels = {};
heat_times = [];
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

    % ── Per-trial CF-end times (variable trial duration in evaluation) ───
    % cf_end_ms{c}(k) = latency (ms, relative to cue onset at t=0) of the
    % 897/898/899 outcome event for the k-th trial of class c, i.e. the end
    % of that trial's continuous-feedback window. Falls back to
    % heat_time_range(2) (no masking) if no outcome event is found.
    cf_end_ms = cell(1, 2);
    for c = 1:2
        trials_c = idx_cls{c};
        cf_end_ms{c} = repmat(heat_time_range(2), 1, numel(trials_c));
        for k = 1:numel(trials_c)
            ev_t = EEG_ep.epoch(trials_c(k)).eventtype;
            ev_l = [EEG_ep.epoch(trials_c(k)).eventlatency{:}];
            if ~iscell(ev_t), ev_t = {ev_t}; end
            for e = 1:numel(ev_t)
                if ev_l(e) > 0 && ev_is(ev_t{e}, OUTCOME_EVENT_TYPES)
                    cf_end_ms{c}(k) = min(ev_l(e), heat_time_range(2));
                    break;
                end
            end
        end
    end

    % ── Channel x time ERD/ERS heatmaps (per band of this file's paradigm)
    p_idx_file = find(strcmp(paradigms, file_paradigm{f}), 1);
    bnd_f = par_bands{p_idx_file};
    heat_pow_f = cell(size(bnd_f,1), 1);
    for idx_b = 1:size(bnd_f,1)
        erd_pct = compute_band_power_pct(EEG_ep, bnd_f(idx_b,:), baseline_ms, ERD_SMOOTH_MS);
        for c = 1:2
            trials_c = idx_cls{c};
            e = erd_pct(:,:,trials_c);
            for k = 1:numel(trials_c)
                e(:, EEG_ep.times > cf_end_ms{c}(k), k) = NaN;
            end
            heat_pow_f{idx_b}.(sprintf('c%d', c)) = squeeze(mean(e, 3, 'omitnan'));
        end
    end
    all_heat_pow{f} = heat_pow_f;
    if f == 1
        heat_times   = EEG_ep.times;
        chanlocs_ref = EEG_ep.chanlocs;
    end
end

%% --- 6. GROUP FILES BY PARADIGM + PER-PARADIGM AVERAGES --------------------
par_mi_mask   = cell(1, n_par);   % {p}: union MI CSP mask across that paradigm's files ([] if none)
par_cvsa_mask = cell(1, n_par);   % {p}: union CVSA CSP mask across that paradigm's files ([] if none)
par_lims      = cell(1, n_par);   % {p}: band_lims cell (one per band of par_bands{p})
par_heat      = cell(1, n_par);   % {p}: {nbands} struct .c1/.c2 [nchan x ntime] grand average
par_heat_mask = cell(1, n_par);   % {p}: channel mask for "sel" heatmaps ([] = all channels)

% Time mask for colour-limit / topoplot-limit computation: CF window only
% (CF onset -> max CF end), same NaN-masked per-trial band power used for
% the heatmaps (per-trial CF-end-aware).
cf_idx = heat_times >= CUE_DURATION_S*1000 & heat_times <= heat_time_range(2);

for p = 1:n_par
    idx_p = find(strcmp(file_paradigm, paradigms{p}));

    par_mi_mask{p}   = union_mask_for(file_csp, idx_p, 'MI',   numel(ref_labels));
    par_cvsa_mask{p} = union_mask_for(file_csp, idx_p, 'CVSA', numel(ref_labels));

    bnd    = par_bands{p};
    origin = par_band_origin{p};

    % ── Channel x time heatmap: grand average across this paradigm's files
    par_heat{p} = cell(size(bnd,1), 1);
    for idx_b = 1:size(bnd,1)
        c1_list = cell(1, numel(idx_p)); c2_list = cell(1, numel(idx_p));
        for ii = 1:numel(idx_p)
            c1_list{ii} = all_heat_pow{idx_p(ii)}{idx_b}.c1;
            c2_list{ii} = all_heat_pow{idx_p(ii)}{idx_b}.c2;
        end
        par_heat{p}{idx_b}.c1 = mean(cat(3, c1_list{:}), 3, 'omitnan');
        par_heat{p}{idx_b}.c2 = mean(cat(3, c2_list{:}), 3, 'omitnan');
    end

    % ── Topoplot colour limits, from the same NaN-masked band power
    %    (over cue->max-CF-end), so topoplots use exactly the data the
    %    per-interval averages below are drawn from.
    lims = cell(size(bnd,1), 1);
    for idx_b = 1:size(bnd,1)
        c1 = par_heat{p}{idx_b}.c1(:, cf_idx);
        c2 = par_heat{p}{idx_b}.c2(:, cf_idx);
        if strcmpi(origin{idx_b}, 'MI')
            neg = min(cat(1, c1, c2), 0);
            mx  = max(-neg(:));
            if mx == 0 || isnan(mx), mx = 1; end
            lims{idx_b} = [-mx 0];
        else   % CVSA: lateralization index = class1 - class2
            d  = c1 - c2;
            mx = max(abs(d(:)));
            if mx == 0 || isnan(mx), mx = 1; end
            lims{idx_b} = [-mx mx];
        end
    end
    par_lims{p} = lims;

    % ── "sel" channel mask for heatmaps: own CSP mask (MI/CVSA), or the
    %    union of MI+CVSA CSP masks for hybrid (used for ALL bands)
    switch paradigms{p}
        case 'mi',   par_heat_mask{p} = par_mi_mask{p};
        case 'cvsa', par_heat_mask{p} = par_cvsa_mask{p};
        case 'hybrid'
            m1 = par_mi_mask{p};   if isempty(m1), m1 = false(numel(ref_labels),1); end
            m2 = par_cvsa_mask{p}; if isempty(m2), m2 = false(numel(ref_labels),1); end
            u = m1 | m2;
            if any(u), par_heat_mask{p} = u; else, par_heat_mask{p} = []; end
        otherwise
            par_heat_mask{p} = [];
    end

    fprintf('Paradigm "%s": %d file(s), MI CSP mask: %s, CVSA CSP mask: %s\n', paradigms{p}, numel(idx_p), ...
            mask_summary(par_mi_mask{p}), mask_summary(par_cvsa_mask{p}));
end

%% --- 7. GENERATE TOPOPLOT GRIDS ---------------------------------------------
fprintf('\n--- GENERATING TOPOPLOTS ---\n');

erd_summary = struct([]);

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
    plot_topo_grid(par_heat{p}, chanlocs_ref, heat_times, bnd, bnames, origin, ...
                   lims, intervals, titles_cols, evt, masks_none, par_dir, sprintf('all_%s', paradigms{p}), fig_vis);
    if any(~cellfun(@isempty, masks_sel))
        plot_topo_grid(par_heat{p}, chanlocs_ref, heat_times, bnd, bnames, origin, ...
                       lims, intervals, titles_cols, evt, masks_sel, par_dir, sprintf('sel_%s', paradigms{p}), fig_vis);
    end

    % ── Channel x time ERD/ERS heatmaps: grand average ───────────────────
    chan_labels = {chanlocs_ref.labels};
    plot_heatmap_grid(par_heat{p}, chan_labels, heat_times, heat_time_range, ...
                       CUE_DURATION_S*1000, bnames, evt, [], par_dir, sprintf('all_%s', paradigms{p}), fig_vis);
    if ~isempty(par_heat_mask{p})
        plot_heatmap_grid(par_heat{p}, chan_labels, heat_times, heat_time_range, ...
                           CUE_DURATION_S*1000, bnames, evt, par_heat_mask{p}, par_dir, sprintf('sel_%s', paradigms{p}), fig_vis);
    end

    % ── ERD/ERS magnitude vs. CSP channel weight (grand average) ─────────
    idx_p = find(strcmp(file_paradigm, paradigms{p}));
    plot_csp_correlation(par_heat{p}, chan_labels, heat_times, heat_time_range, CUE_DURATION_S*1000, ...
                          bnd, bnames, origin, file_csp_info, idx_p, ref_labels, par_dir, paradigms{p}, fig_vis);

    % ── Per-band scalar summary for cross-subject use (mean ERD per class,
    %    discrimination, CSP-weight correlation) ──────────────────────────
    for idx_b = 1:size(bnd,1)
        if strcmpi(origin{idx_b}, 'MI'), mask_b = par_mi_mask{p}; else, mask_b = par_cvsa_mask{p}; end
        row = compute_erd_summary_row(par_heat{p}{idx_b}, heat_times, heat_time_range, CUE_DURATION_S*1000, ...
                                       bnd(idx_b,:), bnames{idx_b}, origin{idx_b}, file_csp_info, idx_p, ref_labels, mask_b);
        row.paradigm = paradigms{p};
        erd_summary(end+1) = row; %#ok<AGROW>
    end

    % ── ERD/ERS lateralization time-course, by CSP-channel hemisphere ROI ─
    plot_csp_roi_timecourse(par_heat{p}, ref_labels, heat_times, heat_time_range, CUE_DURATION_S*1000, ...
                             bnd, bnames, origin, evt, par_mi_mask{p}, par_cvsa_mask{p}, par_dir, paradigms{p}, fig_vis);

    % ── Per-file (same colour scale + bands as this paradigm) ────────────
    for f = idx_p
        plot_topo_grid(all_heat_pow{f}, chanlocs_ref, heat_times, bnd, bnames, origin, ...
                       lims, intervals, titles_cols, evt, masks_none, par_dir_files, sprintf('all_%s', file_label{f}), fig_vis);

        file_mi_mask   = union_mask_for(file_csp, f, 'MI',   numel(ref_labels));
        file_cvsa_mask = union_mask_for(file_csp, f, 'CVSA', numel(ref_labels));
        file_masks_sel = band_masks_for(origin, file_mi_mask, file_cvsa_mask);
        if any(~cellfun(@isempty, file_masks_sel))
            plot_topo_grid(all_heat_pow{f}, chanlocs_ref, heat_times, bnd, bnames, origin, ...
                           lims, intervals, titles_cols, evt, file_masks_sel, par_dir_files, sprintf('sel_%s', file_label{f}), fig_vis);
        end

        plot_csp_roi_timecourse(all_heat_pow{f}, ref_labels, heat_times, heat_time_range, CUE_DURATION_S*1000, ...
                                 bnd, bnames, origin, evt, file_mi_mask, file_cvsa_mask, par_dir_files, file_label{f}, fig_vis);

        plot_heatmap_grid(all_heat_pow{f}, chan_labels, heat_times, heat_time_range, ...
                           CUE_DURATION_S*1000, bnames, evt, [], par_dir_files, sprintf('all_%s', file_label{f}), fig_vis);
        switch paradigms{p}
            case 'mi',   file_heat_mask = file_mi_mask;
            case 'cvsa', file_heat_mask = file_cvsa_mask;
            case 'hybrid'
                m1 = file_mi_mask;   if isempty(m1), m1 = false(numel(ref_labels),1); end
                m2 = file_cvsa_mask; if isempty(m2), m2 = false(numel(ref_labels),1); end
                u = m1 | m2;
                if any(u), file_heat_mask = u; else, file_heat_mask = []; end
            otherwise
                file_heat_mask = [];
        end
        if ~isempty(file_heat_mask)
            plot_heatmap_grid(all_heat_pow{f}, chan_labels, heat_times, heat_time_range, ...
                               CUE_DURATION_S*1000, bnames, evt, file_heat_mask, par_dir_files, sprintf('sel_%s', file_label{f}), fig_vis);
        end
    end
end

save(fullfile(output_dir, 'topo_erders_summary.mat'), 'erd_summary');
fprintf('Saved topo_erders_summary.mat to %s\n', output_dir);

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

function plot_topo_grid(heat_bands, chanlocs, heat_times, bands, band_names, band_origin, ...
                         band_lims, intervals, titles_cols, event_types_p, band_masks, out_dir, prefix, fig_vis)
% PLOT_TOPO_GRID  One figure per band: cols = time interval.
%   heat_bands    : {N} struct .c1/.c2 [nchan x ntime] percent ERD(-)/ERS(+),
%                   same per-trial CF-end-aware NaN-masked band power as the
%                   heatmaps (heat_times in ms, t=0 at cue onset)
%   bands         : [Nx2] frequency bands (Hz), used for figure titles only
%   band_origin   : {N} 'MI' (ERD-only, 2 rows = cues) | 'CVSA' (lateralization
%                   index cue1-cue2, 1 row)
%   band_lims     : {N} colour-axis limits per band
%   intervals     : [Mx2] time windows (s, relative to cue onset)
%   event_types_p : {2} cue/event codes for this paradigm, used as row labels
%   band_masks    : {N} of [] or logical [nchan x 1] — channels not in mask
%                   are forced to 0

    num_intervals = size(intervals, 1);

    for idx_b = 1:size(bands,1)
        low_f = bands(idx_b,1); high_f = bands(idx_b,2);

        c1   = heat_bands{idx_b}.c1;   % [nchan x ntime]
        c2   = heat_bands{idx_b}.c2;
        mask  = band_masks{idx_b};
        is_mi = strcmpi(band_origin{idx_b}, 'MI');

        if is_mi
            num_rows = 2;
            row_labels = {sprintf('Cue %s', event_types_p{1}), sprintf('Cue %s', event_types_p{2})};
        else
            num_rows = 1;
            row_labels = {sprintf('Cue %s - Cue %s', event_types_p{1}, event_types_p{2})};
        end

        h = figure('Name', sprintf('%s — %s %g-%g Hz', prefix, band_origin{idx_b}, low_f, high_f), ...
                    'Color','w', 'NumberTitle','off', 'Visible',fig_vis);
        set(h, 'Units','normalized', 'OuterPosition', [0 0 1 1]);
        tiledlayout(num_rows, num_intervals, 'TileSpacing','compact', 'Padding','tight');

        for r = 1:num_rows
            for c = 1:num_intervals
                nexttile;
                t_idx = heat_times >= intervals(c,1)*1000 & heat_times < intervals(c,2)*1000;

                if is_mi
                    src  = c1; if r == 2, src = c2; end
                    data = mean(src(:, t_idx), 2, 'omitnan');
                    data(data > 0) = 0;   % ERD only
                else
                    d1 = mean(c1(:, t_idx), 2, 'omitnan');
                    d2 = mean(c2(:, t_idx), 2, 'omitnan');
                    data = d1 - d2;
                end
                % No trial of this class has CF still running in this
                % interval for some channel -> mean is NaN; show as 0
                % (= no ERD/ERS), consistent with the heatmap background.
                data(isnan(data)) = 0;
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

        save_path = fullfile(out_dir, sprintf('topo_%s_%s.svg', prefix, band_names{idx_b}));
        saveas(h, save_path, 'svg');
        if ~strcmp(fig_vis, 'on'), close(h); end
        fprintf('  Saved: %s\n', save_path);
    end
end

function tf = ev_is(t, codes)
% EV_IS  True if event type t (char or numeric) matches any of the cellstr codes.
    if isnumeric(t), t = num2str(t); end
    tf = ismember(t, codes);
end

function erd_pct = compute_band_power_pct(EEG_ep, band, baseline_ms, smooth_ms)
% COMPUTE_BAND_POWER_PCT  Per-trial, per-channel band-power percent ERD(-)/ERS(+)
%   relative to that trial's own pre-cue baseline (Pfurtscheller ERD%):
%   bandpass -> instantaneous power -> moving-average smoothing -> percent
%   change vs. mean power in baseline_ms.
%   Returns [nchan x npts x ntrials].
    fs = EEG_ep.srate;
    [nch, npts, ntr] = size(EEG_ep.data);
    X  = reshape(permute(EEG_ep.data, [2 1 3]), npts, nch*ntr);
    Xf = bandpass(X, band, fs, 'ImpulseResponse', 'iir', 'Steepness', 0.5);
    bp = permute(reshape(Xf, npts, nch, ntr), [2 1 3]);

    pw = bp .^ 2;
    smooth_n = max(1, round(smooth_ms/1000 * fs));
    pw = movmean(pw, smooth_n, 2);

    base_idx = EEG_ep.times >= baseline_ms(1) & EEG_ep.times <= baseline_ms(2);
    base_pow = mean(pw(:, base_idx, :), 2);
    erd_pct  = (pw - base_pow) ./ base_pow * 100;
end

function cmap = erd_colormap(n)
% ERD_COLORMAP  Diverging blue (ERD, negative) - white (0) - red (ERS, positive).
    if nargin < 1, n = 256; end
    half = floor(n/2);
    blue_white = [linspace(0,1,half)',   linspace(0,1,half)',   ones(half,1)];
    white_red  = [ones(n-half,1), linspace(1,0,n-half)', linspace(1,0,n-half)'];
    cmap = [blue_white; white_red];
end

function plot_heatmap_grid(heat_bands, chan_labels, heat_times, heat_time_range, ...
                            cue_duration_ms, band_names, event_types_p, sel_mask, out_dir, prefix, fig_vis)
% PLOT_HEATMAP_GRID  One figure per band: 2 rows (one per cue/class), each a
%   channel x time heatmap of % ERD(-)/ERS(+), from cue onset (t=0) to the
%   max CF end (heat_time_range(2)). Diverging colormap centred at 0, with
%   symmetric limits = max(abs(value)) over the data shown in that figure.
%   sel_mask: [] -> all channels; logical mask -> only those channels.

    t_idx = heat_times >= heat_time_range(1) & heat_times <= heat_time_range(2);
    t_s   = heat_times(t_idx) / 1000;

    if isempty(sel_mask)
        ch_idx = 1:numel(chan_labels);
    else
        ch_idx = find(sel_mask);
        if isempty(ch_idx), return; end
    end

    cmap = erd_colormap(256);

    for idx_b = 1:numel(heat_bands)
        c1 = heat_bands{idx_b}.c1(ch_idx, t_idx);
        c2 = heat_bands{idx_b}.c2(ch_idx, t_idx);
        mx = max(abs([c1(:); c2(:)]));
        if mx == 0 || isnan(mx), mx = 1; end

        h = figure('Name', sprintf('%s — %s (heatmap)', prefix, band_names{idx_b}), ...
                    'Color','w', 'NumberTitle','off', 'Visible',fig_vis);
        set(h, 'Units','normalized', 'OuterPosition', [0 0 0.5 1]);
        tiledlayout(2, 1, 'TileSpacing','compact', 'Padding','tight');

        data   = {c1, c2};
        row_labels = {sprintf('Cue %s', event_types_p{1}), sprintf('Cue %s', event_types_p{2})};
        for r = 1:2
            ax = nexttile;
            im = imagesc(ax, t_s, 1:numel(ch_idx), data{r}, [-mx mx]);
            set(im, 'AlphaData', ~isnan(data{r}));
            set(ax, 'Color', [0.85 0.85 0.85], 'YDir','normal', ...
                    'YTick', 1:numel(ch_idx), 'YTickLabel', chan_labels(ch_idx));
            colormap(ax, cmap);
            colorbar(ax);
            hold(ax, 'on');
            xline(ax, cue_duration_ms/1000, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');
            xlabel(ax, 't [s]  (0 = cue onset, dashed = CF onset)');
            ylabel(ax, 'Channel');
            title(ax, row_labels{r}, 'FontWeight','bold');
        end

        sgtitle(sprintf('%s — %s (ERD/ERS heatmap)', prefix, band_names{idx_b}), ...
                'Interpreter','none', 'FontWeight','bold');

        save_path = fullfile(out_dir, sprintf('heat_%s_%s.svg', prefix, band_names{idx_b}));
        saveas(h, save_path, 'svg');
        if ~strcmp(fig_vis, 'on'), close(h); end
        fprintf('  Saved: %s\n', save_path);
    end
end

function imp = csp_importance_for_band(file_csp_info, idx_files, name, band, ref_labels)
% CSP_IMPORTANCE_FOR_BAND  Per-channel CSP "importance" for the given band:
%   sum(|csp_matrices{b}|, 1) over CSP components, scattered into the full
%   ref_labels channel set (0 for channels not selected by the CSP).
%   name : 'MI'|'CVSA' — which CSP block to look up.
%   band : [1x2] frequency band (Hz) to match against file_csp_info(.).bands.
%   Returns [] if no file has CSP info for this (name, band) combination.
    imp = [];
    for f = idx_files
        for v = 1:numel(file_csp_info{f})
            if ~strcmpi(file_csp_info{f}(v).name, name), continue; end
            bnd_v = file_csp_info{f}(v).bands;
            row = find(all(abs(bnd_v - band) < 1e-6, 2), 1);
            if isempty(row), continue; end

            w = sum(abs(file_csp_info{f}(v).csp_matrices{row}), 1);   % [1 x n_sel]
            sel_lbl = lower(strtrim(file_csp_info{f}(v).channels));
            all_lbl = lower(strtrim(ref_labels));
            imp = zeros(numel(ref_labels), 1);
            for s = 1:numel(sel_lbl)
                ch_i = find(strcmp(all_lbl, sel_lbl{s}), 1);
                if ~isempty(ch_i), imp(ch_i) = w(s); end
            end
            return;
        end
    end
end

function row = compute_erd_summary_row(heat_band, heat_times, heat_time_range, cf_start_ms, ...
                                        band, band_name, origin, file_csp_info, idx_files, ref_labels, mask)
% COMPUTE_ERD_SUMMARY_ROW  Per-band scalar summary for cross-subject pooling:
%   mean ERD/ERS per class (averaged over the CSP-selected channels of this
%   band's origin, during the CF window), their discrimination, and the
%   Pearson r between per-channel ERD/ERS discrimination and CSP weight
%   (same r reported in the corr_<paradigm> figure panel for this band).
%   Returns NaN fields if no CSP info is found for this (origin, band).
    cf_idx_l = heat_times >= cf_start_ms & heat_times <= heat_time_range(2);
    row = struct('band_name', band_name, 'band_origin', origin, ...
                 'freq_lo', band(1), 'freq_hi', band(2), ...
                 'mean_erd_c1', NaN, 'mean_erd_c2', NaN, 'discrimination', NaN, ...
                 'csp_r', NaN, 'n_csp_channels', 0);

    imp = csp_importance_for_band(file_csp_info, idx_files, origin, band, ref_labels);
    if isempty(imp), return; end

    d1_full = mean(heat_band.c1(:, cf_idx_l), 2, 'omitnan');
    d2_full = mean(heat_band.c2(:, cf_idx_l), 2, 'omitnan');
    mag = abs(d1_full - d2_full);

    good = ~isnan(mag) & ~isnan(imp);
    if sum(good) >= 2
        cc = corrcoef(mag(good), imp(good));
        row.csp_r = cc(1,2);
    end

    if ~isempty(mask) && any(mask)
        row.mean_erd_c1    = mean(d1_full(mask), 'omitnan');
        row.mean_erd_c2    = mean(d2_full(mask), 'omitnan');
        row.discrimination = abs(row.mean_erd_c1 - row.mean_erd_c2);
        row.n_csp_channels = sum(mask);
    end
end

function plot_csp_correlation(heat_bands, chan_labels, heat_times, heat_time_range, cf_start_ms, ...
                               bands, band_names, band_origin, file_csp_info, idx_files, ref_labels, out_dir, prefix, fig_vis)
% PLOT_CSP_CORRELATION  One figure, one panel per band with CSP info: scatter
%   of per-channel ERD/ERS class discrimination |c1-c2| during continuous
%   feedback (cf_start_ms -> heat_time_range(2), from the same NaN-masked
%   band power as the heatmaps) vs. that channel's CSP weight (sum |w| over
%   CSP components, MI-origin bands use the MI CSP, CVSA-origin bands use
%   the CVSA CSP). |c1-c2| (rather than the per-class average) is used for
%   both MI and CVSA bands: it captures the contralateral ERD/ERS pattern
%   (e.g. C4 ERD for class1 / C3 ERD for class2 both show up as large
%   |c1-c2|, whereas averaging the two classes would cancel them out) and
%   matches the lateralization index already shown for CVSA topoplots.
%   Bands without CSP info (e.g. DEFAULT_BANDS fallback) are skipped.

    cf_idx = heat_times >= cf_start_ms & heat_times <= heat_time_range(2);

    nb    = size(bands,1);
    valid = false(nb,1);
    imps  = cell(nb,1); mags = cell(nb,1); rvals = nan(nb,1);

    for idx_b = 1:nb
        imp = csp_importance_for_band(file_csp_info, idx_files, band_origin{idx_b}, bands(idx_b,:), ref_labels);
        if isempty(imp), continue; end

        d1  = mean(heat_bands{idx_b}.c1(:, cf_idx), 2, 'omitnan');
        d2  = mean(heat_bands{idx_b}.c2(:, cf_idx), 2, 'omitnan');
        mag = abs(d1 - d2);

        good = ~isnan(mag) & ~isnan(imp);
        if sum(good) < 2, continue; end

        imps{idx_b}  = imp;
        mags{idx_b}  = mag;
        cc = corrcoef(mag(good), imp(good));
        rvals(idx_b) = cc(1,2);
        valid(idx_b) = true;
    end

    if ~any(valid), return; end

    h = figure('Name', sprintf('%s — ERD/ERS discrimination vs CSP weight', prefix), ...
                'Color','w', 'NumberTitle','off', 'Visible',fig_vis);
    set(h, 'Units','normalized', 'OuterPosition', [0 0 1 1]);
    tiledlayout(1, sum(valid), 'TileSpacing','compact', 'Padding','tight');

    for idx_b = find(valid)'
        nexttile;
        scatter(mags{idx_b}, imps{idx_b}, 36, 'filled');
        text(mags{idx_b}, imps{idx_b}, chan_labels(:), 'FontSize', 7, ...
             'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'Interpreter','none');
        xlabel('|ERD/ERS class1 - class2| during CF [%]');
        ylabel(sprintf('%s CSP weight (sum |w|)', band_origin{idx_b}));
        title(sprintf('%s (r=%.2f)', band_names{idx_b}, rvals(idx_b)), 'Interpreter','none');
        grid on;
    end

    sgtitle(sprintf('%s — ERD/ERS class discrimination vs. CSP channel weight (CF window)', prefix), ...
            'Interpreter','none', 'FontWeight','bold');

    save_path = fullfile(out_dir, sprintf('corr_%s.svg', prefix));
    saveas(h, save_path, 'svg');
    if ~strcmp(fig_vis, 'on'), close(h); end
    fprintf('  Saved: %s\n', save_path);
end

function roi = hemisphere_of(label)
% HEMISPHERE_OF  'Left'/'Right'/'Midline' from a 10-20 channel label, based
%   on its trailing digit (odd = left, even = right) or trailing 'z'
%   (midline). 'Other' if the label ends in neither.
    lbl  = strtrim(label);
    last = lbl(end);
    if strcmpi(last, 'z')
        roi = 'Midline';
    elseif ~isnan(str2double(last))
        if mod(str2double(last), 2) == 1
            roi = 'Left';
        else
            roi = 'Right';
        end
    else
        roi = 'Other';
    end
end

function plot_csp_roi_timecourse(heat_bands, ref_labels, heat_times, heat_time_range, cue_duration_ms, ...
                                  bands, band_names, band_origin, event_types_p, mi_mask, cvsa_mask, out_dir, prefix, fig_vis)
% PLOT_CSP_ROI_TIMECOURSE  One figure per band: ERD/ERS time-course
%   (cue onset -> max CF end, from the same NaN-masked band power as the
%   heatmaps/topoplots) averaged over the CSP-selected channels of that
%   band's origin, grouped into hemisphere ROIs via HEMISPHERE_OF
%   (Left/Right/Midline; 'Other' channels are dropped).
%   Left column: per-class traces (c1, c2); right column: lateralization
%   (c1-c2). MI-origin bands use mi_mask, CVSA-origin bands use cvsa_mask;
%   bands with no (or empty) CSP mask for their origin are skipped.

    t_idx = heat_times >= heat_time_range(1) & heat_times <= heat_time_range(2);
    t_s   = heat_times(t_idx) / 1000;

    roi_names = {'Left', 'Right', 'Midline'};
    color_c1  = [0.85 0.33 0.10];
    color_c2  = [0.00 0.45 0.74];
    color_lat = [0.49 0.18 0.56];

    for idx_b = 1:size(bands,1)
        if strcmpi(band_origin{idx_b}, 'MI'), mask = mi_mask; else, mask = cvsa_mask; end
        if isempty(mask) || ~any(mask), continue; end

        ch_idx = find(mask);
        hemis  = cellfun(@hemisphere_of, ref_labels(ch_idx), 'UniformOutput', false);

        roi_idx = struct();
        present = {};
        for r = 1:numel(roi_names)
            sel = ch_idx(strcmp(hemis, roi_names{r}));
            if ~isempty(sel)
                roi_idx.(roi_names{r}) = sel;
                present{end+1} = roi_names{r}; %#ok<AGROW>
            end
        end
        if isempty(present), continue; end

        c1 = heat_bands{idx_b}.c1;   % [nchan x ntime]
        c2 = heat_bands{idx_b}.c2;

        h = figure('Name', sprintf('%s — %s ROI time-course', prefix, band_names{idx_b}), ...
                    'Color','w', 'NumberTitle','off', 'Visible',fig_vis);
        set(h, 'Units','normalized', 'OuterPosition', [0 0 0.6 1]);
        tiledlayout(numel(present), 2, 'TileSpacing','compact', 'Padding','tight');

        for r = 1:numel(present)
            roi    = present{r};
            idx    = roi_idx.(roi);
            roi_c1 = mean(c1(idx, t_idx), 1, 'omitnan');
            roi_c2 = mean(c2(idx, t_idx), 1, 'omitnan');
            roi_lat = roi_c1 - roi_c2;

            ax1 = nexttile;
            plot(ax1, t_s, roi_c1, 'Color', color_c1, 'LineWidth', 1.5); hold(ax1, 'on');
            plot(ax1, t_s, roi_c2, 'Color', color_c2, 'LineWidth', 1.5);
            yline(ax1, 0, 'k:', 'HandleVisibility','off');
            xline(ax1, cue_duration_ms/1000, 'k--', 'HandleVisibility','off');
            ylabel(ax1, sprintf('%s (%d ch)\nERD/ERS [%%]', roi, numel(idx)), 'FontWeight','bold');
            if r == 1
                legend(ax1, {sprintf('Cue %s', event_types_p{1}), sprintf('Cue %s', event_types_p{2})}, 'Location','best');
                title(ax1, 'ERD/ERS per class');
            end
            if r == numel(present), xlabel(ax1, 't [s] (0 = cue onset, dashed = CF onset)'); end
            grid(ax1, 'on');

            ax2 = nexttile;
            plot(ax2, t_s, roi_lat, 'Color', color_lat, 'LineWidth', 1.5);
            yline(ax2, 0, 'k:', 'HandleVisibility','off');
            xline(ax2, cue_duration_ms/1000, 'k--', 'HandleVisibility','off');
            ylabel(ax2, 'c1 - c2 [%]');
            if r == 1, title(ax2, 'Lateralization (class1 - class2)'); end
            if r == numel(present), xlabel(ax2, 't [s] (0 = cue onset, dashed = CF onset)'); end
            grid(ax2, 'on');
        end

        sgtitle(sprintf('%s — %s ROI time-course (CSP-selected channels, %s CSP)', ...
                prefix, band_names{idx_b}, band_origin{idx_b}), 'Interpreter','none', 'FontWeight','bold');

        save_path = fullfile(out_dir, sprintf('roi_%s_%s.svg', prefix, band_names{idx_b}));
        saveas(h, save_path, 'svg');
        if ~strcmp(fig_vis, 'on'), close(h); end
        fprintf('  Saved: %s\n', save_path);
    end
end
