%% MAIN_SESSION_OVERVIEW  Quick accuracy overview of BCI evaluation sessions.
%
%   Loads one or more GDF evaluation files. Paradigm is inferred from the
%   filename (first keyword found: 'hybrid' > 'cvsa' > 'mi').
%   Labels are assigned as: mi_1, mi_2, cvsa_1, cvsa_2, hybrid_1, hybrid_2 ...
%
%   Console output, per file: per-trial outcome (HIT/MISS/TIMEOUT) + time
%   for that event, then a session summary with per-file means and
%   per-paradigm experiment totals (pooled across all files of that
%   paradigm). Set VERBOSE_DIAGNOSTIC=true for the per-trial P/argmax
%   class-ordering diagnostic (off by default — very verbose).
%
%   Produces three figures (event-based metrics + offline simulation):
%
%     Fig 1 — Trial accuracy from real GDF outcomes (897 / 898 / 899)
%              Subplot 1: 897 / (897+898+899)
%              Subplot 2: 897 / (897+898)  [no timeout]
%
%     Fig 2 — Time metrics from real GDF event positions
%              Subplot 1: Time to HIT (897 trials)  — mean ± std + individual dots
%              Subplot 2: Time to MISS (898 trials) — mean ± std + individual dots
%
%     Fig 3 — Sample accuracy from offline simulation (MATLAB ≈ ROS, MAE < 1e-3)
%              "At each CF frame, is argmax(P_sLDA) the correct target class?"
%              Subplot 1: MI classifier  — full CF
%              Subplot 2: CVSA classifier — first 3 s only
%              Subplot 3: Hybrid fused   — full CF  (hybrid files only)
%
%   Colors:   MI = orange   CVSA = green   Hybrid = blue

clear; clc; close all;

% --- Display options -------------------------------------------------------
SHOW_FIGURES = false;   % true: figures pop up on screen; false: created hidden (export only)
if SHOW_FIGURES, fig_vis = 'on'; else, fig_vis = 'off'; end

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir,'io'), fullfile(this_dir,'processing'), ...
        fullfile(this_dir,'artifacts'), fullfile(this_dir,'classifier'), ...
        fullfile(this_dir,'integrator'), fullfile(this_dir,'utils'));

HIT_EV = 897;  MISS_EV = 898;  TO_EV = 899;  CF_EV = 781;
CVSA_MAX_S = 3.0;   % seconds used for CVSA sample accuracy

% Set to true to print, for every trial, the per-frame P/argmax diagnostic
% used to debug class-ordering issues. Very verbose — leave false for a
% normal session overview.
VERBOSE_DIAGNOSTIC = false;

COL = struct('mi', [0.85 0.30 0.10], ...
             'cvsa',   [0.10 0.60 0.30], ...
             'hybrid', [0.18 0.45 0.75], ...
             'unknown',[0.55 0.55 0.55]);

% ── File picker ───────────────────────────────────────────────────────────────
default_dir = '/home/paolo/bci_vr_ws/recordings';
if ~isfolder(default_dir), default_dir = this_dir; end
[gdf_names, gdf_dir] = uigetfile({'*.gdf','GDF files (*.gdf)'}, ...
    'Select evaluation GDF file(s)', default_dir, 'MultiSelect','on');
if isequal(gdf_names,0), error('main_session_overview:cancel','No file selected.'); end
if ischar(gdf_names), gdf_names = {gdf_names}; end
n_files = numel(gdf_names);

par_count = struct('mi',0,'cvsa',0,'hybrid',0,'unknown',0);
RES       = {};   % cell array of result structs

%% ═════════════════════════════════════════════════════════════════════════════
%  PER-FILE LOOP
%  ═════════════════════════════════════════════════════════════════════════════
for fi = 1:n_files
    gdf_path = fullfile(gdf_dir, gdf_names{fi});
    [~, basename] = fileparts(gdf_path);
    fprintf('\n[%d/%d]  %s\n', fi, n_files, basename);

    % ── Paradigm detection ────────────────────────────────────────────────────
    paradigm = detect_paradigm(basename);
    par_count.(paradigm) = par_count.(paradigm) + 1;
    file_label = sprintf('%s_%d', paradigm, par_count.(paradigm));
    col = COL.(paradigm);

    % ── Load GDF (raw events only — no YAML needed here) ─────────────────────
    [signal, header, ~] = load_gdf(gdf_path);
    fs  = header.SampleRate;
    POS = header.EVENT.POS;
    TYP = header.EVENT.TYP;

    % ── Event-based trial accuracy + TTH / time-to-miss ──────────────────────
    cf_pos  = POS(TYP == CF_EV);
    n_cf_ev = numel(cf_pos);
    outcomes    = zeros(1, n_cf_ev);
    dt_vals     = nan(1, n_cf_ev);
    tth_vals    = [];
    t_miss_vals = [];
    to_vals     = [];

    for t = 1:n_cf_ev
        after = POS > cf_pos(t);
        idx   = find((TYP==HIT_EV | TYP==MISS_EV | TYP==TO_EV) & after, 1);
        if isempty(idx), continue; end
        outcomes(t) = TYP(idx);
        dt = (POS(idx) - cf_pos(t)) / fs;
        dt_vals(t) = dt;
        if outcomes(t) == HIT_EV
            tth_vals(end+1)    = dt; %#ok<AGROW>
        elseif outcomes(t) == MISS_EV
            t_miss_vals(end+1) = dt; %#ok<AGROW>
        elseif outcomes(t) == TO_EV
            to_vals(end+1)     = dt; %#ok<AGROW>
        end
    end

    % ── Per-trial outcomes ────────────────────────────────────────────────────
    fprintf('  Per-trial outcomes:\n');
    for t = 1:n_cf_ev
        switch outcomes(t)
            case HIT_EV,  out_str = 'HIT';
            case MISS_EV, out_str = 'MISS';
            case TO_EV,   out_str = 'TIMEOUT';
            otherwise,    out_str = '?';
        end
        if isnan(dt_vals(t))
            fprintf('    trial %2d: %-7s\n', t, out_str);
        else
            fprintf('    trial %2d: %-7s  t=%.2fs\n', t, out_str, dt_vals(t));
        end
    end

    n_hit   = sum(outcomes == HIT_EV);
    n_miss  = sum(outcomes == MISS_EV);
    n_to    = sum(outcomes == TO_EV);
    n_tot   = max(1, n_hit + n_miss + n_to);
    hit_tot  = n_hit / n_tot;
    hit_no_to = n_hit / max(1, n_hit + n_miss);
    fprintf('  %s | HIT=%d MISS=%d TO=%d  acc=%.0f%%  no-to=%.0f%%  TTH=%.1fs  Tmiss=%.1fs\n', ...
            file_label, n_hit, n_miss, n_to, hit_tot*100, hit_no_to*100, ...
            msafe(tth_vals), msafe(t_miss_vals));

    % ── Offline simulation for sample accuracy ────────────────────────────────
    sa_mi = NaN;  sa_cvsa = NaN;  sa_fused = NaN;
    sa_mi_hit = NaN; sa_mi_miss = NaN;
    sa_cvsa_hit = NaN; sa_cvsa_miss = NaN;
    sa_fused_hit = NaN; sa_fused_miss = NaN;
    sa_mi_cls = NaN; sa_cvsa_cls = NaN; sa_fused_cls = NaN;
    r_csp_mi        = [];   % populated inside try if MI CSP is available
    r_csp_cvsa      = [];   % populated inside try if CVSA CSP is available
    r_slda_mi       = [];   % sLDA-weighted spatial analysis, MI
    r_slda_cvsa     = [];   % sLDA-weighted spatial analysis, CVSA
    r_erd_mi        = [];   % ERD/ERS epochs, MI CSP channels
    r_erd_cvsa      = [];   % ERD/ERS epochs, CVSA CSP channels

    try
        [params, ~] = load_params_yaml(gdf_path);

        % Optional: set USE_LAUNCH_PARAMS = true (before this loop) to override
        % YAML model paths with the current evaluation.launch values.
        % Useful when the YAML for a specific file references wrong/old models.
        if exist('USE_LAUNCH_PARAMS','var') && USE_LAUNCH_PARAMS
            EVAL_LAUNCH_PATH = fullfile(this_dir,'..','..','launchers_bci','launch','evaluation.launch');
            if exist(EVAL_LAUNCH_PATH,'file')
                params = load_eval_launch_params(params, EVAL_LAUNCH_PATH);
            end
        end

        framerate  = double(params.acquisition.framerate);
        fs_yaml    = double(params.acquisition.samplerate);
        chunk_size = round(fs_yaml / framerate);
        if abs(fs_yaml - fs) > 1, fs_yaml = fs; chunk_size = round(fs/framerate); end

        bufsize_proc = double(params.RingBufferCfg.params.size);
        bufsize_art  = double(params.RingBufferCfgArtifact.params.size);
        eog_names    = to_strcell(params.CarCfg.params.EOG_ch_names);

        do_car_mi   = true;
        do_car_cvsa = true;
        if isfield(params,'processing_fbcsp_mi'),   do_car_mi   = logical(params.processing_fbcsp_mi.do_car);   end
        if isfield(params,'processing_fbcsp_cvsa'), do_car_cvsa = logical(params.processing_fbcsp_cvsa.do_car); end

        use_mi   = ismember(paradigm,{'mi','hybrid'});
        use_cvsa = ismember(paradigm,{'cvsa','hybrid'});

        csp_mi=[]; slda_mi=[]; csp_cvsa=[]; slda_cvsa=[];
        if use_mi,   csp_mi   = load_csp(params,'mi');   slda_mi   = load_slda(params,'mi');   end
        if use_cvsa, csp_cvsa = load_csp(params,'cvsa'); slda_cvsa = load_slda(params,'cvsa'); end

        proc_base = struct('samplerate',fs_yaml,'chunk_size',chunk_size, ...
                           'bufsize',bufsize_proc,'filter_order',4);
        feat_mi=[]; header_mi=[]; feat_cv=[]; header_cv=[];
        if use_mi
            cfg = proc_base; cfg.do_car=do_car_mi; cfg.eog_names=eog_names;
            [feat_mi, header_mi] = apply_processing(signal, header, csp_mi, cfg);
        end
        if use_cvsa
            cfg = proc_base; cfg.do_car=do_car_cvsa; cfg.eog_names=eog_names;
            [feat_cv, header_cv] = apply_processing(signal, header, csp_cvsa, cfg);
        end

        art_cfg = params.ArtifactCfg.params;
        art_cfg.EOG_ch_names = to_strcell(art_cfg.EOG_ch_names);
        [art_flags, ~] = detect_artifacts(signal, header, art_cfg, ...
            struct('samplerate',fs_yaml,'chunk_size',chunk_size,'bufsize_artifact',bufsize_art));

        p_mi_al=[]; if use_mi,   p_mi_al = apply_slda(feat_mi, slda_mi,   csp_mi.bands);   end
        p_cv_al=[]; if use_cvsa, p_cv_al = apply_slda(feat_cv, slda_cvsa, csp_cvsa.bands); end

        int_cfg = params.integrator;
        if ~isfield(int_cfg,'increment'),         int_cfg.increment = 1; end
        if ~isfield(int_cfg,'thresholds_rejection'), int_cfg.thresholds_rejection = []; end
        if ~isfield(int_cfg,'cvsa_influence'),    int_cfg.cvsa_influence = 3.0; end
        if ~isfield(int_cfg,'thresholds') || isempty(int_cfg.thresholds)
            int_cfg.thresholds = params.training_node.thresholds;
        end

        header_chunks = header_mi; if ~use_mi, header_chunks = header_cv; end
        header_chunks.framerate = framerate;

        trials = integrate_signal(p_mi_al, p_cv_al, art_flags, header_chunks, int_cfg, paradigm);
        n3s   = round(CVSA_MAX_S * framerate);
        n_cls = numel(to_vec(int_cfg.classes));

        % ── Class ordering verification ───────────────────────────────────────
        % Classes are always ordered ascending (lower code = class 1):
        %   MI: [769,770]   CVSA: [730,731]   Hybrid: [750,751]
        % sLDA output column 1 = P(class with lower code) = P(class 1)
        classes_sim = sort(to_vec(int_cfg.classes));
        fprintf('  [sim params] classes=[%s]  bufsize=%d  k_gain=%.2f  thresholds=[%s]\n', ...
                num2str(classes_sim(:)','%d '), int_cfg.buffer_size, int_cfg.k_gain, ...
                strjoin(cellfun(@(x)sprintf('%.2f',x),num2cell(to_vec(int_cfg.thresholds)),'un',0),','));

        % Outcome masks for HIT / MISS+TIMEOUT (aligned with trials by CF order)
        n_sim = numel(trials);
        n_ev  = min(n_sim, numel(outcomes));
        hit_mask  = false(1, n_sim);
        miss_mask = false(1, n_sim);
        hit_mask(1:n_ev)  = outcomes(1:n_ev) == HIT_EV;
        miss_mask(1:n_ev) = outcomes(1:n_ev) == MISS_EV | outcomes(1:n_ev) == TO_EV;

        % ── DIAGNOSTIC: inspect ALL trials (HIT and MISS) ────────────────────
        % Prints: onset event code from GDF, target_class, mean P per class, argmax.
        % KEY: if argmax_agrees is LOW on HIT trials → target_class mapping is wrong
        %      (onset event code not matching int_cfg.classes order in YAML).
        if VERBOSE_DIAGNOSTIC
        fprintf('  [DIAGNOSTIC — all trials: onset code | target_class | P | argmax_agrees | outcome]\n');
        POS_h = header_chunks.EVENT.POS;
        TYP_h = header_chunks.EVENT.TYP;
        classes_v = to_vec(int_cfg.classes);
        for t_d = 1:n_sim
            c_d  = trials(t_d).target_class;
            np_d = trials(t_d).n_pre;
            % Find actual onset event code from the chunk-level header
            sc_d = trials(t_d).start_chunk;
            onset_before = NaN;
            for j_d = find(POS_h < sc_d, 1, 'last'):-1:1
                if ismember(TYP_h(j_d), classes_v)
                    onset_before = TYP_h(j_d); break;
                end
            end
            rv_d  = trials(t_d).raw(np_d+1:end, :);
            % Also show first 5 frames (0.25s) separately to catch early HIT dynamics
            n_early_d = min(5, size(rv_d,1));
            rv_early = rv_d(1:n_early_d,:);
            vld_e = ~any(isnan(rv_early),2);
            vld_d = ~any(isnan(rv_d),2);
            if ~any(vld_d) || isnan(c_d)
                fprintf('    trial %2d: onset=%d  target_cls=%s  [no valid frames]\n', ...
                        t_d, onset_before, num2str(c_d));
                continue;
            end
            rv_vld   = rv_d(vld_d,:);
            agrees   = 100*mean(rv_vld(:,c_d)==max(rv_vld,[],2));
            out_str  = '?';
            if t_d <= n_ev
                if outcomes(t_d)==HIT_EV,  out_str='HIT';
                elseif outcomes(t_d)==MISS_EV, out_str='MISS';
                elseif outcomes(t_d)==TO_EV,   out_str='TO';
                end
            end
            % Early frames (first 0.25s): key to detect early CVSA-driven HITs
            p_early_str = '';
            if any(vld_e)
                re = rv_early(vld_e,:);
                p_early_str = sprintf('  early_P=[%.2f,%.2f]', mean(re(:,1)), mean(re(:,2)));
            end
            % Flag anomalies: HIT trial where argmax_agrees < 50% is paradoxical
            flag = '';
            if t_d <= n_ev && outcomes(t_d)==HIT_EV && agrees < 50
                flag = ' ← check early frames!';
            end
            fprintf('    trial %2d: onset=%d  cls=%d  P_all=[%.2f,%.2f]%s  argmax=%.0f%%  %s%s\n', ...
                    t_d, onset_before, c_d, mean(rv_vld(:,1)), mean(rv_vld(:,2)), ...
                    p_early_str, agrees, out_str, flag);
        end
        fprintf('  int_cfg.classes = [%s]  (class 1 = lower code = first column of sLDA)\n', ...
                num2str(classes_v(:)','%d '));
        end

        % ── CSP channel importance ─────────────────────────────────────────────
        %   Alternate component ordering: rows 1,3,5,... maximise class-1 variance;
        %   rows 2,4,6,... maximise class-2 variance (MNE component_order='alternate').
        classes_sorted = sort(to_vec(int_cfg.classes));
        if use_mi && ~isempty(csp_mi)
            [w1, w2, bw1, bw2, sel] = csp_channel_weights(csp_mi);
            r_csp_mi = struct('channels',    {csp_mi.selected_channels}, ...
                              'bands',       csp_mi.bands, ...
                              'n_bands',     csp_mi.n_bands, ...
                              'n_components',csp_mi.n_components, ...
                              'csp_matrices',{csp_mi.csp_matrices}, ...
                              'w_c1',        w1,  'w_c2',      w2, ...
                              'band_w_c1',   bw1, 'band_w_c2', bw2, ...
                              'selectivity', sel, ...
                              'class_codes', classes_sorted(:)');
            t3c1 = top_channels(csp_mi.selected_channels, w1, 3);
            t3c2 = top_channels(csp_mi.selected_channels, w2, 3);
            fprintf('  [CSP-MI]   top ch  c%d: %s  |  c%d: %s\n', ...
                    classes_sorted(1), strjoin(t3c1, '>'), ...
                    classes_sorted(2), strjoin(t3c2, '>'));
        end
        if use_cvsa && ~isempty(csp_cvsa)
            [w1, w2, bw1, bw2, sel] = csp_channel_weights(csp_cvsa);
            r_csp_cvsa = struct('channels',    {csp_cvsa.selected_channels}, ...
                                'bands',       csp_cvsa.bands, ...
                                'n_bands',     csp_cvsa.n_bands, ...
                                'n_components',csp_cvsa.n_components, ...
                                'csp_matrices',{csp_cvsa.csp_matrices}, ...
                                'w_c1',        w1,  'w_c2',      w2, ...
                                'band_w_c1',   bw1, 'band_w_c2', bw2, ...
                                'selectivity', sel, ...
                                'class_codes', classes_sorted(:)');
            t3c1 = top_channels(csp_cvsa.selected_channels, w1, 3);
            t3c2 = top_channels(csp_cvsa.selected_channels, w2, 3);
            fprintf('  [CSP-CVSA] top ch  c%d: %s  |  c%d: %s\n', ...
                    classes_sorted(1), strjoin(t3c1, '>'), ...
                    classes_sorted(2), strjoin(t3c2, '>'));
        end

        % ── sLDA-weighted spatial analysis ─────────────────────────────────────
        %   Per-channel importance = sum_k |sLDA_coef_k| × |CSP_filter_row_k(ch)|
        %   over all sLDA-selected features k. Reveals which channels and bands
        %   the final classifier actually relies on and how strongly.
        if use_mi && ~isempty(csp_mi) && ~isempty(slda_mi)
            r_slda_mi = slda_spatial_weights(csp_mi, slda_mi);
            r_slda_mi.channels    = csp_mi.selected_channels;
            r_slda_mi.bands       = csp_mi.bands;
            r_slda_mi.n_bands     = csp_mi.n_bands;
            r_slda_mi.n_components = csp_mi.n_components;
            r_slda_mi.class_codes = classes_sorted(:)';
            fprintf('  [sLDA-MI]  %d/%d features selected by sLDA\n', ...
                    r_slda_mi.n_selected, csp_mi.n_components * csp_mi.n_bands);
        end
        if use_cvsa && ~isempty(csp_cvsa) && ~isempty(slda_cvsa)
            r_slda_cvsa = slda_spatial_weights(csp_cvsa, slda_cvsa);
            r_slda_cvsa.channels    = csp_cvsa.selected_channels;
            r_slda_cvsa.bands       = csp_cvsa.bands;
            r_slda_cvsa.n_bands     = csp_cvsa.n_bands;
            r_slda_cvsa.n_components = csp_cvsa.n_components;
            r_slda_cvsa.class_codes = classes_sorted(:)';
            fprintf('  [sLDA-CVSA] %d/%d features selected by sLDA\n', ...
                    r_slda_cvsa.n_selected, csp_cvsa.n_components * csp_cvsa.n_bands);
        end

        % ── Sample accuracy: per trial → mean over trials ──────────────────────
        % Formula: for each trial, (# CF frames where argmax(P)==target) / (# valid frames)
        %          then mean over all trials / HIT-only trials / MISS+TO trials.
        sa_mi = NaN; sa_cvsa = NaN; sa_fused = NaN;
        sa_mi_hit = NaN; sa_mi_miss = NaN;
        sa_cvsa_hit = NaN; sa_cvsa_miss = NaN;
        sa_fused_hit = NaN; sa_fused_miss = NaN;
        sa_mi_cls = nan(1,n_cls); sa_cvsa_cls = nan(1,n_cls); sa_fused_cls = nan(1,n_cls);

        if use_mi && ~use_cvsa        % pure MI
            [sa_mi,    sa_mi_cls]   = frame_acc_cls(trials,'raw', Inf, n_cls, []);
            sa_mi_hit  = frame_acc_cls(trials,'raw', Inf, n_cls, hit_mask);
            sa_mi_miss = frame_acc_cls(trials,'raw', Inf, n_cls, miss_mask);
        elseif use_cvsa && ~use_mi    % pure CVSA
            [sa_cvsa,    sa_cvsa_cls]   = frame_acc_cls(trials,'raw', n3s, n_cls, []);
            sa_cvsa_hit  = frame_acc_cls(trials,'raw', n3s, n_cls, hit_mask);
            sa_cvsa_miss = frame_acc_cls(trials,'raw', n3s, n_cls, miss_mask);
        else                          % hybrid
            [sa_mi,    sa_mi_cls]    = frame_acc_cls(trials,'mi',   Inf, n_cls, []);
            [sa_cvsa,  sa_cvsa_cls]  = frame_acc_cls(trials,'cvsa', n3s, n_cls, []);
            [sa_fused, sa_fused_cls] = frame_acc_cls(trials,'raw',  Inf, n_cls, []);
            sa_mi_hit    = frame_acc_cls(trials,'mi',   Inf, n_cls, hit_mask);
            sa_mi_miss   = frame_acc_cls(trials,'mi',   Inf, n_cls, miss_mask);
            sa_cvsa_hit  = frame_acc_cls(trials,'cvsa', n3s, n_cls, hit_mask);
            sa_cvsa_miss = frame_acc_cls(trials,'cvsa', n3s, n_cls, miss_mask);
            sa_fused_hit = frame_acc_cls(trials,'raw',  Inf, n_cls, hit_mask);
            sa_fused_miss= frame_acc_cls(trials,'raw',  Inf, n_cls, miss_mask);
        end

        fprintf('  SA all:     MI=%.1f%%  CVSA(%.0fs)=%.1f%%  Fused=%.1f%%\n', ...
                sa_mi*100, CVSA_MAX_S, sa_cvsa*100, sa_fused*100);
        fprintf('  SA HIT:     MI=%.1f%%  CVSA=%.1f%%  Fused=%.1f%%\n', ...
                sa_mi_hit*100, sa_cvsa_hit*100, sa_fused_hit*100);
        fprintf('  SA MISS/TO: MI=%.1f%%  CVSA=%.1f%%  Fused=%.1f%%\n', ...
                sa_mi_miss*100, sa_cvsa_miss*100, sa_fused_miss*100);
        for cc = 1:n_cls
            fprintf('  SA cls%d:    MI=%.1f%%  CVSA=%.1f%%  Fused=%.1f%%\n', cc, ...
                    sa_mi_cls(cc)*100, sa_cvsa_cls(cc)*100, sa_fused_cls(cc)*100);
        end

        % ── ERD/ERS (Pfurtscheller-style) for spatial validation figures ─────────
        %   Placed AFTER sample-accuracy so any failure cannot abort SA stats.
        %   Baseline = cue period [CF-1.5s, CF]. Result stored for Figs 6 & 7.
        classes_sorted_erd = sort(to_vec(int_cfg.classes));
        if use_mi && ~isempty(csp_mi)
            try
                r_erd_mi = compute_erd_ers(signal, header, csp_mi, ...
                    classes_sorted_erd, eog_names, fs);
                fprintf('  [ERD-MI]   c1=%d c2=%d trials\n', ...
                        r_erd_mi.n_trials_cls(1), r_erd_mi.n_trials_cls(2));
            catch ME2
                r_erd_mi = [];
                fprintf('  [ERD-MI]   skipped [%s] %s\n', ME2.identifier, ME2.message);
            end
        end
        if use_cvsa && ~isempty(csp_cvsa)
            try
                r_erd_cvsa = compute_erd_ers(signal, header, csp_cvsa, ...
                    classes_sorted_erd, eog_names, fs);
                fprintf('  [ERD-CVSA] c1=%d c2=%d trials\n', ...
                        r_erd_cvsa.n_trials_cls(1), r_erd_cvsa.n_trials_cls(2));
            catch ME2
                r_erd_cvsa = [];
                fprintf('  [ERD-CVSA] skipped [%s] %s\n', ME2.identifier, ME2.message);
            end
        end

    catch ME
        fprintf('  [warn] simulation error: %s\n', ME.message);
    end

    % ── Collect results ───────────────────────────────────────────────────────
    r = struct('file_label',file_label,'paradigm',paradigm,'col',col, ...
               'n_hit',n_hit,'n_miss',n_miss,'n_to',n_to, ...
               'hit_tot',hit_tot,'hit_no_to',hit_no_to, ...
               'tth_vals',tth_vals,'t_miss_vals',t_miss_vals,'to_vals',to_vals, ...
               'tth_mean',msafe(tth_vals),'t_miss_mean',msafe(t_miss_vals),'to_mean',msafe(to_vals), ...
               'sa_mi',sa_mi,'sa_cvsa',sa_cvsa,'sa_fused',sa_fused, ...
               'sa_mi_hit',sa_mi_hit,'sa_mi_miss',sa_mi_miss, ...
               'sa_cvsa_hit',sa_cvsa_hit,'sa_cvsa_miss',sa_cvsa_miss, ...
               'sa_fused_hit',sa_fused_hit,'sa_fused_miss',sa_fused_miss, ...
               'sa_mi_cls',sa_mi_cls,'sa_cvsa_cls',sa_cvsa_cls,'sa_fused_cls',sa_fused_cls);
    r.csp_mi          = r_csp_mi;
    r.csp_cvsa        = r_csp_cvsa;
    r.slda_mi_weights   = r_slda_mi;
    r.slda_cvsa_weights = r_slda_cvsa;
    r.erd_mi            = r_erd_mi;
    r.erd_cvsa          = r_erd_cvsa;
    RES{end+1} = r;
end

% ── Sort results by paradigm: MI → CVSA → Hybrid → Unknown ───────────────────
par_order = struct('mi',1,'cvsa',2,'hybrid',3,'unknown',4);
sort_keys = cellfun(@(r) par_order.(r.paradigm), RES);
[~, sort_idx] = sort(sort_keys);
RES   = RES(sort_idx);
n_f   = numel(RES);
xlbls = cellfun(@(r) r.file_label, RES, 'UniformOutput',false);

% Paradigm group boundaries for separator lines and group labels
par_names_sorted = cellfun(@(r) r.paradigm, RES, 'UniformOutput',false);
par_seq = unique(par_names_sorted, 'stable');
grp_start = zeros(1, numel(par_seq));
grp_end   = zeros(1, numel(par_seq));
for g = 1:numel(par_seq)
    idx_g = find(strcmp(par_names_sorted, par_seq{g}));
    grp_start(g) = min(idx_g);
    grp_end(g)   = max(idx_g);
end

%% ═════════════════════════════════════════════════════════════════════════════
%  SESSION SUMMARY — per-file means + per-paradigm experiment totals
%  ═════════════════════════════════════════════════════════════════════════════
fprintf('\n══════════════════ Session summary ══════════════════\n');
for g = 1:numel(par_seq)
    idx_g = grp_start(g):grp_end(g);
    fprintf('  %s:\n', upper(par_seq{g}));
    for i = idx_g
        r = RES{i};
        n_tot_r = max(1, r.n_hit + r.n_miss + r.n_to);
        fprintf('    %-12s HIT=%d MISS=%d TO=%d  acc=%.0f%%  no-to=%.0f%%  TTH=%.2fs  Tmiss=%.2fs  Tto=%.2fs\n', ...
                r.file_label, r.n_hit, r.n_miss, r.n_to, ...
                100*r.n_hit/n_tot_r, 100*r.hit_no_to, r.tth_mean, r.t_miss_mean, r.to_mean);
    end

    n_hit_tot  = sum(cellfun(@(r) r.n_hit,  RES(idx_g)));
    n_miss_tot = sum(cellfun(@(r) r.n_miss, RES(idx_g)));
    n_to_tot   = sum(cellfun(@(r) r.n_to,   RES(idx_g)));
    n_tot_g    = max(1, n_hit_tot + n_miss_tot + n_to_tot);

    tth_all = []; tmiss_all = []; to_all = [];
    for i = idx_g
        tth_all   = [tth_all,   RES{i}.tth_vals];   %#ok<AGROW>
        tmiss_all = [tmiss_all, RES{i}.t_miss_vals]; %#ok<AGROW>
        to_all    = [to_all,    RES{i}.to_vals];     %#ok<AGROW>
    end

    fprintf('    %-12s HIT=%d MISS=%d TO=%d  acc=%.0f%%  no-to=%.0f%%  TTH=%.2fs  Tmiss=%.2fs  Tto=%.2fs  (n=%d files)\n', ...
            'TOTAL', n_hit_tot, n_miss_tot, n_to_tot, ...
            100*n_hit_tot/n_tot_g, 100*n_hit_tot/max(1,n_hit_tot+n_miss_tot), ...
            msafe(tth_all), msafe(tmiss_all), msafe(to_all), numel(idx_g));
end
fprintf('═══════════════════════════════════════════════════════\n');

%% ═════════════════════════════════════════════════════════════════════════════
%  FIG 1 — TRIAL ACCURACY  (real GDF events)
%  ═════════════════════════════════════════════════════════════════════════════
fig1 = figure('Name','Trial Accuracy','Color','w','NumberTitle','off', ...
              'Position',[50 50 max(900,120*n_f) 480], 'Visible', fig_vis);
set(fig1, 'Units','normalized', 'OuterPosition',[0 0 1 1]);

metric_title = {'897 / (897+898+899)', '897 / (897+898)  [no timeout]'};
fields1 = {'hit_tot','hit_no_to'};

for sp = 1:2
    ax = subplot(1,2,sp); hold(ax,'on');
    add_group_decorations(ax, grp_start, grp_end, par_seq, COL, 115);

    % Per-file bars
    for i = 1:n_f
        r = RES{i};
        val = r.(fields1{sp}) * 100;
        bar(ax, i, val, 0.55, 'FaceColor',r.col,'EdgeColor','k','LineWidth',1.2,'HandleVisibility','off');
        text(ax, i, val+2, sprintf('%.0f%%\nH%d M%d T%d', val,r.n_hit,r.n_miss,r.n_to), ...
             'HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
    end

    % Per-paradigm mean line
    for g = 1:numel(par_seq)
        idx_g = grp_start(g):grp_end(g);
        vals_g = cellfun(@(r) r.(fields1{sp})*100, RES(idx_g));
        m = mean(vals_g,'omitnan');
        plot(ax,[idx_g(1)-0.4, idx_g(end)+0.4],[m,m],'--','Color',COL.(par_seq{g})*0.6, ...
             'LineWidth',2,'HandleVisibility','off');
        text(ax, idx_g(end)+0.55, m, sprintf('\\mu=%.0f%%',m), ...
             'FontSize',8,'Color',COL.(par_seq{g})*0.6,'FontWeight','bold');
    end

    % Legend (only sp==1)
    if sp==1
        bar(ax,NaN,NaN,'FaceColor',COL.mi,    'EdgeColor','k','DisplayName','MI');
        bar(ax,NaN,NaN,'FaceColor',COL.cvsa,  'EdgeColor','k','DisplayName','CVSA');
        bar(ax,NaN,NaN,'FaceColor',COL.hybrid,'EdgeColor','k','DisplayName','Hybrid');
        legend(ax,'MI','CVSA','Hybrid','FontSize',9,'Location','south');
    end

    set(ax,'XTick',1:n_f,'XTickLabel',xlbls,'XTickLabelRotation',25,'YLim',[0,115],'XLim',[0.3,n_f+0.7]);
    yline(ax,50,'--k','chance','FontSize',7,'HandleVisibility','off');
    grid(ax,'on'); ylabel('%');
    title(ax, ['Hit rate   ' metric_title{sp}], 'FontSize',11);
end
sgtitle(fig1,'Trial accuracy — real GDF outcomes (897/898/899)','FontSize',13);

%% ═════════════════════════════════════════════════════════════════════════════
%  FIG 2 — TIME TO HIT / TIME TO MISS  (real GDF event positions)
%  ═════════════════════════════════════════════════════════════════════════════
fig2 = figure('Name','Time Metrics','Color','w','NumberTitle','off', ...
              'Position',[60 60 max(900,120*n_f) 480], 'Visible', fig_vis);
set(fig2, 'Units','normalized', 'OuterPosition',[0 0 1 1]);

time_titles = {'Time to HIT  (897 trials)','Time to MISS  (898 trials)'};
time_fields  = {'tth_vals','t_miss_vals'};

for sp = 1:2
    ax = subplot(1,2,sp); hold(ax,'on');
    y_max = 0;
    for i = 1:n_f
        v = RES{i}.(time_fields{sp});
        if ~isempty(v), y_max = max(y_max, max(v)); end
    end
    y_top = max(y_max * 1.25, 1);
    add_group_decorations(ax, grp_start, grp_end, par_seq, COL, y_top);

    for i = 1:n_f
        r  = RES{i};
        vals = r.(time_fields{sp});
        if isempty(vals), continue; end
        m = mean(vals); s = std(vals);
        bar(ax, i, m, 0.5,'FaceColor',r.col,'EdgeColor','k','LineWidth',1.2,'HandleVisibility','off');
        errorbar(ax,i,m,s,'k.','LineWidth',1.5,'CapSize',8,'HandleVisibility','off');
        scatter(ax,i+randn(1,numel(vals))*0.07, vals, 28, r.col,'filled', ...
                'MarkerEdgeColor','k','LineWidth',0.5,'MarkerFaceAlpha',0.8,'HandleVisibility','off');
        text(ax,i,m+s+y_top*0.03, sprintf('%.1f±%.1fs\n(n=%d)',m,s,numel(vals)), ...
             'HorizontalAlignment','center','FontSize',7.5,'FontWeight','bold');
    end

    % Per-paradigm mean
    for g = 1:numel(par_seq)
        idx_g = grp_start(g):grp_end(g);
        vals_g = [];
        for i = idx_g, vals_g = [vals_g, RES{i}.(time_fields{sp})]; end %#ok<AGROW>
        if isempty(vals_g), continue; end
        m = mean(vals_g);
        plot(ax,[idx_g(1)-0.4,idx_g(end)+0.4],[m,m],'--','Color',COL.(par_seq{g})*0.6,'LineWidth',2,'HandleVisibility','off');
        text(ax,idx_g(end)+0.55,m,sprintf('\\mu=%.1fs',m),'FontSize',8,'Color',COL.(par_seq{g})*0.6,'FontWeight','bold');
    end

    set(ax,'XTick',1:n_f,'XTickLabel',xlbls,'XTickLabelRotation',25,'YLim',[0,y_top],'XLim',[0.3,n_f+0.7]);
    grid(ax,'on'); ylabel('seconds'); title(ax,time_titles{sp},'FontSize',11);
end
sgtitle(fig2,'Time metrics — real GDF event positions','FontSize',13);

%% ═════════════════════════════════════════════════════════════════════════════
%  FIG 3 — SAMPLE ACCURACY  (offline simulation, MATLAB ≈ ROS)
%  3 rows × 3 cols:  rows = MI / CVSA / Fused;  cols = ALL / HIT / MISS+TO
%  ═════════════════════════════════════════════════════════════════════════════
sa_row_titles = {sprintf('MI  (full CF)'), ...
                 sprintf('CVSA  (first %.0f s)',CVSA_MAX_S), ...
                 'Hybrid fused  (full CF)'};
sa_all_f   = {'sa_mi',      'sa_cvsa',      'sa_fused'};
sa_hit_f   = {'sa_mi_hit',  'sa_cvsa_hit',  'sa_fused_hit'};
sa_miss_f  = {'sa_mi_miss', 'sa_cvsa_miss', 'sa_fused_miss'};
col_hit    = [0.15 0.65 0.15];
col_miss   = [0.80 0.20 0.10];

fig3 = figure('Name','Sample Accuracy','Color','w','NumberTitle','off', ...
              'Position',[70 70 max(1300,150*n_f) 680], 'Visible', fig_vis);
set(fig3, 'Units','normalized', 'OuterPosition',[0 0 1 1]);

for row = 1:3
    for col_sp = 1:3   % col: 1=all 2=HIT 3=MISS/TO
        ax = subplot(3,3,(row-1)*3+col_sp); hold(ax,'on');
        add_group_decorations(ax, grp_start, grp_end, par_seq, COL, 110);

        if col_sp == 1,    fld = sa_all_f{row};  lab = 'all trials';
        elseif col_sp==2,  fld = sa_hit_f{row};  lab = 'HIT (897)';
        else,              fld = sa_miss_f{row};  lab = 'MISS+TO (898/899)';
        end
        bar_col_edge = [0 0 0];
        if col_sp==2, bar_col_edge = col_hit; end
        if col_sp==3, bar_col_edge = col_miss; end

        for i = 1:n_f
            r   = RES{i};
            val = r.(fld);
            if isnan(val), continue; end
            bface = r.col * (col_sp==1) + col_hit*(col_sp==2) + col_miss*(col_sp==3);
            b_alpha = 1.0 * (col_sp==1) + 0.85*(col_sp==2) + 0.70*(col_sp==3);
            bh = bar(ax,i,val*100,0.55,'FaceColor',bface,'FaceAlpha',b_alpha, ...
                     'EdgeColor',bar_col_edge,'LineWidth',1.2,'HandleVisibility','off');
            text(ax,i,val*100+2.5,sprintf('%.0f%%',val*100), ...
                 'HorizontalAlignment','center','FontSize',7.5,'FontWeight','bold');
        end

        % Per-paradigm mean line
        for g = 1:numel(par_seq)
            idx_g = grp_start(g):grp_end(g);
            vg = cellfun(@(r) r.(fld)*100, RES(idx_g));
            vg = vg(~isnan(vg));
            if isempty(vg), continue; end
            m = mean(vg);
            plot(ax,[idx_g(1)-0.4,idx_g(end)+0.4],[m m],'--', ...
                 'Color',COL.(par_seq{g})*0.6,'LineWidth',1.8,'HandleVisibility','off');
        end

        set(ax,'XTick',1:n_f,'XTickLabel',xlbls,'XTickLabelRotation',25, ...
               'YLim',[0,110],'XLim',[0.3,n_f+0.7]);
        yline(ax,50,'--k','FontSize',7,'HandleVisibility','off');
        grid(ax,'on'); ylabel('%');
        if row==1, title(ax,lab,'FontSize',10,'FontWeight','bold'); end
        if col_sp==1, text(ax,-0.1,0.5,sa_row_titles{row},'Units','normalized', ...
                           'HorizontalAlignment','right','FontSize',9,'FontWeight','bold', ...
                           'Rotation',0,'Color',[0.2 0.2 0.2]); end
    end
end
sgtitle(fig3,sprintf(['Sample accuracy — offline simulation  (MATLAB ≈ ROS)\n' ...
    'per-trial: (frames where argmax(P_{sLDA})==target) / total valid frames,  then mean over trials']), ...
    'FontSize',11);

%% --- Save figures -------------------------------------------------------
out_dir = fullfile(gdf_dir, 'analysis_results', 'overview');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
saveas(fig1, fullfile(out_dir, 'overview_trial_accuracy.svg'),  'svg');
saveas(fig2, fullfile(out_dir, 'overview_time_metrics.svg'),    'svg');
saveas(fig3, fullfile(out_dir, 'overview_sample_accuracy.svg'), 'svg');
fprintf('Saved figures to %s\n', out_dir);
if ~SHOW_FIGURES
    close(fig1); close(fig2); close(fig3);
end

%% ═════════════════════════════════════════════════════════════════════════════
%  FIG 4 — CSP CHANNEL IMPORTANCE
%
%   Spatial filter weights show which scalp regions and frequency bands drive
%   each class, using MNE alternate component ordering:
%     odd filter rows  (1,3,5,...) → class-1-dominant components
%     even filter rows (2,4,6,...) → class-2-dominant components
%
%   Panel layout (one row per CSP type: MI / CVSA):
%     col 1 — class-1 channel importance topoplot (parula, % of total weight)
%     col 2 — class-2 channel importance topoplot (parula, % of total weight)
%     col 3 — per-channel selectivity: (w1-w2)/(w1+w2), red=class1, blue=class2
%     col 4 — per-band weight contribution (grouped bars: class1 vs class2)
%  ═════════════════════════════════════════════════════════════════════════════
csp_type_seen  = {};
csp_data_rows  = {};
csp_row_prefix = {};
for fi2 = 1:numel(RES)
    r2 = RES{fi2};
    if ~isempty(r2.csp_mi) && ~ismember('MI', csp_type_seen)
        csp_type_seen{end+1}  = 'MI';
        csp_data_rows{end+1}  = r2.csp_mi;
        csp_row_prefix{end+1} = sprintf('MI CSP (from %s)', r2.file_label);
    end
    if ~isempty(r2.csp_cvsa) && ~ismember('CVSA', csp_type_seen)
        csp_type_seen{end+1}  = 'CVSA';
        csp_data_rows{end+1}  = r2.csp_cvsa;
        csp_row_prefix{end+1} = sprintf('CVSA CSP (from %s)', r2.file_label);
    end
end
n_csp_f = numel(csp_data_rows);

if n_csp_f > 0
    fig4 = figure('Name','CSP Channel Importance','Color','w', ...
                  'NumberTitle','off','Visible',fig_vis);
    set(fig4,'Units','normalized','Position',[0 0 1 1]);

    % Blue-white-red diverging colormap for selectivity panels
    n2  = 64;
    bwr = [linspace(0,1,n2)', linspace(0,1,n2)', ones(n2,1); ...   % blue→white
           ones(n2,1), linspace(1,0,n2)', linspace(1,0,n2)'];      % white→red

    for ci = 1:n_csp_f
        csp_d  = csp_data_rows{ci};
        cls    = csp_d.class_codes(:)';
        ch     = csp_d.channels;
        prefix = csp_row_prefix{ci};
        nb     = csp_d.n_bands;
        band_lbl = arrayfun(@(b) sprintf('%.0f-%.0f', ...
            csp_d.bands(b,1), csp_d.bands(b,2)), 1:nb, 'UniformOutput', false);

        % col 1: class-1 channel weight — full-scalp topoplot
        %   topo_map bg_zero=true fills all standard 10-20 positions with 0,
        %   giving smooth full-scalp interpolation even when only a subset was
        %   selected. Black-bold labels = selected; gray labels = standard rest.
        ax41 = subplot(n_csp_f, 4, (ci-1)*4 + 1);
        v1 = csp_d.w_c1 * 100;
        mx1 = max(v1(:)); if mx1 < 1e-6, mx1 = 1; end
        topo_map(ch, v1, [0, mx1], ax41, ...
            sprintf('%s | class %d — channel weight (%%)', prefix, cls(1)), true, true);
        colormap(ax41, parula);
        caxis(ax41, [0, mx1]);

        % col 2: class-2 channel weight — full-scalp topoplot
        ax42 = subplot(n_csp_f, 4, (ci-1)*4 + 2);
        v2 = csp_d.w_c2 * 100;
        mx2 = max(v2(:)); if mx2 < 1e-6, mx2 = 1; end
        topo_map(ch, v2, [0, mx2], ax42, ...
            sprintf('%s | class %d — channel weight (%%)', prefix, cls(2)), true, true);
        colormap(ax42, parula);
        caxis(ax42, [0, mx2]);

        % col 3: per-channel class selectivity — full-scalp topoplot
        %   selectivity = (w_c1 - w_c2) / (w_c1 + w_c2) ∈ [-1, +1]
        %   +1 = channel exclusively used by class-1 components
        %   -1 = channel exclusively used by class-2 components
        %    0 = equally shared
        ax43 = subplot(n_csp_f, 4, (ci-1)*4 + 3);
        topo_map(ch, csp_d.selectivity, [-1 1], ax43, ...
            sprintf('%s | selectivity (w_1-w_2)/(w_1+w_2)\n+1=c%d only, -1=c%d only', ...
                    prefix, cls(1), cls(2)), true, true);
        colormap(ax43, bwr);
        caxis(ax43, [-1 1]);

        % col 4: per-band importance — stacked bar chart
        %   Bar height = total |CSP filter| energy in that band (% of all bands).
        %   Stack split = how much comes from class-1-dominant vs class-2-dominant
        %   CSP components (odd rows → class 1, even rows → class 2).
        %   Tall bar = important frequency range; color split = which class drives it.
        ax44 = subplot(n_csp_f, 4, (ci-1)*4 + 4);
        hb   = bar(ax44, 1:nb, [csp_d.band_w_c1(:)*100, csp_d.band_w_c2(:)*100], 'stacked');
        hb(1).FaceColor = [0.20 0.45 0.80];   % blue   = class-1-dominant components
        hb(2).FaceColor = [0.85 0.45 0.10];   % orange = class-2-dominant components
        set(ax44, 'XTick', 1:nb, 'XTickLabel', band_lbl, ...
            'XTickLabelRotation', 45, 'FontSize', 7);
        ylabel(ax44, '% of total |CSP filter| weight', 'FontSize', 8);
        xlabel(ax44, 'Frequency band (Hz)', 'FontSize', 8);
        title(ax44, sprintf('%s\nband importance\nbar height = %% energy; color = class dominance', ...
              prefix), 'FontSize', 8, 'Interpreter', 'none');
        legend(ax44, sprintf('c%d components', cls(1)), sprintf('c%d components', cls(2)), ...
               'Location', 'northeast', 'FontSize', 7);
        box(ax44, 'off');
    end

    sgtitle(fig4, ['CSP spatial filter analysis — channel and band importance' newline ...
        'component order: alternate (odd rows = class-1, even rows = class-2)' newline ...
        'black labels = CSP-selected channels;  gray labels = all other standard 10-20'], ...
        'FontSize', 10);
    saveas(fig4, fullfile(out_dir, 'overview_csp_importance.svg'), 'svg');
    fprintf('Saved %s\n', fullfile(out_dir, 'overview_csp_importance.svg'));
    if ~SHOW_FIGURES, close(fig4); end
end

%% ═════════════════════════════════════════════════════════════════════════════
%  FIG 5 — sLDA-WEIGHTED CHANNEL AND BAND IMPORTANCE
%
%   While Fig 4 shows ALL CSP filter weights, Fig 5 is restricted to the
%   (band, component) features the sLDA actually selected and weighted by
%   |sLDA coefficient| — the direct spatial signature of the classifier's
%   decision function. Answers: "where on the scalp and in which frequency
%   band does the final classifier look, and how strongly?"
%
%   Panel layout (one row per CSP/sLDA type: MI / CVSA):
%     col 1 — channel importance topoplot: sum_k |coef_k| × |CSP_filter_k(ch)|
%              black-bold = CSP-selected channels; gray = rest of 10-20 grid
%     col 2 — band importance bar chart: same sum aggregated per band
%     col 3 — feature selection heatmap (comp × band):
%              color = |sLDA coef|; gray = not selected by sLDA
%  ═════════════════════════════════════════════════════════════════════════════
slda_type_seen  = {};
slda_data_rows  = {};
slda_row_prefix = {};
for fi2 = 1:numel(RES)
    r2 = RES{fi2};
    if ~isempty(r2.slda_mi_weights) && ~ismember('MI', slda_type_seen)
        slda_type_seen{end+1}  = 'MI';
        slda_data_rows{end+1}  = r2.slda_mi_weights;
        slda_row_prefix{end+1} = sprintf('MI sLDA (from %s)', r2.file_label);
    end
    if ~isempty(r2.slda_cvsa_weights) && ~ismember('CVSA', slda_type_seen)
        slda_type_seen{end+1}  = 'CVSA';
        slda_data_rows{end+1}  = r2.slda_cvsa_weights;
        slda_row_prefix{end+1} = sprintf('CVSA sLDA (from %s)', r2.file_label);
    end
end
n_slda_f = numel(slda_data_rows);

if n_slda_f > 0
    fig5 = figure('Name','sLDA Feature Importance','Color','w', ...
                  'NumberTitle','off','Visible',fig_vis);
    set(fig5,'Units','normalized','Position',[0 0 1 1]);

    for ci = 1:n_slda_f
        sd     = slda_data_rows{ci};
        ch     = sd.channels;
        cls    = sd.class_codes(:)';
        prefix = slda_row_prefix{ci};
        n_comp = sd.n_components;
        nb     = sd.n_bands;
        band_lbl = arrayfun(@(b) sprintf('%.0f-%.0f', ...
            sd.bands(b,1), sd.bands(b,2)), 1:nb, 'UniformOutput', false);

        % ── col 1: sLDA-weighted channel importance topoplot ─────────────────
        ax51 = subplot(n_slda_f, 3, (ci-1)*3 + 1);
        v_ch = sd.w_ch * 100;
        mx_ch = max(v_ch(:)); if mx_ch < 1e-6, mx_ch = 1; end
        topo_map(ch, v_ch, [0, mx_ch], ax51, ...
            sprintf('%s\nchannel importance (%%)', prefix), true, true);
        colormap(ax51, parula);
        caxis(ax51, [0, mx_ch]);

        % ── col 2: sLDA-weighted band importance bar chart ──────────────────
        ax52 = subplot(n_slda_f, 3, (ci-1)*3 + 2);
        bar(ax52, 1:nb, sd.band_w * 100, 'FaceColor', [0.25 0.50 0.75]);
        set(ax52, 'XTick', 1:nb, 'XTickLabel', band_lbl, ...
            'XTickLabelRotation', 45, 'FontSize', 7);
        ylabel(ax52, '% of total sLDA importance', 'FontSize', 8);
        xlabel(ax52, 'Frequency band (Hz)', 'FontSize', 8);
        title(ax52, sprintf('%s\nband importance\n(sum_k |coef_k|·|CSP_k|, per band)', prefix), ...
              'FontSize', 8, 'Interpreter', 'none');
        box(ax52, 'off');

        % ── col 3: feature selection heatmap (comp × band) ──────────────────
        %   Each cell (comp, band): color = |sLDA coef| if selected; gray if not.
        %   Odd rows = class-1-dominant; even rows = class-2-dominant (alternate order).
        ax53 = subplot(n_slda_f, 3, (ci-1)*3 + 3);
        fm      = sd.feat_mat;            % [n_comp × n_bands], NaN for unselected
        fm_disp = fm; fm_disp(isnan(fm_disp)) = 0;  % 0 for colormap (masked by AlphaData)
        im = imagesc(ax53, 1:nb, 1:n_comp, fm_disp);
        set(im, 'AlphaData', double(~isnan(fm)));
        set(ax53, 'Color', [0.88 0.88 0.88]);
        colormap(ax53, parula);
        colorbar(ax53, 'FontSize', 6);
        set(ax53, 'XTick', 1:nb, 'XTickLabel', band_lbl, 'XTickLabelRotation', 45, 'FontSize', 6);
        comp_lbl = cell(1, n_comp);
        for cc = 1:n_comp
            cl = cls(1 + mod(cc+1, 2));  % odd cc → cls(1), even cc → cls(2)
            comp_lbl{cc} = sprintf('c%d (↑cl%d)', cc, cl);
        end
        set(ax53, 'YTick', 1:n_comp, 'YTickLabel', comp_lbl, 'FontSize', 6);
        xlabel(ax53, 'Frequency band (Hz)', 'FontSize', 8);
        ylabel(ax53, 'CSP component', 'FontSize', 8);
        title(ax53, sprintf('%s\nfeature selection: %d / %d used\ngray = excluded by sLDA, color = |coef|', ...
              prefix, sd.n_selected, n_comp * nb), ...
              'FontSize', 8, 'Interpreter', 'none');
    end

    sgtitle(fig5, ['sLDA feature selection — spatial and frequency importance' newline ...
        'channel importance = Σ_k |coef_k| · |CSP filter_k(channel)| over sLDA-selected features' newline ...
        'black labels = CSP-selected channels;  gray labels = rest of standard 10-20'], ...
        'FontSize', 10);
    saveas(fig5, fullfile(out_dir, 'overview_slda_importance.svg'), 'svg');
    fprintf('Saved %s\n', fullfile(out_dir, 'overview_slda_importance.svg'));
    if ~SHOW_FIGURES, close(fig5); end
end

%% ═════════════════════════════════════════════════════════════════════════════
%  FIG 6 — ERD/ERS CLASS DISCRIMINATION vs CSP FILTER WEIGHT (CORRELATION)
%
%   Claim: channels weighted heavily by the CSP spatial filter also show
%   stronger ERD/ERS class discrimination during continuous feedback.
%
%   Per panel (one per CSP type × frequency band):
%     x — CSP filter weight: sum(|filter rows|) for this band, per channel
%     y — ERD/ERS class discrimination: |mean_ERD(c1) - mean_ERD(c2)| during CF
%     r — Pearson correlation coefficient across CSP-selected channels
%
%   Uses first file per CSP type (consistent with Figs 4–5).
%  ═════════════════════════════════════════════════════════════════════════════
erd6_type  = {};   % 'MI' / 'CVSA'
erd6_erd   = {};   % erd struct
erd6_csp   = {};   % csp struct

for fi2 = 1:numel(RES)
    r2 = RES{fi2};
    if ~isempty(r2.erd_mi)   && ~isempty(r2.csp_mi)   && ~ismember('MI',   erd6_type)
        erd6_type{end+1} = 'MI';  erd6_erd{end+1} = r2.erd_mi;  erd6_csp{end+1} = r2.csp_mi;
    end
    if ~isempty(r2.erd_cvsa) && ~isempty(r2.csp_cvsa) && ~ismember('CVSA', erd6_type)
        erd6_type{end+1} = 'CVSA'; erd6_erd{end+1} = r2.erd_cvsa; erd6_csp{end+1} = r2.csp_cvsa;
    end
end
n_erd6 = numel(erd6_erd);

if n_erd6 > 0
    nb_max6 = max(cellfun(@(c) c.n_bands, erd6_csp));

    fig6 = figure('Name','ERD/ERS vs CSP Weight','Color','w', ...
                  'NumberTitle','off','Visible',fig_vis);
    set(fig6,'Units','normalized','OuterPosition',[0 0 1 1]);

    for ri = 1:n_erd6
        ed  = erd6_erd{ri};   % erd struct
        cd  = erd6_csp{ri};   % csp struct
        ch  = ed.ch_names;
        nb  = cd.n_bands;
        typ = erd6_type{ri};

        for b = 1:nb
            ax = subplot(n_erd6, nb_max6, (ri-1)*nb_max6 + b);
            hold(ax, 'on');

            % CSP weight per channel: sum |W| over all components, normalised
            W     = cd.csp_matrices{b};          % [n_comp x n_sel]
            csp_w = (sum(abs(W), 1))';           % [n_sel x 1]
            csp_w = csp_w / max(csp_w(:) + eps);

            % ERD/ERS discrimination: |mean_c1 - mean_c2| during CF
            % mean_erd_cf: [n_sel x n_bands x 2]
            e1      = ed.mean_erd_cf(:, b, 1);  % [n_sel x 1]
            e2      = ed.mean_erd_cf(:, b, 2);  % [n_sel x 1]
            discrim = abs(e1 - e2);

            ok = ~isnan(csp_w) & ~isnan(discrim) & isfinite(csp_w) & isfinite(discrim);
            r_val = NaN;
            if sum(ok) >= 3
                C = corrcoef(csp_w(ok), discrim(ok));
                r_val = C(1,2);
            end

            scatter(ax, csp_w(ok), discrim(ok), 40, [0.15 0.40 0.75], 'filled', ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
            % Linear fit line
            if sum(ok) >= 3
                p_fit = polyfit(csp_w(ok), discrim(ok), 1);
                xfit = linspace(0, 1, 50);
                plot(ax, xfit, polyval(p_fit, xfit), '--', 'Color', [0.5 0.5 0.5], ...
                     'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
            % Channel labels
            ch_ok = ch(ok);
            csp_ok = csp_w(ok); disc_ok = discrim(ok);
            for ci = 1:numel(ch_ok)
                text(ax, csp_ok(ci), disc_ok(ci), ['  ' strtrim(ch_ok{ci})], ...
                     'FontSize', 5.5, 'Interpreter', 'none', 'Color', [0.2 0.2 0.2]);
            end

            band_str = sprintf('%.0f-%.0f Hz', cd.bands(b,1), cd.bands(b,2));
            if isnan(r_val)
                r_str = 'r=n/a';
            else
                r_str = sprintf('r = %.2f  (n=%d ch)', r_val, sum(ok));
            end
            title(ax, sprintf('%s | %s\n%s', typ, band_str, r_str), ...
                  'FontSize', 8, 'Interpreter', 'none');
            xlabel(ax, 'CSP filter weight (norm.)', 'FontSize', 7);
            if b == 1
                ylabel(ax, '|ERD(c1) – ERD(c2)|  (%)', 'FontSize', 7);
            end
            box(ax, 'off');  grid(ax, 'on');
        end
    end

    sgtitle(fig6, ['ERD/ERS class discrimination vs CSP spatial filter weight' newline ...
                   'Claim: channels with high CSP weight show stronger ERD/ERS class separation' newline ...
                   '(Pfurtscheller ERD/ERS during CF; Pearson r over CSP-selected channels)'], ...
            'FontSize', 10);
    saveas(fig6, fullfile(out_dir, 'overview_erd_csp_corr.svg'), 'svg');
    fprintf('Saved %s\n', fullfile(out_dir, 'overview_erd_csp_corr.svg'));
    if ~SHOW_FIGURES, close(fig6); end
end

%% ═════════════════════════════════════════════════════════════════════════════
%  FIG 7 — HEMISPHERE ROI ERD/ERS TIME-COURSE
%
%   CSP-selected channels grouped into Left (odd-digit), Right (even-digit),
%   Midline (trailing z) ROIs. ERD/ERS time-course averaged over ROI channels
%   and over all frequency bands, from cue onset (-1.5s) to max CF end (5s).
%
%   Left col  — per-class traces c1(t) and c2(t) with ± across-channel SEM
%   Right col — lateralization c1(t) – c2(t)
%
%   Vertical lines: solid gray = CF onset (t=0), dashed blue = CVSA influence
%   window boundary (3s). Uses first file per CSP type.
%  ═════════════════════════════════════════════════════════════════════════════
erd7_type = {};
erd7_erd  = {};

for fi2 = 1:numel(RES)
    r2 = RES{fi2};
    if ~isempty(r2.erd_mi)   && ~ismember('MI',   erd7_type)
        erd7_type{end+1} = 'MI';   erd7_erd{end+1} = r2.erd_mi;
    end
    if ~isempty(r2.erd_cvsa) && ~ismember('CVSA', erd7_type)
        erd7_type{end+1} = 'CVSA'; erd7_erd{end+1} = r2.erd_cvsa;
    end
end
n_erd7 = numel(erd7_erd);

ROI_NAMES  = {'left', 'right', 'midline'};
ROI_LABELS = {'Left (odd)', 'Right (even)', 'Midline (z)'};
COL_CLS7   = {[0.15 0.50 0.85],  [0.85 0.35 0.10]};   % c1=blue  c2=orange
CVSA_INF_S = 3.0;   % CVSA influence window

if n_erd7 > 0
    % Count non-empty ROI × CSP-type combinations
    has_roi = false(n_erd7, numel(ROI_NAMES));
    roi_ch_cells = cell(n_erd7, numel(ROI_NAMES));  % channel names per ROI
    for ri = 1:n_erd7
        ch = erd7_erd{ri}.ch_names;
        for k = 1:numel(ROI_NAMES)
            roi_mask = strcmp(cellfun(@hemisphere_of, ch(:)', 'UniformOutput', false), ROI_NAMES{k});
            has_roi(ri,k) = any(roi_mask);
            roi_ch_cells{ri,k} = ch(roi_mask);
        end
    end
    n_rows7 = sum(has_roi(:));

    if n_rows7 == 0
        fprintf('  [ROI] no hemisphere ROI channels found, skipping Fig 7\n');
    else
        fig7 = figure('Name','ROI ERD/ERS Time-course','Color','w', ...
                      'NumberTitle','off','Visible',fig_vis);
        set(fig7,'Units','normalized','OuterPosition',[0 0 1 1]);

        row7 = 0;
        for ri = 1:n_erd7
            ed    = erd7_erd{ri};
            tv    = ed.time_vec;              % [1 x T]
            T7    = numel(tv);
            n_cue7 = ed.n_cue;
            ctype = erd7_type{ri};
            ch    = ed.ch_names;

            % Band-average: erd_time [T x n_sel x n_bands x n_cls] → [T x n_sel x n_cls]
            n_cls7  = size(ed.erd_time, 4);
            n_sel7  = size(ed.erd_time, 2);
            et_avg  = reshape(nanmean(ed.erd_time, 3), T7, n_sel7, n_cls7);

            for k = 1:numel(ROI_NAMES)
                if ~has_roi(ri,k), continue; end
                row7 = row7 + 1;
                roi_mask = strcmp(cellfun(@hemisphere_of, ch(:)', 'UniformOutput', false), ROI_NAMES{k});
                n_roi = sum(roi_mask);

                ax1 = subplot(n_rows7, 2, (row7-1)*2 + 1);  hold(ax1,'on');
                ax2 = subplot(n_rows7, 2, (row7-1)*2 + 2);  hold(ax2,'on');

                lat_mean = nan(T7, n_cls7);
                for cls = 1:n_cls7
                    roi_data = et_avg(:, roi_mask, cls);  % [T x n_roi]
                    mn  = nanmean(roi_data, 2);            % [T x 1]
                    sem = nanstd(roi_data, 0, 2) / sqrt(max(n_roi, 1));
                    lat_mean(:, cls) = mn;

                    c_col = COL_CLS7{min(cls, numel(COL_CLS7))};
                    n_tr  = ed.n_trials_cls(min(cls, numel(ed.n_trials_cls)));
                    fill(ax1, [tv, fliplr(tv)], [(mn+sem)', fliplr((mn-sem)')], ...
                         c_col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility','off');
                    plot(ax1, tv, mn', '-', 'Color', c_col, 'LineWidth', 1.5, ...
                         'DisplayName', sprintf('class %d  (n=%d trials)', cls, n_tr));
                end

                % Lateralization c1 – c2
                if n_cls7 >= 2
                    lat = lat_mean(:,1) - lat_mean(:,2);
                    plot(ax2, tv, lat', 'k-', 'LineWidth', 1.5);
                    yline(ax2, 0, '--k', 'LineWidth', 0.8, 'HandleVisibility','off');
                end

                % Common decorations
                for ax_k = [ax1, ax2]
                    xline(ax_k, 0,            '-',  'Color', [0.55 0.55 0.55], 'LineWidth', 1.2, 'HandleVisibility','off');
                    xline(ax_k, CVSA_INF_S,   '--', 'Color', [0.30 0.45 0.80], 'LineWidth', 0.8, 'HandleVisibility','off');
                    xline(ax_k, -1.5,          ':',  'Color', [0.65 0.65 0.65], 'LineWidth', 0.8, 'HandleVisibility','off');
                    xlim(ax_k, [tv(1), tv(end)]);
                    xlabel(ax_k, 'time from CF onset (s)', 'FontSize', 7);
                    grid(ax_k, 'on'); box(ax_k, 'off');
                end
                ylabel(ax1, 'ERD/ERS (%)',           'FontSize', 7);
                ylabel(ax2, 'lateralization c1–c2 (%)','FontSize', 7);

                roi_ch_str = strjoin(roi_ch_cells{ri,k}, ' ');
                title(ax1, sprintf('%s | %s   [%s]  (%d ch, band-avg)', ...
                      ctype, ROI_LABELS{k}, roi_ch_str, n_roi), ...
                      'FontSize', 7.5, 'Interpreter', 'none');
                title(ax2, 'lateralization  (blue dashed = CVSA influence window)', ...
                      'FontSize', 7.5, 'Interpreter', 'none');
                legend(ax1, 'Location', 'best', 'FontSize', 6.5);
            end
        end

        sgtitle(fig7, ['ERD/ERS hemisphere ROI time-courses  (band-averaged, Pfurtscheller)' newline ...
                       'Baseline = cue period [–1.5s, 0s] | gray line = CF onset | dotted = cue onset | blue dashed = CVSA window (3s)'], ...
                'FontSize', 10);
        saveas(fig7, fullfile(out_dir, 'overview_roi_timecourse.svg'), 'svg');
        fprintf('Saved %s\n', fullfile(out_dir, 'overview_roi_timecourse.svg'));
        if ~SHOW_FIGURES, close(fig7); end
    end
end

fprintf('\nDone.\n');

%% ── Local functions ──────────────────────────────────────────────────────────

function par = detect_paradigm(fname)
    s = lower(fname);
    if contains(s,'hybrid'),   par = 'hybrid';
    elseif contains(s,'cvsa'), par = 'cvsa';
    elseif contains(s,'mi'),   par = 'mi';
    else,                      par = 'unknown';
    end
end

function [sa, sa_cls] = frame_acc_cls(trials, modality, n_max, n_cls, outcome_mask)
% FRAME_ACC_CLS  Per-trial sample accuracy, total and per-class.
%   For each trial: (# CF frames where argmax(P_sLDA)==target_class) / (# valid frames).
%   Then mean over selected trials.
%
%   modality:     'raw' (native/fused) | 'mi' | 'cvsa'
%   n_max:        max CF frames to consider (Inf = full trial; use round(3*fr) for CVSA)
%   outcome_mask: logical vector selecting which trials to include ([] = all)
%
%   Class 1 = lower class code (769/730/750), class 2 = higher code (770/731/751).
%   Column 1 of P_sLDA = P(class 1), column 2 = P(class 2). Consistent with apply_slda.

    n_tr = numel(trials);
    if nargin < 5 || isempty(outcome_mask)
        outcome_mask = true(1, n_tr);
    end
    outcome_mask = logical(outcome_mask(:)');
    if numel(outcome_mask) < n_tr
        outcome_mask(end+1:n_tr) = false;
    end

    accs   = nan(n_tr, 1);
    cls_ok = nan(n_tr, 1);

    for t = 1:n_tr
        if ~outcome_mask(t), continue; end
        c  = trials(t).target_class;
        np = trials(t).n_pre;
        if isnan(c) || c < 1, continue; end
        n_use = min(trials(t).n_cf, n_max);
        if n_use < 1, continue; end
        switch modality
            case 'raw',  raw = trials(t).raw(np+1:np+n_use, :);
            case 'mi',   raw = trials(t).p_mi(np+1:np+n_use, :);
            case 'cvsa', raw = trials(t).p_cvsa(np+1:np+n_use, :);
            otherwise,   raw = trials(t).raw(np+1:np+n_use, :);
        end
        vld = ~any(isnan(raw), 2);
        if ~any(vld), continue; end
        % Correct frame: target class column has the highest probability
        accs(t)   = mean(raw(vld, c) == max(raw(vld, :), [], 2));
        cls_ok(t) = c;
    end

    sa     = mean(accs, 'omitnan');
    sa_cls = nan(1, n_cls);
    for c = 1:n_cls
        idx_c = cls_ok == c;
        if any(idx_c), sa_cls(c) = mean(accs(idx_c), 'omitnan'); end
    end
end

function v = msafe(x)
    if isempty(x) || all(isnan(x(:))), v = NaN; else, v = mean(x(:),'omitnan'); end
end

%% ═════════════════════════════════════════════════════════════════════════════
%  SHARED HELPER: add separators and group labels to an axes
%  ═════════════════════════════════════════════════════════════════════════════
function add_group_decorations(ax, grp_start, grp_end, par_seq, COL, ylim_top)
    % Shade background per group
    for g = 1:numel(par_seq)
        x0 = grp_start(g) - 0.5;  x1 = grp_end(g) + 0.5;
        bg = COL.(par_seq{g}) * 0.12 + 0.88;   % very light tint
        patch(ax,[x0,x1,x1,x0],[0,0,ylim_top,ylim_top], bg, ...
              'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');
        % Group label at top
        text(ax, (x0+x1)/2, ylim_top*0.97, upper(par_seq{g}), ...
             'HorizontalAlignment','center','FontSize',10,'FontWeight','bold', ...
             'Color', COL.(par_seq{g})*0.7, 'VerticalAlignment','top');
    end
    % Separator lines between groups
    for g = 1:numel(par_seq)-1
        xline(ax, grp_end(g)+0.5, '--','Color',[0.5 0.5 0.5],'LineWidth',1, ...
              'HandleVisibility','off');
    end
end

function [w_c1, w_c2, band_w_c1, band_w_c2, selectivity] = csp_channel_weights(csp)
%CSP_CHANNEL_WEIGHTS  Per-channel and per-band importance from CSP filter matrices.
%   MNE component_order='alternate': odd filter rows (1,3,...) maximise class-1
%   variance; even rows (2,4,...) maximise class-2 variance.
%   w_c1 and w_c2 are each normalised independently to sum to 1.
%   band_w_c1/band_w_c2 are jointly normalised (sum over both = 1).
    n_bands = csp.n_bands;
    n_sel   = numel(csp.selected_channels);
    w1_raw = zeros(1, n_sel);
    w2_raw = zeros(1, n_sel);
    bw1    = zeros(1, n_bands);
    bw2    = zeros(1, n_bands);
    for b = 1:n_bands
        W      = csp.csp_matrices{b};       % [n_comp x n_sel]
        n_comp = size(W, 1);
        idx1   = 1:2:n_comp;                % class-1 dominant (alternate)
        idx2   = 2:2:n_comp;                % class-2 dominant
        a1     = sum(abs(W(idx1, :)), 1);
        a2     = sum(abs(W(idx2, :)), 1);
        w1_raw = w1_raw + a1;
        w2_raw = w2_raw + a2;
        bw1(b) = sum(a1);
        bw2(b) = sum(a2);
    end
    w_c1       = w1_raw / max(sum(w1_raw), eps);
    w_c2       = w2_raw / max(sum(w2_raw), eps);
    tot        = sum(bw1) + sum(bw2);
    band_w_c1  = bw1 / max(tot, eps);
    band_w_c2  = bw2 / max(tot, eps);
    selectivity = (w_c1 - w_c2) ./ max(w_c1 + w_c2, eps);
end

function top = top_channels(ch_names, weights, n)
%TOP_CHANNELS  Return the n channel names with the highest weights.
    [~, ord] = sort(weights(:)', 'descend');
    top = ch_names(ord(1:min(n, numel(ch_names))));
end

function erd = compute_erd_ers(signal, header, csp, class_codes, eog_names, fs)
%COMPUTE_ERD_ERS  Pfurtscheller-style ERD/ERS for CSP-selected channels.
%
%   Epochs each CF trial from (CF - 1.5s) to (CF + 5s).
%   Baseline = cue period [CF-1.5s, CF] (= 1.5s immediately preceding CF onset).
%   ERD/ERS = (power - baseline) / baseline * 100.
%   Power = 200ms moving average of instantaneous squared signal.
%   Filtering matches apply_processing.m: causal LP then HP Butterworth.
%
%   Returns struct:
%     .time_vec      [1 x T]                   s relative to CF onset
%     .mean_erd_cf   [n_sel x n_bands x n_cls] mean ERD during CF per class
%     .erd_time      [T x n_sel x n_bands x n_cls] time-course (NaN = no data)
%     .ch_names      {1 x n_sel}               channel names
%     .n_cue         scalar                    # cue samples
%     .n_trials_cls  [1 x n_cls]               trials per class

    CUE_S    = 1.5;
    CF_MAX_S = 5.0;
    MA_S     = 0.2;
    CF_EV    = 781;
    OUT_EVS  = [897, 898, 899];

    [n_samp, n_ch_all] = size(signal);
    all_labels = to_strcell(header.Label);
    n_cls = numel(class_codes);

    % Resolve CSP-selected channel indices (safe: skip unresolved)
    sel_names = csp.selected_channels;
    n_all = numel(sel_names);
    sel_idx = zeros(1, n_all);
    for i = 1:n_all
        m = find(strcmpi(all_labels, strtrim(sel_names{i})), 1);
        if ~isempty(m), sel_idx(i) = m; end
    end
    valid    = sel_idx > 0;
    sel_idx  = sel_idx(valid);
    sel_names = sel_names(valid);
    n_sel    = numel(sel_idx);
    if n_sel == 0
        wanted = to_strcell(csp.selected_channels);
        fprintf('  [ERD] no CSP channels matched GDF labels — wanted: %s\n', ...
                strjoin(wanted(:)', ', '));
        fprintf('  [ERD] GDF has: %s\n', strjoin(all_labels(:)', ', '));
        erd = []; return;
    end

    % EOG indices for CAR (matches apply_processing.m)
    eog_idx = [];
    for i = 1:numel(eog_names)
        m = find(strcmpi(all_labels, strtrim(eog_names{i})), 1);
        if ~isempty(m), eog_idx(end+1) = m; end  %#ok<AGROW>
    end
    non_eog = setdiff(1:n_ch_all, eog_idx);

    % CAR → keep only CSP-selected channels (processed one-shot, not chunk-based)
    car_mean = mean(signal(:, non_eog), 2);
    sig_sel  = signal(:, sel_idx) - car_mean;   % [n_samp x n_sel]

    % Trial discovery: CF events + preceding class events
    POS = header.EVENT.POS;
    TYP = header.EVENT.TYP;
    cf_all  = POS(TYP == CF_EV);
    n_cf_ev = numel(cf_all);
    n_cue   = round(CUE_S * fs);
    n_cf    = round(CF_MAX_S * fs);
    T       = n_cue + n_cf;

    trial_cls    = nan(1, n_cf_ev);
    trial_cf_dur = nan(1, n_cf_ev);
    for t = 1:n_cf_ev
        p0 = cf_all(t);
        if p0 - n_cue < 1 || p0 + n_cf - 1 > n_samp, continue; end
        pre = find(POS < p0 & ismember(TYP, class_codes), 1, 'last');
        if isempty(pre), continue; end
        ci = find(class_codes == TYP(pre), 1);
        if isempty(ci), continue; end
        trial_cls(t) = ci;
        post = find(POS > p0 & ismember(TYP, OUT_EVS), 1);
        if ~isempty(post)
            trial_cf_dur(t) = min(POS(post) - p0, n_cf);
        else
            trial_cf_dur(t) = n_cf;
        end
    end
    n_trials_cls = arrayfun(@(c) sum(trial_cls == c), 1:n_cls);

    % Running sum/count (one band at a time to limit peak memory)
    n_bands = csp.n_bands;
    ep_sum  = zeros(T, n_sel, n_bands, n_cls);
    ep_cnt  = zeros(T, n_cls);

    % Count valid time points per trial (same across bands)
    for t = 1:n_cf_ev
        if isnan(trial_cls(t)), continue; end
        vd = n_cue + trial_cf_dur(t);
        ci = trial_cls(t);
        ep_cnt(1:vd, ci) = ep_cnt(1:vd, ci) + 1;
    end

    % Filter + accumulate per band
    nyq  = fs / 2;
    n_ma = max(1, round(MA_S * fs));
    for b = 1:n_bands
        lo = csp.bands(b,1);  hi = csp.bands(b,2);
        [b_lp, a_lp] = butter(4, hi/nyq, 'low');
        [b_hp, a_hp] = butter(4, lo/nyq, 'high');
        tmp   = filter(b_lp, a_lp, sig_sel);
        tmp   = filter(b_hp, a_hp, tmp);
        pwr_b = movmean(tmp .^ 2, n_ma, 1);   % [n_samp x n_sel], 200ms MA
        clear tmp;

        for t = 1:n_cf_ev
            if isnan(trial_cls(t)), continue; end
            p0  = cf_all(t);
            ep0 = p0 - n_cue;
            ep_pwr = pwr_b(ep0 : p0 + n_cf - 1, :);  % [T x n_sel]
            base   = mean(ep_pwr(1:n_cue, :), 1);      % [1 x n_sel]
            erd_ep = (ep_pwr - base) ./ (abs(base) + eps) * 100;
            vd  = n_cue + trial_cf_dur(t);
            ci  = trial_cls(t);
            ep_sum(1:vd, :, b, ci) = ep_sum(1:vd, :, b, ci) + erd_ep(1:vd, :);
        end
        clear pwr_b;
    end
    clear sig_sel;

    % Compute mean (NaN where count = 0)
    erd_time = nan(T, n_sel, n_bands, n_cls);
    for cls = 1:n_cls
        cnt = max(ep_cnt(:, cls), 1);   % avoid /0; zero rows set to NaN below
        for b = 1:n_bands
            tmp = ep_sum(:, :, b, cls) ./ cnt;
            tmp(ep_cnt(:, cls) == 0, :) = NaN;
            erd_time(:, :, b, cls) = tmp;
        end
    end

    % Mean during CF window (t > 0)
    cf_samp     = (n_cue+1):T;
    mean_erd_cf = reshape(nanmean(erd_time(cf_samp,:,:,:), 1), n_sel, n_bands, n_cls);
    time_vec    = ((-n_cue):(n_cf-1)) / fs;

    erd = struct('time_vec',     time_vec, ...
                 'mean_erd_cf',  mean_erd_cf, ...
                 'erd_time',     erd_time, ...
                 'ch_names',     {sel_names(:)'}, ...
                 'n_cue',        n_cue, ...
                 'n_trials_cls', n_trials_cls);
end

function roi = hemisphere_of(ch)
%HEMISPHERE_OF  Assign a 10-20 channel to a hemisphere ROI.
%   Returns 'left' (odd trailing digit), 'right' (even trailing digit),
%   'midline' (trailing z/Z), or '' (unknown/non-standard).
    ch   = strtrim(lower(ch));
    last = ch(end);
    if last == 'z'
        roi = 'midline';
    elseif last >= '1' && last <= '9'
        if mod(str2double(last), 2) == 1, roi = 'left'; else, roi = 'right'; end
    else
        roi = '';
    end
end

function r = slda_spatial_weights(csp, slda)
%SLDA_SPATIAL_WEIGHTS  Channel and band importance weighted by |sLDA coefficient|.
%
%   For each sLDA-selected feature k (= a specific CSP component in a specific
%   frequency band), contributes |coef_k| × |CSP_filter_row_k| to the per-channel
%   importance map. This gives the spatial signature of the classifier's decision.
%
%   Band-major feature vector: index k (1-based) → band = ceil(k/n_comp),
%   component = mod(k-1, n_comp)+1. If feature selection is active, slda.weights
%   is already indexed in selection order (w(i) belongs to selected_feature_indices(i)).
    n_comp  = csp.n_components;
    n_bands = csp.n_bands;
    n_sel_ch = numel(csp.selected_channels);

    if isempty(slda.selected_feature_indices)
        sel_idx = (1:(n_comp * n_bands))';   % no FS: all features used
    else
        sel_idx = slda.selected_feature_indices(:);  % 1-based, FS active
    end

    w_ch     = zeros(1, n_sel_ch);
    band_w   = zeros(1, n_bands);
    feat_mat = nan(n_comp, n_bands);   % NaN = not selected

    for fi = 1:numel(sel_idx)
        k    = sel_idx(fi);
        band = ceil(k / n_comp);
        comp = mod(k-1, n_comp) + 1;
        if band < 1 || band > n_bands || comp < 1 || comp > n_comp, continue; end
        lda_w = abs(slda.weights(fi));
        filt  = abs(csp.csp_matrices{band}(comp, :));
        w_ch            = w_ch + lda_w * filt;
        band_w(band)    = band_w(band) + lda_w * sum(filt);
        feat_mat(comp, band) = lda_w;
    end

    sw = sum(w_ch);   if sw > 0, w_ch   = w_ch   / sw;  end
    sb = sum(band_w); if sb > 0, band_w = band_w / sb;   end

    r = struct('w_ch', w_ch, 'band_w', band_w, 'feat_mat', feat_mat, ...
               'n_selected', numel(sel_idx));
end
