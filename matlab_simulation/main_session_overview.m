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
              'Position',[50 50 max(900,120*n_f) 480]);

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
              'Position',[60 60 max(900,120*n_f) 480]);

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
              'Position',[70 70 max(1300,150*n_f) 680]);

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
out_dir = fullfile(gdf_dir, 'analysis_results');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
exportgraphics(fig1, fullfile(out_dir, 'overview_trial_accuracy.png'), 'Resolution', 150);
exportgraphics(fig2, fullfile(out_dir, 'overview_time_metrics.png'),   'Resolution', 150);
exportgraphics(fig3, fullfile(out_dir, 'overview_sample_accuracy.png'),'Resolution', 150);
fprintf('Saved figures to %s\n', out_dir);

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
