%% MAIN_EVALUATE_SINGLE  Per-file detailed evaluation: trials + metrics + topoplots.
%
%   Picks one or more GDF recording(s) via GUI. For each file:
%     1. Runs the full pipeline (identical to main_simulate).
%     2. Shows all CF trials in one figure (plot_trials).
%     3. Computes and plots a metric summary:
%          trial accuracy (simulated / event-based 897-898 / no-artifact-reject),
%          per-class hit rate, frame-level sample accuracy (total + mean/trial),
%          mean time-to-hit per class, mean time of peak output on miss trials,
%          artifact rate.
%     4. Plots ERD/ERS band-power time courses (selected channels × CF time,
%          averaged across trials, one row per band, class 1 vs class 2).
%     5. Plots toposcatter maps of mean band power per class × band.
%     6. Saves figures 1-4 as SVG to <gdf_dir>/results/<basename>/.
%     7. Saves eval_single_<paradigm>_<basename>.mat (used by main_compare_sessions).
%
%   Event codes assumed: 897 = hit, 898 = miss/timeout.

clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir,'io'), fullfile(this_dir,'processing'), ...
        fullfile(this_dir,'artifacts'), fullfile(this_dir,'classifier'), ...
        fullfile(this_dir,'integrator'), fullfile(this_dir,'plotting'), ...
        fullfile(this_dir,'utils'));

HIT_CODE      = 897;
MISS_CODE     = 898;
TIMEOUT_CODE  = 899;
CF_CODE       = 781;
ART_REJECT_THR = 0.30;   % trial rejected when > 30 % of CF frames have artifact

% ── File picker ──────────────────────────────────────────────────────────────
default_dir = '/home/paolo/bci_vr_ws/recordings';
if ~isfolder(default_dir), default_dir = fileparts(this_dir); end
[gdf_names, gdf_dir] = uigetfile({'*.gdf','GDF recordings (*.gdf)'}, ...
    'Select GDF file(s)', default_dir, 'MultiSelect','on');
if isequal(gdf_names,0), error('main_evaluate_single:cancel','No file selected.'); end
if ischar(gdf_names), gdf_names = {gdf_names}; end
n_files = numel(gdf_names);

% ── Per-file loop ─────────────────────────────────────────────────────────────
for fi = 1:n_files
    gdf_path = fullfile(gdf_dir, gdf_names{fi});
    [~, basename] = fileparts(gdf_path);
    fprintf('\n[%d/%d]  %s\n', fi, n_files, basename);

    res_dir = fullfile(gdf_dir, 'results', basename);
    if ~isfolder(res_dir), mkdir(res_dir); end

    % ── Load ─────────────────────────────────────────────────────────────────
    [signal, header, ~] = load_gdf(gdf_path);
    [params, ~]          = load_params_yaml(gdf_path);
    header_orig          = header;   % sample-level events for 897/898 lookup

    paradigm   = params.integrator.paradigm;
    fs         = double(params.acquisition.samplerate);
    framerate  = double(params.acquisition.framerate);
    chunk_size = round(fs / framerate);
    if abs(fs - header.SampleRate) > 1e-3
        fs = header.SampleRate;  chunk_size = round(fs / framerate);
    end
    bufsize_proc = double(params.RingBufferCfg.params.size);
    bufsize_art  = double(params.RingBufferCfgArtifact.params.size);
    eog_names    = to_strcell(params.CarCfg.params.EOG_ch_names);

    do_car_mi   = true; if isfield(params,'processing_fbcsp_mi'),   do_car_mi   = logical(params.processing_fbcsp_mi.do_car);   end
    do_car_cvsa = true; if isfield(params,'processing_fbcsp_cvsa'), do_car_cvsa = logical(params.processing_fbcsp_cvsa.do_car); end

    use_mi   = ismember(paradigm, {'mi','hybrid'});
    use_cvsa = ismember(paradigm, {'cvsa','hybrid'});
    csp_mi=[]; slda_mi=[]; csp_cvsa=[]; slda_cvsa=[];
    if use_mi,   csp_mi   = load_csp(params,'mi');   slda_mi   = load_slda(params,'mi');   end
    if use_cvsa, csp_cvsa = load_csp(params,'cvsa'); slda_cvsa = load_slda(params,'cvsa'); end

    % ── Pipeline (identical to main_simulate, but 4th output kept) ───────────
    feat_mi=[]; header_mi=[]; feat_pre_mi=[];
    feat_cv=[]; header_cv=[]; feat_pre_cv=[]; info_proc=[];
    proc_base = struct('samplerate',fs,'chunk_size',chunk_size,'bufsize',bufsize_proc,'filter_order',4);
    if use_mi
        cfg = proc_base; cfg.do_car = do_car_mi; cfg.eog_names = eog_names;
        [feat_mi, header_mi, info_proc, feat_pre_mi] = apply_processing(signal, header, csp_mi, cfg);
    end
    if use_cvsa
        cfg = proc_base; cfg.do_car = do_car_cvsa; cfg.eog_names = eog_names;
        [feat_cv, header_cv, info2, feat_pre_cv]  = apply_processing(signal, header, csp_cvsa, cfg);
        if isempty(info_proc), info_proc = info2; end
    end

    art_cfg = params.ArtifactCfg.params;
    art_cfg.EOG_ch_names = to_strcell(art_cfg.EOG_ch_names);
    [art_flags, ~] = detect_artifacts(signal, header, art_cfg, ...
        struct('samplerate',fs,'chunk_size',chunk_size,'bufsize_artifact',bufsize_art));

    p_mi_al=[]; if use_mi,   p_mi_al   = apply_slda(feat_mi, slda_mi,   csp_mi.bands);   end
    p_cv_al=[]; if use_cvsa, p_cv_al   = apply_slda(feat_cv, slda_cvsa, csp_cvsa.bands); end

    int_cfg = params.integrator;
    if ~isfield(int_cfg,'increment'),            int_cfg.increment = 1; end
    if ~isfield(int_cfg,'thresholds_rejection'), int_cfg.thresholds_rejection = []; end
    if ~isfield(int_cfg,'cvsa_influence'),       int_cfg.cvsa_influence = 2.5; end
    if ~isfield(int_cfg,'thresholds') || isempty(int_cfg.thresholds)
        int_cfg.thresholds = params.training_node.thresholds;
    end

    header_chunks = header_mi; if ~use_mi, header_chunks = header_cv; end
    header_chunks.framerate = framerate;

    trials = integrate_signal(p_mi_al, p_cv_al, art_flags, header_chunks, int_cfg, paradigm);
    n_trials = numel(trials);
    if n_trials == 0
        fprintf('  [warn] no trials found — skipping.\n'); continue;
    end

    % ── Metrics ──────────────────────────────────────────────────────────────
    thresholds = to_vec(int_cfg.thresholds);
    classes    = to_vec(int_cfg.classes);
    n_cls      = numel(classes);
    init_val   = to_vec(int_cfg.init_val);
    p_rest     = init_val(1);

    art_rate_per_trial = arrayfun(@(tr) mean(tr.artifact(tr.n_pre+1:end)), trials);
    art_rate           = mean(art_rate_per_trial);

    % Event-based accuracy from GDF outcome events
    POS_o = header_orig.EVENT.POS;  TYP_o = header_orig.EVENT.TYP;
    cf_sample_pos = POS_o(TYP_o == CF_CODE);
    trial_outcome_ev = zeros(1,n_trials);
    for t = 1:min(n_trials, numel(cf_sample_pos))
        after = POS_o > cf_sample_pos(t);
        idx   = find((TYP_o==HIT_CODE | TYP_o==MISS_CODE | TYP_o==TIMEOUT_CODE) & after, 1);
        if ~isempty(idx), trial_outcome_ev(t) = TYP_o(idx); end
    end
    n_hits_ev    = sum(trial_outcome_ev == HIT_CODE);
    n_misses_ev  = sum(trial_outcome_ev == MISS_CODE);
    n_timeout_ev = sum(trial_outcome_ev == TIMEOUT_CODE);
    % acc including timeouts: 897 / (897+898+899)
    trial_acc_total      = n_hits_ev / max(1, n_hits_ev + n_misses_ev + n_timeout_ev);
    % acc excluding timeouts: 897 / (897+898)
    trial_acc_no_timeout = n_hits_ev / max(1, n_hits_ev + n_misses_ev);

    % Per-class metrics
    tth_per_class       = nan(1,n_cls);
    t_miss_peak_cls     = nan(1,n_cls);
    hit_rate_cls        = nan(1,n_cls);
    confidence_cls      = nan(1,n_cls);
    peak_norm_cls       = nan(1,n_cls);

    for c = 1:n_cls
        mask_c = arrayfun(@(tr) ~isnan(tr.target_class) && tr.target_class==c, trials);
        if ~any(mask_c), continue; end
        trs_c     = trials(mask_c);
        pass_c    = [trs_c.pass];
        hit_rate_cls(c) = mean(pass_c);

        tth_v = []; t_miss_v = []; conf_v = []; peak_v = [];
        for t = 1:numel(trs_c)
            np = trs_c(t).n_pre;
            int_cf  = trs_c(t).integrated(np+1:end, c);
            raw_cf  = trs_c(t).raw(np+1:end, c);
            norm_cf = trs_c(t).normalized(np+1:end, c);
            n_cf_t  = trs_c(t).n_cf;
            t_axis  = (0:n_cf_t-1) / framerate;

            if trs_c(t).pass
                hf = find(int_cf >= thresholds(c), 1);
                if ~isempty(hf), tth_v(end+1) = t_axis(hf); end %#ok<AGROW>
            else
                [~, pf] = max(int_cf);
                if ~isempty(pf), t_miss_v(end+1) = t_axis(pf); end %#ok<AGROW>
            end
            vld = ~isnan(raw_cf);
            if any(vld)
                conf_v(end+1) = mean(raw_cf(vld)); %#ok<AGROW>
                peak_v(end+1) = max(norm_cf(vld)); %#ok<AGROW>
            end
        end
        if ~isempty(tth_v),    tth_per_class(c)   = mean(tth_v);    end
        if ~isempty(t_miss_v), t_miss_peak_cls(c)  = mean(t_miss_v); end
        if ~isempty(conf_v),   confidence_cls(c)   = mean(conf_v);   end
        if ~isempty(peak_v),   peak_norm_cls(c)    = mean(peak_v);   end
    end

    % Frame-level sample accuracy per class
    sample_acc_cls = nan(1, n_cls);
    for c = 1:n_cls
        cls_correct = [];
        for t = 1:n_trials
            if isnan(trials(t).target_class) || trials(t).target_class ~= c, continue; end
            np  = trials(t).n_pre;
            raw = trials(t).raw(np+1:end,:);
            vld = ~any(isnan(raw),2);
            if ~any(vld), continue; end
            ok = raw(vld,c) == max(raw(vld,:),[],2);
            cls_correct = [cls_correct; ok]; %#ok<AGROW>
        end
        if ~isempty(cls_correct), sample_acc_cls(c) = mean(cls_correct); end
    end

    cls_names = arrayfun(@(x) num2str(x), classes, 'UniformOutput', false);
    fprintf('  paradigm=%-8s  trials=%d  acc_total=%.0f%%  acc_no_timeout=%.0f%%  art=%.0f%%\n', ...
            paradigm, n_trials, 100*trial_acc_total, 100*trial_acc_no_timeout, 100*art_rate);
    for c = 1:n_cls
        fprintf('  class %s: hit=%.0f%%  sa=%.0f%%  tth=%.2fs  t_miss=%.2fs\n', ...
                cls_names{c}, 100*hit_rate_cls(c), 100*sample_acc_cls(c), ...
                tth_per_class(c), t_miss_peak_cls(c));
    end

    % ── Fig 1: metric summary (2×2) ───────────────────────────────────────────
    fig1 = figure('Name',sprintf('Metrics — %s',basename),'Color','w', ...
                  'NumberTitle','off','Position',[50 50 900 580]);

    subplot(2,2,1);
    bh = bar([trial_acc_total, trial_acc_no_timeout]*100, 'FaceColor','flat');
    bh.CData = [.22 .50 .78; .28 .68 .30];
    set(gca,'XTickLabel',{'897/(897+898+899)','897/(897+898)'}, ...
            'XTickLabelRotation',15,'YLim',[0,105]); grid on;
    ylabel('%'); title('Trial accuracy (event-based)','FontSize',9);
    yline(50,'--k','chance','FontSize',7,'LabelHorizontalAlignment','left');
    text(0.98,0.95,sprintf('hit=%d  miss=%d  timeout=%d', n_hits_ev,n_misses_ev,n_timeout_ev), ...
         'Units','normalized','HorizontalAlignment','right','FontSize',7);

    subplot(2,2,2);
    bar(1:n_cls, sample_acc_cls*100, 'FaceColor',[.45 .35 .75]);
    set(gca,'XTick',1:n_cls,'XTickLabel',cls_names,'YLim',[0,105]); grid on;
    ylabel('%'); title('sLDA frame accuracy per class','FontSize',9);
    yline(50,'--k');

    subplot(2,2,3);
    bar(1:n_cls, tth_per_class, 'FaceColor',[.3 .6 .9]);
    set(gca,'XTick',1:n_cls,'XTickLabel',cls_names); grid on;
    ylabel('s'); title('Mean time to hit','FontSize',9);

    subplot(2,2,4);
    bar(1:n_cls, t_miss_peak_cls, 'FaceColor',[.9 .5 .25]);
    set(gca,'XTick',1:n_cls,'XTickLabel',cls_names); grid on;
    ylabel('s'); title('Mean time to peak (miss trials)','FontSize',9);

    sgtitle(sprintf('%s  [%s]  %d trials', basename, upper(paradigm), n_trials), ...
            'FontSize',10,'Interpreter','none');
    saveas(fig1, fullfile(res_dir,'01_metrics.svg'),'svg');

    % ── ERD/ERS data (primary stream, ERD% vs 1s pre-CF baseline) ───────────
    if use_mi
        feat_pre = feat_pre_mi;  csp_topo = csp_mi;  sel_ch = csp_mi.selected_channels;
    else
        feat_pre = feat_pre_cv;  csp_topo = csp_cvsa; sel_ch = csp_cvsa.selected_channels;
    end
    n_sel_ch = numel(sel_ch);
    n_bands  = csp_topo.n_bands;
    bands    = csp_topo.bands;
    feat3d   = permute(reshape(feat_pre.', n_sel_ch, n_bands, size(feat_pre,1)), [3,1,2]);

    n_cf_min = min(arrayfun(@(tr) tr.n_cf, trials));
    n_bl     = round(framerate);   % 1 s pre-CF baseline
    erd_sum  = zeros(n_cls, n_cf_min, n_sel_ch, n_bands);
    erd_cnt  = zeros(n_cls, 1);

    for t = 1:n_trials
        c = trials(t).target_class;
        if isnan(c), continue; end
        sc    = trials(t).start_chunk;
        n_use = min(trials(t).n_cf, n_cf_min);
        bl_r  = max(1,sc-n_bl) : sc-1;
        cf_r  = sc : sc + n_use - 1;
        bl_r  = bl_r(bl_r>=1 & bl_r<=size(feat3d,1));
        cf_r  = cf_r(cf_r>=1 & cf_r<=size(feat3d,1));
        if isempty(bl_r) || numel(cf_r)<n_use, continue; end
        pbl = mean(feat3d(bl_r,:,:), 1);
        pcf = feat3d(cf_r,:,:);
        if any(isnan(pbl(:))) || any(isnan(pcf(:))), continue; end
        erd_t = (pcf - pbl) ./ max(pbl,eps) * 100;
        erd_sum(c,1:n_use,:,:) = erd_sum(c,1:n_use,:,:) + ...
            reshape(erd_t, 1, n_use, n_sel_ch, n_bands);
        erd_cnt(c) = erd_cnt(c) + 1;
    end
    erd_mean = erd_sum;
    for c = 1:n_cls
        if erd_cnt(c)>0, erd_mean(c,:,:,:) = erd_sum(c,:,:,:) / erd_cnt(c); end
    end
    t_ax_cf = (0:n_cf_min-1) / framerate;

    % Static ERD: mean over CF window, sum over bands → [n_cls × n_sel_ch]
    erd_static = nan(n_cls, n_sel_ch);
    for c = 1:n_cls
        if erd_cnt(c)>0
            tmp = squeeze(mean(erd_mean(c,:,:,:), 2));   % [n_sel_ch × n_bands]
            tmp = reshape(tmp, n_sel_ch, n_bands);
            erd_static(c,:) = sum(tmp, 2);
        end
    end

    % ── Figs 2–(n_cls+1): ERD/ERS heatmap, one figure per class ─────────────
    cmap_erd = redblue_cmap();
    for c = 1:n_cls
        all_d    = squeeze(erd_mean(c,:,:,:));
        clim_erd = max(1, max(abs(all_d(:))));
        fig_erd  = figure('Name',sprintf('ERD/ERS %s — %s',cls_names{c},basename), ...
                          'Color','w','NumberTitle','off', ...
                          'Position',[60+c*30 30 700 130*n_bands+60]);
        for b = 1:n_bands
            ax   = subplot(n_bands,1,b);
            data = squeeze(erd_mean(c,:,:,b));   % [n_cf_min × n_sel_ch]
            if ~any(isnan(data(:)))
                imagesc(ax, t_ax_cf, 1:n_sel_ch, data.');
                colormap(ax, cmap_erd);  caxis(ax, [-clim_erd, clim_erd]);
                colorbar(ax,'FontSize',6);
                set(ax,'YTick',1:n_sel_ch,'YTickLabel',sel_ch,'FontSize',6);
                if b==n_bands, xlabel(ax,'Time from CF onset (s)','FontSize',7); end
                title(ax,sprintf('[%.0f–%.0f Hz]',bands(b,1),bands(b,2)),'FontSize',8);
                xline(ax,0,'k--');
            else
                axis(ax,'off');
                text(ax,0.5,0.5,'no data','HorizontalAlignment','center','Units','normalized');
            end
        end
        sgtitle(fig_erd, sprintf('ERD/ERS [%%Δ vs 1s baseline]  class %s — %s [%s]', ...
                cls_names{c},basename,upper(paradigm)),'FontSize',10,'Interpreter','none');
        saveas(fig_erd, fullfile(res_dir,sprintf('%02d_erd_class%s.svg',c+1,cls_names{c})),'svg');
    end

    % ── Fig (n_cls+2): spatial ERD/ERS — static + time snapshots ─────────────
    n_tpts   = min(4, n_cf_min);
    t_idx    = round(linspace(1, n_cf_min, n_tpts));
    n_rows   = 1 + n_tpts;
    all_sv   = erd_static(:);
    all_tv   = zeros(n_cls*n_cf_min*n_sel_ch,1);
    for c=1:n_cls
        if erd_cnt(c)>0
            tmp = sum(reshape(erd_mean(c,:,:,:), n_cf_min, n_sel_ch, n_bands), 3);
            all_tv = [all_tv; tmp(:)]; %#ok<AGROW>
        end
    end
    clim_s = sym_clim([all_sv; all_tv]);

    fig_topo = figure('Name',sprintf('Spatial ERD/ERS — %s',basename),'Color','w', ...
                      'NumberTitle','off','Position',[80 50 310*n_cls 260*n_rows+60]);
    for c = 1:n_cls
        ax = subplot(n_rows, n_cls, c);
        if erd_cnt(c)>0
            topo_scatter(sel_ch(:), erd_static(c,:).', clim_s, ax, ...
                sprintf('class %s — mean CF (Σbands)', cls_names{c}));
        else
            axis(ax,'off'); text(ax,0.5,0.5,'no data','HorizontalAlignment','center','Units','normalized');
        end
        for ti = 1:n_tpts
            ax = subplot(n_rows, n_cls, n_cls + (ti-1)*n_cls + c);
            if erd_cnt(c)>0
                snap = sum(squeeze(erd_mean(c,t_idx(ti),:,:)), 2);
                topo_scatter(sel_ch(:), snap, clim_s, ax, ...
                    sprintf('t = %.1f s', t_ax_cf(t_idx(ti))));
            else
                axis(ax,'off');
            end
        end
    end
    sgtitle(fig_topo, sprintf('Spatial ERD/ERS — %s [%s]  (blue=ERD↓  red=ERS↑)', ...
            basename,upper(paradigm)),'FontSize',10,'Interpreter','none');
    saveas(fig_topo, fullfile(res_dir,sprintf('%02d_spatial_erd.svg',n_cls+2)),'svg');

    % ── Save mat ─────────────────────────────────────────────────────────────
    results = struct( ...
        'subject',              subject_from_basename(basename), ...
        'basename',             basename, ...
        'gdf_path',             gdf_path, ...
        'paradigm',             paradigm, ...
        'classes',              classes, ...
        'n_trials',             n_trials, ...
        'trial_acc_total',      trial_acc_total, ...
        'trial_acc_no_timeout', trial_acc_no_timeout, ...
        'n_hits_events',        n_hits_ev, ...
        'n_misses_events',      n_misses_ev, ...
        'n_timeout_events',     n_timeout_ev, ...
        'tth_per_class',        tth_per_class, ...
        't_miss_peak_cls',      t_miss_peak_cls, ...
        'sample_acc_cls',       sample_acc_cls, ...
        'art_rate',             art_rate, ...
        'hit_rate_cls',         hit_rate_cls, ...
        'confidence_cls',       confidence_cls, ...
        'peak_norm_cls',        peak_norm_cls, ...
        'erd_mean',             erd_mean, ...     % [n_cls × n_cf_min × n_sel_ch × n_bands] ERD%
        'erd_static',           erd_static, ...   % [n_cls × n_sel_ch] mean-CF Σbands
        't_ax_cf',              t_ax_cf, ...
        'selected_channels',    {sel_ch}, ...
        'bands',                bands, ...
        'framerate',            framerate, ...
        'chunk_size',           chunk_size ...
    );

    mat_path = fullfile(gdf_dir, sprintf('eval_single_%s_%s.mat', paradigm, basename));
    save(mat_path, 'results', 'trials', 'int_cfg');
    fprintf('  saved  %s\n         %s\n', res_dir, mat_path);
end

% ── Local helpers ─────────────────────────────────────────────────────────────
function s = subject_from_basename(bn)
    parts = strsplit(bn, '.');
    if numel(parts) >= 1, s = parts{1}; else, s = bn; end
end

function cm = redblue_cmap(n)
    if nargin < 1, n = 256; end
    h = floor(n/2);
    r = [linspace(0,1,h); linspace(1,1,h)];
    g = [linspace(0,1,h); linspace(1,0,h)];
    b = [linspace(1,1,h); linspace(1,0,h)];
    cm = [r(:), g(:), b(:)];
end

function cl = sym_clim(vals)
    v = max(abs(vals(isfinite(vals))));
    if isempty(v) || v == 0, v = 1; end
    cl = [-v, v];
end
