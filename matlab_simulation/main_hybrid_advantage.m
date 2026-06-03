%% MAIN_HYBRID_ADVANTAGE  Offline evidence for the hybrid MI+CVSA BCI benefit.
%
%   Loads one or more GDF files (paradigm = hybrid, mi, or cvsa).
%
%   For each HYBRID GDF the SAME EEG is replayed in three conditions:
%     • Hybrid  — LOP fusion with cosine-annealed CVSA prior
%     • MI-only — leaky integrator driven by P_MI alone
%     • CVSA-only — leaky integrator driven by P_CVSA alone
%   Trial matching is exact because all three conditions share the same
%   event 781 positions and the same art_flags.
%
%   For MI-only or CVSA-only GDFs: baseline metrics only (no comparison).
%
%   Figures produced (per hybrid file):
%     Fig 1  — Temporal dynamics: P_CVSA vs P_MI over CF time
%              (Is CVSA informative early while MI is still weak?)
%     Fig 2  — Buffer trajectories: target-class buffer mean ± std
%              for HIT trials and MISS trials, all three conditions
%     Fig 3  — Hit rate + time-to-hit: barplot + CDF comparison
%     Fig 4  — Saved trials: MI fails → hybrid succeeds
%              (Direct mechanistic evidence for the advantage)
%     Fig 5  — Classifier correlation ρ(P_MI, P_CVSA)
%              (Low ρ → LOP adds real information; high ρ → overcounting)
%   Plus Fig 6 (multi-file overview) when > 1 file is loaded.
%
%   Saves figures as SVG to <gdf_dir>/hybrid_advantage/
%   and advantage_<paradigm>_<basename>.mat per file.

clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir,'io'), fullfile(this_dir,'processing'), ...
        fullfile(this_dir,'artifacts'), fullfile(this_dir,'classifier'), ...
        fullfile(this_dir,'integrator'), fullfile(this_dir,'plotting'), ...
        fullfile(this_dir,'utils'));

N_EARLY_S = 0.5;   % seconds considered "trial onset" for saved-trial analysis

% ── File picker ───────────────────────────────────────────────────────────────
default_dir = '/home/paolo/bci_vr_ws/recordings';
if ~isfolder(default_dir), default_dir = fileparts(this_dir); end
[gdf_names, gdf_dir] = uigetfile({'*.gdf','GDF files (*.gdf)'}, ...
    'Select GDF file(s)', default_dir, 'MultiSelect','on');
if isequal(gdf_names,0), error('main_hybrid_advantage:cancel','No file selected.'); end
if ischar(gdf_names), gdf_names = {gdf_names}; end
n_files = numel(gdf_names);

out_dir = fullfile(gdf_dir, 'hybrid_advantage');
if ~isfolder(out_dir), mkdir(out_dir); end

all_res_cell = {};   % cell array — avoids dissimilar-struct issues
par_count    = struct('mi',0,'cvsa',0,'hybrid',0);  % counter for short labels

%% ═════════════════════════════════════════════════════════════════════════════
%  PER-FILE PROCESSING LOOP
%  ═════════════════════════════════════════════════════════════════════════════
for fi = 1:n_files
    gdf_path = fullfile(gdf_dir, gdf_names{fi});
    [~, basename] = fileparts(gdf_path);
    fprintf('\n[%d/%d]  %s\n', fi, n_files, basename);

    % ── Load GDF + YAML ───────────────────────────────────────────────────────
    [signal, header, ~] = load_gdf(gdf_path);
    [params, ~]         = load_params_yaml(gdf_path);

    paradigm   = params.integrator.paradigm;
    is_hybrid  = strcmp(paradigm,'hybrid');
    fs         = double(params.acquisition.samplerate);
    framerate  = double(params.acquisition.framerate);
    chunk_size = round(fs / framerate);
    if abs(fs - header.SampleRate) > 1e-3
        fs = header.SampleRate; chunk_size = round(fs / framerate);
    end

    bufsize_proc = double(params.RingBufferCfg.params.size);
    bufsize_art  = double(params.RingBufferCfgArtifact.params.size);
    eog_names    = to_strcell(params.CarCfg.params.EOG_ch_names);

    do_car_mi   = true;  if isfield(params,'processing_fbcsp_mi'),   do_car_mi   = logical(params.processing_fbcsp_mi.do_car);   end
    do_car_cvsa = true;  if isfield(params,'processing_fbcsp_cvsa'), do_car_cvsa = logical(params.processing_fbcsp_cvsa.do_car); end

    use_mi   = ismember(paradigm,{'mi','hybrid'});
    use_cvsa = ismember(paradigm,{'cvsa','hybrid'});

    csp_mi=[]; slda_mi=[]; csp_cvsa=[]; slda_cvsa=[];
    if use_mi,   csp_mi   = load_csp(params,'mi');   slda_mi   = load_slda(params,'mi');   end
    if use_cvsa, csp_cvsa = load_csp(params,'cvsa'); slda_cvsa = load_slda(params,'cvsa'); end

    proc_base = struct('samplerate',fs,'chunk_size',chunk_size, ...
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
        struct('samplerate',fs,'chunk_size',chunk_size,'bufsize_artifact',bufsize_art));

    p_mi_al=[]; if use_mi,   p_mi_al = apply_slda(feat_mi, slda_mi,   csp_mi.bands);   end
    p_cv_al=[]; if use_cvsa, p_cv_al = apply_slda(feat_cv, slda_cvsa, csp_cvsa.bands); end

    int_cfg = params.integrator;
    if ~isfield(int_cfg,'increment'),            int_cfg.increment = 1; end
    if ~isfield(int_cfg,'thresholds_rejection'), int_cfg.thresholds_rejection = []; end
    if ~isfield(int_cfg,'cvsa_influence'),       int_cfg.cvsa_influence = 3.0; end
    if ~isfield(int_cfg,'thresholds') || isempty(int_cfg.thresholds)
        int_cfg.thresholds = params.training_node.thresholds;
    end

    thresholds = to_vec(int_cfg.thresholds);
    classes    = to_vec(int_cfg.classes);
    n_cls      = numel(classes);
    init_val   = to_vec(int_cfg.init_val);
    p_rest     = init_val(1);
    cls_names  = arrayfun(@(x) num2str(x), classes, 'UniformOutput',false);

    header_chunks = header_mi; if ~use_mi, header_chunks = header_cv; end
    header_chunks.framerate = framerate;

    % ── Integrate: native paradigm + available single-modality simulations ───────
    trials_hyb  = integrate_signal(p_mi_al, p_cv_al, art_flags, header_chunks, int_cfg, paradigm);
    trials_mi   = struct([]);
    trials_cvsa = struct([]);
    if ~isempty(p_mi_al)
        trials_mi   = integrate_signal(p_mi_al, [],      art_flags, header_chunks, int_cfg, 'mi');
    end
    if ~isempty(p_cv_al)
        trials_cvsa = integrate_signal([],      p_cv_al, art_flags, header_chunks, int_cfg, 'cvsa');
    end

    n_trials = numel(trials_hyb);
    if n_trials == 0, fprintf('  [warn] no trials\n'); continue; end

    % ── Short label (mi_1, mi_2, cvsa_1, hybrid_1 …) ─────────────────────────
    par_count.(paradigm) = par_count.(paradigm) + 1;
    file_label = sprintf('%s_%d', paradigm, par_count.(paradigm));

    % ── Event-based accuracy from actual GDF outcomes (897/898/899) ───────────
    HIT_CODE     = 897;  MISS_CODE = 898;  TIMEOUT_CODE = 899;  CF_CODE = 781;
    POS_ev = header.EVENT.POS;  TYP_ev = header.EVENT.TYP;
    cf_samp = POS_ev(TYP_ev == CF_CODE);
    trial_outcome_ev = zeros(1, n_trials);
    for t = 1:min(n_trials, numel(cf_samp))
        after = POS_ev > cf_samp(t);
        idx   = find((TYP_ev==HIT_CODE|TYP_ev==MISS_CODE|TYP_ev==TIMEOUT_CODE)&after, 1);
        if ~isempty(idx), trial_outcome_ev(t) = TYP_ev(idx); end
    end
    n_hits_ev    = sum(trial_outcome_ev == HIT_CODE);
    n_misses_ev  = sum(trial_outcome_ev == MISS_CODE);
    n_timeout_ev = sum(trial_outcome_ev == TIMEOUT_CODE);
    n_total_ev   = max(1, n_hits_ev + n_misses_ev + n_timeout_ev);
    hit_rate_ev       = n_hits_ev / n_total_ev;
    hit_rate_no_to    = n_hits_ev / max(1, n_hits_ev + n_misses_ev);
    timeout_rate      = n_timeout_ev / n_total_ev;

    % ── Real TTH from GDF: time from CF onset (781) to first 897 event ───────
    % Only computed for trials that actually produced event 897 online.
    real_tth_cls      = nan(1, max(2, numel(to_vec(params.integrator.classes))));
    real_tth_cls_vals = cell(size(real_tth_cls));
    for c_tmp = 1:numel(real_tth_cls)
        tth_v_real = [];
        for t = 1:min(n_trials, numel(cf_samp))
            if trial_outcome_ev(t) ~= HIT_CODE, continue; end
            tr_cls = trials_hyb(t).target_class;
            if isnan(tr_cls) || tr_cls ~= c_tmp, continue; end
            % First 897 after this CF onset
            mask897 = POS_ev > cf_samp(t) & TYP_ev == HIT_CODE;
            if any(mask897)
                tth_v_real(end+1) = (min(POS_ev(mask897)) - cf_samp(t)) / fs; %#ok<AGROW>
            end
        end
        real_tth_cls_vals{c_tmp} = tth_v_real;
        if ~isempty(tth_v_real), real_tth_cls(c_tmp) = mean(tth_v_real); end
    end
    real_tth_all = [real_tth_cls_vals{:}];  % all hits, pooled across classes

    % Condition labels and handles (only show conditions available)
    cond_trials  = {};
    cond_labels  = {};
    cond_colors  = [];
    if is_hybrid
        cond_trials = {trials_hyb, trials_mi, trials_cvsa};
        cond_labels = {'Hybrid','MI-only','CVSA-only'};
        cond_colors = [0.18 0.45 0.75;  0.85 0.30 0.10;  0.10 0.60 0.30];
    elseif strcmp(paradigm,'mi')
        cond_trials = {trials_mi};
        cond_labels = {'MI-only'};
        cond_colors = [0.85 0.30 0.10];
    else
        cond_trials = {trials_cvsa};
        cond_labels = {'CVSA-only'};
        cond_colors = [0.10 0.60 0.30];
    end
    n_cond = numel(cond_trials);

    fprintf('  paradigm=%-7s  trials=%d', paradigm, n_trials);
    for ic = 1:n_cond
        fprintf('  %s=%.0f%%', cond_labels{ic}, 100*mean([cond_trials{ic}.pass]));
    end
    fprintf('\n');

    n_cf_min = min(arrayfun(@(tr) tr.n_cf, trials_hyb));  % trials_hyb = native paradigm output
    t_ax     = (0:n_cf_min-1) / framerate;
    n_early  = max(1, round(N_EARLY_S * framerate));

    % ── Per-condition metrics ─────────────────────────────────────────────────
    hit_rate  = nan(1,n_cond);
    tth_cls   = cell(n_cls,n_cond);   % TTH values per class per condition

    for ic = 1:n_cond
        tr = cond_trials{ic};
        hit_rate(ic) = mean([tr.pass]);
        for c = 1:n_cls
            tth_v = [];
            for t = 1:n_trials
                if isnan(tr(t).target_class)||tr(t).target_class~=c||~tr(t).pass, continue; end
                np    = tr(t).n_pre;
                iv    = tr(t).integrated(np+1:end, c);
                t_ax_t = (0:tr(t).n_cf-1)/framerate;
                hf = find(iv >= thresholds(c), 1);
                if ~isempty(hf), tth_v(end+1) = t_ax_t(hf); end %#ok<AGROW>
            end
            tth_cls{c,ic} = tth_v;
        end
    end

    % ── Buffer trajectory: target-class buffer averaged across HIT/MISS ───────
    % For each trial t: buffer(target_class, t) over CF time
    buf_hit  = cell(n_cond,1);   % mean buffer for HIT trials [n_cf_min x 1]
    buf_miss = cell(n_cond,1);
    buf_hit_sd  = cell(n_cond,1);
    buf_miss_sd = cell(n_cond,1);

    for ic = 1:n_cond
        tr = cond_trials{ic};
        acc_h = zeros(n_cf_min,1); cnt_h = 0;
        acc_m = zeros(n_cf_min,1); cnt_m = 0;
        sq_h  = zeros(n_cf_min,1); sq_m  = zeros(n_cf_min,1);

        for t = 1:n_trials
            c = tr(t).target_class; if isnan(c), continue; end
            np = tr(t).n_pre;
            iv = tr(t).integrated(np+1:end, c);  % target-class buffer
            nu = min(numel(iv), n_cf_min);
            if any(isnan(iv(1:nu))), continue; end
            if tr(t).pass
                acc_h(1:nu) = acc_h(1:nu) + iv(1:nu);
                sq_h(1:nu)  = sq_h(1:nu)  + iv(1:nu).^2;
                cnt_h = cnt_h + 1;
            else
                acc_m(1:nu) = acc_m(1:nu) + iv(1:nu);
                sq_m(1:nu)  = sq_m(1:nu)  + iv(1:nu).^2;
                cnt_m = cnt_m + 1;
            end
        end
        if cnt_h > 0
            buf_hit{ic} = acc_h/cnt_h;
            buf_hit_sd{ic} = sqrt(max(0, sq_h/cnt_h - buf_hit{ic}.^2));
        end
        if cnt_m > 0
            buf_miss{ic} = acc_m/cnt_m;
            buf_miss_sd{ic} = sqrt(max(0, sq_m/cnt_m - buf_miss{ic}.^2));
        end
    end

    % ── Temporal dynamics: mean classifier output over CF time per class ─────────
    % For hybrid: P_MI and P_CVSA separately (from .p_mi / .p_cvsa fields).
    % For single-modality: the native classifier output from .raw.
    p_mi_mean  = nan(n_cls,n_cf_min);
    p_cv_mean  = nan(n_cls,n_cf_min);
    p_mi_std   = nan(n_cls,n_cf_min);
    p_cv_std   = nan(n_cls,n_cf_min);
    p_nat_mean = nan(n_cls,n_cf_min);  % native paradigm output (always available)
    p_nat_std  = nan(n_cls,n_cf_min);

    for c = 1:n_cls
        acc_nat=zeros(n_cf_min,1); sq_nat=zeros(n_cf_min,1);
        acc_mi =zeros(n_cf_min,1); sq_mi =zeros(n_cf_min,1);
        acc_cv =zeros(n_cf_min,1); sq_cv =zeros(n_cf_min,1);
        cnt_nat=0; cnt_mi=0; cnt_cv=0;
        for t = 1:n_trials
            if isnan(trials_hyb(t).target_class)||trials_hyb(t).target_class~=c, continue; end
            np = trials_hyb(t).n_pre;
            nu = min(trials_hyb(t).n_cf, n_cf_min);

            % Native output (works for any paradigm)
            rv = trials_hyb(t).raw(np+1:np+nu, c);
            vn = ~isnan(rv);
            if any(vn)
                rv(~vn)=0;
                acc_nat(1:nu)=acc_nat(1:nu)+rv; sq_nat(1:nu)=sq_nat(1:nu)+rv.^2;
                cnt_nat=cnt_nat+1;
            end

            % MI and CVSA individual outputs (hybrid only — else stays NaN)
            if is_hybrid
                pm  = trials_hyb(t).p_mi(np+1:np+nu, c);
                pc  = trials_hyb(t).p_cvsa(np+1:np+nu, c);
                vm  = ~isnan(pm); vc = ~isnan(pc);
                if any(vm)
                    pm(~vm)=0;
                    acc_mi(1:nu)=acc_mi(1:nu)+pm; sq_mi(1:nu)=sq_mi(1:nu)+pm.^2;
                    cnt_mi=cnt_mi+1;
                end
                if any(vc)
                    pc(~vc)=0;
                    acc_cv(1:nu)=acc_cv(1:nu)+pc; sq_cv(1:nu)=sq_cv(1:nu)+pc.^2;
                    cnt_cv=cnt_cv+1;
                end
            end
        end
        if cnt_nat>0, p_nat_mean(c,:)=acc_nat'/cnt_nat;
            p_nat_std(c,:)=sqrt(max(0,sq_nat'/cnt_nat-p_nat_mean(c,:).^2)); end
        if cnt_mi>0, p_mi_mean(c,:)=acc_mi'/cnt_mi;
            p_mi_std(c,:)=sqrt(max(0,sq_mi'/cnt_mi-p_mi_mean(c,:).^2)); end
        if cnt_cv>0, p_cv_mean(c,:)=acc_cv'/cnt_cv;
            p_cv_std(c,:)=sqrt(max(0,sq_cv'/cnt_cv-p_cv_mean(c,:).^2)); end
    end

    % ── Saved trials, early values, correlation (hybrid only) ─────────────────
    saved_mask = false(1,n_trials);
    early_cv_saved = []; early_mi_saved = [];
    early_cv_other = []; early_mi_other = [];
    rho_trial  = nan(n_trials,1);
    % Per-trial mean output for HIT vs MISS (all paradigms)
    pt_mean_hit  = nan(n_trials,1);
    pt_mean_miss = nan(n_trials,1);
    for t = 1:n_trials
        c = trials_hyb(t).target_class; if isnan(c), continue; end
        np = trials_hyb(t).n_pre;
        rv = trials_hyb(t).raw(np+1:end, c);
        rv = rv(~isnan(rv));
        if trials_hyb(t).pass && ~isempty(rv), pt_mean_hit(t)  = mean(rv); end
        if ~trials_hyb(t).pass && ~isempty(rv), pt_mean_miss(t) = mean(rv); end
    end

    if is_hybrid
        saved_mask = ~[trials_mi.pass] & [trials_hyb.pass];
        for t = 1:n_trials
            c = trials_hyb(t).target_class; if isnan(c), continue; end
        end

        % Saved trials and early classifier values
        saved_mask = ~[trials_mi.pass] & [trials_hyb.pass];
        for t = 1:n_trials
            c = trials_hyb(t).target_class; if isnan(c), continue; end
            np = trials_hyb(t).n_pre;
            ne = min(n_early, trials_hyb(t).n_cf);
            pm_e = trials_hyb(t).p_mi(np+1:np+ne, c);
            pc_e = trials_hyb(t).p_cvsa(np+1:np+ne, c);
            pm_e = pm_e(~isnan(pm_e));  pc_e = pc_e(~isnan(pc_e));
            if saved_mask(t)
                if ~isempty(pm_e), early_mi_saved(end+1)  = mean(pm_e); end %#ok<AGROW>
                if ~isempty(pc_e), early_cv_saved(end+1) = mean(pc_e); end %#ok<AGROW>
            else
                if ~isempty(pm_e), early_mi_other(end+1)  = mean(pm_e); end %#ok<AGROW>
                if ~isempty(pc_e), early_cv_other(end+1) = mean(pc_e); end %#ok<AGROW>
            end

            % Correlation per trial: P_MI(c1) vs P_CVSA(c1) across CF frames
            pm_cf = trials_hyb(t).p_mi(np+1:end, 1);
            pc_cf = trials_hyb(t).p_cvsa(np+1:end, 1);
            vld   = ~isnan(pm_cf) & ~isnan(pc_cf);
            if sum(vld) > 4
                rho_trial(t) = corr(pm_cf(vld), pc_cf(vld));
            end
        end
        fprintf('  Saved trials (MI-fail→hybrid-HIT): %d / %d  mean ρ(P_MI,P_CVSA)=%.3f\n', ...
                sum(saved_mask), n_trials, mean(rho_trial,'omitnan'));
    end

    %% ────────────────────────────────────────────────────────────────────────
    %  FIGURES
    %  ────────────────────────────────────────────────────────────────────────

    % ── Fig 1: Temporal dynamics — all paradigms ─────────────────────────────
    % Hybrid: shows P_MI (red) and P_CVSA (blue) separately → reveals complementarity.
    % Single-modality: shows the native classifier output (P_MI or P_CVSA).
    fig1 = figure('Name',sprintf('Temporal Dynamics — %s',basename), ...
                  'Color','w','NumberTitle','off', ...
                  'Position',[40 40 max(600,280*n_cls*2) 380]);
    for c = 1:n_cls
        ax = subplot(1,n_cls,c);
        hold(ax,'on');

        % Native output (available for all paradigms — fused for hybrid, raw for single)
        t_n = t_ax(~isnan(p_nat_mean(c,:)));
        nat_m = p_nat_mean(c,~isnan(p_nat_mean(c,:)));
        nat_s = p_nat_std(c,~isnan(p_nat_mean(c,:)));
        if ~isempty(nat_m)
            lbl_nat = upper(paradigm);  % 'HYBRID', 'MI', 'CVSA'
            col_nat = cond_colors(1,:);
            fill(ax,[t_n,fliplr(t_n)],[nat_m+nat_s,fliplr(nat_m-nat_s)], ...
                 col_nat,'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
            plot(ax,t_n,nat_m,'Color',col_nat,'LineWidth',2.5,'DisplayName',lbl_nat);
        end

        if is_hybrid
            % Additionally show individual classifiers
            t_cv = t_ax(~isnan(p_cv_mean(c,:))); cv_m=p_cv_mean(c,~isnan(p_cv_mean(c,:))); cv_s=p_cv_std(c,~isnan(p_cv_mean(c,:)));
            t_mi = t_ax(~isnan(p_mi_mean(c,:))); mi_m=p_mi_mean(c,~isnan(p_mi_mean(c,:))); mi_s=p_mi_std(c,~isnan(p_mi_mean(c,:)));
            if ~isempty(cv_m)
                fill(ax,[t_cv,fliplr(t_cv)],[cv_m+cv_s,fliplr(cv_m-cv_s)], ...
                     [0.4 0.6 0.95],'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
                plot(ax,t_cv,cv_m,'Color',[0.15 0.35 0.85],'LineWidth',2.0,'LineStyle','--','DisplayName','P_{CVSA}');
            end
            if ~isempty(mi_m)
                fill(ax,[t_mi,fliplr(t_mi)],[mi_m+mi_s,fliplr(mi_m-mi_s)], ...
                     [0.95 0.55 0.35],'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
                plot(ax,t_mi,mi_m,'Color',[0.80 0.20 0.05],'LineWidth',2.0,'LineStyle','--','DisplayName','P_{MI}');
            end
            % Mark fusion schedule
            xline(ax,int_cfg.cvsa_influence/2,'--','Color',[0.5 0.5 0.5], ...
                  'Label','α=0.5','LabelVerticalAlignment','top','HandleVisibility','off');
            xline(ax,int_cfg.cvsa_influence,'--','Color',[0.3 0.3 0.3], ...
                  'Label','α=0','LabelVerticalAlignment','top','HandleVisibility','off');
        end

        yline(ax,0.5,'k:','LineWidth',1,'HandleVisibility','off');
        legend(ax,'FontSize',9,'Location','best');
        set(ax,'XLim',[0,t_ax(end)],'YLim',[0.35,1.0]); grid(ax,'on');
        xlabel(ax,'Time from CF onset (s)','FontSize',9);
        ylabel(ax,'P(target class)','FontSize',9);
        title(ax,sprintf('Class %s  —  all trials (mean ± std)',cls_names{c}),'FontSize',10);
    end
    if is_hybrid
        sgtitle(fig1,sprintf('Temporal dynamics: fused / P_{MI} / P_{CVSA} — %s [HYBRID]',basename), ...
                'FontSize',11,'Interpreter','none');
    else
        sgtitle(fig1,sprintf('Classifier output over CF time — %s [%s]',basename,upper(paradigm)), ...
                'FontSize',11,'Interpreter','none');
    end
    saveas(fig1, fullfile(out_dir,sprintf('%s_01_temporal_dynamics.svg',basename)),'svg');

    % ── Fig 2: Buffer trajectories ────────────────────────────────────────────
    fig2 = figure('Name',sprintf('Buffer Trajectories — %s',basename), ...
                  'Color','w','NumberTitle','off','Position',[50 50 900 400]);
    ax_h = subplot(1,2,1); hold(ax_h,'on');
    ax_m = subplot(1,2,2); hold(ax_m,'on');

    for ic = 1:n_cond
        col = cond_colors(ic,:);
        if ~isempty(buf_hit{ic})
            bh = buf_hit{ic}; sd = buf_hit_sd{ic};
            fill(ax_h,[t_ax,fliplr(t_ax)],[bh+sd;flipud(bh-sd)]', ...
                 col,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
            plot(ax_h,t_ax,bh,'Color',col,'LineWidth',2.5,'DisplayName',cond_labels{ic});
        end
        if ~isempty(buf_miss{ic})
            bm = buf_miss{ic}; sd = buf_miss_sd{ic};
            fill(ax_m,[t_ax,fliplr(t_ax)],[bm+sd;flipud(bm-sd)]', ...
                 col,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
            plot(ax_m,t_ax,bm,'Color',col,'LineWidth',2.5,'DisplayName',cond_labels{ic});
        end
    end
    thr_ref = max(thresholds);
    for ax = [ax_h, ax_m]
        yline(ax,thr_ref,'k--','LineWidth',1.2,'Label','threshold','FontSize',7);
        yline(ax,p_rest,'k:','LineWidth',0.8);
        set(ax,'XLim',[0,t_ax(end)],'YLim',[max(0,p_rest-0.05),1.0]);
        xlabel(ax,'Time from CF onset (s)','FontSize',9);
        ylabel(ax,'Buffer — target class','FontSize',9);
        legend(ax,'FontSize',8,'Location','best'); grid(ax,'on');
    end
    title(ax_h,'HIT trials  (target-class buffer mean ± std)','FontSize',10);
    title(ax_m,'MISS trials  (target-class buffer mean ± std)','FontSize',10);
    sgtitle(fig2,sprintf('Buffer trajectories — %s  [%s]',basename,upper(paradigm)), ...
            'FontSize',12,'Interpreter','none');
    saveas(fig2, fullfile(out_dir,sprintf('%s_02_buffer_trajectories.svg',basename)),'svg');

    % ── Fig 3: Hit rate (real + sim) + TTH + CDF ─────────────────────────────
    fig3 = figure('Name',sprintf('Performance — %s',basename), ...
                  'Color','w','NumberTitle','off','Position',[60 60 1100 420]);

    ax1 = subplot(1,3,1); hold(ax1,'on');

    % Real event-based bars (solid, thick outline)
    real_labels = {'897/(897+898+899)','897/(897+898)'};
    real_vals   = [hit_rate_ev, hit_rate_no_to]*100;
    real_cols   = [0.18 0.45 0.75; 0.20 0.65 0.25];
    for k = 1:2
        bar(ax1, k, real_vals(k), 0.45, 'FaceColor',real_cols(k,:), ...
            'EdgeColor',[0 0 0],'LineWidth',1.5,'DisplayName',real_labels{k});
        text(ax1, k, real_vals(k)+2, sprintf('%.0f%%\n(H=%d M=%d T=%d)', ...
             real_vals(k), n_hits_ev, n_misses_ev, n_timeout_ev), ...
             'HorizontalAlignment','center','FontSize',7.5,'FontWeight','bold');
    end

    % ── Hybrid: add MI-sim and CVSA-sim bars (NOT Hybrid-sim) ────────────────
    % Non-hybrid files: only the two real bars above.
    if is_hybrid
        % ic=2 → MI-sim, ic=3 → CVSA-sim  (ic=1 = hybrid is the native → use real values)
        sim_ic   = 2:n_cond;   % indices of non-native simulations
        sim_cols_h = cond_colors(sim_ic,:);
        sim_lbl_h  = cond_labels(sim_ic);   % {'MI-only','CVSA-only'}
        for k = 1:numel(sim_ic)
            ic  = sim_ic(k);
            val = hit_rate(ic)*100;
            bar(ax1, 2+k, val, 0.45, ...
                'FaceColor',sim_cols_h(k,:),'FaceAlpha',0.55, ...
                'EdgeColor',sim_cols_h(k,:),'LineWidth',1.2,'LineStyle','--', ...
                'DisplayName',[sim_lbl_h{k} ' (sim)']);
            text(ax1, 2+k, val+2, sprintf('%.0f%%',val), ...
                 'HorizontalAlignment','center','FontSize',7.5,'Color',[0.4 0.4 0.4]);
        end
        xline(ax1, 2.5, '--', 'Color',[0.6 0.6 0.6],'LineWidth',0.8, ...
              'Label','real | sim','LabelHorizontalAlignment','center', ...
              'FontSize',7,'HandleVisibility','off');
        x_ticks  = 1:(2+numel(sim_ic));
        x_labels = [real_labels, cellfun(@(l)[l ' (sim)'],sim_lbl_h,'un',0)];
        title(ax1,{'Hit rate'; 'solid = real (GDF events)   dashed = offline sim'},'FontSize',9);
    else
        x_ticks  = 1:2;
        x_labels = real_labels;
        title(ax1,'Hit rate — real GDF events (897/898/899)','FontSize',10);
    end

    legend(ax1,'FontSize',7,'Location','north');
    yline(ax1,50,'--k','FontSize',7,'HandleVisibility','off');
    set(ax1,'XTick',x_ticks,'XTickLabel',x_labels,'XTickLabelRotation',20,'YLim',[0,112]);
    ylabel('%'); grid(ax1,'on');

    % ── TTH: real (always) + MI-sim/CVSA-sim (hybrid only) ──────────────────
    ax2 = subplot(1,3,2); hold(ax2,'on');
    tth_sim_mat = cellfun(@(v) nanmean_safe(v), tth_cls);   % [n_cls x n_cond]
    tth_real    = real_tth_cls(1:n_cls);

    if is_hybrid
        sim_ic_tth = 2:n_cond;   % MI-sim and CVSA-sim only
        n_bars_grp = 1 + numel(sim_ic_tth);
    else
        sim_ic_tth = [];
        n_bars_grp = 1;
    end
    grp_w  = 0.7;
    bar_w  = grp_w / n_bars_grp;
    offs   = linspace(-grp_w/2+bar_w/2, grp_w/2-bar_w/2, n_bars_grp);

    for c = 1:n_cls
        % Real bar (black solid)
        if ~isnan(tth_real(c))
            bar(ax2, c+offs(1), tth_real(c), bar_w*0.9, ...
                'FaceColor',[0.15 0.15 0.15],'EdgeColor','k','LineWidth',1.5, ...
                'DisplayName','Real (GDF 897)');
            text(ax2, c+offs(1), tth_real(c)+0.12, sprintf('%.1fs',tth_real(c)), ...
                 'HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
        end
        % Simulated bars (hybrid only: MI-sim, CVSA-sim)
        for k = 1:numel(sim_ic_tth)
            ic  = sim_ic_tth(k);
            val = tth_sim_mat(c, ic);
            if ~isnan(val) && val > 0
                bar(ax2, c+offs(1+k), val, bar_w*0.9, ...
                    'FaceColor',cond_colors(ic,:),'FaceAlpha',0.55, ...
                    'EdgeColor',cond_colors(ic,:),'LineWidth',1.2,'LineStyle','--', ...
                    'DisplayName',[cond_labels{ic} ' (sim)']);
                text(ax2, c+offs(1+k), val+0.12, sprintf('%.1fs',val), ...
                     'HorizontalAlignment','center','FontSize',6.5,'Color',[0.4 0.4 0.4]);
            end
        end
    end
    if is_hybrid
        lbl_tth = [{'Real (GDF 897)'}, cellfun(@(l)[l ' (sim)'],cond_labels(2:n_cond),'un',0)];
        title(ax2,{'Mean TTH — HIT trials only'; 'black = real online   dashed = offline sim'},'FontSize',9);
    else
        lbl_tth = {'Real (GDF 897)'};
        title(ax2,'Mean TTH — HIT trials only (real GDF events)','FontSize',10);
    end
    legend(ax2, lbl_tth,'FontSize',7,'Location','north');
    set(ax2,'XTick',1:n_cls,'XTickLabel',cls_names,'XTickLabelRotation',15,'YLim',[0,inf]);
    ylabel('s'); grid(ax2,'on');

    % ── CDF of TTH ─────────────────────────────────────────────────────────
    ax3 = subplot(1,3,3); hold(ax3,'on');
    % Real TTH: thick black (always)
    if ~isempty(real_tth_all)
        xs=sort(real_tth_all); ys=(1:numel(xs))/numel(xs)*100;
        stairs(ax3,xs,ys,'k-','LineWidth',3,'DisplayName','Real (GDF 897)');
    end
    % Simulated: MI-sim and CVSA-sim (hybrid only)
    for k = 1:numel(sim_ic_tth)
        ic = sim_ic_tth(k);
        all_tth = [];
        for c = 1:n_cls
            if ~isempty(tth_cls{c,ic}), all_tth=[all_tth,tth_cls{c,ic}]; end %#ok<AGROW>
        end
        if ~isempty(all_tth)
            xs=sort(all_tth); ys=(1:numel(xs))/numel(xs)*100;
            stairs(ax3,xs,ys,'Color',cond_colors(ic,:),'LineWidth',2,'LineStyle','--', ...
                   'DisplayName',[cond_labels{ic} ' (sim)']);
        end
    end
    legend(ax3,'FontSize',8,'Location','best'); grid(ax3,'on');
    xlabel(ax3,'Time to hit (s)'); ylabel(ax3,'Cumulative % HIT trials');
    if is_hybrid
        title(ax3,{'CDF — time to hit'; 'black = real online   dashed = offline sim'},'FontSize',9);
    else
        title(ax3,'CDF — time to hit (real GDF events)','FontSize',10);
    end

    if is_hybrid
        sgtitle(fig3,sprintf('Performance — %s  [%s]   |  solid border = real GDF events  |  dashed = offline simulation', ...
                basename,upper(paradigm)),'FontSize',10,'Interpreter','none');
    else
        sgtitle(fig3,sprintf('Performance — %s  [%s]   |  accuracy from real GDF events (897/898/899)', ...
                basename,upper(paradigm)),'FontSize',10,'Interpreter','none');
    end
    saveas(fig3, fullfile(out_dir,sprintf('%s_03_performance.svg',basename)),'svg');

    % ── Fig 4: Saved trials (hybrid) or Hit/Miss output distribution (single) ───
    fig4 = figure('Name',sprintf('Trial Analysis — %s',basename), ...
                  'Color','w','NumberTitle','off','Position',[70 70 1000 420]);

    if is_hybrid
        % ── Hybrid: show saved-trial breakdown + early classifier values ──────
        n_saved       = sum(saved_mask);
        n_both_miss   = sum(~[trials_mi.pass] & ~[trials_hyb.pass]);
        n_both_hit    = sum( [trials_mi.pass] &  [trials_hyb.pass]);
        n_mi_only_hit = sum( [trials_mi.pass] & ~[trials_hyb.pass]);
        fprintf('  Both HIT=%d  Saved=%d  Both MISS=%d  HybridWorse=%d\n', ...
                n_both_hit, n_saved, n_both_miss, n_mi_only_hit);

        ax1 = subplot(1,3,1);
        vals_bar = [n_both_hit, n_saved, n_mi_only_hit, n_both_miss];
        bh4 = bar(ax1, vals_bar,'FaceColor','flat');
        bh4.CData = [0.2 0.7 0.2; 0.2 0.5 0.8; 0.9 0.5 0.1; 0.7 0.2 0.2];
        set(ax1,'XTickLabel',{'Both HIT','Saved','MI only','Both MISS'},'XTickLabelRotation',15);
        ylabel('trials'); grid(ax1,'on');
        title(ax1,sprintf('Trial classification  (n=%d)',n_trials),'FontSize',10);
        for k=1:numel(vals_bar)
            text(ax1,k,vals_bar(k)+0.3,sprintf('%d (%.0f%%)',vals_bar(k),100*vals_bar(k)/n_trials), ...
                 'HorizontalAlignment','center','FontSize',8);
        end

        ax2 = subplot(1,3,2);
        if ~isempty(early_cv_saved) && ~isempty(early_cv_other)
            boxplot(ax2,[early_cv_saved(:);early_cv_other(:)], ...
                    [ones(numel(early_cv_saved),1);2*ones(numel(early_cv_other),1)], ...
                    'Labels',{'Saved','Other'},'Colors','bk','Symbol','o');
        end
        hold(ax2,'on'); yline(ax2,0.5,'k--','FontSize',8);
        ylabel(sprintf('P_{CVSA}(target) at onset [%.0fs]',N_EARLY_S),'FontSize',9);
        title(ax2,'CVSA at trial onset','FontSize',10); grid(ax2,'on');
        text(ax2,0.5,0.02,sprintf('Saved:%.2f±%.2f  Other:%.2f±%.2f', ...
             mean_safe(early_cv_saved),std_safe(early_cv_saved), ...
             mean_safe(early_cv_other),std_safe(early_cv_other)), ...
             'Units','normalized','HorizontalAlignment','center','FontSize',8);

        ax3 = subplot(1,3,3);
        if ~isempty(early_mi_saved) && ~isempty(early_mi_other)
            boxplot(ax3,[early_mi_saved(:);early_mi_other(:)], ...
                    [ones(numel(early_mi_saved),1);2*ones(numel(early_mi_other),1)], ...
                    'Labels',{'Saved','Other'},'Colors','rk','Symbol','o');
        end
        hold(ax3,'on'); yline(ax3,0.5,'k--','FontSize',8);
        ylabel(sprintf('P_{MI}(target) at onset [%.0fs]',N_EARLY_S),'FontSize',9);
        title(ax3,'MI at trial onset','FontSize',10); grid(ax3,'on');
        text(ax3,0.5,0.02,sprintf('Saved:%.2f±%.2f  Other:%.2f±%.2f', ...
             mean_safe(early_mi_saved),std_safe(early_mi_saved), ...
             mean_safe(early_mi_other),std_safe(early_mi_other)), ...
             'Units','normalized','HorizontalAlignment','center','FontSize',8);

        sgtitle(fig4,sprintf('Saved trials: MI fails → hybrid succeeds — %s',basename), ...
                'FontSize',11,'Interpreter','none');

    else
        % ── Single-modality: show classifier output distribution HIT vs MISS ─
        hit_vals  = pt_mean_hit(~isnan(pt_mean_hit));
        miss_vals = pt_mean_miss(~isnan(pt_mean_miss));

        ax1 = subplot(1,3,1);
        if ~isempty(hit_vals) && ~isempty(miss_vals)
            boxplot(ax1,[hit_vals(:);miss_vals(:)], ...
                    [ones(numel(hit_vals),1);2*ones(numel(miss_vals),1)], ...
                    'Labels',{'HIT','MISS'},'Colors','gr','Symbol','o');
        end
        hold(ax1,'on'); yline(ax1,0.5,'k--'); grid(ax1,'on');
        ylabel('Mean P(target class) during CF','FontSize',9);
        title(ax1,'Classifier output: HIT vs MISS','FontSize',10);
        text(ax1,0.5,0.02,sprintf('HIT:%.2f±%.2f  MISS:%.2f±%.2f', ...
             mean_safe(hit_vals),std_safe(hit_vals), ...
             mean_safe(miss_vals),std_safe(miss_vals)), ...
             'Units','normalized','HorizontalAlignment','center','FontSize',8);

        ax2 = subplot(1,3,2);
        hold(ax2,'on');
        hit_idx  = find(~isnan(pt_mean_hit));
        miss_idx = find(~isnan(pt_mean_miss));
        scatter(ax2,hit_idx,  pt_mean_hit(hit_idx),  45,[0.1 0.6 0.1],'filled','DisplayName','HIT');
        scatter(ax2,miss_idx, pt_mean_miss(miss_idx), 45,[0.8 0.1 0.1],'s','DisplayName','MISS');
        yline(ax2,0.5,'k:'); legend(ax2,'FontSize',8); grid(ax2,'on');
        xlabel(ax2,'trial index'); ylabel(ax2,'mean P(target class)');
        title(ax2,'Per-trial output progression','FontSize',10);

        ax3 = subplot(1,3,3);
        all_early = arrayfun(@(t) mean_early_raw(trials_hyb(t), n_early, n_trials), 1:n_trials);
        vld_e = ~isnan(all_early);
        pass_v = [trials_hyb.pass];
        if any(vld_e)
            grp = 2 - double(pass_v(vld_e));  % 1=HIT, 2=MISS
            boxplot(ax3, all_early(vld_e), grp, ...
                    'Labels',{'HIT','MISS'},'Colors','gr','Symbol','o');
        end
        hold(ax3,'on'); yline(ax3,0.5,'k--'); grid(ax3,'on');
        ylabel(sprintf('P(target class) at onset [first %.0fs]',N_EARLY_S),'FontSize',9);
        title(ax3,sprintf('%s output at trial onset',upper(paradigm)),'FontSize',10);

        sgtitle(fig4,sprintf('Classifier output analysis — %s [%s]',basename,upper(paradigm)), ...
                'FontSize',11,'Interpreter','none');
    end
    saveas(fig4,fullfile(out_dir,sprintf('%s_04_trial_analysis.svg',basename)),'svg');

    % ── Fig 5: Classifier correlation (hybrid only) ───────────────────────────
    if is_hybrid
        rho_valid = rho_trial(~isnan(rho_trial));
        fig5 = figure('Name',sprintf('Classifier Correlation — %s',basename), ...
                      'Color','w','NumberTitle','off','Position',[80 80 1000 380]);

        ax1 = subplot(1,3,1);
        histogram(ax1,rho_valid,15,'FaceColor',[0.5 0.65 0.88],'EdgeColor','w','Normalization','probability');
        hold(ax1,'on');
        xline(ax1,mean(rho_valid),'r-','LineWidth',2.5, ...
              'Label',sprintf(' mean=%.2f',mean(rho_valid)),'FontSize',8);
        xline(ax1,0,'k--','LineWidth',1.2);
        xline(ax1,0.3,'--','Color',[0.8 0.5 0.1],'LineWidth',1, ...
              'Label',' ρ=0.3','LabelVerticalAlignment','bottom','FontSize',7);
        xlabel(ax1,'ρ(P_{MI}, P_{CVSA}) per trial'); ylabel('Fraction of trials');
        title(ax1,'Per-trial correlation','FontSize',10); grid(ax1,'on');
        text(ax1,0.05,0.93,'ρ < 0.3  →  LOP adds info','Units','normalized','FontSize',8,'Color',[0.3 0.5 0.1]);
        text(ax1,0.05,0.82,'ρ > 0.5  →  risk of overcounting','Units','normalized','FontSize',8,'Color',[0.7 0.2 0.1]);

        ax2 = subplot(1,3,2);
        % CDF of ρ
        xs = sort(rho_valid); ys = (1:numel(xs))/numel(xs)*100;
        plot(ax2,xs,ys,'b-','LineWidth',2);
        hold(ax2,'on');
        xline(ax2,0,'k--'); xline(ax2,0.3,'--','Color',[0.8 0.5 0.1]);
        xlabel(ax2,'ρ per trial'); ylabel(ax2,'Cumulative %');
        title(ax2,'CDF of per-trial ρ','FontSize',10); grid(ax2,'on');
        text(ax2,0.05,0.85,sprintf('%.0f%% trials have ρ < 0.3', ...
             100*mean(rho_valid<0.3)),'Units','normalized','FontSize',9);

        ax3 = subplot(1,3,3);
        % Scatter P_MI vs P_CVSA for all CF frames (subsample for speed)
        all_pm = []; all_pc = [];
        for t = 1:n_trials
            np = trials_hyb(t).n_pre;
            pm = trials_hyb(t).p_mi(np+1:end,1);
            pc = trials_hyb(t).p_cvsa(np+1:end,1);
            vld = ~isnan(pm)&~isnan(pc);
            all_pm = [all_pm; pm(vld)]; all_pc = [all_pc; pc(vld)]; %#ok<AGROW>
        end
        n_plot = min(1500, numel(all_pm));
        idx = randperm(numel(all_pm), n_plot);
        scatter(ax3,all_pm(idx),all_pc(idx),6,'filled', ...
                'MarkerFaceColor',[0.5 0.5 0.8],'MarkerFaceAlpha',0.3);
        hold(ax3,'on');
        xline(ax3,0.5,'k:'); yline(ax3,0.5,'k:');
        if numel(all_pm)>1
            r_all = corr(all_pm,all_pc);
            title(ax3,sprintf('P_{MI} vs P_{CVSA}  (all CF frames)\nρ_{global}=%.3f',r_all),'FontSize',10);
        end
        xlabel(ax3,'P_{MI}(c1)'); ylabel(ax3,'P_{CVSA}(c1)');
        axis(ax3,'square'); grid(ax3,'on');

        sgtitle(fig5,sprintf('Classifier correlation — %s  [LOP valid if ρ low]',basename), ...
                'FontSize',11,'Interpreter','none');
        saveas(fig5,fullfile(out_dir,sprintf('%s_05_correlation.svg',basename)),'svg');
    end

    % ── Print console summary ─────────────────────────────────────────────────
    fprintf('\n  ┌─ Summary ─────────────────────────────────────────\n');
    for ic = 1:n_cond
        all_tth_ic = [];
        for c=1:n_cls; all_tth_ic=[all_tth_ic,tth_cls{c,ic}]; end %#ok<AGROW>
        fprintf('  │  %-12s  hit=%.0f%%  tth=%.2fs±%.2fs\n', cond_labels{ic}, ...
                100*hit_rate(ic), mean_safe(all_tth_ic), std_safe(all_tth_ic));
    end
    if is_hybrid
        fprintf('  │  Saved trials: %d / %d  (%.0f%%)\n', ...
                sum(saved_mask), n_trials, 100*sum(saved_mask)/n_trials);
        fprintf('  │  ρ(P_MI,P_CVSA): mean=%.3f  std=%.3f\n', ...
                mean(rho_trial,'omitnan'), std(rho_trial,'omitnan'));
    end
    fprintf('  └───────────────────────────────────────────────────\n');

    % ── Save per-file mat ─────────────────────────────────────────────────────
    res = struct( ...
        'file_label',     file_label, ...
        'basename',       basename, ...
        'paradigm',       paradigm, ...
        'n_trials',       n_trials, ...
        'hit_rate',       hit_rate, ...          % simulated (offline integrator)
        'hit_rate_ev',    hit_rate_ev, ...       % event-based: 897/(897+898+899)
        'hit_rate_no_to', hit_rate_no_to, ...    % event-based: 897/(897+898)
        'timeout_rate',   timeout_rate, ...      % event-based: 899/(897+898+899)
        'n_hits_ev',      n_hits_ev, ...
        'n_misses_ev',    n_misses_ev, ...
        'n_timeout_ev',   n_timeout_ev, ...
        'real_tth_cls',   real_tth_cls, ...      % real TTH per class (from 897 events)
        'real_tth_all',   real_tth_all, ...      % pooled real TTH (all HIT trials)
        'cond_labels',    {cond_labels}, ...
        'tth_cls',        {tth_cls}, ...
        'n_saved',        sum(saved_mask), ...
        'rho_mean',       mean(rho_trial,'omitnan'), ...
        'rho_std',        std(rho_trial,'omitnan'), ...
        'p_mi_mean',      p_mi_mean, ...
        'p_cv_mean',      p_cv_mean, ...
        'buf_hit',        {buf_hit}, ...
        'buf_miss',       {buf_miss}, ...
        'early_cv_saved', early_cv_saved, ...
        'early_mi_saved', early_mi_saved, ...
        'saved_mask',     saved_mask, ...
        'rho_trial',      rho_trial, ...
        'thresholds',     thresholds, ...
        'framerate',      framerate, ...
        't_ax',           t_ax ...
    );
    mat_path = fullfile(gdf_dir, sprintf('advantage_%s_%s.mat',paradigm,basename));
    save(mat_path,'res','trials_hyb','trials_mi','trials_cvsa');
    fprintf('  Saved: %s\n\n', mat_path);

    all_res_cell{end+1} = res;
end

%% ═════════════════════════════════════════════════════════════════════════════
%  FIG 6: Multi-file overview (if more than one file)
%  ═════════════════════════════════════════════════════════════════════════════
if numel(all_res_cell) > 1
    n_f        = numel(all_res_cell);
    file_lbls  = cellfun(@(r) r.file_label, all_res_cell, 'UniformOutput',false);

    fig6 = figure('Name','Multi-File Overview','Color','w','NumberTitle','off', ...
                  'Position',[90 90 1400 480]);

    % ── Panel 1: event-based accuracy per file (3 bars) ──────────────────────
    % Data from actual GDF outcomes (897/898/899), not from offline simulation.
    ax1 = subplot(1,3,1);
    ev_mat = nan(n_f,3);  % [hit_rate_ev, hit_rate_no_to, timeout_rate] per file
    for i = 1:n_f
        r = all_res_cell{i};
        ev_mat(i,:) = [r.hit_rate_ev, r.hit_rate_no_to, r.timeout_rate];
    end
    ev_plot = ev_mat*100; ev_plot(isnan(ev_plot)) = 0;
    bev = bar(ax1, ev_plot, 'grouped');
    ev_colors = [0.18 0.45 0.75; 0.20 0.65 0.25; 0.85 0.30 0.10];
    for k = 1:3
        if isa(bev(k),'matlab.graphics.chart.primitive.Bar')
            bev(k).FaceColor = ev_colors(k,:);
        end
    end
    legend(ax1,{'897/(897+898+899)','897/(897+898)','timeout 899/total'}, ...
           'FontSize',7,'Location','north');
    set(ax1,'XTick',1:n_f,'XTickLabel',file_lbls,'XTickLabelRotation',25,'YLim',[0,108]);
    yline(ax1,50,'--k','FontSize',7); ylabel('%'); grid(ax1,'on');
    % Annotate hit count above each bar group
    for i = 1:n_f
        r = all_res_cell{i};
        text(ax1,i,max(ev_plot(i,:))+3, ...
             sprintf('H=%d M=%d T=%d',r.n_hits_ev,r.n_misses_ev,r.n_timeout_ev), ...
             'HorizontalAlignment','center','FontSize',6.5,'Color',[0.3 0.3 0.3]);
    end
    title(ax1,'Event-based accuracy (actual GDF outcomes)','FontSize',10);

    % ── Panel 2: per-paradigm summary (mean ± std across files) ──────────────
    ax2 = subplot(1,3,2);
    paradigms_found = unique(cellfun(@(r) r.paradigm, all_res_cell, 'UniformOutput',false));
    par_colors = struct('mi',[0.85 0.30 0.10],'cvsa',[0.10 0.60 0.30],'hybrid',[0.18 0.45 0.75]);
    hold(ax2,'on');
    for pi = 1:numel(paradigms_found)
        par = paradigms_found{pi};
        vals = cellfun(@(r) r.hit_rate_ev*100, all_res_cell(strcmp(cellfun(@(r)r.paradigm,all_res_cell,'un',0),par)));
        col_p = par_colors.(par);
        xp = pi;
        bar(ax2, xp, mean(vals,'omitnan'), 0.5, 'FaceColor',col_p,'FaceAlpha',0.65,'EdgeColor',col_p);
        errorbar(ax2, xp, mean(vals,'omitnan'), std(vals,'omitnan'), 'k.','LineWidth',1.5,'CapSize',8);
        scatter(ax2, xp + randn(numel(vals),1)*0.07, vals, 40, col_p, 'filled', ...
                'MarkerEdgeColor','k','LineWidth',0.5);
    end
    set(ax2,'XTick',1:numel(paradigms_found),'XTickLabel',upper(paradigms_found), ...
            'XTickLabelRotation',15,'YLim',[0,108]);
    yline(ax2,50,'--k','FontSize',7); ylabel('%'); grid(ax2,'on');
    title(ax2,'Hit rate by paradigm  (mean ± std + dots)','FontSize',10);

    % ── Panel 3: within-hybrid advantage (LOP vs MI-sim, hybrid files only) ──
    ax3 = subplot(1,3,3);
    hold(ax3,'on');
    hyb_idx = find(strcmp(cellfun(@(r)r.paradigm,all_res_cell,'un',0),'hybrid'));
    if ~isempty(hyb_idx)
        delta_h = nan(numel(hyb_idx),1);
        lbl_h   = cell(numel(hyb_idx),1);
        for k = 1:numel(hyb_idx)
            r = all_res_cell{hyb_idx(k)};
            % hit_rate: [Hybrid, MI-sim, CVSA-sim]
            if numel(r.hit_rate) >= 2
                delta_h(k) = (r.hit_rate(1) - r.hit_rate(2)) * 100;
            end
            lbl_h{k} = r.file_label;
        end
        bh3 = bar(ax3, delta_h, 'FaceColor','flat');
        if isa(bh3,'matlab.graphics.chart.primitive.Bar')
            for k = 1:numel(hyb_idx)
                if isnan(delta_h(k)),   bh3.CData(k,:) = [0.7 0.7 0.7];
                elseif delta_h(k) >= 0, bh3.CData(k,:) = [0.18 0.55 0.25];
                else,                    bh3.CData(k,:) = [0.80 0.20 0.10];
                end
            end
        end
        set(ax3,'XTick',1:numel(hyb_idx),'XTickLabel',lbl_h,'XTickLabelRotation',20);
        yline(ax3,0,'k-','LineWidth',1.2);
        for k = 1:numel(hyb_idx)
            if ~isnan(delta_h(k))
                text(ax3,k,delta_h(k)+sign(delta_h(k))*1.5,sprintf('%.0f%%',delta_h(k)), ...
                     'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
            end
        end
    else
        text(ax3,0.5,0.5,'No hybrid files loaded','Units','normalized', ...
             'HorizontalAlignment','center','FontSize',10,'Color',[0.5 0.5 0.5]);
        axis(ax3,'off');
    end
    grid(ax3,'on'); ylabel('Δ hit rate  hybrid − MI-only (%)','FontSize',9);
    title(ax3,'Within-hybrid advantage  (same EEG, LOP vs MI alone)','FontSize',10);
    subtitle(ax3,'green = fusion helped, red = fusion hurt','FontSize',8);

    sgtitle(fig6,'Multi-file overview — hybrid BCI advantage','FontSize',13);
    saveas(fig6, fullfile(out_dir,'00_multifile_overview.svg'),'svg');
end

fprintf('\nAll done. Figures saved to:\n  %s\n', out_dir);

%% ── Local helpers ────────────────────────────────────────────────────────────
function v = mean_safe(x)
    if isempty(x), v = NaN; else, v = mean(x(:),'omitnan'); end
end

function v = std_safe(x)
    if numel(x)<2, v = NaN; else, v = std(x(:),'omitnan'); end
end

function v = nanmean_safe(x)
    if isempty(x), v = NaN; else, v = mean(x(:),'omitnan'); end
end

function v = mean_early_raw(tr, n_early, ~)
    c  = tr.target_class;
    np = tr.n_pre;
    if isnan(c)
        v = NaN; return
    end
    ne = min(n_early, tr.n_cf);
    rv = tr.raw(np+1:np+ne, c);
    rv = rv(~isnan(rv));
    if isempty(rv), v = NaN; else, v = mean(rv); end
end
