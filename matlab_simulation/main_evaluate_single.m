%% MAIN_EVALUATE_SINGLE  Per-file detailed evaluation: metrics + EEG analysis.
%
%   For each selected GDF:
%     Fig 1 — metrics summary: trial acc, sLDA frame acc (MI/CVSA/fused for hybrid),
%              time metrics (TTH+miss peak), mean classifier output on miss/timeout
%     Fig 2,3 — ERD/ERS per class: channel time-courses + heatmaps
%     Fig 4   — temporal topoplot (3 rows: cls1 / cls2 / diff × CUE + 0.5s CF windows)
%     Fig 5   — FBCSP log-variance features over CF time + sLDA discriminant topoplot
%     Fig 6   — EEG-behavior correlation (sLDA output distribution + integrator curves)
%   Saves SVG to <gdf_dir>/results/<basename>/ and eval_single_<paradigm>_<basename>.mat.

clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir,'io'), fullfile(this_dir,'processing'), ...
        fullfile(this_dir,'artifacts'), fullfile(this_dir,'classifier'), ...
        fullfile(this_dir,'integrator'), fullfile(this_dir,'plotting'), ...
        fullfile(this_dir,'utils'));

HIT_CODE     = 897;
MISS_CODE    = 898;
TIMEOUT_CODE = 899;
CF_CODE      = 781;

% ── File picker ───────────────────────────────────────────────────────────────
default_dir = '/home/paolo/bci_vr_ws/recordings';
if ~isfolder(default_dir), default_dir = fileparts(this_dir); end
[gdf_names, gdf_dir] = uigetfile({'*.gdf','GDF recordings (*.gdf)'}, ...
    'Select GDF file(s)', default_dir, 'MultiSelect','on');
if isequal(gdf_names,0), error('main_evaluate_single:cancel','No file selected.'); end
if ischar(gdf_names), gdf_names = {gdf_names}; end

for fi = 1:numel(gdf_names)
    gdf_path = fullfile(gdf_dir, gdf_names{fi});
    [~, basename] = fileparts(gdf_path);
    fprintf('\n[%d/%d]  %s\n', fi, numel(gdf_names), basename);

    res_dir = fullfile(gdf_dir, 'results', basename);
    if ~isfolder(res_dir), mkdir(res_dir); end

    % ── Load ─────────────────────────────────────────────────────────────────
    [signal, header, ~] = load_gdf(gdf_path);
    [params, ~]         = load_params_yaml(gdf_path);
    header_orig         = header;

    paradigm   = params.integrator.paradigm;
    fs         = double(params.acquisition.samplerate);
    framerate  = double(params.acquisition.framerate);
    chunk_size = round(fs / framerate);
    if abs(fs - header.SampleRate) > 1e-3
        fs = header.SampleRate; chunk_size = round(fs / framerate);
    end
    bufsize_proc = double(params.RingBufferCfg.params.size);
    bufsize_art  = double(params.RingBufferCfgArtifact.params.size);
    eog_names    = to_strcell(params.CarCfg.params.EOG_ch_names);

    do_car_mi   = true;
    do_car_cvsa = true;
    if isfield(params,'processing_fbcsp_mi'),   do_car_mi   = logical(params.processing_fbcsp_mi.do_car);   end
    if isfield(params,'processing_fbcsp_cvsa'), do_car_cvsa = logical(params.processing_fbcsp_cvsa.do_car); end

    use_mi   = ismember(paradigm, {'mi','hybrid'});
    use_cvsa = ismember(paradigm, {'cvsa','hybrid'});
    csp_mi=[]; slda_mi=[]; csp_cvsa=[]; slda_cvsa=[];
    if use_mi,   csp_mi   = load_csp(params,'mi');   slda_mi   = load_slda(params,'mi');   end
    if use_cvsa, csp_cvsa = load_csp(params,'cvsa'); slda_cvsa = load_slda(params,'cvsa'); end

    feat_mi=[]; header_mi=[]; feat_pre_mi=[];
    feat_cv=[]; header_cv=[]; feat_pre_cv=[]; info_proc=[];
    proc_base = struct('samplerate',fs,'chunk_size',chunk_size, ...
                       'bufsize',bufsize_proc,'filter_order',4);
    if use_mi
        cfg = proc_base; cfg.do_car = do_car_mi; cfg.eog_names = eog_names;
        [feat_mi, header_mi, info_proc, feat_pre_mi] = apply_processing(signal, header, csp_mi, cfg);
    end
    if use_cvsa
        cfg = proc_base; cfg.do_car = do_car_cvsa; cfg.eog_names = eog_names;
        [feat_cv, header_cv, info2, feat_pre_cv] = apply_processing(signal, header, csp_cvsa, cfg);
        if isempty(info_proc), info_proc = info2; end
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
    if ~isfield(int_cfg,'cvsa_influence'),       int_cfg.cvsa_influence = 2.5; end
    if ~isfield(int_cfg,'thresholds') || isempty(int_cfg.thresholds)
        int_cfg.thresholds = params.training_node.thresholds;
    end

    header_chunks = header_mi; if ~use_mi, header_chunks = header_cv; end
    header_chunks.framerate = framerate;

    trials   = integrate_signal(p_mi_al, p_cv_al, art_flags, header_chunks, int_cfg, paradigm);
    n_trials = numel(trials);
    if n_trials == 0
        fprintf('  [warn] no trials — skipping.\n'); continue
    end

    thresholds = to_vec(int_cfg.thresholds);
    classes    = to_vec(int_cfg.classes);
    n_cls      = numel(classes);
    init_val   = to_vec(int_cfg.init_val);
    p_rest     = init_val(1);
    cls_names  = arrayfun(@(x) num2str(x), classes, 'UniformOutput', false);

    % ── Event-based accuracy ──────────────────────────────────────────────────
    POS_o = header_orig.EVENT.POS; TYP_o = header_orig.EVENT.TYP;
    cf_pos_samp = POS_o(TYP_o == CF_CODE);
    trial_outcome_ev = zeros(1, n_trials);
    for t = 1:min(n_trials, numel(cf_pos_samp))
        after = POS_o > cf_pos_samp(t);
        idx   = find((TYP_o==HIT_CODE|TYP_o==MISS_CODE|TYP_o==TIMEOUT_CODE) & after, 1);
        if ~isempty(idx), trial_outcome_ev(t) = TYP_o(idx); end
    end
    n_hits_ev    = sum(trial_outcome_ev == HIT_CODE);
    n_misses_ev  = sum(trial_outcome_ev == MISS_CODE);
    n_timeout_ev = sum(trial_outcome_ev == TIMEOUT_CODE);
    trial_acc_events     = n_hits_ev / max(1, n_hits_ev + n_misses_ev + n_timeout_ev);
    trial_acc_no_timeout = n_hits_ev / max(1, n_hits_ev + n_misses_ev);
    trial_acc            = mean([trials.pass]);

    % ── Artifact rate ─────────────────────────────────────────────────────────
    art_rate_per_trial = arrayfun(@(tr) mean(tr.artifact(tr.n_pre+1:end)), trials);
    art_rate           = mean(art_rate_per_trial);
    trial_acc_no_rej   = mean([trials(art_rate_per_trial <= 0.30).pass]);

    % ── Per-class metrics ─────────────────────────────────────────────────────
    tth_per_class    = nan(1,n_cls);
    t_miss_peak_cls  = nan(1,n_cls);
    hit_rate_cls     = nan(1,n_cls);
    miss_max_fused   = nan(1,n_cls);
    miss_mean_mi     = nan(1,n_cls);
    miss_mean_cvsa   = nan(1,n_cls);
    sa_miss_fused    = nan(1,n_cls);
    sa_timeout_fused = nan(1,n_cls);
    sa_miss_mi       = nan(1,n_cls);
    sa_timeout_mi    = nan(1,n_cls);
    sa_miss_cvsa     = nan(1,n_cls);
    sa_timeout_cvsa  = nan(1,n_cls);

    for c = 1:n_cls
        mask_c_idx = find(arrayfun(@(tr) ~isnan(tr.target_class)&&tr.target_class==c, trials));
        if isempty(mask_c_idx), continue; end
        trs_c = trials(mask_c_idx);
        hit_rate_cls(c) = mean([trs_c.pass]);
        tth_v=[]; miss_v=[]; maxint_v=[]; miv=[]; cvv=[];
        fmiss=[]; ftout=[]; mi_miss=[]; mi_tout=[]; cv_miss=[]; cv_tout=[];
        for t = 1:numel(trs_c)
            np   = trs_c(t).n_pre;
            ic   = trs_c(t).integrated(np+1:end, c);
            t_ax = (0:trs_c(t).n_cf-1) / framerate;
            if trs_c(t).pass
                hf = find(ic >= thresholds(c), 1);
                if ~isempty(hf), tth_v(end+1) = t_ax(hf); end %#ok<AGROW>
            else
                [~,pf] = max(ic);
                miss_v(end+1) = t_ax(pf); %#ok<AGROW>
                maxint_v(end+1) = max(ic(~isnan(ic))); %#ok<AGROW>
                if strcmp(paradigm,'hybrid')
                    pm = trs_c(t).p_mi(np+1:end,c);
                    pv = trs_c(t).p_cvsa(np+1:end,c);
                    miv(end+1)  = mean(pm(~isnan(pm)));  %#ok<AGROW>
                    cvv(end+1)  = mean(pv(~isnan(pv)));  %#ok<AGROW>
                end
            end
            % sample accuracy for miss/timeout trials
            oc = trial_outcome_ev(mask_c_idx(t));
            if oc == MISS_CODE || oc == TIMEOUT_CODE
                raw_cf = trs_c(t).raw(np+1:end,:);
                vld = ~any(isnan(raw_cf),2);
                if any(vld)
                    sa_f = mean(raw_cf(vld,c) == max(raw_cf(vld,:),[],2));
                    if oc == MISS_CODE
                        fmiss(end+1) = sa_f; %#ok<AGROW>
                    else
                        ftout(end+1) = sa_f; %#ok<AGROW>
                    end
                    if strcmp(paradigm,'hybrid')
                        pm2 = trs_c(t).p_mi(np+1:end,:);   vm2 = ~any(isnan(pm2),2);
                        pv2 = trs_c(t).p_cvsa(np+1:end,:); vv2 = ~any(isnan(pv2),2);
                        if any(vm2)
                            sa_mi2 = mean(pm2(vm2,c)==max(pm2(vm2,:),[],2));
                            if oc==MISS_CODE, mi_miss(end+1)=sa_mi2; else, mi_tout(end+1)=sa_mi2; end %#ok<AGROW>
                        end
                        if any(vv2)
                            sa_cv2 = mean(pv2(vv2,c)==max(pv2(vv2,:),[],2));
                            if oc==MISS_CODE, cv_miss(end+1)=sa_cv2; else, cv_tout(end+1)=sa_cv2; end %#ok<AGROW>
                        end
                    end
                end
            end
        end
        if ~isempty(tth_v),     tth_per_class(c)    = mean(tth_v);     end
        if ~isempty(miss_v),    t_miss_peak_cls(c)   = mean(miss_v);   end
        if ~isempty(maxint_v),  miss_max_fused(c)    = mean(maxint_v); end
        if ~isempty(miv),       miss_mean_mi(c)      = mean(miv);      end
        if ~isempty(cvv),       miss_mean_cvsa(c)    = mean(cvv);      end
        if ~isempty(fmiss),     sa_miss_fused(c)     = mean(fmiss);    end
        if ~isempty(ftout),     sa_timeout_fused(c)  = mean(ftout);    end
        if ~isempty(mi_miss),   sa_miss_mi(c)        = mean(mi_miss);  end
        if ~isempty(mi_tout),   sa_timeout_mi(c)     = mean(mi_tout);  end
        if ~isempty(cv_miss),   sa_miss_cvsa(c)      = mean(cv_miss);  end
        if ~isempty(cv_tout),   sa_timeout_cvsa(c)   = mean(cv_tout);  end
    end

    % ── Sample accuracy per class per modality ────────────────────────────────
    sample_acc_fused  = nan(1,n_cls);
    sample_acc_mi_cls = nan(1,n_cls);
    sample_acc_cv_cls = nan(1,n_cls);
    trial_sa_list     = nan(n_trials,1);

    for t = 1:n_trials
        c = trials(t).target_class; if isnan(c), continue; end
        np  = trials(t).n_pre;
        raw = trials(t).raw(np+1:end,:);
        vld = ~any(isnan(raw),2);
        if any(vld), trial_sa_list(t) = mean(raw(vld,c)==max(raw(vld,:),[],2)); end
    end
    for c = 1:n_cls
        cf_f=[]; cf_mi=[]; cf_cv=[];
        for t = 1:n_trials
            if isnan(trials(t).target_class)||trials(t).target_class~=c, continue; end
            np = trials(t).n_pre;
            rf = trials(t).raw(np+1:end,:);
            vf = ~any(isnan(rf),2);
            if any(vf), cf_f = [cf_f; rf(vf,c)==max(rf(vf,:),[],2)]; end %#ok<AGROW>
            if strcmp(paradigm,'hybrid')
                rm = trials(t).p_mi(np+1:end,:);
                vm = ~any(isnan(rm),2);
                if any(vm), cf_mi = [cf_mi; rm(vm,c)==max(rm(vm,:),[],2)]; end %#ok<AGROW>
                rc = trials(t).p_cvsa(np+1:end,:);
                vc = ~any(isnan(rc),2);
                if any(vc), cf_cv = [cf_cv; rc(vc,c)==max(rc(vc,:),[],2)]; end %#ok<AGROW>
            end
        end
        if ~isempty(cf_f),  sample_acc_fused(c)   = mean(cf_f);  end
        if ~isempty(cf_mi), sample_acc_mi_cls(c)  = mean(cf_mi); end
        if ~isempty(cf_cv), sample_acc_cv_cls(c)  = mean(cf_cv); end
    end
    sample_acc_total      = mean(sample_acc_fused, 'omitnan');
    sample_acc_mean_trial = mean(trial_sa_list,    'omitnan');

    fprintf('  paradigm=%-8s  trials=%d  acc_ev=%.0f%%  acc_no_tout=%.0f%%  art=%.0f%%\n', ...
            paradigm, n_trials, 100*trial_acc_events, 100*trial_acc_no_timeout, 100*art_rate);

    % ── Primary modality for EEG features ────────────────────────────────────
    if use_mi
        feat_main = feat_mi; feat_pre = feat_pre_mi;
        csp_main  = csp_mi;  slda_main = slda_mi; mod_label = 'MI';
    else
        feat_main = feat_cv; feat_pre = feat_pre_cv;
        csp_main  = csp_cvsa; slda_main = slda_cvsa; mod_label = 'CVSA';
    end
    sel_ch   = csp_main.selected_channels;
    n_sel_ch = numel(sel_ch);
    n_bands  = csp_main.n_bands;
    bands    = csp_main.bands;
    n_comp   = csp_main.n_components;

    % feat3d: [n_chunks × n_sel_ch × n_bands]
    feat3d = permute(reshape(feat_pre.', n_sel_ch, n_bands, size(feat_pre,1)), [3,1,2]);

    % ── ERD/ERS (1s pre-CF baseline) ─────────────────────────────────────────
    n_cf_min = min(arrayfun(@(tr) tr.n_cf, trials));
    n_bl     = round(framerate);
    erd_sum  = zeros(n_cls, n_cf_min, n_sel_ch, n_bands);
    erd_cnt  = zeros(n_cls,1);

    for t = 1:n_trials
        c = trials(t).target_class; if isnan(c), continue; end
        sc    = trials(t).start_chunk;
        n_use = min(trials(t).n_cf, n_cf_min);
        bl_r  = max(1,sc-n_bl):sc-1;
        cf_r  = sc:sc+n_use-1;
        bl_r  = bl_r(bl_r>=1 & bl_r<=size(feat3d,1));
        cf_r  = cf_r(cf_r>=1 & cf_r<=size(feat3d,1));
        if isempty(bl_r)||numel(cf_r)<n_use, continue; end
        pbl = mean(feat3d(bl_r,:,:),1);
        pcf = feat3d(cf_r,:,:);
        if any(isnan(pbl(:)))||any(isnan(pcf(:))), continue; end
        erd_t = (pcf-pbl)./max(pbl,eps)*100;
        erd_sum(c,1:n_use,:,:) = erd_sum(c,1:n_use,:,:) + reshape(erd_t,1,n_use,n_sel_ch,n_bands);
        erd_cnt(c) = erd_cnt(c)+1;
    end
    erd_mean = erd_sum;
    for c = 1:n_cls
        if erd_cnt(c)>0, erd_mean(c,:,:,:) = erd_sum(c,:,:,:)/erd_cnt(c); end
    end
    t_ax_cf = (0:n_cf_min-1)/framerate;

    erd_static = nan(n_cls,n_sel_ch);
    for c = 1:n_cls
        if erd_cnt(c)>0
            tmp = squeeze(mean(erd_mean(c,:,:,:),2));
            erd_static(c,:) = sum(reshape(tmp,n_sel_ch,n_bands),2).';
        end
    end

    %% ── Fig 1: Metrics summary ───────────────────────────────────────────────
    fig1 = figure('Name',sprintf('Metrics — %s',basename),'Color','w', ...
                  'NumberTitle','off','Position',[40 40 1200 660]);

    subplot(2,3,1);
    bh = bar([trial_acc_events, trial_acc_no_timeout]*100);
    bh.FaceColor='flat'; bh.CData=[.22 .50 .78;.28 .68 .30];
    set(gca,'XTick',1:2,'XTickLabel',{'897/(897+898+899)','897/(897+898)'},...
        'XTickLabelRotation',15,'YLim',[0,105]); grid on;
    ylabel('%'); title('Trial accuracy (event-based)','FontSize',9);
    yline(50,'--k','chance','FontSize',6,'LabelHorizontalAlignment','left');
    text(0.98,0.95,sprintf('H=%d  M=%d  TO=%d',n_hits_ev,n_misses_ev,n_timeout_ev),...
         'Units','normalized','HorizontalAlignment','right','FontSize',7);

    subplot(2,3,2);
    if strcmp(paradigm,'hybrid')
        bar([sample_acc_fused;sample_acc_mi_cls;sample_acc_cv_cls].'*100,'grouped');
        legend({'Fused','MI','CVSA'},'FontSize',7,'Location','south');
    else
        bar(sample_acc_fused*100,'FaceColor',[.45 .35 .75]);
    end
    set(gca,'XTick',1:n_cls,'XTickLabel',cls_names,'YLim',[0,105]); grid on;
    ylabel('%'); title('sLDA frame accuracy per class','FontSize',9);
    yline(50,'--k');

    subplot(2,3,3);
    cnt_mat = zeros(n_cls,3);
    for c = 1:n_cls
        mc = arrayfun(@(tr)~isnan(tr.target_class)&&tr.target_class==c, trials);
        ov = trial_outcome_ev(mc);
        cnt_mat(c,:) = [sum(ov==HIT_CODE), sum(ov==MISS_CODE), sum(ov==TIMEOUT_CODE)];
    end
    bar(cnt_mat,'stacked');
    legend({'Hit (897)','Miss (898)','Timeout (899)'},'FontSize',7,'Location','north');
    set(gca,'XTick',1:n_cls,'XTickLabel',cls_names); grid on;
    ylabel('count'); title('Trial outcomes per class','FontSize',9);

    subplot(2,3,4);
    t_data = [tth_per_class; t_miss_peak_cls];
    bh2 = bar(t_data.','grouped');
    bh2(1).FaceColor=[.30 .60 .90]; bh2(2).FaceColor=[.90 .50 .25];
    legend({'Time to hit','Time to miss peak'},'FontSize',7,'Location','north');
    set(gca,'XTick',1:n_cls,'XTickLabel',cls_names); grid on;
    ylabel('s'); title('Time metrics per class','FontSize',9);

    ax5 = subplot(2,3,5);
    if strcmp(paradigm,'hybrid')
        d5 = [sa_miss_fused; sa_timeout_fused; sa_miss_mi; sa_timeout_mi; sa_miss_cvsa; sa_timeout_cvsa].';
        bar(ax5, d5, 'grouped');
        legend(ax5,{'fused-miss','fused-to','MI-miss','MI-to','CVSA-miss','CVSA-to'},...
               'FontSize',5,'Location','southoutside','Orientation','horizontal');
    else
        d5 = [sa_miss_fused; sa_timeout_fused].';
        bar(ax5, d5, 'grouped');
        legend(ax5,{'miss','timeout'},'FontSize',7,'Location','south');
    end
    hold(ax5,'on');
    yline(ax5,0.5,'--k','chance','FontSize',6,'LabelHorizontalAlignment','left');
    set(ax5,'XTick',1:n_cls,'XTickLabel',cls_names,'YLim',[0,1.05]); grid(ax5,'on');
    ylabel(ax5,'fraction correct frames'); title(ax5,'Sample accuracy on miss/timeout','FontSize',9);

    ax6 = subplot(2,3,6);
    scatter(ax6, 1:n_trials, art_rate_per_trial*100, 38, [.20 .55 .20], 'filled');
    hold(ax6,'on');
    yline(ax6, 100*art_rate, '--r', sprintf('mean=%.0f%%',100*art_rate), ...
          'FontSize',7, 'LabelHorizontalAlignment','right', 'LineWidth',1.2);
    yline(ax6, 30, ':k', '30%', 'FontSize',6, 'LabelHorizontalAlignment','left');
    xlabel(ax6,'trial index'); ylabel(ax6,'artifact frames (%)');
    title(ax6,sprintf('Artifact rate per trial  (mean=%.0f%%)',100*art_rate),'FontSize',9);
    set(ax6,'XLim',[0,n_trials+1],'YLim',[0,max(105,max(art_rate_per_trial)*100+10)]);
    grid(ax6,'on');

    sgtitle(sprintf('%s  [%s]  %d trials', basename, upper(paradigm), n_trials), ...
            'FontSize',11,'Interpreter','none');
    saveas(fig1, fullfile(res_dir,'01_metrics.svg'),'svg');

    %% ── Fig 2,3: ERD/ERS time course + heatmap, one figure per class ────────
    cmap_rd = rdbu_cmap();
    for c = 1:n_cls
        fig_erd = figure('Name',sprintf('ERD/ERS cls %s — %s',cls_names{c},basename), ...
                         'Color','w','NumberTitle','off', ...
                         'Position',[60+c*20 30 900 max(320,160*n_bands+80)]);
        clim_e = max(1, max(abs(squeeze(erd_mean(c,:,:,:))), [], 'all'));

        for b = 1:n_bands
            band_lbl = sprintf('[%.0f–%.0f Hz]',bands(b,1),bands(b,2));
            data_b   = squeeze(erd_mean(c,:,:,b));  % [n_cf_min × n_sel_ch]

            axL = subplot(n_bands,2,(b-1)*2+1);
            if ~any(isnan(data_b(:)))
                plot(axL, t_ax_cf, data_b, 'Color',[.75 .75 .75],'LineWidth',.7);
                hold(axL,'on');
                plot(axL, t_ax_cf, mean(data_b,2),'k-','LineWidth',2.0,'DisplayName','mean');
                yline(axL,0,'k:','LineWidth',.8);
                set(axL,'XLim',[t_ax_cf(1),t_ax_cf(end)]);
                ylim_v = [-clim_e*0.6, clim_e*0.6]; if diff(ylim_v)>0, ylim(axL,ylim_v); end
                if b==1, legend(axL,{'channels','mean'},'FontSize',6,'Location','best'); end
            else; axis(axL,'off'); end
            if b==n_bands, xlabel(axL,'time from CF onset (s)','FontSize',7); end
            ylabel(axL,'ERD% vs baseline','FontSize',7);
            title(axL,band_lbl,'FontSize',8); grid(axL,'on');

            axR = subplot(n_bands,2,(b-1)*2+2);
            if ~any(isnan(data_b(:)))
                imagesc(axR, t_ax_cf, 1:n_sel_ch, data_b.');
                colormap(axR,cmap_rd); caxis(axR,[-clim_e,clim_e]);
                colorbar(axR,'FontSize',5);
                set(axR,'YTick',1:n_sel_ch,'YTickLabel',sel_ch,'FontSize',5.5);
                if b==n_bands, xlabel(axR,'time from CF onset (s)','FontSize',7); end
                xline(axR,0,'w--','LineWidth',.8);
            else; axis(axR,'off'); end
            title(axR,band_lbl,'FontSize',8);
        end
        sgtitle(fig_erd, sprintf('ERD/ERS [%%\\Delta vs 1s baseline]  class %s (%s)  —  %s  [%s]', ...
                cls_names{c},mod_label,basename,upper(paradigm)),'FontSize',9,'Interpreter','none');
        saveas(fig_erd, fullfile(res_dir,sprintf('%02d_erd_cls%s.svg',c+1,cls_names{c})),'svg');
    end

    %% ── Fig 4: Temporal topoplot evolution ───────────────────────────────────
    n_step   = max(1, round(0.5*framerate));   % columns every 0.5 s
    n_bl_cue = round(framerate);               % 1 s baseline before cue onset

    onset_chunks = nan(1,n_trials);
    for t = 1:n_trials
        oc = trials(t).onset_code; sc = trials(t).start_chunk;
        if isnan(oc), continue; end
        cand = find((header_chunks.EVENT.TYP==oc)&(header_chunks.EVENT.POS<sc),1,'last');
        if ~isempty(cand)
            onset_chunks(t) = header_chunks.EVENT.POS(cand);
        else
            onset_chunks(t) = sc - round(1.5*framerate);
        end
    end

    n_cf_win    = floor(n_cf_min / n_step);
    n_win_total = 1 + n_cf_win;
    win_lbl     = [{'Cue'}, arrayfun(@(w) sprintf('%.1f-%.1fs',(w-1)*.5,w*.5),1:n_cf_win,'un',0)];
    row_lbl     = [cls_names, {'diff cls1-cls2'}];
    n_row_topo  = min(3, n_cls + 1);

    % One figure per modality (MI and/or CVSA depending on paradigm)
    topo_mod_list = {};
    if use_mi
        ns_mi = numel(csp_mi.selected_channels); nb_mi = csp_mi.n_bands;
        fd_mi = permute(reshape(feat_pre_mi.', ns_mi, nb_mi, size(feat_pre_mi,1)), [3,1,2]);
        topo_mod_list{end+1} = struct('fd', fd_mi, 'sel_ch', {csp_mi.selected_channels}, ...
                                       'n_sel_ch', ns_mi, 'n_bands', nb_mi, 'label', 'MI');
    end
    if use_cvsa
        ns_cv = numel(csp_cvsa.selected_channels); nb_cv = csp_cvsa.n_bands;
        fd_cv = permute(reshape(feat_pre_cv.', ns_cv, nb_cv, size(feat_pre_cv,1)), [3,1,2]);
        topo_mod_list{end+1} = struct('fd', fd_cv, 'sel_ch', {csp_cvsa.selected_channels}, ...
                                       'n_sel_ch', ns_cv, 'n_bands', nb_cv, 'label', 'CVSA');
    end

    topo_mean = [];   % first modality result — saved in results struct below
    for mi_idx = 1:numel(topo_mod_list)
        m    = topo_mod_list{mi_idx};
        fd_m = m.fd;

        topo_sum_m = zeros(n_cls, n_win_total, m.n_sel_ch, m.n_bands);
        topo_cnt_m = zeros(n_cls,1);
        for t = 1:n_trials
            c = trials(t).target_class; if isnan(c), continue; end
            sc = trials(t).start_chunk; oc = onset_chunks(t);
            if isnan(oc), continue; end
            bl_r = max(1,oc-n_bl_cue):oc-1;
            bl_r = bl_r(bl_r>=1 & bl_r<=size(fd_m,1));
            if numel(bl_r)<2, continue; end
            pbl_m = mean(fd_m(bl_r,:,:),1);
            if any(isnan(pbl_m(:))), continue; end
            cue_r = oc:sc-1;
            cue_r = cue_r(cue_r>=1 & cue_r<=size(fd_m,1));
            if isempty(cue_r), continue; end
            pcue_m = mean(fd_m(cue_r,:,:),1);
            topo_sum_m(c,1,:,:) = topo_sum_m(c,1,:,:) + ...
                reshape((pcue_m-pbl_m)./max(pbl_m,eps)*100, 1, 1, m.n_sel_ch, m.n_bands);
            for w = 1:n_cf_win
                rs = sc+(w-1)*n_step; re = min(rs+n_step-1,size(fd_m,1));
                if rs>size(fd_m,1), break; end
                pw_m = mean(fd_m(rs:re,:,:),1);
                topo_sum_m(c,1+w,:,:) = topo_sum_m(c,1+w,:,:) + ...
                    reshape((pw_m-pbl_m)./max(pbl_m,eps)*100, 1, 1, m.n_sel_ch, m.n_bands);
            end
            topo_cnt_m(c) = topo_cnt_m(c)+1;
        end
        topo_mean_m = topo_sum_m;
        for c=1:n_cls
            if topo_cnt_m(c)>0, topo_mean_m(c,:,:,:) = topo_sum_m(c,:,:,:)/topo_cnt_m(c); end
        end
        topo_bnd_m = reshape(sum(topo_mean_m,4), n_cls, n_win_total, m.n_sel_ch);
        if isempty(topo_mean), topo_mean = topo_mean_m; end

        % Tight manual layout: no colorbar to save space
        fig4 = figure('Name',sprintf('Topo %s — %s',m.label,basename),'Color','w', ...
                      'NumberTitle','off', ...
                      'Position',[100 40 max(260*n_win_total,680) 280*n_row_topo+100]);
        lm4=0.03; rm4=0.01; tmf=0.11; bmf=0.04; hg4=0.003; vg4=0.018;
        cw4 = (1-lm4-rm4-(n_win_total-1)*hg4)/n_win_total;
        ch4 = (1-tmf-bmf-(n_row_topo-1)*vg4)/n_row_topo;

        for row = 1:n_row_topo
            bot4 = 1 - tmf - row*(ch4+vg4) + vg4;
            for col = 1:n_win_total
                lft4 = lm4 + (col-1)*(cw4+hg4);
                ax4  = axes('Parent',fig4,'Position',[lft4, bot4, cw4, ch4]); %#ok<LAXES>
                if row <= n_cls
                    if topo_cnt_m(row)==0, axis(ax4,'off'); continue; end
                    vals4 = squeeze(topo_bnd_m(row,col,:));
                else
                    if n_cls<2||any(topo_cnt_m(1:2)==0), axis(ax4,'off'); continue; end
                    vals4 = squeeze(topo_bnd_m(1,col,:)) - squeeze(topo_bnd_m(2,col,:));
                end
                row_str4=''; if col==1, row_str4=[row_lbl{row} '  ']; end
                topo_map(m.sel_ch(:), vals4(:), [], ax4, [row_str4 win_lbl{col}], col==n_win_total, true);
            end
        end
        sgtitle(fig4, sprintf('ERD/ERS topoplots (%s)  --  %s  [%s]  (blue=ERD  red=ERS)', ...
                m.label,basename,upper(paradigm)),'FontSize',9,'Interpreter','none');
        saveas(fig4, fullfile(res_dir,sprintf('%02d_topo_%s.svg',n_cls+1+mi_idx,m.label)),'svg');
    end
    if isempty(topo_mean), topo_mean = zeros(n_cls,n_win_total,n_sel_ch,n_bands); end

    %% ── Fig 5: FBCSP log-variance features + sLDA discriminant ──────────────
    % Common trial vectors (fused raw for hybrid) — kept for results struct
    trial_pass_v   = [trials.pass]';
    trial_cls_v    = [trials.target_class]';
    trial_mean_raw = nan(n_trials,n_cls);
    trial_max_int  = nan(n_trials,n_cls);
    trial_mean_erd = nan(n_trials,1);
    for t = 1:n_trials
        np  = trials(t).n_pre;
        raw = trials(t).raw(np+1:end,:);
        int = trials(t).integrated(np+1:end,:);
        vld = ~any(isnan(raw),2);
        if any(vld)
            trial_mean_raw(t,:) = mean(raw(vld,:),1);
            trial_max_int(t,:)  = max(int(vld,:),[],1);
        end
        sc   = trials(t).start_chunk;
        bl_r = max(1,sc-n_bl):sc-1;
        cf_r = sc:min(sc+n_cf_min-1,size(feat3d,1));
        bl_r = bl_r(bl_r>=1 & bl_r<=size(feat3d,1));
        if ~isempty(bl_r) && ~isempty(cf_r)
            pbl_t = mean(feat3d(bl_r,:,:),1);
            pcf_t = mean(feat3d(cf_r,:,:),1);
            if ~any(isnan(pbl_t(:)))&&~any(isnan(pcf_t(:)))
                trial_mean_erd(t) = mean((pcf_t-pbl_t)./max(pbl_t,eps)*100,'all');
            end
        end
    end

    % Build per-modality struct list for Fig 5 + Fig 6
    fig_mod_list = {};
    if use_mi
        ns5 = numel(csp_mi.selected_channels); nb5 = csp_mi.n_bands;
        fd5 = permute(reshape(feat_pre_mi.', ns5, nb5, size(feat_pre_mi,1)), [3,1,2]);
        src5 = 'raw'; if strcmp(paradigm,'hybrid'), src5 = 'p_mi'; end
        fig_mod_list{end+1} = struct('feat', feat_mi, 'feat3d', fd5, ...
            'csp', csp_mi, 'slda', slda_mi, 'sel_ch', {csp_mi.selected_channels}, ...
            'n_sel_ch', ns5, 'n_bands', nb5, 'n_comp', csp_mi.n_components, ...
            'label', 'MI', 'src', src5);
    end
    if use_cvsa
        ns5 = numel(csp_cvsa.selected_channels); nb5 = csp_cvsa.n_bands;
        fd5 = permute(reshape(feat_pre_cv.', ns5, nb5, size(feat_pre_cv,1)), [3,1,2]);
        src5 = 'raw'; if strcmp(paradigm,'hybrid'), src5 = 'p_cvsa'; end
        fig_mod_list{end+1} = struct('feat', feat_cv, 'feat3d', fd5, ...
            'csp', csp_cvsa, 'slda', slda_cvsa, 'sel_ch', {csp_cvsa.selected_channels}, ...
            'n_sel_ch', ns5, 'n_bands', nb5, 'n_comp', csp_cvsa.n_components, ...
            'label', 'CVSA', 'src', src5);
    end

    for fmi = 1:numel(fig_mod_list)
        fm = fig_mod_list{fmi};

        %% Fig 5: FBCSP log-variance features + sLDA discriminant
        n_feat_raw5 = fm.n_comp * fm.n_bands;
        feat_cf_mean5 = nan(n_cls, n_cf_min, n_feat_raw5);
        for c = 1:n_cls
            acc_f5 = zeros(n_cf_min,n_feat_raw5); cnt_f5 = 0;
            for t = 1:n_trials
                if isnan(trials(t).target_class)||trials(t).target_class~=c, continue; end
                sc    = trials(t).start_chunk;
                n_use = min(trials(t).n_cf, n_cf_min);
                idx   = sc:sc+n_use-1;
                idx   = idx(idx>=1 & idx<=size(fm.feat,1));
                if numel(idx)<n_use, continue; end
                v5 = log(max(fm.feat(idx,:),eps));
                if any(isnan(v5(:))), continue; end
                acc_f5(1:n_use,:) = acc_f5(1:n_use,:) + v5; cnt_f5 = cnt_f5+1;
            end
            if cnt_f5>0, feat_cf_mean5(c,:,:) = acc_f5/cnt_f5; end
        end

        n_bm5 = size(fm.slda.bands,1); n_feat_f5 = fm.slda.n_components * n_bm5;
        full_w5 = zeros(1,n_feat_f5);
        if isempty(fm.slda.selected_feature_indices) && numel(fm.slda.weights)==n_feat_f5
            full_w5 = fm.slda.weights;
        elseif ~isempty(fm.slda.selected_feature_indices)
            full_w5(fm.slda.selected_feature_indices) = fm.slda.weights;
        else
            nu5 = min(numel(fm.slda.weights),n_feat_f5);
            full_w5(1:nu5) = fm.slda.weights(1:nu5);
        end
        W_slda5 = reshape(full_w5, fm.slda.n_components, n_bm5).';
        ch_discr5 = zeros(1,fm.n_sel_ch);
        for bs = 1:n_bm5
            d5 = max(abs(fm.csp.bands - fm.slda.bands(bs,:)),[],2);
            bc5 = find(d5<1e-3,1);
            if isempty(bc5), continue; end
            ch_discr5 = ch_discr5 + W_slda5(bs,:) * fm.csp.csp_matrices{bc5};
        end

        p_mean_cls5 = nan(n_cls,n_cf_min);
        for c = 1:n_cls
            acc_p5=zeros(n_cf_min,1); cnt_p5=0;
            for t = 1:n_trials
                if isnan(trials(t).target_class)||trials(t).target_class~=c, continue; end
                np = trials(t).n_pre;
                rv5 = trials(t).(fm.src)(np+1:end,c);
                n_u5 = min(numel(rv5),n_cf_min);
                if any(isnan(rv5(1:n_u5))), continue; end
                acc_p5(1:n_u5) = acc_p5(1:n_u5)+rv5(1:n_u5); cnt_p5=cnt_p5+1;
            end
            if cnt_p5>0, p_mean_cls5(c,:) = acc_p5'/cnt_p5; end
        end

        n_col5 = n_cls+1;
        fig5 = figure('Name',sprintf('Features & sLDA %s — %s',fm.label,basename),...
                      'Color','w','NumberTitle','off','Position',[120 40 1100 600]);
        cols_c = lines(n_cls);
        for c = 1:n_cls
            ax5c = subplot(2, n_col5, c);
            fc5 = squeeze(feat_cf_mean5(c,:,:));
            if ~any(isnan(fc5(:)))
                bl_f5 = mean(fc5(1:min(3,n_cf_min),:),1);
                fc_r5 = fc5 - repmat(bl_f5,n_cf_min,1);
                imagesc(ax5c, t_ax_cf, 1:n_feat_raw5, fc_r5.');
                colormap(ax5c, rdbu_cmap()); caxis(ax5c, sym_clim(fc_r5(:)));
                colorbar(ax5c,'FontSize',5);
                for b = 0:fm.n_bands, yline(ax5c,b*fm.n_comp+0.5,'k-','LineWidth',.5); end
                set(ax5c,'YTick',ceil(fm.n_comp/2):fm.n_comp:n_feat_raw5, ...
                    'YTickLabel',arrayfun(@(b)sprintf('b%d',b),1:fm.n_bands,'un',0),'FontSize',7);
                xlabel(ax5c,'time (s)','FontSize',7);
            else; axis(ax5c,'off'); end
            title(ax5c,sprintf('Log-var features  cls %s',cls_names{c}),'FontSize',9);
        end
        ax_topo5 = subplot(2, n_col5, n_col5);
        topo_map(fm.sel_ch(:), ch_discr5(:), [], ax_topo5, ...
                 sprintf('sLDA discriminant\n(+=%s  -=%s)',cls_names{min(2,n_cls)},cls_names{1}), true, true);
        ax_prob5 = subplot(2, n_col5, n_col5+1 : 2*n_col5);
        hold(ax_prob5,'on');
        for c = 1:n_cls
            if ~any(isnan(p_mean_cls5(c,:)))
                plot(ax_prob5, t_ax_cf, p_mean_cls5(c,:), 'Color',cols_c(c,:), ...
                     'LineWidth',2,'DisplayName',sprintf('cls %s',cls_names{c}));
            end
        end
        yline(ax_prob5,0.5,'k--','FontSize',7);
        for c=1:n_cls, yline(ax_prob5,thresholds(c),':','Color',cols_c(c,:),'LineWidth',1); end
        legend(ax_prob5,'FontSize',8,'Location','best');
        xlabel(ax_prob5,'time from CF onset (s)','FontSize',8);
        ylabel(ax_prob5,'P(class)','FontSize',8);
        title(ax_prob5,sprintf('Mean raw sLDA output per class (%s)',fm.label),'FontSize',9);
        set(ax_prob5,'XLim',[t_ax_cf(1),t_ax_cf(end)],'YLim',[0,1]); grid(ax_prob5,'on');
        sgtitle(fig5,sprintf('FBCSP features & sLDA (%s) — %s [%s]',fm.label,basename,upper(paradigm)),...
                'FontSize',10,'Interpreter','none');
        saveas(fig5, fullfile(res_dir,sprintf('%02d_features_slda_%s.svg',n_cls+2+fmi,fm.label)),'svg');

        %% Fig 6: EEG-behavior correlation
        trial_mean_raw6 = nan(n_trials,n_cls);
        trial_mean_erd6 = nan(n_trials,1);
        for t = 1:n_trials
            np6 = trials(t).n_pre;
            raw6 = trials(t).(fm.src)(np6+1:end,:);
            vld6 = ~any(isnan(raw6),2);
            if any(vld6), trial_mean_raw6(t,:) = mean(raw6(vld6,:),1); end
            sc6  = trials(t).start_chunk;
            bl6  = max(1,sc6-n_bl):sc6-1;
            cf6  = sc6:min(sc6+n_cf_min-1,size(fm.feat3d,1));
            bl6  = bl6(bl6>=1 & bl6<=size(fm.feat3d,1));
            if ~isempty(bl6) && ~isempty(cf6)
                pb6 = mean(fm.feat3d(bl6,:,:),1); pc6 = mean(fm.feat3d(cf6,:,:),1);
                if ~any(isnan(pb6(:)))&&~any(isnan(pc6(:)))
                    trial_mean_erd6(t) = mean((pc6-pb6)./max(pb6,eps)*100,'all');
                end
            end
        end

        pass_col = [.15 .65 .15]; miss_col = [.85 .15 .15]; cols_c6 = lines(max(n_cls,2));
        fig6 = figure('Name',sprintf('EEG-Behavior %s — %s',fm.label,basename),...
                      'Color','w','NumberTitle','off','Position',[140 40 1000 500]);

        ax6 = subplot(2,2,1); bar_d6=nan(n_cls,2); err_d6=nan(n_cls,2);
        for c = 1:n_cls
            mh6 = trial_cls_v==c & trial_pass_v; mm6 = trial_cls_v==c & ~trial_pass_v;
            vh6 = trial_mean_raw6(mh6,c); vm6 = trial_mean_raw6(mm6,c);
            if ~isempty(vh6), bar_d6(c,1)=mean(vh6); err_d6(c,1)=std(vh6); end
            if ~isempty(vm6), bar_d6(c,2)=mean(vm6); err_d6(c,2)=std(vm6); end
        end
        bh6 = bar(ax6, bar_d6,'grouped');
        bh6(1).FaceColor=pass_col; if numel(bh6)>1, bh6(2).FaceColor=miss_col; end
        hold(ax6,'on');
        for c=1:n_cls
            for j=1:2
                if ~isnan(err_d6(c,j))
                    xp6 = c + (j==1)*(-0.15) + (j==2)*(0.15);
                    errorbar(ax6,xp6,bar_d6(c,j),err_d6(c,j),'k.','LineWidth',1.2);
                end
            end
        end
        yline(ax6,0.5,'k--');
        set(ax6,'XTick',1:n_cls,'XTickLabel',cls_names,'YLim',[0,1]); grid on;
        legend(ax6,{'hit','miss'},'FontSize',7);
        ylabel(ax6,sprintf('mean sLDA %s (target class)',fm.label));
        title(ax6,sprintf('Classifier confidence: hit vs miss (%s)',fm.label),'FontSize',9);

        ax6 = subplot(2,2,2); hold(ax6,'on');
        for c = 1:n_cls
            mc6 = trial_cls_v==c & ~isnan(trial_cls_v);
            xv6 = trial_mean_erd6(mc6); yv6 = trial_mean_raw6(mc6,c); pv6 = trial_pass_v(mc6);
            scatter(ax6,xv6(pv6), yv6(pv6), 45, pass_col,'filled','DisplayName',sprintf('%s hit',cls_names{c}));
            scatter(ax6,xv6(~pv6),yv6(~pv6),45, miss_col,'s','DisplayName',sprintf('%s miss',cls_names{c}));
        end
        xlabel(ax6,sprintf('mean ERD%% (%s channels)',fm.label));
        ylabel(ax6,sprintf('mean sLDA %s output',fm.label));
        title(ax6,sprintf('ERD vs classifier output (%s)',fm.label),'FontSize',9);
        legend(ax6,'FontSize',7,'Location','best'); grid(ax6,'on'); yline(ax6,0.5,'k--');

        ax6 = subplot(2,2,3); hold(ax6,'on');
        for c = 1:n_cls
            int_h6=zeros(n_cf_min,1); nh6=0; int_m6=zeros(n_cf_min,1); nm6=0;
            for t = 1:n_trials
                if isnan(trial_cls_v(t))||trial_cls_v(t)~=c, continue; end
                np6=trials(t).n_pre; iv6=trials(t).integrated(np6+1:end,c);
                n_u6=min(numel(iv6),n_cf_min);
                if any(isnan(iv6(1:n_u6))), continue; end
                if trials(t).pass, int_h6(1:n_u6)=int_h6(1:n_u6)+iv6(1:n_u6); nh6=nh6+1;
                else,              int_m6(1:n_u6)=int_m6(1:n_u6)+iv6(1:n_u6); nm6=nm6+1; end
            end
            if nh6>0, plot(ax6,t_ax_cf,int_h6/nh6,'Color',pass_col,'LineWidth',2,...
                          'DisplayName',sprintf('%s hit (n=%d)',cls_names{c},nh6)); end
            if nm6>0, plot(ax6,t_ax_cf,int_m6/nm6,'Color',miss_col,'LineWidth',2,'LineStyle','--',...
                          'DisplayName',sprintf('%s miss (n=%d)',cls_names{c},nm6)); end
            if isfinite(thresholds(c)), yline(ax6,thresholds(c),':','Color',cols_c6(c,:),'LineWidth',1); end
        end
        yline(ax6,p_rest,'k--');
        xlabel(ax6,'time from CF onset (s)'); ylabel(ax6,'integrator value');
        int_lbl = 'Mean integrator: hit vs miss';
        if strcmp(paradigm,'hybrid'), int_lbl = [int_lbl ' (fused)']; end
        title(ax6,int_lbl,'FontSize',9);
        legend(ax6,'FontSize',7,'Location','best'); grid(ax6,'on');
        set(ax6,'XLim',[t_ax_cf(1),t_ax_cf(end)],'YLim',[0,1]);

        ax6 = subplot(2,2,4); hold(ax6,'on');
        for t = 1:n_trials
            c = trial_cls_v(t); if isnan(c)||c>n_cls, continue; end
            col6 = pass_col*trial_pass_v(t) + miss_col*(1-trial_pass_v(t));
            scatter(ax6,t,trial_mean_raw6(t,c),28,col6,'filled');
        end
        xlabel(ax6,'trial index'); ylabel(ax6,sprintf('mean sLDA %s (target class)',fm.label));
        title(ax6,'Per-trial output (green=hit  red=miss)','FontSize',9);
        yline(ax6,0.5,'k--'); grid(ax6,'on');

        sgtitle(fig6,sprintf('EEG-Behavior [%s]  %s  [%s]',fm.label,basename,upper(paradigm)),...
                'FontSize',10,'Interpreter','none');
        saveas(fig6, fullfile(res_dir,sprintf('%02d_eeg_behavior_%s.svg',n_cls+2+numel(fig_mod_list)+fmi,fm.label)),'svg');
    end

    %% ── Save mat ─────────────────────────────────────────────────────────────
    results = struct( ...
        'subject',               subject_from_basename(basename), ...
        'basename',              basename, ...
        'gdf_path',              gdf_path, ...
        'paradigm',              paradigm, ...
        'classes',               classes, ...
        'n_trials',              n_trials, ...
        'trial_acc',             trial_acc, ...
        'trial_acc_events',      trial_acc_events, ...
        'trial_acc_no_timeout',  trial_acc_no_timeout, ...
        'trial_acc_no_rej',      trial_acc_no_rej, ...
        'n_hits_events',         n_hits_ev, ...
        'n_misses_events',       n_misses_ev, ...
        'n_timeout_events',      n_timeout_ev, ...
        'tth_per_class',         tth_per_class, ...
        't_miss_peak_cls',       t_miss_peak_cls, ...
        'miss_max_fused',        miss_max_fused, ...
        'miss_mean_mi',          miss_mean_mi, ...
        'miss_mean_cvsa',        miss_mean_cvsa, ...
        'sample_acc_fused',      sample_acc_fused, ...
        'sample_acc_mi_cls',     sample_acc_mi_cls, ...
        'sample_acc_cv_cls',     sample_acc_cv_cls, ...
        'sample_acc_total',      sample_acc_total, ...
        'sample_acc_mean_trial', sample_acc_mean_trial, ...
        'hit_rate_cls',          hit_rate_cls, ...
        'art_rate',              art_rate, ...
        'erd_mean',              erd_mean, ...
        'erd_static',            erd_static, ...
        'topo_mean',             topo_mean, ...
        't_ax_cf',               t_ax_cf, ...
        'selected_channels',     {sel_ch}, ...
        'bands',                 bands, ...
        'framerate',             framerate, ...
        'chunk_size',            chunk_size, ...
        'trial_mean_raw',        trial_mean_raw, ...
        'trial_max_int',         trial_max_int, ...
        'trial_mean_erd',        trial_mean_erd, ...
        'trial_pass_v',          trial_pass_v, ...
        'trial_cls_v',           trial_cls_v ...
    );
    mat_path = fullfile(gdf_dir, sprintf('eval_single_%s_%s.mat', paradigm, basename));
    save(mat_path, 'results', 'trials', 'int_cfg');
    fprintf('  saved  %s\n         %s\n', res_dir, mat_path);
end

% ── Local helpers ──────────────────────────────────────────────────────────────
function s = subject_from_basename(bn)
    p = strsplit(bn,'.'); if numel(p)>=1, s=p{1}; else, s=bn; end
end

function cm = rdbu_cmap(n)
    if nargin<1, n=256; end
    h = floor(n/2);
    r = [linspace(0.17,1,h), linspace(1,0.70,h)];
    g = [linspace(0.51,1,h), linspace(1,0.09,h)];
    b = [linspace(0.73,1,h), linspace(1,0.07,h)];
    cm = [r(:), g(:), b(:)];
end

function cl = sym_clim(vals)
    v = max(abs(vals(isfinite(vals))));
    if isempty(v)||v==0, v=1; end
    cl = [-v, v];
end
