%% MAIN_BATCH_EVALUATE  Multi-file offline evaluation with metrics and visualisation.
%
%   Replays the full BCI pipeline (same as main_simulate) over one or more
%   GDF recordings, computes performance metrics per file and their average,
%   and optionally plots all CF trials grouped by class.
%
%   Results are saved to a .mat file for cross-paradigm comparison with
%   main_compare_paradigms.m.
%
%   USAGE
%     Set show_all_trials and ARTIFACT_REJECT_THR below, then run.
%     A multi-select dialog picks the GDFs; each must have a sibling YAML.

clear; clc; close all;

% ── Configuration ─────────────────────────────────────────────────────────────
show_all_trials     = true;   % plot all CF trials grouped by class (one figure per class)
ARTIFACT_REJECT_THR = 0.30;  % CF frames fraction with artifact -> trial "rejected" for no-reject metric

% ── Paths ─────────────────────────────────────────────────────────────────────
this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
addpath(fullfile(this_dir, 'io'));
addpath(fullfile(this_dir, 'processing'));
addpath(fullfile(this_dir, 'artifacts'));
addpath(fullfile(this_dir, 'classifier'));
addpath(fullfile(this_dir, 'integrator'));
addpath(fullfile(this_dir, 'plotting'));
addpath(fullfile(this_dir, 'utils'));

% ── File selection (multi-select) ─────────────────────────────────────────────
default_dir = '/home/paolo/bci_vr_ws/recordings';
if ~isfolder(default_dir), default_dir = pwd; end

[gdf_names, gdf_dir] = uigetfile({'*.gdf','GDF recordings (*.gdf)'}, ...
    'Select one or more GDF recordings (Ctrl/Cmd = multi-select)', ...
    default_dir, 'MultiSelect', 'on');
if isequal(gdf_names, 0)
    error('main_batch_evaluate:cancel', 'No files selected.');
end
if ischar(gdf_names), gdf_names = {gdf_names}; end   % single file -> cell
n_files   = numel(gdf_names);
gdf_paths = cellfun(@(n) fullfile(gdf_dir, n), gdf_names, 'UniformOutput', false);

fprintf('\n=== BATCH EVALUATION: %d file(s) ===\n\n', n_files);

% ── Per-file processing ────────────────────────────────────────────────────────
per_file(n_files) = struct();

for f = 1:n_files
    gdf_path = gdf_paths{f};
    fprintf('--- File %d/%d: %s ---\n', f, n_files, gdf_names{f});

    % 1) Load GDF + companion YAML
    [signal, header, basename] = load_gdf(gdf_path);
    [params, ~]                = load_params_yaml(gdf_path);

    paradigm   = params.integrator.paradigm;
    fs         = double(params.acquisition.samplerate);
    framerate  = double(params.acquisition.framerate);
    chunk_size = round(fs / framerate);
    if abs(fs - header.SampleRate) > 1e-3
        log_step('main_batch_evaluate: YAML fs=%.1f != GDF fs=%.1f -> using GDF', ...
                 fs, header.SampleRate);
        fs         = header.SampleRate;
        chunk_size = round(fs / framerate);
    end

    bufsize_proc = double(params.RingBufferCfg.params.size);
    bufsize_art  = double(params.RingBufferCfgArtifact.params.size);
    eog_names    = to_strcell(params.CarCfg.params.EOG_ch_names);

    do_car_mi   = true;
    do_car_cvsa = true;
    if isfield(params, 'processing_fbcsp_mi'),   do_car_mi   = logical(params.processing_fbcsp_mi.do_car);   end
    if isfield(params, 'processing_fbcsp_cvsa'), do_car_cvsa = logical(params.processing_fbcsp_cvsa.do_car); end

    use_mi   = ismember(paradigm, {'mi','hybrid'});
    use_cvsa = ismember(paradigm, {'cvsa','hybrid'});

    % 2) CSP + sLDA models
    csp_mi = []; slda_mi = []; csp_cvsa = []; slda_cvsa = [];
    if use_mi,   csp_mi   = load_csp(params,'mi');   slda_mi   = load_slda(params,'mi');   end
    if use_cvsa, csp_cvsa = load_csp(params,'cvsa'); slda_cvsa = load_slda(params,'cvsa'); end

    % 3) Processing streams
    features_mi = []; header_mi = []; features_cv = []; header_cv = []; info_proc = [];
    if use_mi
        proc_cfg = struct('samplerate',fs,'chunk_size',chunk_size,'bufsize',bufsize_proc, ...
                          'filter_order',4,'do_car',do_car_mi,'eog_names',{eog_names});
        [features_mi, header_mi, info_proc] = apply_processing(signal, header, csp_mi, proc_cfg);
    end
    if use_cvsa
        proc_cfg = struct('samplerate',fs,'chunk_size',chunk_size,'bufsize',bufsize_proc, ...
                          'filter_order',4,'do_car',do_car_cvsa,'eog_names',{eog_names});
        [features_cv, header_cv, info2] = apply_processing(signal, header, csp_cvsa, proc_cfg);
        if isempty(info_proc), info_proc = info2; end
    end

    % 4) Artifact detection
    art_cfg = params.ArtifactCfg.params;
    art_cfg.EOG_ch_names = to_strcell(art_cfg.EOG_ch_names);
    cfg_art = struct('samplerate',fs,'chunk_size',chunk_size,'bufsize_artifact',bufsize_art);
    [art_flags, ~] = detect_artifacts(signal, header, art_cfg, cfg_art);

    % 5) sLDA classification
    p_mi_a = []; if use_mi,   p_mi_a = apply_slda(features_mi, slda_mi,   csp_mi.bands);   end
    p_cv_a = []; if use_cvsa, p_cv_a = apply_slda(features_cv, slda_cvsa, csp_cvsa.bands); end

    % 6) Trial integrator
    int_cfg = params.integrator;
    if ~isfield(int_cfg,'increment'),            int_cfg.increment            = 1;   end
    if ~isfield(int_cfg,'thresholds_rejection'), int_cfg.thresholds_rejection = [];  end
    if ~isfield(int_cfg,'cvsa_influence'),       int_cfg.cvsa_influence       = 2.5; end
    if use_mi, hdr_ch = header_mi; else, hdr_ch = header_cv; end
    hdr_ch.framerate = framerate;
    trials = integrate_signal(p_mi_a, p_cv_a, art_flags, hdr_ch, int_cfg, paradigm);

    % 7) Compute metrics
    metrics = compute_metrics(trials, int_cfg, framerate, ARTIFACT_REJECT_THR);

    % 8) Store
    per_file(f).basename  = basename;
    per_file(f).gdf_path  = gdf_path;
    per_file(f).paradigm  = paradigm;
    per_file(f).classes   = to_vec(int_cfg.classes);
    per_file(f).trials    = trials;
    per_file(f).metrics   = metrics;
    per_file(f).int_cfg   = int_cfg;
    per_file(f).framerate = framerate;

    print_file_report(metrics, basename, paradigm);

    % 9) Optional per-class trial visualisation
    if show_all_trials
        plot_trials_by_class(trials, int_cfg, framerate, paradigm, basename);
    end
end

% ── Aggregate across files ─────────────────────────────────────────────────────
paradigm_label = per_file(1).paradigm;
aggregate      = compute_aggregate(per_file);

print_aggregate_report(aggregate, paradigm_label, n_files);

% ── Save ──────────────────────────────────────────────────────────────────────
[~, first_base] = fileparts(gdf_paths{1});
save_path = fullfile(gdf_dir, sprintf('eval_%s_%s.mat', paradigm_label, first_base));
save(save_path, 'per_file', 'aggregate', 'paradigm_label', 'gdf_paths');
fprintf('Results saved to:\n  %s\n\n', save_path);


%% ═══════════════════════════════════════════════════════════════════════════
%   LOCAL FUNCTIONS
%% ═══════════════════════════════════════════════════════════════════════════

% ── compute_metrics ───────────────────────────────────────────────────────────
function m = compute_metrics(trials, int_cfg, framerate, art_reject_thr)
% Returns a struct of scalar and per-class metrics for one session.
%
%   Metrics computed:
%     trial_acc            overall hit rate (%)
%     trial_acc_no_reject  hit rate excluding high-artifact trials (%)
%     mean_tth             mean time-to-hit over HIT trials (s)
%     mean_sample_acc      mean frame-level classification accuracy (%)
%     mean_confidence      mean P(target class) over CF frames (%)
%     mean_art_rate        mean fraction of CF frames with artifact (%)
%     mean_peak_norm       mean maximum normalised target probability
%     n_trials / n_valid / n_rejected
%
%   Per-class (m.per_class(c)):
%     hit_rate / recall    = TP / (TP+FN)  (%)
%     precision            = TP / (TP+FP)  (%)
%       where FP = HIT trials whose target was a *different* class
%     mean_tth, sample_acc, confidence, peak_norm

    n_trials  = numel(trials);
    classes   = to_vec(int_cfg.classes);
    n_cls     = numel(classes);

    is_hit         = false(n_trials, 1);
    tth            = NaN(n_trials, 1);
    sample_acc     = NaN(n_trials, 1);
    confidence     = NaN(n_trials, 1);
    art_rate       = NaN(n_trials, 1);
    peak_norm      = NaN(n_trials, 1);
    is_rejected    = false(n_trials, 1);
    target_cls_arr = NaN(n_trials, 1);

    for t = 1:n_trials
        tr       = trials(t);
        cf_range = (tr.n_pre + 1) : size(tr.normalized, 1);
        if isempty(cf_range), continue; end

        is_hit(t)         = tr.pass;
        target_cls_arr(t) = tr.target_class;
        art_rate(t)       = mean(tr.artifact(cf_range));
        is_rejected(t)    = art_rate(t) >= art_reject_thr;

        if isnan(tr.target_class), continue; end

        % Time-to-hit: first frame where normalized(target) >= 1.0
        if tr.pass
            hit_f = find(tr.normalized(cf_range, tr.target_class) >= 1.0, 1);
            if ~isempty(hit_f), tth(t) = (hit_f - 1) / framerate; end
        end

        % Sample accuracy: argmax(raw P) per frame vs target class
        valid = ~any(isnan(tr.raw(cf_range,:)), 2);
        if any(valid)
            [~, pred]    = max(tr.raw(cf_range(valid), :), [], 2);
            sample_acc(t) = mean(pred == tr.target_class) * 100;
        end

        % Confidence: mean P(target class) over valid CF frames
        p_tgt = tr.raw(cf_range, tr.target_class);
        confidence(t) = mean(p_tgt(~isnan(p_tgt))) * 100;

        % Peak normalised value for target class
        peak_norm(t) = max(tr.normalized(cf_range, tr.target_class));
    end

    valid_t = ~isnan(target_cls_arr);
    n_valid = sum(valid_t);

    m.n_trials              = n_trials;
    m.n_valid               = n_valid;
    m.n_rejected            = sum(is_rejected & valid_t);
    m.trial_acc             = safe_mean(is_hit(valid_t)) * 100;
    m.mean_tth              = mean(tth(is_hit & valid_t), 'omitnan');
    m.mean_sample_acc       = mean(sample_acc(valid_t),   'omitnan');
    m.mean_confidence       = mean(confidence(valid_t),   'omitnan');
    m.mean_art_rate         = mean(art_rate(valid_t),     'omitnan') * 100;
    m.mean_peak_norm        = mean(peak_norm(valid_t),    'omitnan');

    non_rej                 = valid_t & ~is_rejected;
    m.trial_acc_no_reject   = safe_mean(is_hit(non_rej)) * 100;

    % Per-class breakdown
    m.per_class = struct([]);
    for c = 1:n_cls
        cls_mask = valid_t & (target_cls_arr == c);
        oth_mask = valid_t & (target_cls_arr ~= c);
        n_c      = sum(cls_mask);
        tp       = sum(is_hit(cls_mask));
        fn       = sum(~is_hit(cls_mask));
        fp       = sum(is_hit(oth_mask));

        m.per_class(c).class_code  = classes(c);
        m.per_class(c).n_trials    = n_c;
        m.per_class(c).hit_rate    = safe_div(tp, tp+fn) * 100;  % = recall
        m.per_class(c).recall      = safe_div(tp, tp+fn) * 100;
        m.per_class(c).precision   = safe_div(tp, tp+fp) * 100;
        m.per_class(c).mean_tth    = mean(tth(cls_mask & is_hit), 'omitnan');
        m.per_class(c).sample_acc  = mean(sample_acc(cls_mask),   'omitnan');
        m.per_class(c).confidence  = mean(confidence(cls_mask),   'omitnan');
        m.per_class(c).peak_norm   = mean(peak_norm(cls_mask),    'omitnan');
    end

    % Raw arrays (used by compute_aggregate)
    m.raw.is_hit       = is_hit;
    m.raw.tth          = tth;
    m.raw.sample_acc   = sample_acc;
    m.raw.confidence   = confidence;
    m.raw.art_rate     = art_rate;
    m.raw.peak_norm    = peak_norm;
    m.raw.is_rejected  = is_rejected;
    m.raw.target_cls   = target_cls_arr;
end

% ── compute_aggregate ─────────────────────────────────────────────────────────
function agg = compute_aggregate(per_file)
    n_files   = numel(per_file);
    flds      = {'trial_acc','trial_acc_no_reject','mean_tth','mean_sample_acc', ...
                 'mean_confidence','mean_art_rate','mean_peak_norm'};
    agg       = struct();
    for fi = 1:numel(flds)
        fld  = flds{fi};
        vals = arrayfun(@(pf) pf.metrics.(fld), per_file);
        agg.(fld).mean = mean(vals, 'omitnan');
        agg.(fld).std  = std(vals, 0, 'omitnan');
        agg.(fld).vals = vals;
    end

    n_cls    = numel(per_file(1).metrics.per_class);
    pc_flds  = {'hit_rate','recall','precision','mean_tth','sample_acc','confidence','peak_norm'};
    agg.per_class = struct([]);
    for c = 1:n_cls
        agg.per_class(c).class_code = per_file(1).metrics.per_class(c).class_code;
        for fi2 = 1:numel(pc_flds)
            fld  = pc_flds{fi2};
            vals = arrayfun(@(pf) pf.metrics.per_class(c).(fld), per_file);
            agg.per_class(c).(fld).mean = mean(vals, 'omitnan');
            agg.per_class(c).(fld).std  = std(vals, 0, 'omitnan');
            agg.per_class(c).(fld).vals = vals;
        end
    end
    agg.n_files = n_files;
end

% ── print_file_report ─────────────────────────────────────────────────────────
function print_file_report(m, basename, paradigm)
    fprintf('\n  [%s | %s]\n', upper(paradigm), basename);
    fprintf('  Trials: %d  |  Valid: %d  |  Rejected (art >= %.0f%%): %d\n', ...
            m.n_trials, m.n_valid, 30, m.n_rejected);
    fprintf('  Trial accuracy         : %5.1f%%\n', m.trial_acc);
    fprintf('  Trial acc (no-reject)  : %5.1f%%\n', m.trial_acc_no_reject);
    fprintf('  Mean time-to-hit       : %5.2f s\n', m.mean_tth);
    fprintf('  Sample accuracy        : %5.1f%%\n', m.mean_sample_acc);
    fprintf('  Mean sLDA confidence   : %5.1f%%\n', m.mean_confidence);
    fprintf('  Artifact rate          : %5.1f%%\n', m.mean_art_rate);
    fprintf('  Mean peak norm target  : %5.3f\n',   m.mean_peak_norm);
    for c = 1:numel(m.per_class)
        pc = m.per_class(c);
        tth_str = 'N/A';
        if ~isnan(pc.mean_tth), tth_str = sprintf('%.2f s', pc.mean_tth); end
        fprintf('  Class %d (code %d)  n=%d  hit=%.0f%%  prec=%.0f%%  recall=%.0f%%  tth=%s\n', ...
                c, pc.class_code, pc.n_trials, pc.hit_rate, pc.precision, pc.recall, tth_str);
    end
end

% ── print_aggregate_report ────────────────────────────────────────────────────
function print_aggregate_report(agg, paradigm, n_files)
    fprintf('\n%s\n', repmat('=', 1, 68));
    fprintf('  AGGREGATE  [%s]  %d file(s)  (mean +/- std)\n', upper(paradigm), n_files);
    fprintf('%s\n', repmat('=', 1, 68));
    flds   = {'trial_acc','trial_acc_no_reject','mean_tth','mean_sample_acc', ...
              'mean_confidence','mean_art_rate','mean_peak_norm'};
    labels = {'Trial accuracy (%)', 'Trial acc no-reject (%)', 'Mean TTH (s)', ...
              'Sample accuracy (%)', 'Mean confidence (%)', 'Artifact rate (%)', ...
              'Mean peak norm'};
    for i = 1:numel(flds)
        v = agg.(flds{i});
        fprintf('  %-26s : %6.2f +/- %.2f\n', labels{i}, v.mean, v.std);
    end
    for c = 1:numel(agg.per_class)
        pc = agg.per_class(c);
        fprintf('  Class %d (code %d)  hit=%.1f+/-%.1f%%  prec=%.1f+/-%.1f%%  tth=%.2f+/-%.2fs\n', ...
                c, pc.class_code, ...
                pc.hit_rate.mean, pc.hit_rate.std, ...
                pc.precision.mean, pc.precision.std, ...
                pc.mean_tth.mean, pc.mean_tth.std);
    end
    fprintf('%s\n\n', repmat('=', 1, 68));
end

% ── plot_trials_by_class ──────────────────────────────────────────────────────
function plot_trials_by_class(trials, int_cfg, framerate, paradigm, basename)
% Two figures: one per class. Each subplot = one trial targeting that class.
% Same P(c1) view as plot_trials.m.

    classes    = to_vec(int_cfg.classes);
    init_val   = to_vec(int_cfg.init_val);
    p_rest     = init_val(1);
    thresholds = to_vec(int_cfg.thresholds);
    thr_up     = thresholds(1);
    thr_dn     = 1 - thresholds(2);

    for c = 1:numel(classes)
        cls_code = classes(c);
        mask     = [trials.onset_code] == cls_code;
        cls_tr   = trials(mask);
        if isempty(cls_tr), continue; end

        n_tr  = numel(cls_tr);
        n_col = min(4, n_tr);
        n_row = ceil(n_tr / n_col);

        fig_name = sprintf('[%s] %s — class %d (code %d)', ...
                           upper(paradigm), basename, c, cls_code);
        figure('Name', fig_name, 'Color', 'w', 'NumberTitle', 'off');

        for t = 1:n_tr
            tr      = cls_tr(t);
            ax      = subplot(n_row, n_col, t);
            hold(ax, 'on'); grid(ax, 'on');
            nframes  = size(tr.integrated, 1);
            time_s   = ((0:nframes-1) - tr.n_pre) / framerate;
            cf_end_t = (tr.n_cf - 1) / framerate;

            % Artifact shading
            if any(tr.artifact)
                yl = [-0.05, 1.10];
                in_seg = false; seg_t0 = 0;
                for k = 1:nframes
                    if tr.artifact(k) && ~in_seg
                        seg_t0 = time_s(k); in_seg = true;
                    elseif (~tr.artifact(k) || k == nframes) && in_seg
                        patch(ax, [seg_t0 time_s(k) time_s(k) seg_t0], ...
                              [yl(1) yl(1) yl(2) yl(2)], [.92 .92 .92], ...
                              'EdgeColor','none','HandleVisibility','off');
                        in_seg = false;
                    end
                end
            end

            yline(ax, thr_up, 'k--', 'HandleVisibility','off');
            yline(ax, thr_dn, 'k--', 'HandleVisibility','off');
            yline(ax, p_rest, ':',   'Color',[.4 .4 .4],'HandleVisibility','off');
            xline(ax, 0,        ':', 'Color',[.4 .4 .4],'HandleVisibility','off');
            xline(ax, cf_end_t, ':', 'Color',[.4 .4 .4],'HandleVisibility','off');

            scatter(ax, time_s, tr.raw(:,1), 18, [0.85 0.40 0.10], 'filled', ...
                    'MarkerFaceAlpha',.55,'MarkerEdgeColor','none');
            plot(ax, time_s, tr.integrated(:,1),   '-',  'LineWidth', 2.0, ...
                 'Color', [0.05 0.35 0.85]);
            plot(ax, time_s, tr.normalized_pc1,    '--', 'LineWidth', 1.6, ...
                 'Color', [0.10 0.65 0.30]);

            ylim(ax, [-0.05, 1.10]);
            xlim(ax, [min(time_s), max(time_s)]);
            xlabel(ax, 't [s]');
            ylabel(ax, sprintf('P(c%d)', classes(1)));
            result_str  = ternary(tr.pass, 'PASS', 'miss');
            result_clr  = ternary(tr.pass, [0 .5 0], [.7 0 0]);
            title(ax, sprintf('trial %d  %s', t, result_str), ...
                  'Color', result_clr, 'FontSize', 8);

            if t == 1
                legend(ax, {'sLDA P(c1)','integrated','norm'}, ...
                       'Location','best','Box','off','FontSize',7);
            end
        end

        sgtitle(sprintf('[%s] %s — class %d (code %d)  %d/%d PASS', ...
                        upper(paradigm), basename, c, cls_code, ...
                        sum([cls_tr.pass]), n_tr), 'FontSize', 11);
    end
end

% ── Tiny helpers ──────────────────────────────────────────────────────────────
function v = safe_mean(x)
    if isempty(x), v = NaN; else, v = mean(x, 'omitnan'); end
end
function v = safe_div(a, b)
    if b == 0, v = NaN; else, v = a / b; end
end
function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
