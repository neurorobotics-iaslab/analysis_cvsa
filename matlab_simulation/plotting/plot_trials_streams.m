function fig = plot_trials_streams(trials_hyb, trials_mi, trials_cvsa, int_cfg, framerate, ...
                                     basename, trial_outcome_real, oc, show_figure)
% PLOT_TRIALS_STREAMS  One panel per 781 trial, comparing the leaky-integrator
%   "control signal" P(class 1) for three counterfactual input streams that
%   all share the same buffer/threshold parameters (int_cfg):
%     - Hybrid   : trials_hyb  — fused MI+CVSA (as actually run)
%     - MI-only  : trials_mi   — same buffer fed with raw MI sLDA output
%     - CVSA-only: trials_cvsa — same buffer fed with raw CVSA sLDA output
%
%   Horizontal lines mark the thresholds (upper = thr(c1), lower = 1-thr(c2))
%   and p_rest. Light-grey shading = artifact gate (shared across streams,
%   since the artifact flag does not depend on paradigm).
%
%   Background colour = REAL outcome of the hybrid recording (HIT/MISS/TIMEOUT
%   from GDF events 897/898/899). The title shows the *simulated* outcome and
%   time-to-event for each of the three streams.
%
%   trial_outcome_real : [n_trials] real GDF outcome codes (897/898/899/0)
%   oc                 : struct with fields .hyb/.mi/.cvsa, each [n_trials x 2]
%                        = [outcome_code, t_event_s]; outcome_code:
%                        0=unknown, 1=HIT, 2=MISS, 3=TIMEOUT

    if nargin < 7 || isempty(trial_outcome_real)
        trial_outcome_real = zeros(1, numel(trials_hyb));
    end
    if nargin < 9 || isempty(show_figure)
        show_figure = true;
    end
    HIT_EV = 897;  MISS_EV = 898;  TO_EV = 899;
    OC_STR = {'?', 'HIT', 'MISS', 'TO'};   % index = outcome_code + 1

    if isempty(trials_hyb)
        log_step('plot_trials_streams: no trials -> nothing to plot');
        fig = gobjects(0);
        return;
    end

    n_trials = numel(trials_hyb);
    n_cols   = min(4, n_trials);
    n_rows   = ceil(n_trials / n_cols);

    classes    = to_vec(int_cfg.classes);
    init_val   = to_vec(int_cfg.init_val);
    p_rest     = init_val(1);
    thresholds = to_vec(int_cfg.thresholds);
    thr_up     = thresholds(1);
    thr_dn     = 1 - thresholds(2);

    COL_HYB = [0.05 0.35 0.85];
    COL_MI  = [0.85 0.40 0.10];
    COL_CV  = [0.49 0.18 0.56];

    fig_name = sprintf('Hybrid vs MI-only vs CVSA-only [control signal] %s', basename);
    if show_figure, vis = 'on'; else, vis = 'off'; end
    fig = figure('Name', fig_name, 'Color', 'w', 'NumberTitle', 'off', 'Visible', vis);
    set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(fig, 'InvertHardcopy', 'off');   % preserve per-trial axes background colours on export

    for t = 1:n_trials
        tr_h = trials_hyb(t); tr_m = trials_mi(t); tr_c = trials_cvsa(t);
        ax = subplot(n_rows, n_cols, t); hold(ax, 'on'); grid(ax, 'on');

        nframes  = size(tr_h.integrated, 1);
        n_pre    = tr_h.n_pre;
        time_s   = ((0:nframes-1) - n_pre) / framerate;
        cf_end_t = (tr_h.n_cf - 1) / framerate;

        if nframes < 2
            % Recording ended before this trial's CF started (n_cf == 0 in
            % integrate_signal) -- nothing to plot, avoid the degenerate
            % xlim([x,x]) call below.
            text(ax, 0.5, 0.5, sprintf('Trial %d\nno CF data\n(recording ended early)', t), ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 8);
            set(ax, 'XTick', [], 'YTick', []);
            continue;
        end

        % --- artifact shading (shared across streams) ---
        if any(tr_h.artifact)
            yl = [-0.05, 1.10];
            in_seg = false; seg_t0 = 0;
            for k = 1:nframes
                if tr_h.artifact(k) && ~in_seg
                    seg_t0 = time_s(k); in_seg = true;
                elseif (~tr_h.artifact(k) || k == nframes) && in_seg
                    patch(ax, [seg_t0, time_s(k), time_s(k), seg_t0], ...
                              [yl(1), yl(1), yl(2), yl(2)], ...
                              [0.92 0.92 0.92], 'EdgeColor', 'none', ...
                              'HandleVisibility', 'off');
                    in_seg = false;
                end
            end
        end

        % --- threshold + reference lines ---
        yline(ax, thr_up,  'k--', 'HandleVisibility', 'off');
        yline(ax, thr_dn,  'k--', 'HandleVisibility', 'off');
        yline(ax, p_rest,  ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
        xline(ax, 0,        ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
        xline(ax, cf_end_t, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');

        % --- the three control-signal traces (integrator buffer, P(c1)) ---
        plot(ax, time_s, tr_h.integrated(:, 1), '-', 'LineWidth', 2.0, 'Color', COL_HYB);
        plot(ax, time_s, tr_m.integrated(:, 1), '-', 'LineWidth', 1.6, 'Color', COL_MI);
        plot(ax, time_s, tr_c.integrated(:, 1), '-', 'LineWidth', 1.6, 'Color', COL_CV);

        ylim(ax, [-0.05, 1.10]);
        xlim(ax, [min(time_s), max(time_s)]);
        xlabel(ax, 't [s]');
        ylabel(ax, sprintf('P(class %d)', classes(1)));

        if isnan(tr_h.target_class)
            title(ax, sprintf('Trial %d  —  no onset event', t), 'FontSize', 8);
        else
            real_ev = trial_outcome_real(t);
            switch real_ev
                case HIT_EV,  bg_col = [0.88 1.00 0.88];
                case MISS_EV, bg_col = [1.00 0.88 0.88];
                case TO_EV,   bg_col = [1.00 0.96 0.82];
                otherwise,    bg_col = [0.96 0.96 0.96];
            end
            set(ax, 'Color', bg_col);

            str_h = sprintf('HYB:%s(%.1fs)',  OC_STR{oc.hyb(t,1)+1},  oc.hyb(t,2));
            str_m = sprintf('MI:%s(%.1fs)',   OC_STR{oc.mi(t,1)+1},   oc.mi(t,2));
            str_c = sprintf('CVSA:%s(%.1fs)', OC_STR{oc.cvsa(t,1)+1}, oc.cvsa(t,2));
            title(ax, sprintf('Trial %d  |  class %d (%d)\n%s   %s   %s', ...
                              t, tr_h.target_class, classes(tr_h.target_class), str_h, str_m, str_c), ...
                  'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'none');
        end

        if t == 1
            legend(ax, {'Hybrid (fused)', 'MI-only', 'CVSA-only'}, ...
                   'Location', 'best', 'Box', 'off', 'FontSize', 7);
        end
    end

    sgtitle(sprintf('%s — control signal P(c1): Hybrid vs MI-only vs CVSA-only (same buffer/thresholds, background = real hybrid outcome)', ...
                    basename), 'Interpreter', 'none');
end
