function fig = plot_trials(trials, int_cfg, framerate, paradigm, basename, trial_outcome_real, show_figure)
% PLOT_TRIALS  One panel per 781 trial. All curves are in the P(class 1)
%   reference frame (= the first class returned by sLDA, which is the only
%   probability the integrator consumes). Each panel shows:
%
%     - scatter of the raw sLDA / fused-sLDA P(c1) per frame
%     - line of integrated raw P(c1)
%     - line of integrated normalized P(c1)
%     - horizontal threshold lines:
%         upper = thr(c1)              (class 1 wins when raw integrated >= upper)
%         lower = 1 - thr(c2)          (class 2 wins when raw integrated <= lower)
%     - dotted line at p_rest = init_val(1)
%     - light-grey shading where the artifact gate fired
%
%   trial_outcome_real (optional): vector of actual GDF outcomes per trial
%     897=HIT, 898=MISS, 899=TIMEOUT, 0=unknown.
%     When provided, titles show REAL outcome; simulated pass/fail shown in brackets.
%
%   Title: REAL outcome (HIT/MISS/TIMEOUT) from GDF + [sim PASS/miss] in brackets.

    if nargin < 6 || isempty(trial_outcome_real)
        trial_outcome_real = zeros(1, numel(trials));
    end
    if nargin < 7 || isempty(show_figure)
        show_figure = true;
    end
    HIT_EV = 897;  MISS_EV = 898;  TO_EV = 899;

    if isempty(trials)
        log_step('plot_trials: no trials -> nothing to plot');
        fig = gobjects(0);
        return;
    end

    n_trials   = numel(trials);
    n_cols     = min(4, n_trials);
    n_rows     = ceil(n_trials / n_cols);

    classes    = to_vec(int_cfg.classes);
    init_val   = to_vec(int_cfg.init_val);
    p_rest     = init_val(1);
    thresholds = to_vec(int_cfg.thresholds);
    thr_up     = thresholds(1);
    thr_dn     = 1 - thresholds(2);

    % Count real outcomes for the super-title
    n_hit_r  = sum(trial_outcome_real == HIT_EV);
    n_miss_r = sum(trial_outcome_real == MISS_EV);
    n_to_r   = sum(trial_outcome_real == TO_EV);
    has_real = any(trial_outcome_real ~= 0);

    fig_name = sprintf('Continuous feedback [%s] %s', paradigm, basename);
    if show_figure, vis = 'on'; else, vis = 'off'; end
    fig = figure('Name', fig_name, 'Color', 'w', 'NumberTitle', 'off', 'Visible', vis);
    set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(fig, 'InvertHardcopy', 'off');   % preserve per-trial axes background colours (HIT/MISS/TIMEOUT) on export
    log_step('plot_trials: %d trials in a %dx%d grid (P(c1) view, thr_up=%.2f, thr_dn=%.2f)', ...
             n_trials, n_rows, n_cols, thr_up, thr_dn);

    for t = 1:n_trials
        tr      = trials(t);
        ax      = subplot(n_rows, n_cols, t); hold(ax, 'on'); grid(ax, 'on');
        nframes = size(tr.integrated, 1);
        n_pre   = tr.n_pre;
        % t = 0 is the first CF chunk (POS); the reset publish sits at
        % t = -1/framerate; the last CF chunk (POS + DUR) sits at
        % t = (n_cf - 1) / framerate.
        time_s   = ((0:nframes-1) - n_pre) / framerate;
        cf_end_t = (tr.n_cf - 1) / framerate;

        if nframes < 2
            % Recording ended before this trial's CF started (n_cf == 0 in
            % integrate_signal) -- nothing to plot, avoid the degenerate
            % xlim([x,x]) call below.
            text(ax, 0.5, 0.5, sprintf('Trial %d\nno CF data\n(recording ended early)', t), ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 8);
            set(ax, 'XTick', [], 'YTick', []);
            continue;
        end

        % --- artifact shading ---
        if any(tr.artifact)
            yl = [-0.05, 1.10];
            in_seg   = false;
            seg_t0   = 0;
            for k = 1:nframes
                if tr.artifact(k) && ~in_seg
                    seg_t0 = time_s(k); in_seg = true;
                elseif (~tr.artifact(k) || k == nframes) && in_seg
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
        % CF-window boundaries: dotted vertical lines at t=0 (reset/CF start)
        % and t = cf_end_t (last CF chunk). Anything outside is the pre/post
        % visualization context.
        xline(ax, 0,        ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
        xline(ax, cf_end_t, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');

        % --- traces (all P(class 1) only) ---
        is_hybrid = strcmp(paradigm, 'hybrid');
        has_sub = is_hybrid && isfield(tr, 'p_mi') && ~all(isnan(tr.p_mi(:, 1)));
        if has_sub
            scatter(ax, time_s, tr.p_mi(:, 1),   12, [0.20 0.60 0.90], 'filled', ...
                    'MarkerFaceAlpha', 0.50, 'MarkerEdgeColor', 'none');
            scatter(ax, time_s, tr.p_cvsa(:, 1), 12, [0.85 0.20 0.70], 'filled', ...
                    'MarkerFaceAlpha', 0.50, 'MarkerEdgeColor', 'none');
        end
        scatter(ax, time_s, tr.raw(:, 1),       18, [0.85 0.40 0.10], 'filled', ...
                'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none');
        plot   (ax, time_s, tr.integrated(:, 1),   '-',  'LineWidth', 2.0, ...
                'Color', [0.05 0.35 0.85]);
        plot   (ax, time_s, tr.normalized_pc1,     '--', 'LineWidth', 1.6, ...
                'Color', [0.10 0.65 0.30]);

        ylim(ax, [-0.05, 1.10]);
        xlim(ax, [min(time_s), max(time_s)]);
        xlabel(ax, 't [s]');
        ylabel(ax, sprintf('P(class %d)', classes(1)));

        if isnan(tr.target_class)
            title(ax, sprintf('Trial %d  —  no onset event', t), 'FontSize',8);
        else
            real_ev = trial_outcome_real(t);
            if real_ev == HIT_EV
                outcome_str = 'HIT';    bg_col = [0.88 1.00 0.88];
            elseif real_ev == MISS_EV
                outcome_str = 'MISS';   bg_col = [1.00 0.88 0.88];
            elseif real_ev == TO_EV
                outcome_str = 'TIMEOUT'; bg_col = [1.00 0.96 0.82];
            else
                outcome_str = '???';    bg_col = [0.96 0.96 0.96];
            end
            set(ax,'Color',bg_col);
            title(ax, sprintf('Trial %d  |  class %d (%d)  |  %s', ...
                              t, tr.target_class, classes(tr.target_class), outcome_str), ...
                  'FontSize', 8, 'FontWeight','bold', 'Interpreter','none');
        end

        if t == 1
            if has_sub
                legend(ax, {'MI P(c1)', 'CVSA P(c1)', 'fused P(c1)', 'integrated raw', 'integrated norm'}, ...
                       'Location', 'best', 'Box', 'off');
            else
                legend(ax, {'sLDA P(c1)', 'integrated raw', 'integrated norm'}, ...
                       'Location', 'best', 'Box', 'off');
            end
        end
    end

    if has_real
        sgtitle(sprintf('%s — %s  |  REAL: %d HIT  %d MISS  %d TIMEOUT  (sim: %d PASS / %d)', ...
                        paradigm, basename, n_hit_r, n_miss_r, n_to_r, ...
                        sum([trials.pass]), n_trials), 'Interpreter','none');
    else
        sgtitle(sprintf('%s — %s  |  sim: %d PASS / %d trials', ...
                        paradigm, basename, sum([trials.pass]), n_trials), 'Interpreter','none');
    end
end


function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
