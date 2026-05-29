function plot_trials(trials, int_cfg, framerate, paradigm, basename)
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
%   Title: PASS if integrated/raw[target_class] >= thresholds[target_class] within CF window
%          (matches Training.cpp evaluation mode: is_target_hit checks raw >= threshold).

    if isempty(trials)
        log_step('plot_trials: no trials -> nothing to plot');
        return;
    end

    n_trials   = numel(trials);
    n_cols     = min(4, n_trials);
    n_rows     = ceil(n_trials / n_cols);

    classes    = to_vec(int_cfg.classes);
    init_val   = to_vec(int_cfg.init_val);
    p_rest     = init_val(1);
    thresholds = to_vec(int_cfg.thresholds);
    thr_up     = thresholds(1);            % class-1 win threshold in P(c1) space
    thr_dn     = 1 - thresholds(2);        % class-2 win threshold in P(c1) space

    fig_name = sprintf('Continuous feedback [%s] %s', paradigm, basename);
    figure('Name', fig_name, 'Color', 'w', 'NumberTitle', 'off');
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
            title(ax, sprintf('trial %d  (no onset)', t));
        else
            title(ax, sprintf('trial %d  target=c%d (%d)  %s', ...
                              t, tr.target_class, classes(tr.target_class), ...
                              ternary(tr.pass, 'PASS', 'miss')));
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

    sgtitle(sprintf('%s — continuous-feedback trials  (%d PASS / %d trials)', ...
                    paradigm, sum([trials.pass]), n_trials));
end


function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
