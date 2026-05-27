function trials = integrate_signal(p_mi, p_cvsa, art_flags, header_chunks, int_cfg, paradigm)
% INTEGRATE_SIGNAL  Per-trial binary leaky integrator. Matches the validated
%   reference in src/test_pipeline/src/test_full_pipeline.m
%   (init_integrator_state + step_integrator + compute_fusion + normalize_probs).
%
%   For every event 781 in header_chunks (POS/DUR in chunk units after
%   apply_processing rescaled the header), we run a fresh integrator state
%   for DUR chunks:
%     - Frame 1 (= reset frame): output = init_val, no step. This is what
%       the ROS resetIntegrator() publishes at the 781 chunk.
%     - Frames 2..DUR: leaky binary integrator
%           p_new = p_prev + sign(p_in - 0.5) * step
%           step  = min(|p_max - 0.5| * 2 * k_gain, 1) / buffer_size
%       artifact-flagged chunks freeze p_prev (p_new = p_prev).
%     - Hybrid: per-frame Bayesian fusion of MI and CVSA, with the CVSA
%       prior decaying via alpha(t) over half_life = 2.5 s.
%
% Output: trials struct array, one entry per 781 event:
%   .start_chunk, .dur_chunks   trial extent (chunk units)
%   .onset_code                 onset event code preceding 781
%   .target_class               1-based index into classes
%   .raw         [DUR x 2]      per-frame fused/raw sLDA P(c)
%   .integrated  [DUR x 2]      buffer state (= init_val at frame 1)
%   .normalized  [DUR x 2]      per-class linear stretch into [0, 1]
%   .artifact    [DUR x 1]      logical
%   .pass        bool           target hit normalized >= 1.0 within trial

    CF_CODE   = 781;
    HALF_LIFE = 2.5;    % default; overridden by int_cfg.cvsa_influence if present
    if isfield(int_cfg, 'cvsa_influence')
        HALF_LIFE = double(int_cfg.cvsa_influence);
    end

    classes   = to_vec(int_cfg.classes);
    n_cls     = numel(classes);
    init_val  = to_vec(int_cfg.init_val);
    if isempty(init_val), init_val = repmat(1/n_cls, 1, n_cls); end
    p_rest    = init_val(1);                 % binary scalar (= P(class 1) prior)
    thresholds = to_vec(int_cfg.thresholds);
    bsize     = double(int_cfg.buffer_size);
    k_gain    = double(int_cfg.k_gain);

    framerate = header_chunks.framerate;     % set by caller (main_simulate)

    POS = header_chunks.EVENT.POS;
    TYP = header_chunks.EVENT.TYP;
    DUR = header_chunks.EVENT.DUR;
    n_chunks = max([size(p_mi, 1), size(p_cvsa, 1), numel(art_flags)]);

    cf_idx = find(TYP == CF_CODE);
    if isempty(cf_idx)
        error('integrate_signal:no_cf', 'No event 781 in the GDF.');
    end
    log_step('integrate_signal: %d CF trials (paradigm=%s)', numel(cf_idx), paradigm);

    trials = struct([]);

    N_PRE = 1;   % one frame BEFORE the CF: the reset publish (p_rest).
                 % No post-CF frame: integration runs strictly inside
                 % [POS, POS+DUR] (inclusive both ends).

    for t = 1:numel(cf_idx)
        i_cf        = cf_idx(t);
        start_chunk = POS(i_cf);
        dur_offset  = max(0, DUR(i_cf));              % offset to the LAST CF chunk
        end_chunk   = min(n_chunks, start_chunk + dur_offset);
        n_cf        = end_chunk - start_chunk + 1;    % number of CF chunks = DUR+1

        onset_code   = find_onset_before(TYP, i_cf, classes);
        target_class = find(classes == onset_code, 1);
        if isempty(target_class), target_class = NaN; end

        % Visualisation frames = N_PRE (the reset publish) + n_cf (the CF
        % chunks). The CF spans [POS, POS+DUR] inclusive; we run the
        % integrator for those n_cf chunks. The last chunk may not exist
        % if the recording ends mid-trial; n_cf already accounts for that.
        n_total     = N_PRE + n_cf;
        raw_trial   = NaN(n_total, 2);
        integrated  = zeros(n_total, 2);
        art_trial   = false(n_total, 1);

        % --- N_PRE reset frame(s): the value ROS publishes when 781 fires,
        %     before any integration. No input consumed.
        for k = 1:N_PRE
            integrated(k, :) = [p_rest, 1 - p_rest];
            % raw_trial(k, :) stays NaN, art_trial(k) stays false
        end

        % Reset state (matches ROS resetIntegrator() before the next chunk).
        p_prev      = p_rest;
        frame_count = 0;

        for j = 1:n_cf
            k = j + N_PRE;
            c = start_chunk + (j - 1);

            art = false;
            if c <= numel(art_flags), art = art_flags(c); end
            art_trial(k) = art;

            % --- per-frame raw / fused sLDA input ------------------------
            switch paradigm
                case 'mi'
                    p_in = p_mi(c, 1);
                    raw_trial(k, :) = p_mi(c, :);
                case 'cvsa'
                    p_in = p_cvsa(c, 1);
                    raw_trial(k, :) = p_cvsa(c, :);
                case 'hybrid'
                    if any(isnan(p_mi(c, :))) || any(isnan(p_cvsa(c, :)))
                        p_in = NaN;
                        raw_trial(k, :) = [NaN, NaN];
                    else
                        t_sec  = frame_count / framerate;
                        p_fus  = bayesian_fuse(p_mi(c, :), p_cvsa(c, :), t_sec, HALF_LIFE);
                        p_in   = p_fus(1);
                        raw_trial(k, :) = p_fus;
                    end
                otherwise
                    error('integrate_signal:paradigm', 'Unknown paradigm "%s".', paradigm);
            end

            frame_count = frame_count + 1;

            % --- leaky binary integrator step (matches step_integrator) ---
            if ~art && ~isnan(p_in)
                p_max = max(p_in, 1 - p_in);
                vel   = min(abs(p_max - 0.5) * 2 * k_gain, 1);
                step  = vel / bsize;
                p_prev = max(0, min(1, p_prev + sign(p_in - 0.5) * step));
            end
            integrated(k, :) = [p_prev, 1 - p_prev];
        end

        % --- per-class linear stretch normalization (matches ROS
        %     Integrator::normalize_input, applied independently to each class).
        normalized = integrated;
        for c = 1:2
            thr = thresholds(c);
            if thr > p_rest
                slope = (1 - p_rest) / (thr - p_rest);
                normalized(:, c) = max(0, min(1, p_rest + (integrated(:, c) - p_rest) * slope));
            end
        end

        % --- P(c1) "view" of the normalization. ROS normalizes the two
        %     classes independently (normalize_input is per-element with the
        %     class's own threshold). When plotting only P(c1), we still want
        %     BOTH thresholds to shape the curve: above p_rest the curve uses
        %     thr(c1) (heading toward a class-1 decision); below p_rest it
        %     uses thr(c2) (heading toward a class-2 decision, viewed in
        %     P(c1) space as 1 - normalized_c2). The formula is identical to
        %     normalize_input applied symmetrically around p_rest.
        normalized_pc1 = integrated(:, 1);
        upper = integrated(:, 1) >= p_rest;
        normalized_pc1( upper) = normalized( upper, 1);
        normalized_pc1(~upper) = 1 - normalized(~upper, 2);

        % PASS check is only over the CF window (frames N_PRE+1 .. n_total),
        % not the reset frame.
        cf_range = (N_PRE + 1):n_total;
        pass = false;
        if ~isnan(target_class) && ~isempty(cf_range)
            pass = any(normalized(cf_range, target_class) >= 1.0);
        end

        trials(end+1).start_chunk    = start_chunk; %#ok<AGROW>
        trials(end).n_cf             = n_cf;             % number of CF chunks (POS..POS+DUR)
        trials(end).n_pre            = N_PRE;            % reset frames before CF
        trials(end).onset_code       = onset_code;
        trials(end).target_class     = target_class;
        trials(end).raw              = raw_trial;
        trials(end).integrated       = integrated;
        trials(end).normalized       = normalized;       % per-class (ROS style)
        trials(end).normalized_pc1   = normalized_pc1;   % P(c1) view (both thresholds)
        trials(end).artifact         = art_trial;
        trials(end).pass             = pass;

        if ~isnan(target_class)
            log_step('integrate_signal: trial %d/%d  onset=%d  target=c%d  %s', ...
                     t, numel(cf_idx), onset_code, target_class, ternary(pass, 'PASS', 'miss'));
        else
            log_step('integrate_signal: trial %d/%d  onset=? (no onset class before 781)', ...
                     t, numel(cf_idx));
        end
    end
end


function code = find_onset_before(TYP, i_cf, classes)
% Most-recent onset-class event strictly before the i_cf-th event.
    code = NaN;
    for j = i_cf-1:-1:1
        if ismember(TYP(j), classes)
            code = TYP(j);
            return;
        end
    end
end


function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
