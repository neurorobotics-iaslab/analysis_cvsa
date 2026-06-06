function [art_flags, info] = detect_artifacts(signal, header, art_cfg, cfg)
% DETECT_ARTIFACTS  Faithful port of artifacts_bci, identical to the
%   validated reference in src/artifacts_bci/test/test_artifacts.m
%   (MAE vs ROS ~ 1e-7).
%
%   Per chunk:
%     1) CAR on the chunk using non-EOG channels
%     2) EOG path:   chunk_car -> LP(freq_low_EOG) -> HP(freq_high_EOG)
%     3) peaks path: chunk_car -> HP(freq_high_peaks)
%     4) push both into circular ring buffers of size bufsize
%     5) when the ring is full:
%          hEOG = max(abs(buf_eog(:,c1) - buf_eog(:,c2)))
%          vEOG = if >=3 EOG refs: max(abs((c1+c2)/2 - c3)),
%                 else max(abs((c1+c2)/2))
%          peaks: max(abs(buf_peaks(:, non_eog)))
%          has_artifact = (hEOG > th_hEOG) || (vEOG > th_vEOG)
%                                          || (peaks > th_peaks)
%
% Inputs
%   signal   [N x C]
%   header   .Label
%   art_cfg  params.ArtifactCfg.params (yamlmatlab struct)
%   cfg      .samplerate, .chunk_size, .bufsize_artifact

    [N, C] = size(signal);
    chunk_size = cfg.chunk_size;
    bufsize    = cfg.bufsize_artifact;
    n_chunks   = floor(N / chunk_size);
    fs         = cfg.samplerate;
    nyq        = fs / 2;

    EOG_ch  = resolve_channels(art_cfg.EOG_ch_names, header.Label);
    non_eog = setdiff(1:C, EOG_ch);
    if numel(EOG_ch) < 2
        error('detect_artifacts:eog', 'Need at least 2 EOG-reference channels.');
    end

    th_hEOG  = double(art_cfg.th_hEOG);
    th_vEOG  = double(art_cfg.th_vEOG);
    th_peaks = double(art_cfg.th_peaks);

    [b_lp,  a_lp]  = butter(art_cfg.filterOrder_EOG,   art_cfg.freq_low_EOG    / nyq, 'low');
    [b_hp,  a_hp]  = butter(art_cfg.filterOrder_EOG,   art_cfg.freq_high_EOG   / nyq, 'high');
    [b_hpk, a_hpk] = butter(art_cfg.filterOrder_peaks, art_cfg.freq_high_peaks / nyq, 'high');
    zi_lp  = []; % zeros(max(length(a_lp),  length(b_lp))  - 1, C);
    zi_hp  = []; % zeros(max(length(a_hp),  length(b_hp))  - 1, C);
    zi_hpk = []; % zeros(max(length(a_hpk), length(b_hpk)) - 1, C);

    buf_eog   = nan(bufsize, C);
    buf_peaks = nan(bufsize, C);

    art_flags = false(n_chunks, 1);
    first_valid_chunk = ceil(bufsize / chunk_size);

    log_step('detect_artifacts: %d chunks, bufsize=%d -> first valid chunk = %d', ...
             n_chunks, bufsize, first_valid_chunk);

    for k = 1:n_chunks
        s0 = (k - 1) * chunk_size + 1;
        s1 = k * chunk_size;
        chunk = signal(s0:s1, :);

        % CAR on non-EOG channels (matches validated test)
        chunk_car = chunk - mean(chunk(:, non_eog), 2);

        % EOG: LP -> HP
        [eog_lp, zi_lp] = filter(b_lp, a_lp, chunk_car, zi_lp);
        [eog_bp, zi_hp] = filter(b_hp, a_hp, eog_lp,    zi_hp);
        % peaks: HP
        [pks_hp, zi_hpk] = filter(b_hpk, a_hpk, chunk_car, zi_hpk, 1);

        % circular ring buffers
        buf_eog(:,:) = [buf_eog(chunk_size+1:end,:); eog_bp];
        buf_peaks(:,:) = [buf_peaks(chunk_size+1:end,:); pks_hp];

        if any(isnan(buf_eog(:))), continue; end

        % EOG check
        heog = buf_eog(:, EOG_ch(1)) - buf_eog(:, EOG_ch(2));
        if numel(EOG_ch) >= 3
            veog = (buf_eog(:, EOG_ch(1)) + buf_eog(:, EOG_ch(2))) / 2 - buf_eog(:, EOG_ch(3));
        else
            veog = (buf_eog(:, EOG_ch(1)) + buf_eog(:, EOG_ch(2))) / 2;
        end
        has_art = (max(abs(heog)) > th_hEOG) || (max(abs(veog)) > th_vEOG);

        % Peaks check (non-EOG channels only)
        if max(abs(buf_peaks(:, non_eog)), [], 'all') > th_peaks
            has_art = true;
        end
        art_flags(k) = has_art;

        if mod(k, 500) == 0
            log_step('detect_artifacts: chunk %d/%d, artifacts so far = %d', ...
                     k, n_chunks, sum(art_flags));
        end
    end

    info = struct('first_valid_chunk', first_valid_chunk, 'n_chunks', n_chunks);
    log_step('detect_artifacts: %d/%d chunks flagged (%.2f%%)', ...
             sum(art_flags), n_chunks, 100*sum(art_flags)/max(1,n_chunks));
end
