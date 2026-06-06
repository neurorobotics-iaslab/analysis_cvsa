function [features, header_out, info, features_pre] = apply_processing(signal, header, csp, cfg)
% APPLY_PROCESSING  FBCSP pipeline, identical to the validated reference in
%   src/slda_bci/test/test_slda.m (GDF mode, MAE vs ROS < 1e-6):
%
%     raw chunk -> CAR (mean of non-EOG channels) -> for each band:
%       chunk_car -> LP(hi) -> HP(lo) stateful causal Butterworth ->
%       push into NaN-init ring buffer (full nchannels x bufsize x nbands).
%     When the buffer is full (no NaN), CSP per band on the selected
%     channels, then mean power = sum(x^2)/bufsize per CSP component.
%     Features are flattened column-major: [comp1_band1..compN_band1,
%     comp1_band2..]. EVENT.POS/.DUR are rescaled to the chunk timeline.
%
% Inputs
%   signal  [N x C]    raw EEG
%   header  struct     uses .Label, .EVENT.POS/.DUR
%   csp     struct     load_csp output
%   cfg     struct     .samplerate, .chunk_size, .bufsize, .filter_order,
%                      .do_car, .eog_names
%
% Outputs
%   features    [n_chunks x (n_components*n_bands)]  raw mean-power features
%               (NaN before the ringbuffer is full)
%   header_out  header with EVENT.POS/.DUR in chunk-timeline units
%   info        .first_valid_chunk, .n_chunks, .chunk_size, .bufsize,
%               .feature_dim, .n_bands, .n_components
%   features_pre  [n_chunks x (n_sel_ch*n_bands)]  pre-CSP per-channel mean
%               power, band-major layout [b0_ch0, b0_ch1, ..., b1_ch0, ...]
%               NaN before ring buffer is full. Optional 4th output.

    [N, C] = size(signal);
    chunk_size = cfg.chunk_size;
    n_chunks   = floor(N / chunk_size);

    n_bands  = csp.n_bands;
    n_comp   = csp.n_components;
    feat_dim = n_comp * n_bands;
    features = NaN(n_chunks, feat_dim);

    csp_ch   = resolve_channels(csp.selected_channels, header.Label);
    n_sel_ch = numel(csp_ch);
    features_pre = NaN(n_chunks, n_sel_ch * n_bands);
    if cfg.do_car
        eog_idx = resolve_channels(cfg.eog_names, header.Label);
    else
        eog_idx = [];
    end
    non_eog = setdiff(1:C, eog_idx);

    fs    = cfg.samplerate;
    order = cfg.filter_order;
    nyq   = fs / 2;
    b_lp = cell(n_bands, 1); a_lp = cell(n_bands, 1);
    b_hp = cell(n_bands, 1); a_hp = cell(n_bands, 1);
    zi_lp = cell(n_bands, 1); zi_hp = cell(n_bands, 1);
    for b = 1:n_bands
        lo = csp.bands(b, 1);
        hi = csp.bands(b, 2);
        [b_lp{b}, a_lp{b}] = butter(order, hi / nyq, 'low');  
        [b_hp{b}, a_hp{b}] = butter(order, lo / nyq, 'high');
        zi_lp{b} = []; % zeros(max(numel(a_lp{b}), numel(b_lp{b})) - 1, C);
        zi_hp{b} = []; % zeros(max(numel(a_hp{b}), numel(b_hp{b})) - 1, C);
    end

    bufsize = cfg.bufsize;
    bufs    = NaN(bufsize, C, n_bands);

    first_valid_chunk = ceil(bufsize / chunk_size);
    log_step('apply_processing: %d chunks, chunk_size=%d, bufsize=%d -> first valid chunk = %d', ...
             n_chunks, chunk_size, bufsize, first_valid_chunk);

    for k = 1:n_chunks
        s0 = (k - 1) * chunk_size + 1;
        s1 = k * chunk_size;
        chunk = signal(s0:s1, :);

        if cfg.do_car && ~isempty(non_eog)
            chunk_car = chunk - mean(chunk(:, non_eog), 2);
        else
            chunk_car = chunk;
        end

        for b = 1:n_bands
            [lp_out, zi_lp{b}] = filter(b_lp{b}, a_lp{b}, chunk_car, zi_lp{b});
            [bp_out, zi_hp{b}] = filter(b_hp{b}, a_hp{b}, lp_out,    zi_hp{b});
            bufs(:,:,b) = [bufs(chunk_size+1:end,:,b); bp_out];
        end

        if any(isnan(bufs(:))), continue; end

        csp_feats     = zeros(n_comp,    n_bands);
        pre_csp_feats = zeros(n_sel_ch,  n_bands);
        for b = 1:n_bands
            buf_sel            = bufs(:, csp_ch, b);
            csp_out            = buf_sel * csp.csp_matrices{b}.';
            csp_feats(:, b)    = sum(csp_out  .^ 2, 1).' / bufsize;
            pre_csp_feats(:,b) = sum(buf_sel  .^ 2, 1).' / bufsize;
        end
        features(k, :)     = reshape(csp_feats,     1, []);
        features_pre(k, :) = reshape(pre_csp_feats, 1, []);  % band-major: rows=ch cols=band, col-major → [ch1_b1..chM_b1, ch1_b2..]

        if mod(k, 500) == 0
            log_step('apply_processing: chunk %d/%d', k, n_chunks);
        end
    end

    header_out = header;
    if ~isempty(header.EVENT.POS)
        header_out.EVENT.POS = ceil(header.EVENT.POS / chunk_size);
        header_out.EVENT.DUR = ceil(header.EVENT.DUR / chunk_size);
    end

    n_valid = sum(~any(isnan(features), 2));
    info = struct('first_valid_chunk', first_valid_chunk, ...
                  'n_chunks', n_chunks, 'chunk_size', chunk_size, ...
                  'bufsize', bufsize, 'feature_dim', feat_dim, ...
                  'n_bands', n_bands, 'n_components', n_comp);
    log_step('apply_processing: produced %d feature vectors (dim=%d)', n_valid, feat_dim);
end
