function S = load_hybrid_streams(gdf_path)
% LOAD_HYBRID_STREAMS  Shared GDF->sLDA-stream loader for hybrid-only analyses.
%   Loads a GDF + companion YAML, runs FBCSP + artifact detection + sLDA for
%   both MI and CVSA, and returns everything main_hybrid_advantage_probs.m and
%   main_hybrid_advantage_integ.m need downstream. If the file's paradigm is
%   not 'hybrid', S.skip = true and the caller should skip to the next file.
%
%   Fields: skip, signal, header, basename, paradigm, fs, framerate,
%   chunk_size, csp_mi, csp_cvsa, slda_mi, slda_cvsa, p_mi_aligned,
%   p_cvsa_aligned, art_flags, int_cfg, header_chunks.

    S = struct('skip', false);
    [S.signal, S.header, S.basename] = load_gdf(gdf_path);
    [params, ~] = load_params_yaml(gdf_path);

    S.paradigm = params.integrator.paradigm;
    if ~strcmp(S.paradigm, 'hybrid')
        S.skip = true;
        return;
    end

    S.fs        = double(params.acquisition.samplerate);
    S.framerate = double(params.acquisition.framerate);
    S.chunk_size = round(S.fs / S.framerate);
    if abs(S.fs - S.header.SampleRate) > 1e-3
        log_step('load_hybrid_streams: YAML samplerate=%.1f != GDF samplerate=%.1f -> using GDF', ...
                 S.fs, S.header.SampleRate);
        S.fs = S.header.SampleRate;
        S.chunk_size = round(S.fs / S.framerate);
    end

    bufsize_proc = double(params.RingBufferCfg.params.size);
    bufsize_art  = double(params.RingBufferCfgArtifact.params.size);
    eog_names    = to_strcell(params.CarCfg.params.EOG_ch_names);

    do_car_mi   = logical(params.processing_fbcsp_mi.do_car);
    do_car_cvsa = logical(params.processing_fbcsp_cvsa.do_car);
    nchannels   = double(params.processing_fbcsp_mi.nchannels);

    S.signal = S.signal(:, 1:nchannels);

    log_step('load_hybrid_streams: paradigm=%s, fs=%g, framerate=%g, chunk=%d, bufproc=%d, bufart=%d, nchannels=%d', ...
             S.paradigm, S.fs, S.framerate, S.chunk_size, bufsize_proc, bufsize_art, nchannels);

    S.csp_mi    = load_csp(params, 'mi');
    S.slda_mi   = load_slda(params, 'mi');
    S.csp_cvsa  = load_csp(params, 'cvsa');
    S.slda_cvsa = load_slda(params, 'cvsa');

    proc_cfg_mi = struct('samplerate', S.fs, 'chunk_size', S.chunk_size, ...
                         'bufsize', bufsize_proc, 'filter_order', 4, ...
                         'do_car', do_car_mi, 'eog_names', {eog_names});
    [features_mi, header_mi, info_proc] = apply_processing(S.signal, S.header, S.csp_mi, proc_cfg_mi);

    proc_cfg_cvsa = struct('samplerate', S.fs, 'chunk_size', S.chunk_size, ...
                           'bufsize', bufsize_proc, 'filter_order', 4, ...
                           'do_car', do_car_cvsa, 'eog_names', {eog_names});
    [features_cv, ~, ~] = apply_processing(S.signal, S.header, S.csp_cvsa, proc_cfg_cvsa);

    art_cfg = params.ArtifactCfg.params;
    art_cfg.EOG_ch_names = to_strcell(art_cfg.EOG_ch_names);
    cfg_art = struct('samplerate', S.fs, 'chunk_size', S.chunk_size, 'bufsize_artifact', bufsize_art);
    [S.art_flags, ~] = detect_artifacts(S.signal, S.header, art_cfg, cfg_art);

    S.p_mi_aligned   = apply_slda(features_mi, S.slda_mi,   S.csp_mi.bands);
    S.p_cvsa_aligned = apply_slda(features_cv, S.slda_cvsa, S.csp_cvsa.bands);

    log_step('load_hybrid_streams: streams aligned (n_chunks=%d, first_valid_proc=%d)', ...
             info_proc.n_chunks, info_proc.first_valid_chunk);

    S.int_cfg = params.integrator;
    if ~isfield(S.int_cfg, 'increment'),            S.int_cfg.increment = 1; end
    if ~isfield(S.int_cfg, 'thresholds_rejection'), S.int_cfg.thresholds_rejection = []; end
    if ~isfield(S.int_cfg, 'cvsa_influence'),       S.int_cfg.cvsa_influence = 2.5; end
    if ~isfield(S.int_cfg, 'thresholds'),           S.int_cfg.thresholds = params.training_node.thresholds; end

    S.header_chunks = header_mi;
    S.header_chunks.framerate = S.framerate;
end
