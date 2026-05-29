%% MAIN_SIMULATE  Offline MATLAB simulator of the ROS BCI pipeline.
%
%   1. Pop up a file picker to choose a GDF recording.
%   2. Load the sibling YAML (same basename) that the bag_bci recorder
%      saved: it has every node's parameters (CAR, artifact, CSP, sLDA,
%      integrator).
%   3. Replay the pipeline chunk-by-chunk like ROS would:
%        raw EEG -> CAR + per-band Butterworth + ringbuf + CSP + power
%                                                              -> features
%        raw EEG -> artifact detector                          -> flags
%        features -> sLDA log + sigmoid                        -> P(c)
%        P(c)    -> leaky-WTA integrator (resets at every 781) -> integrated
%   4. Plot one panel per continuous-feedback trial (event 781..781+DUR-1).
%
%   Works for paradigm = 'mi' | 'cvsa' | 'hybrid' automatically — the
%   paradigm is read from the YAML.

clear; clc; close all;

% --- Make every subfolder visible to the simulator -----------------------
this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
addpath(fullfile(this_dir, 'io'));
addpath(fullfile(this_dir, 'processing'));
addpath(fullfile(this_dir, 'artifacts'));
addpath(fullfile(this_dir, 'classifier'));
addpath(fullfile(this_dir, 'integrator'));
addpath(fullfile(this_dir, 'plotting'));
addpath(fullfile(this_dir, 'utils'));

% --- GUI: pick the GDF; everything else flows from the sibling YAML -----
default_dir = '/home/paolo/bci_vr_ws/test_node_data/test_evaluation/recorded';
if ~isfolder(default_dir)
    default_dir = '/home/paolo/bci_vr_ws/recordings';
end
[gdf_name, gdf_dir] = uigetfile({'*.gdf', 'GDF recordings (*.gdf)'}, ...
                                'Select a GDF recording', default_dir);
if isequal(gdf_name, 0)
    error('main_simulate:cancel', 'No GDF selected.');
end
gdf_path = fullfile(gdf_dir, gdf_name);

%% --- 1) Load GDF --------------------------------------------------------
[signal, header, basename] = load_gdf(gdf_path);

% --- 2) Load companion YAML (rosparam dump) -----------------------------
[params, ~] = load_params_yaml(gdf_path);

paradigm  = params.integrator.paradigm;
fs        = double(params.acquisition.samplerate);
framerate = double(params.acquisition.framerate);
chunk_size = round(fs / framerate);
if abs(fs - header.SampleRate) > 1e-3
    log_step('main_simulate: YAML samplerate=%.1f != GDF samplerate=%.1f -> using GDF', ...
             fs, header.SampleRate);
    fs = header.SampleRate;
    chunk_size = round(fs / framerate);
end

bufsize_proc = double(params.RingBufferCfg.params.size);
bufsize_art  = double(params.RingBufferCfgArtifact.params.size);
eog_names    = to_strcell(params.CarCfg.params.EOG_ch_names);

% Dynamically resolve do_car for each stream
do_car_mi = true;
if isfield(params, 'processing_fbcsp_mi')
    do_car_mi = logical(params.processing_fbcsp_mi.do_car);
end
do_car_cvsa = true;
if isfield(params, 'processing_fbcsp_cvsa')
    do_car_cvsa = logical(params.processing_fbcsp_cvsa.do_car);
end

log_step('main_simulate: paradigm=%s, fs=%g, framerate=%g, chunk=%d, bufproc=%d, bufart=%d', ...
         paradigm, fs, framerate, chunk_size, bufsize_proc, bufsize_art);

% --- 3) Load CSP + sLDA per used paradigm --------------------------------
use_mi   = ismember(paradigm, {'mi', 'hybrid'});
use_cvsa = ismember(paradigm, {'cvsa', 'hybrid'});
csp_mi   = []; slda_mi   = [];
csp_cvsa = []; slda_cvsa = [];
if use_mi
    csp_mi   = load_csp(params, 'mi');
    slda_mi  = load_slda(params, 'mi');
end
if use_cvsa
    csp_cvsa  = load_csp(params, 'cvsa');
    slda_cvsa = load_slda(params, 'cvsa');
end

% --- 4) Apply processing (one stream per paradigm) ----------------------
features_mi = []; header_mi = [];
features_cv = []; header_cv = [];
info_proc   = [];
if use_mi
    proc_cfg_mi = struct('samplerate', fs, 'chunk_size', chunk_size, ...
                         'bufsize', bufsize_proc, 'filter_order', 4, ...
                         'do_car', do_car_mi, 'eog_names', {eog_names});
    [features_mi, header_mi, info_proc] = apply_processing(signal, header, csp_mi, proc_cfg_mi);
end
if use_cvsa
    proc_cfg_cvsa = struct('samplerate', fs, 'chunk_size', chunk_size, ...
                           'bufsize', bufsize_proc, 'filter_order', 4, ...
                           'do_car', do_car_cvsa, 'eog_names', {eog_names});
    [features_cv, header_cv, info2] = apply_processing(signal, header, csp_cvsa, proc_cfg_cvsa);
    if isempty(info_proc), info_proc = info2; end
end

% --- 5) Apply artifact detection on raw signal --------------------------
art_cfg = params.ArtifactCfg.params;
art_cfg.EOG_ch_names = to_strcell(art_cfg.EOG_ch_names);
cfg_art = struct('samplerate', fs, 'chunk_size', chunk_size, 'bufsize_artifact', bufsize_art);
[art_flags, info_art] = detect_artifacts(signal, header, art_cfg, cfg_art);

% --- 6) Apply sLDA over the whole feature stream ------------------------
%   Note on alignment: features (apply_processing) and art_flags
%   (detect_artifacts) both live on the same chunk-index axis (k = 1..n_chunks).
%   The artifact ringbuf fills first (bufsize_art / chunk_size chunks earlier
%   than the processing ringbuf), so art_flags are already "valid" by the
%   time features become valid -- the +bufsize_proc/chunk offset is implicit
%   in starting integration only at the first 781, which is well past both
%   buffers' fill points.
n_chunks    = info_proc.n_chunks;
first_valid = info_proc.first_valid_chunk;
art_lead    = max(0, round((bufsize_proc - bufsize_art) / chunk_size));

p_mi_aligned   = []; if use_mi,   p_mi_aligned   = apply_slda(features_mi, slda_mi,   csp_mi.bands  ); end
p_cvsa_aligned = []; if use_cvsa, p_cvsa_aligned = apply_slda(features_cv, slda_cvsa, csp_cvsa.bands); end

log_step('main_simulate: streams aligned (n_chunks=%d, first_valid_proc=%d, art_lead=%d chunks)', ...
         n_chunks, first_valid, art_lead);

% --- 8) Integrate per trial ---------------------------------------------
int_cfg = params.integrator;
% Pull the dynamic_reconfigure-able fields with safe defaults
if ~isfield(int_cfg, 'increment'),               int_cfg.increment = 1; end
if ~isfield(int_cfg, 'thresholds_rejection'),    int_cfg.thresholds_rejection = []; end
if ~isfield(int_cfg, 'cvsa_influence'),          int_cfg.cvsa_influence = 2.5; end
if ~isfield(int_cfg, 'thresholds'),              int_cfg.thresholds = params.training_node.thresholds; end

% header_chunks for trial extents (use whichever paradigm produced a header)
if use_mi,   header_chunks = header_mi;
else,        header_chunks = header_cv;
end
header_chunks.framerate = framerate;

trials = integrate_signal(p_mi_aligned, p_cvsa_aligned, art_flags, ...
                          header_chunks, int_cfg, paradigm);

% --- 9) Plot per-trial panels -------------------------------------------
plot_trials(trials, int_cfg, framerate, paradigm, basename);

log_step('main_simulate: done -> %d/%d trials PASS', sum([trials.pass]), numel(trials));



% ----------------------------------------------------------------------
%   Local helper: paradigm key used in the YAML processing_fbcsp_<key>.
%   For hybrid we use the MI processing-node block to read do_car (both
%   blocks share the same value at launch time).
% ----------------------------------------------------------------------
function k = first_paradigm_key(paradigm)
    switch paradigm
        case 'cvsa', k = 'cvsa';
        otherwise,   k = 'mi';
    end
end
