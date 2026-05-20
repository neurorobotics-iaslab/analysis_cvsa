%% MAIN_EVALUATE_METRICS  Offline simulator and detailed performance metric evaluator of the ROS BCI pipeline.
%
%   This script:
%   1. Prompts the user to select a GDF recording (or uses a predefined path).
%   2. Automatically locates the companion YAML parameter dump (saved during online session).
%   3. Replays the exact causal BCI pipeline (CAR -> per-band Butterworth filtering ->
%      ring buffering -> CSP spatial filtering -> LDA classification -> Integrator).
%   4. Segment the results by trial (events 781 to 33549).
%   5. Performs extensive sample-by-sample and trial-by-trial metrics analysis:
%       - Single-sample classification accuracy, target confidence, and artifact rates.
%       - Trial Hit/Miss classification and Time-to-Hit (TTH).
%       - Class-specific metrics to diagnose biases.
%   6. Prints a detailed ASCII report to the MATLAB command window.
%   7. Displays a diagnostic plot of the session.
%
%   Author: Antigravity

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

% --- File selection ------------------------------------------------------
default_dir = '/home/paolo/bci_vr_ws/test_node_data/test_evaluation/recorded';
if ~isfolder(default_dir)
    default_dir = '/home/paolo/bci_vr_ws/recordings';
end

[gdf_name, gdf_dir] = uigetfile({'*.gdf', 'GDF recordings (*.gdf)'}, ...
                                'Select a GDF recording for metrics evaluation', default_dir);
if isequal(gdf_name, 0)
    fprintf('Evaluation canceled.\n');
    return;
end
gdf_path = fullfile(gdf_dir, gdf_name);

%% --- 1) Load Data and Parameters ----------------------------------------
fprintf('\n========================================================================\n');
fprintf('  OFFLINE METRIC EVALUATION FOR: %s\n', gdf_name);
fprintf('========================================================================\n');

[signal, header, basename] = load_gdf(gdf_path);
[params, ~] = load_params_yaml(gdf_path);

paradigm  = params.integrator.paradigm;
fs        = double(params.acquisition.samplerate);
framerate = double(params.acquisition.framerate);
chunk_size = round(fs / framerate);

if abs(fs - header.SampleRate) > 1e-3
    fprintf('[Warning] YAML samplerate (%.1f Hz) differs from GDF SampleRate (%.1f Hz). Using GDF.\n', ...
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

% --- 2) Load CSP and sLDA Models -----------------------------------------
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

% --- 3) Run Processing Streams -------------------------------------------
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

% --- 4) Run Causal Artifact Detection ------------------------------------
art_cfg = params.ArtifactCfg.params;
art_cfg.EOG_ch_names = to_strcell(art_cfg.EOG_ch_names);
cfg_art = struct('samplerate', fs, 'chunk_size', chunk_size, 'bufsize_artifact', bufsize_art);
[art_flags, info_art] = detect_artifacts(signal, header, art_cfg, cfg_art);

% --- 5) Apply sLDA Classification ----------------------------------------
p_mi_aligned   = []; if use_mi,   p_mi_aligned   = apply_slda(features_mi, slda_mi,   csp_mi.bands  ); end
p_cvsa_aligned = []; if use_cvsa, p_cvsa_aligned = apply_slda(features_cv, slda_cvsa, csp_cvsa.bands); end

% --- 6) Run Trial Integrator ---------------------------------------------
int_cfg = params.integrator;
if ~isfield(int_cfg, 'increment'),               int_cfg.increment = 1; end
if ~isfield(int_cfg, 'thresholds_rejection'),    int_cfg.thresholds_rejection = []; end

if use_mi,   header_chunks = header_mi;
else,        header_chunks = header_cv;
end
header_chunks.framerate = framerate;

trials = integrate_signal(p_mi_aligned, p_cvsa_aligned, art_flags, ...
                          header_chunks, int_cfg, paradigm);

%% --- 7) Detailed Metrics Computations -----------------------------------
n_trials = numel(trials);
classes = to_vec(int_cfg.classes);

% Metrics storage
t_outcome     = cell(n_trials, 1);    % 'HIT' or 'MISS'
t_tth         = NaN(n_trials, 1);     % Time-to-hit in seconds
t_duration    = NaN(n_trials, 1);     % Trial continuous feedback duration in seconds
t_sample_acc  = NaN(n_trials, 1);     % Sample accuracy percentage
t_mean_conf   = NaN(n_trials, 1);     % Mean confidence (prob of correct class)
t_art_rate    = NaN(n_trials, 1);     % Artifact rate percentage

for t = 1:n_trials
    tr = trials(t);
    cf_range = (tr.n_pre + 1) : size(tr.normalized, 1); % Exclude initial reset frame
    
    t_duration(t) = numel(cf_range) / framerate;
    
    % Hit/Miss outcome
    if tr.pass
        t_outcome{t} = 'HIT';
        % Find first frame in CF window where normalized target reaches 1.0
        hit_frame = find(tr.normalized(cf_range, tr.target_class) >= 1.0, 1);
        if ~isempty(hit_frame)
            t_tth(t) = (hit_frame - 1) / framerate;
        end
    else
        t_outcome{t} = 'MISS';
    end
    
    % Single-sample classification accuracy
    if ~isnan(tr.target_class)
        % For each sample, prediction is the class with max sLDA raw probability
        [~, pred_classes] = max(tr.raw(cf_range, :), [], 2);
        correct_samples = (pred_classes == tr.target_class);
        t_sample_acc(t) = mean(correct_samples) * 100;
        
        % Mean confidence: average sLDA raw probability of the true target class
        t_mean_conf(t)  = mean(tr.raw(cf_range, tr.target_class)) * 100;
    else
        t_sample_acc(t) = NaN;
        t_mean_conf(t)  = NaN;
    end
    
    % Artifact rate
    t_art_rate(t) = mean(tr.artifact(cf_range)) * 100;
end

% --- 8) Print Session Report ---------------------------------------------
fprintf('\n%-8s | %-12s | %-8s | %-9s | %-13s | %-15s | %-14s | %-12s\n', ...
        'Trial #', 'Target Class', 'Outcome', 'TTH (s)', 'Duration (s)', 'Sample Acc (%)', 'Confidence (%)', 'Artifact (%)');
fprintf('%s\n', repmat('-', 1, 107));

for t = 1:n_trials
    target_str = sprintf('Class %d (%d)', trials(t).target_class, trials(t).onset_code);
    if isnan(trials(t).target_class)
        target_str = 'Unknown';
    end
    
    tth_str = 'N/A';
    if ~isnan(t_tth(t))
        tth_str = sprintf('%.2f', t_tth(t));
    end
    
    fprintf('  %-6d | %-12s |   %-5s   |   %-5s   |     %-7.2f |      %-9.1f |       %-8.1f |     %-7.1f\n', ...
            t, target_str, t_outcome{t}, tth_str, t_duration(t), t_sample_acc(t), t_mean_conf(t), t_art_rate(t));
end
fprintf('%s\n', repmat('=', 1, 107));

% --- 9) Compute Session-Level Summary Metrics ----------------------------
hits = strcmp(t_outcome, 'HIT');
total_hits = sum(hits);
hit_rate = (total_hits / n_trials) * 100;
mean_tth = mean(t_tth(~isnan(t_tth)));
mean_sample_acc = mean(t_sample_acc, 'omitnan');
mean_confidence = mean(t_mean_conf, 'omitnan');
mean_artifact = mean(t_art_rate);

% Class-specific summaries
fprintf('  GLOBAL SUMMARY STATISTICS:\n');
fprintf('  - Total Trials Processed: %d\n', n_trials);
fprintf('  - Global Hit Rate:        %.1f%% (%d / %d trials)\n', hit_rate, total_hits, n_trials);
if total_hits > 0
    fprintf('  - Average Time-to-Hit:    %.2f seconds (HIT trials only)\n', mean_tth);
else
    fprintf('  - Average Time-to-Hit:    N/A seconds (0 HIT trials)\n');
end
fprintf('  - Overall Sample-Accuracy: %.1f%%\n', mean_sample_acc);
fprintf('  - Overall sLDA Confidence: %.1f%%\n', mean_confidence);
fprintf('  - Overall Artifact Rate:   %.1f%%\n', mean_artifact);
fprintf('%s\n', repmat('-', 1, 50));

% Loop through each configured class to diagnose potential asymmetry/bias
for c_idx = 1:numel(classes)
    cls_code = classes(c_idx);
    class_mask = [trials.onset_code] == cls_code;
    n_class_trials = sum(class_mask);
    
    if n_class_trials > 0
        cls_hits = hits(class_mask);
        cls_hit_rate = (sum(cls_hits) / n_class_trials) * 100;
        cls_tth = t_tth(class_mask);
        cls_mean_tth = mean(cls_tth(~isnan(cls_tth)));
        cls_sample_acc = mean(t_sample_acc(class_mask), 'omitnan');
        cls_confidence = mean(t_mean_conf(class_mask), 'omitnan');
        
        fprintf('  Class %d (Code: %d) Summary:\n', c_idx, cls_code);
        fprintf('    * Trials Count:     %d\n', n_class_trials);
        fprintf('    * Hit Rate:         %.1f%% (%d / %d trials)\n', cls_hit_rate, sum(cls_hits), n_class_trials);
        if sum(cls_hits) > 0
            fprintf('    * Avg Time-to-Hit:  %.2f seconds\n', cls_mean_tth);
        else
            fprintf('    * Avg Time-to-Hit:  N/A seconds\n');
        end
        fprintf('    * Sample Accuracy:  %.1f%%\n', cls_sample_acc);
        fprintf('    * Mean Confidence:  %.1f%%\n', cls_confidence);
    else
        fprintf('  Class %d (Code: %d): No trials found in session.\n', c_idx, cls_code);
    end
    fprintf('%s\n', repmat('-', 1, 50));
end
fprintf('========================================================================\n\n');

% --- 10) Plot Trials Diagnostic Panel ------------------------------------
plot_trials(trials, int_cfg, framerate, paradigm, basename);



% -------------------------------------------------------------------------
%   Helper function: local paradigm key mapping
% -------------------------------------------------------------------------
function k = first_paradigm_key(paradigm)
    switch paradigm
        case 'cvsa', k = 'cvsa';
        otherwise,   k = 'mi';
    end
end

function s = to_strcell(val)
    if ischar(val)
        s = {val};
    elseif iscell(val)
        s = cellfun(@char, val, 'UniformOutput', false);
    else
        s = {};
    end
end
