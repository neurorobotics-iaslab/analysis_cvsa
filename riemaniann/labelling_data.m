%% prepare the data for the rimannian space
clear all; % close all;

addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_cvsa/equal_ros')
addpath(genpath('/home/paolo/cvsa/ic_cvsa_ws/src/qda_cvsa/test'))
addpath(genpath('/home/paolo/Local/Matlab/yamlmatlab'));

%% Initialization
threshold_gmm_ic = 0.7;
band = [8 14];
band_str = mat2str(band);
signals = [];
raw_signal =[];
artifacts = [];
headers.TYP = [];
headers.POS = [];
headers.DUR = [];
classes = [730 731];      
cf_event = 781;
fix_event = 786;
nchannels = 39;
nclasses = length(classes);
filterOrder = 4;
avg = 1;% 0.75;

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);

datapath = '/home/paolo/cvsa/ic_cvsa_ws/src/';
path_gmm = [datapath ,'gmm_cvsa/cfg/gmm_c7_26112025_102244.yaml'];
save_riemman = [datapath, 'analysis_cvsa/riemaniann/test.mat'];

%% concatenate the files
nFiles = length(filenames);
for idx_file= 1: nFiles
    fullpath_file_shift = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_shift);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;
    sampleRate = header.SampleRate;


    disp('   [proc] power band');

    % for power band using hilbert transformation
    bufferSize = floor(avg*sampleRate);
    chunkSize = 32;
    eog.filterOrder = 4;
    eog.band = [1 7];
    eog.label = {'FP1', 'FP2', 'EOG'};
    eog.h_threshold = 80;
    eog.v_threshold = 80;
    muscle.filterOrder = 4;
    muscle.freq = 1; % remove antneuro problems
    muscle.threshold = 100;
    [signal_processed, header_processed] = processing_onlineROS_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize);
    artifact = artifact_rejection(c_signal, header, nchannels, bufferSize, chunkSize, eog, muscle);

    c_header = headers;
    c_header.sampleRate = header_processed.SampleRate/chunkSize;
    c_header.channels_labels = header_processed.Label;
    if isempty(find(header_processed.EVENT.TYP == 2, 1)) % no eye calibration
        c_header.TYP = cat(1, c_header.TYP, header_processed.EVENT.TYP);
        c_header.DUR = cat(1, c_header.DUR, header_processed.EVENT.DUR);
        c_header.POS = cat(1, c_header.POS, header_processed.EVENT.POS + size(signals, 1));
    else
        k = find(header_processed.EVENT.TYP == 1, 1);
        c_header.TYP = cat(1, c_header.TYP, header_processed.EVENT.TYP(k:end));
        c_header.DUR = cat(1, c_header.DUR, header_processed.EVENT.DUR(k:end));
        c_header.POS = cat(1, c_header.POS, header_processed.EVENT.POS(k:end) + size(signals, 1));
    end
    signals = cat(1, signals, signal_processed(:,:));
    artifacts = cat(1, artifacts, artifact(:,:));
    headers = c_header;

    raw_signal = cat(1, raw_signal, c_signal);
end


%% labels for the data
events = headers;
sampleRate = events.sampleRate;
cuePOS = events.POS(ismember(events.TYP, classes));
cueDUR = events.DUR(ismember(events.TYP, classes));
cueTYP = events.TYP(ismember(events.TYP, classes));

fixPOS = events.POS(events.TYP == 786);
fixDUR = events.DUR(events.TYP == 786);

cfPOS = events.POS(events.TYP == 781);
cfDUR = events.DUR(events.TYP == 781);

minDurCue = min(cueDUR);
minDurFix = min(fixDUR);
ntrial = length(cuePOS);


%% ----------------- GMM -----------------
% load the gmm and apply it
disp(['Loading model from: ', path_gmm]);
try
    modelData = ReadYaml(path_gmm);
    params = modelData.GmmModelCfg.params;
    model_params = modelData.GmmModelCfg.model_params;
catch ME
    disp('Error in the loading of the YAML file. Is the file path correct? Do you have YAMLMatlab installed? Is the path correct?.');
    disp(ME.message);
    return;
end
mu_datal = cellfun(@(x) x(1), params.mu);
sigma_data = cellfun(@(x) x(1), params.sigma);
nfeatures = params.nfeatures;

% model params
K = model_params.K;
ic_index = find(cell2mat(model_params.classes) == 1);
weights = cell2mat(model_params.weights);
means_cell = model_params.means;
means = cell2mat(means_cell);
cov_cell = model_params.covariances; 
covariances = zeros(nfeatures, nfeatures, K);
for k = 1:K
    % Extract the k-th matrix
    matrix_cell = cov_cell{k}; 
    
    matrix_k = cell2mat(matrix_cell);
    covariances(:, :, k) = matrix_k;
end
gmm_model = gmdistribution(means, covariances, weights);

fprintf('GMM model (K=%d, NFeatures=%d) loaded.\n', K, nfeatures);

% features extraction and classification
sparsity = nan(size(signals, 1), nfeatures);
for idx_sample = 1:size(signals,1)
    c_signal = signals(idx_sample,:);
    [sparsity(idx_sample,:), ~, ~, ~, ~, ~, ~, ~] = compute_features_gmm(c_signal);
end

data_standardized = (sparsity - mu_datal) ./ sigma_data;

gmm_prob_matlab = posterior(gmm_model, data_standardized);

disp('Data classified.');


%% extract the chunk for which the gmm says it is IC
mask_gmm = zeros(size(raw_signal,1), 1);
for i = find(gmm_prob_matlab(:,ic_index) >= threshold_gmm_ic)'
    mask_gmm((i-1)*chunkSize+1:i*chunkSize) = 1;
end

mask_artifacts = zeros(size(raw_signal, 1), 1);
for i = 1:size(artifacts, 1)
    if artifacts(i) == 1
        mask_artifacts((i-1)*chunkSize+1:i*chunkSize) = 1;
    end
end

mask_trial = zeros(size(raw_signal,1), 1);
mask_labels = zeros(size(raw_signal,1),1);
for i = 1:ntrial
    mask_trial(cfPOS(i)*chunkSize:cfPOS(i)*chunkSize+cfDUR(i)*chunkSize-1) = 1;
    mask_labels(cfPOS(i)*chunkSize:cfPOS(i)*chunkSize+cfDUR(i)*chunkSize-1) = cueTYP(i);
end

occipital = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'}; [~, ch_occipital] = ismember(occipital, channels_label);

save(save_riemman, 'raw_signal', 'mask_trial', 'mask_gmm', 'ntrial', 'channels_label', 'mask_labels', 'occipital', 'ch_occipital', 'mask_artifacts')
disp(save_riemman)

function [sparsity, label_sparsity, o_l, o_r, frontal, c_l, c_r, excluded_chs] = compute_features_gmm(c_signal) %% in the features must be update for subbands
    % take the one which contribute at the 95% of the energy
    sparsity = nan(3,1);
    label_sparsity = [{'LI'},{'GI'},{'GB'}];
    o_l = sort([29 13 30 37 33 34 17]); o_r = sort([31 15 32 35 36 38 18]);
    frontal = sort([3 4 5 20 21]); c_l = sort([6 22 25 8 11 27]); c_r = sort([7 24 26 10 12 28]);
    excluded_chs = [1,2,19];

    % --- LAP --- Calcola il LAP per tutti i punti nella finestra passata (window_signal)
    % show the occipital lateralization that is strongand present during the CVSA
    P_left_window  = mean(c_signal(o_l));
    P_right_window = mean(c_signal(o_r));
    LAP_history = (P_right_window - P_left_window) ./ (P_right_window + P_left_window + eps);
    sparsity(1) = abs(LAP_history); % LAP_Mean

    % --- Gini Index + Occipital Power ---  -> when high there is a zone stronger, so IC
    % show the focusing is weighted in with the power in the occipital part, in this way ig CVSA then strong value 
    non_zeros_chs = setdiff(1:size(c_signal,2), excluded_chs);
    global_mean = mean(c_signal(non_zeros_chs)); % car filter
    current_signal_normalized = c_signal - global_mean; % remove the global energy
    mean_roi_raw = [mean(current_signal_normalized(frontal)), mean(current_signal_normalized(c_l)), ...
        mean(current_signal_normalized(c_r)), mean(current_signal_normalized(o_l)), ...
        mean(current_signal_normalized(o_r))];
    mean_roi = abs(mean_roi_raw); % make sure the energy is positive--> we are using peak and valli with same significance
    mean_roi_ordered = sort(mean_roi);
    n = length(mean_roi_ordered);
    sum_roi_p = 0;
    for i = 1:n
        sum_roi_p = sum_roi_p + (n+1-i) * mean_roi_ordered(i);
    end
    total_sum = sum(mean_roi_ordered);
    if total_sum > 0
        gi = (1/n) * (n+1-2*sum_roi_p/total_sum);
    else
        gi = 0;
    end
    % compute the weight factor
    mean_roi_raw = [mean(c_signal(frontal)), mean(c_signal(c_l)), ...
                    mean(c_signal(c_r)), mean(c_signal(o_l)), ...
                    mean(c_signal(o_r))];
    pot_occipital = max(mean_roi_raw(4), mean_roi_raw(5));
    pot_total_roi = sum(mean_roi_raw); % Somma di F, CL, CR, OL, OR
    if pot_total_roi > 0
        occipital_power = pot_occipital / pot_total_roi;
    else
        occipital_power = 0; 
    end
    sparsity(2) = occipital_power * gi;

    % --- GB ---
    % return the global power mean, 
    sparsity(3) = global_mean;
end