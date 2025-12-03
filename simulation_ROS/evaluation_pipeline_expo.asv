% load of gmm, qda and for each trial shows the results
clear all; % close all;

addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_cvsa/equal_ros')
addpath(genpath('/home/paolo/Local/Matlab/yamlmatlab'));

%% Initialization
nchannels = 39;
classes = [730 731];   
nclasses = length(classes);

[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
nFiles = length(filenames);

%% load the file and perfo
for idx_file = 1:nFiles
    fullpath_file_gdf = fullfile(pathname, filenames{idx_file});
    disp(['File ' num2str(idx_file) '/' num2str(nFiles)]);
    disp(['   Loading gdf file : ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_gdf);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;
    sampleRate = header.SampleRate;

    disp(['   Loading parameters file: ', filenames{idx_file}(1:end-3) 'yaml'])
    fullpath_file_parameters = [pathname(1:end-4) 'parameters/' filenames{idx_file}(1:end-3) 'yaml'];
    [ringBufferCfg, artifactCfg, processingCfg, gmmCfg, qdaCfg, integratorCfg] = loadParameters(fullpath_file_parameters);

    disp(['   Loading GMM file: ', gmmCfg.file_name])
    gmm_path = [pathname(1:end-4) gmmCfg.file_name];
    [gmmCfg.model, mu_data, sigma_data] = loadGMM(gmm_path);

    disp(['   Loading QDA file: ', qdaCfg.file_name])
    qda_path = [pathname(1:end-4) qdaCfg.file_name];
    qdaCfg.model = loadQDA(qda_path);

    disp('   ARTIFACT')
    bufferSize = ringBufferCfg.size;
    chunkSize = processingCfg.chunkSize;
    eog.filterOrder = artifactCfg.filterOrder_EOG;
    eog.band = [artifactCfg.freq_low_EOG artifactCfg.freq_high_EOG];
    eog.label = channels_label(cell2mat(artifact.EOG_ch));
    eog.h_threshold = artifactCfg.th_hEOG;
    eog.v_threshold = artifactCfg.th_vEOG;
    picks.filterOrder = artifactCfg.filterOrder_peaks;
    picks.freq = artifactCfg.freq_high_peaks; % remove antneuro problems
    picks.threshold = artifactCfg.th_peaks;
    artifact = artifact_rejection(c_signal, header, nchannels, bufferSize, chunkSize, eog, picks);

    disp('   PROCESSING DATA') ------
    [signal_processed, header_processed] = processing_onlineROS_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize);
end


threshold_gmm_ic = 0.7;
band = [8 14];
band_str = mat2str(band);
filterOrder = 4;
avg = 1;% 0.75;
eog_threshold = 500;

%% concatenate the files
for idx_file= 1: nFiles
    disp('   [proc] power band');

    % for power band using hilbert transformation
    
    [signal_processed, header_processed] = processing_onlineROS_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize);
    artifact = artifact_rejection(c_signal, header, nchannels, bufferSize, chunkSize, eog, picks);

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

o_l = cell2mat(params.occipital_left_idx);
o_r = cell2mat(params.occipital_right_idx);
c_l = cell2mat(params.central_left_idx);
c_r = cell2mat(params.central_right_idx);
excl_chs = cell2mat(params.excluded_idx);

type = model_params.type;

% features extraction and classification
sparsity = nan(size(signals, 1), nfeatures);
for idx_sample = 1:size(signals,1)
    c_signal = signals(idx_sample,:);
    [sparsity(idx_sample,:), ~] = compute_features_icnic(c_signal, type, o_l, o_r, c_l, c_r, excl_chs);
end

data_standardized = (sparsity - mu_datal) ./ sigma_data;

gmm_prob_matlab = posterior(gmm_model, data_standardized);

disp('Data classified.');

%% ----------------- QDA -----------------
features = signals;

qdaCfg = loadQDA(yaml_QDA_path);

prob_qda_matlab = apply_qda_matrix(qdaCfg, log(features(:,qdaCfg.idchans{1})));

%% ----------------- plot prob integrated -----------------
alpha = 0.98;
bufferSize = 48;
for idx_trial =1:ntrial
    start_trial = cuePOS(idx_trial);
    
    end_trial = cfPOS(idx_trial) + cfDUR(idx_trial)-1;

    nsamples_trial = end_trial - start_trial;
    c_gmm_prob = gmm_prob_matlab(start_trial:end_trial,:);
    c_qda_prob = prob_qda_matlab(start_trial:end_trial,:);
    c_artifact = artifacts(start_trial:end_trial,:);
    trial_start_cf = cfPOS(idx_trial) - start_trial;
    c_power = log(signals(start_trial:end_trial,:));
    mask = [];

    c_integrated = ones(nsamples_trial, 1)*1/nclasses;

    buffer = ones(1, nsamples_trial)*1/nclasses;

    for idx_sample = 1:nsamples_trial
        if idx_sample >= trial_start_cf
            if c_gmm_prob(idx_sample,ic_index) >= threshold_gmm_ic && c_artifact(idx_sample) == 0
                if c_qda_prob(idx_sample, 1) >= 0.5
                    inc = 1/bufferSize;
                    buffer(idx_sample)=min(buffer(idx_sample-1) + inc,1);
                    tmp = 1;
                else
                    inc = -1/bufferSize;
                    buffer(idx_sample)=max(buffer(idx_sample-1) + inc,0);
                    tmp = 0;
                end
                
                mask = cat(1, mask, 1);
                c_integrated(idx_sample) = c_integrated(idx_sample-1) * alpha + (1-alpha) * tmp;
            else
                mask = cat(1, mask, 0);
                c_integrated(idx_sample) = c_integrated(idx_sample-1);
                buffer(idx_sample)=buffer(idx_sample-1);
            end
        else
            mask = cat(1, mask, 0);
        end
    end

    figure();
    subplot(411)
    imagesc(c_power')
    hold on;
    xline(cueDUR(idx_trial), 'LineStyle','-');
    hold off;
    yticks(1:nchannels); yticklabels(channels_label)
    title('Log band power')

    subplot(412)
    plot(c_gmm_prob(:,ic_index))
    hold on;
    plot(c_artifact);
    yline(threshold_gmm_ic, 'LineStyle','--');
    xline(cueDUR(idx_trial), 'LineStyle','-');
    hold off;
    legend('gmm prob', 'artifact', 'threshold ic');
    title('artifacts and gmm probabilities')


    subplot(413)
    tmp_prob = c_qda_prob;
    tmp_prob(mask == 0,1) = nan;
    plot(c_qda_prob(:,1))
    hold on
    scatter(1:size(c_qda_prob, 1), c_qda_prob(:,1), 15, 'black', 'filled')
    scatter(1:size(c_qda_prob, 1), tmp_prob(:,1), 15, 'green', 'filled')
    yline(0.5, 'LineStyle','--');
    xline(cueDUR(idx_trial), 'LineStyle','-');
    hold off
    ylim([0 1])
    legend('qda prob','qda not used', 'qda prob used')
    title('classifier probability')

    subplot(414)
    plot(c_integrated)
    hold on
    plot(buffer);
    yline(0.7, 'LineStyle','--');
    yline(0.3, 'LineStyle','--');
    xline(cueDUR(idx_trial), 'LineStyle','-');
    hold off
    legend('buffer', 'integrated prob')
    ylim([0 1])
    title('integrated signal')

    sgtitle(['trial ' num2str(idx_trial) ' | class aked ' num2str(cueTYP(idx_trial))])
end
