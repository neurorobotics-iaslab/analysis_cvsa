% load of gmm, qda and for each trial shows the results
clear all; % close all;

addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_cvsa/equal_ros')
addpath(genpath('/home/paolo/Local/Matlab/yamlmatlab'));

%% ----------------- Initialization -----------------
nchannels = 39;
classes = [730 731];   
nclasses = length(classes);

[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
nFiles = length(filenames);

%% reasoning for one file
for idx_file = 1:nFiles
    fullpath_file_gdf = fullfile(pathname, filenames{idx_file});
    disp(['File ' num2str(idx_file) '/' num2str(nFiles)]);
    disp(['   Loading gdf file : ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_gdf);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;

    %% ----------------- load the file -----------------
    disp(['   Loading parameters file: ', filenames{idx_file}(1:end-3) 'yaml'])
    fullpath_file_parameters = [pathname(1:end-4) 'parameters/' filenames{idx_file}(1:end-3) 'yaml'];
    [ringBufferCfg, artifactCfg, processingCfg, gmmCfg, qdaCfg, integratorCfg] = loadParameters(fullpath_file_parameters);

    disp(['   Loading GMM file: ', gmmCfg.file_name])
    gmm_path = [pathname(1:end-4) gmmCfg.file_name];
    [gmmCfg.model, gmmCfg.params] = loadGMM(gmm_path);

    disp(['   Loading QDA file: ', qdaCfg.file_name])
    qda_path = [pathname(1:end-4) qdaCfg.file_name];
    qdaCfg.model = loadQDA(qda_path);

    %% ----------------- Artifact -----------------
    disp('   marking signal for artifact remotion')
    bufferSize = ringBufferCfg.size;
    chunkSize = processingCfg.chunkSize;
    eog.filterOrder = artifactCfg.filterOrder_EOG;
    eog.band = [artifactCfg.freq_low_EOG artifactCfg.freq_high_EOG];
    eog.label = channels_label(cell2mat(artifactCfg.EOG_ch));
    eog.h_threshold = artifactCfg.th_hEOG;
    eog.v_threshold = artifactCfg.th_vEOG;
    picks.filterOrder = artifactCfg.filterOrder_peaks;
    picks.freq = artifactCfg.freq_high_peaks; % remove antneuro problems
    picks.threshold = artifactCfg.th_peaks;
    artifact = artifact_rejection(c_signal, header, nchannels, bufferSize, chunkSize, eog, picks);

    %% ----------------- data processing -----------------
    disp('   processing EEG data') 
    filterOrder = processingCfg.filterOrder;
    band = processingCfg.bands; % i knwo we are using one band
    [signal_processed, header_processed] = processing_onlineROS_CSD_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize);

    %% ----------------- labels for the data -----------------
    disp('   extracting labels trials') 
    events = header_processed.EVENT;
    sampleRate = header_processed.SampleRate;
    cueDUR = events.DUR(ismember(events.TYP, classes));
    cueTYP = events.TYP(ismember(events.TYP, classes));

    fixPOS = events.POS(events.TYP == 786);
    fixDUR = events.DUR(events.TYP == 786);

    cfPOS = events.POS(events.TYP == 781);
    cfDUR = events.DUR(events.TYP == 781);

    boom = events.TYP(ismember(events.TYP, [897, 898, 899]));

    minDurCue = min(cueDUR);
    minDurFix = min(fixDUR);
    ntrial = length(fixDUR);

    %% ----------------- GMM -----------------
    disp('   applying the GMM to all the signal') 
    mu_data = cellfun(@(x) x(1), gmmCfg.params.mu); 
    sigma_data = cellfun(@(x) x(1), gmmCfg.params.sigma);
    o_l = cell2mat(gmmCfg.params.occipital_left_idx);
    o_r = cell2mat(gmmCfg.params.occipital_right_idx);
    c_l = cell2mat(gmmCfg.params.central_left_idx);
    c_r = cell2mat(gmmCfg.params.central_right_idx);
    excl_chs = cell2mat(gmmCfg.params.excluded_idx);

    type = gmmCfg.params.type;
    nfeatures = gmmCfg.params.nfeatures;

    % features extraction and classification
    sparsity = nan(size(signal_processed, 1), nfeatures);
    for idx_sample = 1:size(signal_processed,1)
        tmp_c_signal = signal_processed(idx_sample,:);
        [sparsity(idx_sample,:), ~] = compute_features_icnic(tmp_c_signal, type, o_l, o_r, c_l, c_r, nfeatures);
    end

    data_standardized = (sparsity - mu_data) ./ sigma_data;

    gmm_prob = posterior(gmmCfg.model, data_standardized);

    %% ----------------- QDA -----------------
    disp('   applying QDA to all the signal') 
    qda_prob = apply_qda_matrix(qdaCfg.model, log(signal_processed(:,qdaCfg.model.idchans{1})));

    %% ----------------- integrated prob -----------------
    event_start = 781;
    [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, gmm_prob, qda_prob, events, event_start, cell2mat(gmmCfg.params.classes));

    %% TODO ragiona online e mostra la lateralizzazione trial per trial

end