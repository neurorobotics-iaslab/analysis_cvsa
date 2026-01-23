% load of gmm, qda and for each trial shows the results
clear all; % close all;

addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/equal_ros')
addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/utils')
addpath(genpath('/home/paolo/Local/Matlab/yamlmatlab'));

%% ----------------------------------
%             MY METHOD
%  ----------------------------------
% ----------------- inizialization -----------------
nchannels = 16;
classes = [769 770, 783];   

[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files My method', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
nFiles = length(filenames);
Database = [];

%% reasoning for one file
for idx_file = 1:nFiles
    fullpath_file_gdf = fullfile(pathname, filenames{idx_file});
    disp(['File ' num2str(idx_file) '/' num2str(nFiles)]);
    disp(['   Loading gdf file : ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_gdf);
    c_signal = c_signal(:,1:nchannels);
    channels_label =  {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'Fp1', 'CP1', 'CPz', 'CP2', 'Fp2'};

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
    excl_ch = {'Fp1', 'Fp2', 'EOG'};
    [found, indices] = ismember(excl_ch, channels_label);
    excl_chs = indices(found);
    bands = processingCfg.bands; 
    nbands = size(bands, 1);
    signals = cell(1, nbands);
    for idx_band=1:nbands
        band = bands(idx_band,:);
        [signal_processed, header_processed] = processing_onlineROS_CAR_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize, excl_chs);
        signals{idx_band} = signal_processed;
    end
    

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
    c_l = cell2mat(gmmCfg.params.central_left_idx);
    c_r = cell2mat(gmmCfg.params.central_right_idx);
    c_c = cell2mat(gmmCfg.params.central_idx);
    excl_chs = cell2mat(gmmCfg.params.excluded_idx);

    type = gmmCfg.params.type;
    nfeatures = gmmCfg.params.nfeatures;

    % features extraction and classification
    sparsity = nan(size(signals{1}, 1), nfeatures);
    for idx_sample = 1:size(signals{1},1)
        tmp_sparsity = [];
        for idx_band=1:nbands
            tmp = compute_features_icnic_mi(signals{idx_band}(idx_sample,:), c_l, c_r, c_c);
            tmp_sparsity = [tmp_sparsity, tmp];
        end
        sparsity(idx_sample,:) = tmp_sparsity;
    end

    data_standardized = (sparsity - mu_data) ./ sigma_data;

    gmm_prob = posterior(gmmCfg.model, data_standardized);

    %% ----------------- QDA -----------------
    disp('   applying QDA to all the signal') 
    X = [];
    for idx_band=1:nbands
        for idx_band2=1:nbands
            if all(bands(idx_band,:) == qdaCfg.model.bands(idx_band2,:))
                chs = qdaCfg.model.idchans{idx_band2};
                tmp = signals{idx_band}(:,chs); 
            end
        end
        X = [X, tmp];
    end
    qda_prob = apply_qda_matrix(qdaCfg.model, log(X));

    %% ----------------- integrated prob -----------------
    event_start = 781;
    rejection_look_gmm = false;
    [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, gmm_prob, qda_prob, events, event_start, cell2mat(gmmCfg.params.classes), rejection_look_gmm);

    %% ----------------- plot prob integrated -----------------
    ic_index = find(cell2mat(gmmCfg.params.classes) == integratorCfg.ic_class_label);
    sampleRate_ros = bufferSize/chunkSize;
    do_plot = false;
    r_square_data_all_1 = []; r_square_data_all_2 = []; r_square_label_all = [];
    r_square_data_gmm_1 = []; r_square_data_gmm_2 = []; r_square_label_gmm = [];
    for idx_trial = 1:ntrial
        start_trial = fixPOS(idx_trial);
        end_trial = cfPOS(idx_trial) + cfDUR(idx_trial)-1;
        c_power_1 = log(signals{1}(start_trial:end_trial,:));
        c_power_2 = log(signals{2}(start_trial:end_trial,:));
        c_gmm_prob = gmm_prob(start_trial:end_trial, :);
        c_artifact = artifact(start_trial:end_trial);
        c_qda_prob = qda_prob(start_trial:end_trial,:);
        c_mask = mask(start_trial:end_trial);
        c_integrated = integrated_prob(start_trial:end_trial,:);
        trial_dur = size(c_power_1, 1);

        % for metrics r^2
        if ismember(cueTYP(idx_trial), [769 770])
            r_square_data_all_1 = [r_square_data_all_1; c_power_1];
            r_square_data_all_2 = [r_square_data_all_2; c_power_2];
            r_square_label_all = [r_square_label_all; repmat(cueTYP(idx_trial), trial_dur, 1)];
            tmp_c_power = c_power_1(c_gmm_prob(:,ic_index) > 0.5,:);
            r_square_data_gmm_1 = [r_square_data_gmm_1; tmp_c_power];
            tmp_c_power = c_power_2(c_gmm_prob(:,ic_index) > 0.5,:);
            r_square_data_gmm_2 = [r_square_data_gmm_2; tmp_c_power];
            r_square_label_gmm = [r_square_label_gmm; repmat(cueTYP(idx_trial), sum(c_gmm_prob(:,ic_index) > integratorCfg.ic_threshold), 1)];
        end

        if do_plot
            figure();
            subplot(511)
            imagesc(c_power_1')
            hold on;
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off;
            yticks(1:nchannels); yticklabels(channels_label)
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            title(['Log band power ' num2str(bands(1,1)) '-' num2str(bands(1,2))])

            subplot(512)
            imagesc(c_power_2')
            hold on;
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off;
            yticks(1:nchannels); yticklabels(channels_label)
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            title(['Log band power ' num2str(bands(2,1)) '-' num2str(bands(2,2))])

            subplot(513)
            plot(c_gmm_prob(:,ic_index))
            hold on;
            plot(c_artifact);
            yline(integratorCfg.ic_threshold, 'LineStyle','--');
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off;
            legend('gmm prob', 'artifact', 'threshold ic');
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            title('artifacts and gmm probabilities')


            subplot(514)
            tmp_prob = c_qda_prob;
            tmp_prob(c_mask == 0,1) = nan;
            plot(c_qda_prob(:,1))
            hold on
            scatter(1:size(c_qda_prob, 1), c_qda_prob(:,1), 15, 'black', 'filled')
            scatter(1:size(c_qda_prob, 1), tmp_prob(:,1), 15, 'green', 'filled')
            yline(0.5, 'LineStyle','--');
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off
            ylim([0 1])
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            legend('qda prob','qda not used', 'qda prob used')
            title('classifier probability')

            subplot(515)
            plot(c_integrated(:,1))
            hold on
            yline(integratorCfg.feedbackThs(1), 'LineStyle','--');
            yline(1-integratorCfg.feedbackThs(2), 'LineStyle','--');
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off
            legend('integrated prob')
            ylim([0 1])
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            title('integrated signal')

            if boom(idx_trial) == 897
                strboom =  'HIT';
            elseif boom(idx_trial) == 898
                strboom = 'MISS';
            elseif boom(idx_trial) == 899
                strboom = 'TIMEOUT';
            else
                disp('ERROR')
            end
            sgtitle(['trial ' num2str(idx_trial) ' | class aked ' num2str(cueTYP(idx_trial)) ' | ' strboom])
        end
    end

    %% take accuracy
    [m] = computeMetrics(integratorCfg, artifact, gmm_prob, qda_prob, integrated_prob, events, event_start, cell2mat(gmmCfg.params.classes), [769 770 783]);
    current_method = 'gmm';
    % Crea entry per il database
    entry = struct();
    entry.File = filenames{idx_file};
    entry.Method = current_method; % 'GMM' o 'Traditional' (impostalo tu nel loop)
    
    % --- ACTIVE TASK METRICS ---
    % Trial Accuracy
    entry.Act_Acc_Raw = m.accuracy.trial.active.raw * 100;
    entry.Act_Acc_NoTimeout = m.accuracy.trial.active.no_timeout * 100;
    
    % Sample Accuracy
    entry.Act_SampAcc_All = m.accuracy.sample.active.all * 100;
    entry.Act_SampAcc_NoTimeout = m.accuracy.sample.active.no_timeout * 100;
    
    % Time (Convertito in secondi)
    % Assicurati che sampleRate_ros sia corretto (es. 16 Hz)
    entry.Act_Time_Hit = m.time.active.hit_avg / sampleRate_ros;
    entry.Act_Time_Miss = m.time.active.miss_avg / sampleRate_ros;
    
    % Stability
    entry.Act_Wobble = m.stability.active.wobble_avg;
    entry.Act_WDR = m.stability.active.wdr_avg * 100;
    
    % --- REST TASK METRICS ---
    entry.Rest_Acc = m.accuracy.trial.rest.acc * 100; % Corretti (Timeout)
    entry.Rest_FPR = m.accuracy.trial.rest.fpr * 100; % False Positives
    entry.Rest_Time_FP = m.time.rest.fp_avg / sampleRate_ros;
    entry.Rest_MaxDev = m.stability.rest.max_dev_avg;

    % Aggiungi al database
    % (Inizializza Database = [] fuori dal loop)
    if ~exist('Database','var'), Database = []; end
    Database = [Database; entry];

    [r2_values] = calc_r2_from_data(r_square_data_gmm_1, r_square_label_gmm, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['traditional | QDA data | ' num2str(size(r_square_data_gmm_1,1)) ' | band ' num2str(bands(1,1)) '-' num2str(bands(1,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all_1, r_square_label_all, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all_1,1)) ' | band ' num2str(bands(1,1)) '-' num2str(bands(1,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_gmm_2, r_square_label_gmm, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_gmm_2,1)) ' | band ' num2str(bands(2,1)) '-' num2str(bands(2,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all_2, r_square_label_all, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all_2,1)) ' | band ' num2str(bands(2,1)) '-' num2str(bands(2,2))]);

end








%% -------------------------------------------------------------------------------
%                            TRADITIONAL METHOD
%  -------------------------------------------------------------------------------
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files Traditional method', 'MultiSelect', 'on');
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
    channels_label =  {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'Fp1', 'CP1', 'CPz', 'CP2', 'Fp2'};

    %% ----------------- load the file -----------------
    disp(['   Loading parameters file: ', filenames{idx_file}(1:end-3) 'yaml'])
    fullpath_file_parameters = [pathname(1:end-4) 'parameters/' filenames{idx_file}(1:end-3) 'yaml'];
    [ringBufferCfg, artifactCfg, processingCfg, gmmCfg, qdaCfg, integratorCfg] = loadParameters(fullpath_file_parameters);

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
    excl_ch = {'Fp1', 'Fp2', 'EOG'};
    [found, indices] = ismember(excl_ch, channels_label);
    excl_chs = indices(found);
    bands = processingCfg.bands; 
    nbands = size(bands, 1);
    signals = cell(1, nbands);
    for idx_band=1:nbands
        band = bands(idx_band,:);
        [signal_processed, header_processed] = processing_onlineROS_CAR_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize, excl_chs);
        signals{idx_band} = signal_processed;
    end
    

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
    disp('   applying fake GMM to all the signal') 

    gmm_prob = zeros(size(signals{1},1),2);
    gmm_prob(:,1) = 1;

    gmmCfg.params.classes = [1, 0];

    %% ----------------- QDA -----------------
    disp('   applying QDA to all the signal') 
    X = [];
    for idx_band=1:nbands
        for idx_band2=1:nbands
            if all(bands(idx_band,:) == qdaCfg.model.bands(idx_band2,:))
                chs = qdaCfg.model.idchans{idx_band2};
                tmp = signals{idx_band}(:,chs); 
            end
        end
        X = [X, tmp];
    end
    qda_prob = apply_qda_matrix(qdaCfg.model, log(X));

    %% ----------------- integrated prob -----------------
    event_start = 781;
    rejection_look_gmm = false;
    [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, gmm_prob, qda_prob, events, event_start, gmmCfg.params.classes, rejection_look_gmm);

    %% ----------------- plot prob integrated -----------------
    ic_index = find(gmmCfg.params.classes == integratorCfg.ic_class_label);
    sampleRate_ros = bufferSize/chunkSize;
    do_plot = false;
    r_square_data_all_1 = []; r_square_data_all_2 = []; r_square_label_all = [];
    r_square_data_gmm_1 = []; r_square_data_gmm_2 = []; r_square_label_gmm = [];
    for idx_trial = 1:ntrial
        start_trial = fixPOS(idx_trial);
        end_trial = cfPOS(idx_trial) + cfDUR(idx_trial)-1;
        c_power_1 = log(signals{1}(start_trial:end_trial,:));
        c_power_2 = log(signals{2}(start_trial:end_trial,:));
        c_gmm_prob = gmm_prob(start_trial:end_trial, :);
        c_artifact = artifact(start_trial:end_trial);
        c_qda_prob = qda_prob(start_trial:end_trial,:);
        c_mask = mask(start_trial:end_trial);
        c_integrated = integrated_prob(start_trial:end_trial,:);
        trial_dur = size(c_power_1, 1);

        % for metrics r^2
        if ismember(cueTYP(idx_trial), [769 770])
            r_square_data_all_1 = [r_square_data_all_1; c_power_1];
            r_square_data_all_2 = [r_square_data_all_2; c_power_2];
            r_square_label_all = [r_square_label_all; repmat(cueTYP(idx_trial), trial_dur, 1)];
            tmp_c_power = c_power_1(c_gmm_prob(:,ic_index) > 0.5,:);
            r_square_data_gmm_1 = [r_square_data_gmm_1; tmp_c_power];
            tmp_c_power = c_power_2(c_gmm_prob(:,ic_index) > 0.5,:);
            r_square_data_gmm_2 = [r_square_data_gmm_2; tmp_c_power];
            r_square_label_gmm = [r_square_label_gmm; repmat(cueTYP(idx_trial), sum(c_gmm_prob(:,ic_index) > integratorCfg.ic_threshold), 1)];
        end

        if do_plot
            figure();
            subplot(511)
            imagesc(c_power_1')
            hold on;
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off;
            yticks(1:nchannels); yticklabels(channels_label)
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            title(['Log band power ' num2str(bands(1,1)) '-' num2str(bands(1,2))])

            subplot(512)
            imagesc(c_power_2')
            hold on;
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off;
            yticks(1:nchannels); yticklabels(channels_label)
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            title(['Log band power ' num2str(bands(2,1)) '-' num2str(bands(2,2))])

            subplot(513)
            plot(c_gmm_prob(:,ic_index))
            hold on;
            plot(c_artifact);
            yline(integratorCfg.ic_threshold, 'LineStyle','--');
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off;
            legend('gmm prob', 'artifact', 'threshold ic');
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            title('artifacts and gmm probabilities')


            subplot(514)
            tmp_prob = c_qda_prob;
            tmp_prob(c_mask == 0,1) = nan;
            plot(c_qda_prob(:,1))
            hold on
            scatter(1:size(c_qda_prob, 1), c_qda_prob(:,1), 15, 'black', 'filled')
            scatter(1:size(c_qda_prob, 1), tmp_prob(:,1), 15, 'green', 'filled')
            yline(0.5, 'LineStyle','--');
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off
            ylim([0 1])
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            legend('qda prob','qda not used', 'qda prob used')
            title('classifier probability')

            subplot(515)
            plot(c_integrated(:,1))
            hold on
            yline(integratorCfg.feedbackThs(1), 'LineStyle','--');
            yline(1-integratorCfg.feedbackThs(2), 'LineStyle','--');
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off
            legend('integrated prob')
            ylim([0 1])
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            title('integrated signal')

            if boom(idx_trial) == 897
                strboom =  'HIT';
            elseif boom(idx_trial) == 898
                strboom = 'MISS';
            elseif boom(idx_trial) == 899
                strboom = 'TIMEOUT';
            else
                disp('ERROR')
            end
            sgtitle(['trial ' num2str(idx_trial) ' | class aked ' num2str(cueTYP(idx_trial)) ' | ' strboom])
        end
    end

    %% take accuracy
    [m] = computeMetrics(integratorCfg, artifact, gmm_prob, qda_prob, integrated_prob, events, event_start, gmmCfg.params.classes, [769 770 783]);
    current_method = 'traditional';
    % Crea entry per il database
    entry = struct();
    entry.File = filenames{idx_file};
    entry.Method = current_method; % 'GMM' o 'Traditional' (impostalo tu nel loop)
    
    % --- ACTIVE TASK METRICS ---
    % Trial Accuracy
    entry.Act_Acc_Raw = m.accuracy.trial.active.raw * 100;
    entry.Act_Acc_NoTimeout = m.accuracy.trial.active.no_timeout * 100;
    
    % Sample Accuracy
    entry.Act_SampAcc_All = m.accuracy.sample.active.all * 100;
    entry.Act_SampAcc_NoTimeout = m.accuracy.sample.active.no_timeout * 100;
    
    % Time (Convertito in secondi)
    % Assicurati che sampleRate_ros sia corretto (es. 16 Hz)
    entry.Act_Time_Hit = m.time.active.hit_avg / sampleRate_ros;
    entry.Act_Time_Miss = m.time.active.miss_avg / sampleRate_ros;
    
    % Stability
    entry.Act_Wobble = m.stability.active.wobble_avg;
    entry.Act_WDR = m.stability.active.wdr_avg * 100;
    
    % --- REST TASK METRICS ---
    entry.Rest_Acc = m.accuracy.trial.rest.acc * 100; % Corretti (Timeout)
    entry.Rest_FPR = m.accuracy.trial.rest.fpr * 100; % False Positives
    entry.Rest_Time_FP = m.time.rest.fp_avg / sampleRate_ros;
    entry.Rest_MaxDev = m.stability.rest.max_dev_avg;

    % Aggiungi al database
    % (Inizializza Database = [] fuori dal loop)
    if ~exist('Database','var'), Database = []; end
    Database = [Database; entry];


    [r2_values] = calc_r2_from_data(r_square_data_gmm_1, r_square_label_gmm, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['traditional | QDA data | ' num2str(size(r_square_data_gmm_1,1)) ' | band ' num2str(bands(1,1)) '-' num2str(bands(1,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all_1, r_square_label_all, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all_1,1)) ' | band ' num2str(bands(1,1)) '-' num2str(bands(1,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_gmm_2, r_square_label_gmm, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_gmm_2,1)) ' | band ' num2str(bands(2,1)) '-' num2str(bands(2,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all_2, r_square_label_all, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all_2,1)) ' | band ' num2str(bands(2,1)) '-' num2str(bands(2,2))]);

end






%% VISUALIZZAZIONE RISULTATI: GMM vs TRADITIONAL
% Richiede la struct 'Database' popolata

methods = {Database.Method};
is_gmm = strcmp(methods, 'gmm');
is_trad = strcmp(methods, 'traditional'); % o 'traditional'

% Colori
c_gmm = [0 0.4470 0.7410];
c_trad = [0.6350 0.0780 0.1840];
grp_colors = [c_gmm; c_trad];

%% FIGURA 1: ACTIVE TASK PERFORMANCE (Accuracy)
figure('Name', 'Active Task Accuracy', 'Color', 'w', 'Position', [100 100 1000 400]);

subplot(1, 2, 1);
data_raw = [[Database(is_gmm).Act_Acc_Raw]', [Database(is_trad).Act_Acc_Raw]'];
boxplot(data_raw, 'Labels', {'GMM', 'Traditional'});
ylabel('Accuracy (%)'); title('Trial Accuracy (Raw - inc. Timeouts)');
grid on;

subplot(1, 2, 2);
data_noto = [[Database(is_gmm).Act_Acc_NoTimeout]', [Database(is_trad).Act_Acc_NoTimeout]'];
boxplot(data_noto, 'Labels', {'GMM', 'Traditional'});
ylabel('Accuracy (%)'); title('Trial Accuracy (No Timeouts)');
grid on;

%% FIGURA 2: ACTIVE TASK TIME & SAMPLE ACCURACY
figure('Name', 'Active Time & Sample Acc', 'Color', 'w', 'Position', [100 550 1000 400]);

subplot(1, 3, 1);
% Time to Hit comparison
data_time = [[Database(is_gmm).Act_Time_Hit]', [Database(is_trad).Act_Time_Hit]'];
boxplot(data_time, 'Labels', {'GMM', 'Traditional'});
ylabel('Time (s)'); title('Avg Time to Hit');
grid on;

subplot(1, 3, 2);
% Sample Accuracy All
data_samp = [[Database(is_gmm).Act_SampAcc_All]', [Database(is_trad).Act_SampAcc_All]'];
boxplot(data_samp, 'Labels', {'GMM', 'Traditional'});
ylabel('Sample Acc (%)'); title('Sample Acc (All Trials)');
grid on;

subplot(1, 3, 3);
% Sample Accuracy No Timeout
data_samp_nt = [[Database(is_gmm).Act_SampAcc_NoTimeout]', [Database(is_trad).Act_SampAcc_NoTimeout]'];
boxplot(data_samp_nt, 'Labels', {'GMM', 'Traditional'});
ylabel('Sample Acc (%)'); title('Sample Acc (Hit/Miss Only)');
grid on;

%% FIGURA 3: ACTIVE STABILITY
figure('Name', 'Active Stability', 'Color', 'w', 'Position', [100 100 800 400]);

subplot(1, 2, 1);
data_wobble = [[Database(is_gmm).Act_Wobble]', [Database(is_trad).Act_Wobble]'];
boxplot(data_wobble, 'Labels', {'GMM', 'Traditional'});
ylabel('Path Length'); title('Wobble (Lower is Better)');
grid on;

subplot(1, 2, 2);
data_wdr = [[Database(is_gmm).Act_WDR]', [Database(is_trad).Act_WDR]'];
boxplot(data_wdr, 'Labels', {'GMM', 'Traditional'});
ylabel('% Wrong Direction'); title('WDR (Lower is Better)');
grid on;

%% FIGURA 4: REST TASK ANALYSIS
% Controlla se ci sono dati di rest (non-NaN)
if any(~isnan([Database.Rest_FPR]))
    figure('Name', 'Rest Task Analysis', 'Color', 'w', 'Position', [100 550 1000 400]);

    subplot(1, 3, 1);
    data_rest_acc = [[Database(is_gmm).Rest_Acc]', [Database(is_trad).Rest_Acc]'];
    boxplot(data_rest_acc, 'Labels', {'GMM', 'Traditional'});
    ylabel('Success Rate (%)'); title('Rest Accuracy (Correct Rejection)');
    grid on;

    subplot(1, 3, 2);
    data_fpr = [[Database(is_gmm).Rest_FPR]', [Database(is_trad).Rest_FPR]'];
    boxplot(data_fpr, 'Labels', {'GMM', 'Traditional'});
    ylabel('FP Rate (%)'); title('False Positives (Lower is Better)');
    grid on;

    subplot(1, 3, 3);
    data_dev = [[Database(is_gmm).Rest_MaxDev]', [Database(is_trad).Rest_MaxDev]'];
    boxplot(data_dev, 'Labels', {'GMM', 'Traditional'});
    ylabel('Deviation'); title('Max Deviation from 0.5');
    grid on;
end
