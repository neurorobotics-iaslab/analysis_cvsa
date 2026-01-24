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
    entry = struct();
    entry.File = filenames{idx_file};
    entry.Method = current_method; 
    
    if ~exist('sampleRate_ros', 'var'), sampleRate_ros = 16; end

    % 1. ACTIVE METRICS
    entry.Act_Trial_Acc = m.act.acc.trial * 100;           % %
    entry.Act_Sample_QDA = m.act.acc.sample_qda * 100;     % %
    entry.Act_GMM_Conf = m.act.gmm.high_conf_ratio * 100;  % %
    entry.Act_Acc_NoTo  = m.act.acc.no_timeout * 100;
    entry.Act_Prog_Score = m.act.prog.timeout_score;       % Score (+1/-1)
    
    entry.Act_Time_Hit = m.act.time.hit / sampleRate_ros;  % Secondi
    
    % 2. REST METRICS
    entry.Rest_Acc = m.rest.acc.trial * 100;               % %
    entry.Rest_Safe_Zone = m.rest.stab.safe_time_ratio * 100; % %
    entry.Rest_Time_Err = m.rest.time.err / sampleRate_ros;   % Secondi

    % Aggiungi al database
    if ~exist('Database','var'), Database = []; end
    Database = [Database; entry];

    [r2_values] = calc_r2_from_data(r_square_data_gmm_1, r_square_label_gmm, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | QDA data | ' num2str(size(r_square_data_gmm_1,1)) ' | band ' num2str(bands(1,1)) '-' num2str(bands(1,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all_1, r_square_label_all, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all_1,1)) ' | band ' num2str(bands(1,1)) '-' num2str(bands(1,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_gmm_2, r_square_label_gmm, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_gmm_2,1)) ' | band ' num2str(bands(2,1)) '-' num2str(bands(2,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all_2, r_square_label_all, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all_2,1)) ' | band ' num2str(bands(2,1)) '-' num2str(bands(2,2))]);

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
    entry = struct();
    entry.File = filenames{idx_file};
    entry.Method = current_method; 
    
    if ~exist('sampleRate_ros', 'var'), sampleRate_ros = 16; end

    % 1. ACTIVE METRICS
    entry.Act_Trial_Acc = m.act.acc.trial * 100;           % %
    entry.Act_Sample_QDA = m.act.acc.sample_qda * 100;     % %
    entry.Act_Acc_NoTo  = m.act.acc.no_timeout * 100;
    entry.Act_GMM_Conf = m.act.gmm.high_conf_ratio * 100;  % %
    entry.Act_Prog_Score = m.act.prog.timeout_score;       % Score (+1/-1)
    
    entry.Act_Time_Hit = m.act.time.hit / sampleRate_ros;  % Secondi
    
    % 2. REST METRICS
    entry.Rest_Acc = m.rest.acc.trial * 100;               % %
    entry.Rest_Safe_Zone = m.rest.stab.safe_time_ratio * 100; % %
    entry.Rest_Time_Err = m.rest.time.err / sampleRate_ros;   % Secondi

    % Aggiungi al database
    if ~exist('Database','var'), Database = []; end
    Database = [Database; entry];

    [r2_values] = calc_r2_from_data(r_square_data_gmm_1, r_square_label_gmm, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | QDA data | ' num2str(size(r_square_data_gmm_1,1)) ' | band ' num2str(bands(1,1)) '-' num2str(bands(1,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all_1, r_square_label_all, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all_1,1)) ' | band ' num2str(bands(1,1)) '-' num2str(bands(1,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_gmm_2, r_square_label_gmm, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_gmm_2,1)) ' | band ' num2str(bands(2,1)) '-' num2str(bands(2,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all_2, r_square_label_all, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all_2,1)) ' | band ' num2str(bands(2,1)) '-' num2str(bands(2,2))]);

end





%% ========================================================================
%  VISUALIZZAZIONE FINALE AGGIORNATA
%  (Eseguire SOLO DOPO aver popolato 'Database' con i nuovi campi)
% ========================================================================

% Setup Metodi
if ~exist('Database', 'var') || isempty(Database)
    error('La variabile Database è vuota o non esiste. Esegui prima i cicli di analisi!');
end

methods = {Database.Method};
is_gmm = strcmpi(methods, 'gmm');
is_trad = strcmpi(methods, 'traditional') | strcmpi(methods, 'trad');

% Colori: Blu (GMM) vs Arancio (Traditional)
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];

% Helper Function Handle
get_d = @(field) extract_data_nan_zero(Database, is_gmm, is_trad, field);

% FIGURA 1: ACTIVE TASK - PERFORMANCE
figure('Name', 'Active Task: Performance', 'Color', 'w', 'Position', [50 50 1200 400]);

% 1. Trial Accuracy (Raw - include i Timeout)
subplot(1, 4, 1);
[d, g] = get_d('Act_Trial_Acc');
custom_boxplot(d, g, colors, 'Accuracy (%)', 'Raw Trial Accuracy (inc. Timeouts)', [0 105]);

% 2. Sample Accuracy (Classificatore)
subplot(1, 4, 2);
[d, g] = get_d('Act_Sample_QDA');
custom_boxplot(d, g, colors, 'Accuracy (%)', 'QDA Sample Accuracy', [50 100]);

% 3. NO-TIMEOUT ACCURACY (Hit / (Hit + Miss))
subplot(1, 4, 3);
% Qui usiamo try-catch nel caso avessi dimenticato di aggiornare il Database
try
    [d, g] = get_d('Act_Acc_NoTo');
    custom_boxplot(d, g, colors, 'Accuracy (%)', 'Trial Acc (No Timeouts)', [0 105]);
catch
    title('Dato Act_Acc_NoTo Mancante!');
end

% 4. Timeout Progress Score
subplot(1, 4, 4);
[d, g] = get_d('Act_Prog_Score');
custom_boxplot(d, g, colors, 'Score (+1/-1)', 'Timeout Direction Quality', [-1.1 1.1]);
yline(0, 'k--', 'LineWidth', 1);

% FIGURA 2: ACTIVE TASK - TIME & QUALITY
figure('Name', 'Active Task: Time & Quality', 'Color', 'w', 'Position', [100 200 800 400]);

% 1. Time to Hit
subplot(1, 2, 1);
[d, g] = get_d('Act_Time_Hit');
custom_boxplot(d, g, colors, 'Time (s)', 'Time to Hit', []);

% 2. Timeout Progress Score (Ripetuto per confronto)
subplot(1, 2, 2);
[d, g] = get_d('Act_Prog_Score');
custom_boxplot(d, g, colors, 'Avg Score', 'Timeout Goodness (+1=Close, -1=Wrong)', [-1.1 1.1]);
yline(0, 'k--', 'LineWidth', 1);

% FIGURA 3: REST TASK - SAFETY
figure('Name', 'Rest Task: Safety', 'Color', 'w', 'Position', [150 350 1200 400]);

% 1. Rest Accuracy
subplot(1, 3, 1);
[d, g] = get_d('Rest_Acc');
custom_boxplot(d, g, colors, 'Accuracy (%)', 'Rest Accuracy (Correct Rejection)', [-5 105]);

% 2. Safe Zone Ratio
subplot(1, 3, 2);
[d, g] = get_d('Rest_Safe_Zone');
custom_boxplot(d, g, colors, '% Time', 'Time in Dead Zone (0.4-0.6)', [-5 105]);

% 3. Time to Error
subplot(1, 3, 3);
[d, g] = get_d('Rest_Time_Err');
custom_boxplot(d, g, colors, 'Time (s)', 'Duration of False Positives', []);


% --- FUNZIONI DI SUPPORTO LOCALI (DEVONO ESSERE ALLA FINE) ---

function [data, groups] = extract_data_nan_zero(db, idx_gmm, idx_trad, field)
    % Verifica esistenza campo
    if ~isfield(db, field)
        error(['Il campo "' field '" non esiste nel Database. Rilancia i cicli di analisi!']);
    end

    d_gmm = [db(idx_gmm).(field)]';
    d_trad = [db(idx_trad).(field)]';
    
    % Gestione NaN -> 0
    d_gmm(isnan(d_gmm)) = 0;
    d_trad(isnan(d_trad)) = 0;
    
    if isempty(d_gmm), d_gmm = []; end
    if isempty(d_trad), d_trad = []; end
    
    data = [d_gmm; d_trad];
    groups = [repmat({'GMM'}, length(d_gmm), 1); repmat({'Traditional'}, length(d_trad), 1)];
end

function custom_boxplot(data, groups, colors, y_label, t_title, y_lims)
    if isempty(data)
        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
        title(t_title); axis off; return;
    end
    
    h = boxplot(data, groups, 'Colors', 'k', 'Symbol', 'o', 'Widths', 0.5);
    set(h, 'LineWidth', 1.2);
    ylabel(y_label, 'FontWeight', 'bold');
    title(t_title, 'FontSize', 10, 'FontWeight', 'normal');
    grid on;
    if ~isempty(y_lims), ylim(y_lims); end
    
    % Logica colorazione robusta
    h_box = findobj(gca, 'Tag', 'Box');
    idx = length(h_box);
    if idx >= 1
        try
            % Se ci sono due box, il primo handle è l'ultimo dato (Traditional)
            if idx > 1
                patch(get(h_box(idx),'XData'), get(h_box(idx),'YData'), colors(1,:), 'FaceAlpha', 0.5); % GMM
                patch(get(h_box(idx-1),'XData'), get(h_box(idx-1),'YData'), colors(2,:), 'FaceAlpha', 0.5); % Trad
            else
                % Se c'è un solo box, colora in base al gruppo
                u_gr = unique(groups);
                if strcmpi(u_gr{1}, 'GMM')
                    patch(get(h_box(1),'XData'), get(h_box(1),'YData'), colors(1,:), 'FaceAlpha', 0.5);
                else
                    patch(get(h_box(1),'XData'), get(h_box(1),'YData'), colors(2,:), 'FaceAlpha', 0.5);
                end
            end
        catch
            % Fallback
        end
    end
end





