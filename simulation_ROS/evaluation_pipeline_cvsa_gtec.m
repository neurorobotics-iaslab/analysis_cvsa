% load of gmm, qda and for each trial shows the results
clear all; % close all;

addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/equal_ros')
addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/utils')
addpath(genpath('/home/paolo/Local/Matlab/yamlmatlab'));


% --- inizialization ---
nchannels = 16;
classes = [730 731 783];   

[filenames_my, pathname_my] = uigetfile('*.gdf', 'Select GDF Files My method', 'MultiSelect', 'on');
if ischar(filenames_my)
    filenames_my = {filenames_my};
end
nFiles_my = length(filenames_my);

[filenames_trad, pathname_trad] = uigetfile('*.gdf', 'Select GDF Files Traditional method', 'MultiSelect', 'on');
if ischar(filenames_trad)
    filenames_trad = {filenames_trad};
end
nFiles_trad = length(filenames_trad);

Database = [];

%% -------------------------------------------------------------------------------
%                                    MY METHOD
%  -------------------------------------------------------------------------------
for idx_file = 1:nFiles_my
    fullpath_file_gdf = fullfile(pathname_my, filenames_my{idx_file});
    disp(['File ' num2str(idx_file) '/' num2str(nFiles_my)]);
    disp(['   Loading gdf file : ', filenames_my{idx_file}]);
    [c_signal,header] = sload(fullpath_file_gdf);
    c_signal = c_signal(:,1:nchannels);
    channels_label =  {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'Fp1', 'CP1', 'CPz', 'CP2', 'Fp2'};

    %% ----------------- load the file for the my method -----------------
    disp(['   Loading parameters file: ', filenames_my{idx_file}(1:end-3) 'yaml'])
    fullpath_file_parameters = [pathname_my(1:end-4) 'parameters/' filenames_my{idx_file}(1:end-3) 'yaml'];
    [ringBufferCfg, artifactCfg, processingCfg, gmmCfg, qdaCfg, integratorCfg] = loadParameters(fullpath_file_parameters);

    disp(['   Loading GMM file: ', gmmCfg.file_name])
    gmm_path = [pathname_my(1:end-4) gmmCfg.file_name];
    [gmmCfg.model, gmmCfg.params] = loadGMM(gmm_path);

    disp(['   Loading QDA file (my): ', qdaCfg.file_name])
    qda_path = [pathname_my(1:end-4) qdaCfg.file_name];
    qdaCfg.model = loadQDA(qda_path);

    % --- load the qda for teh traditional method ---
    k_gain_trad = nan(1, nFiles_trad);
disp('   Loading trad files to have the QDA');
    for idx_file_t = 1:nFiles_trad
        fullpath_file_parameters = [pathname_trad(1:end-4) 'parameters/' filenames_trad{idx_file_t}(1:end-3) 'yaml'];
        [~, ~, ~, ~, qdaCfg_trad, integratorCfg_trad] = loadParameters(fullpath_file_parameters);
        k_gain_trad(idx_file_t) = integratorCfg_trad.k_gain;
    end
    integratorCfg_trad.k_gain = mean(k_gain_trad);

    disp(['   Loading QDA file (trad): ', qdaCfg_trad.file_name])
    qda_path = [pathname_trad(1:end-4) qdaCfg_trad.file_name];
    qdaCfg_trad.model = loadQDA(qda_path);

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
    band = processingCfg.bands; 
    [signals, header_processed] = processing_onlineROS_CAR_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize, excl_chs);
    

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
    o_l = cell2mat(gmmCfg.params.occipital_left_idx);
    o_r = cell2mat(gmmCfg.params.occipital_right_idx);
    excl_chs = cell2mat(gmmCfg.params.excluded_idx);

    type = gmmCfg.params.type;
    nfeatures = gmmCfg.params.nfeatures;

    % features extraction and classification
    sparsity = nan(size(signals, 1), nfeatures);
    for idx_sample = 1:size(signals,1)
        sparsity(idx_sample,:) = compute_features_icnic_cvsa(signals(idx_sample,:), type, o_l, o_r, c_l, c_r, nfeatures);
    end

    data_standardized = (sparsity - mu_data) ./ sigma_data;

    gmm_prob = posterior(gmmCfg.model, data_standardized);

    %% ----------------- QDA -----------------
    disp('   applying QDA to all the signal') 
    X = [];
    if all(band == qdaCfg.model.bands)
        chs = qdaCfg.model.idchans{1};
        X = signals(:, chs);
    else
        disp('ERROR in the band!')
    end
    qda_prob = apply_qda_matrix(qdaCfg.model, log(X)); 

    % --- qda trad ---
    X_trad = [];
    if all(band == qdaCfg_trad.model.bands)
        chs = qdaCfg_trad.model.idchans{1};
        X_trad = signals(:, chs);
    else
        disp('ERROR in the band!')
    end
    qda_prob_trad = apply_qda_matrix(qdaCfg_trad.model, log(X_trad));

    %% ----------------- integrated prob -----------------
    event_start = 781;
    rejection_look_gmm = false;
    [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, gmm_prob, qda_prob, events, event_start, cell2mat(gmmCfg.params.classes), rejection_look_gmm);

    % --- integrated prob simulation of traditional ---
    gmm_prob_fake = zeros(size(signals,1),2);
    gmm_prob_fake(:,1) = 1;
    [integrated_prob_fake, mask_fake] = applyIntegration(integratorCfg_trad, artifact, gmm_prob_fake, qda_prob_trad, events, event_start, [1, 0], rejection_look_gmm);
    [curve_rest, curve_act, curve_act_noTimeout, n_act_timeout] = calc_pareto_matrix(integrated_prob_fake, events, classes, event_start);
    
    %% ----------------- plot prob integrated -----------------
    ic_index = find(cell2mat(gmmCfg.params.classes) == integratorCfg.ic_class_label);
    sampleRate_ros = bufferSize/chunkSize;
    do_plot = true;
    r_square_data_all = []; r_square_label_all = [];
    r_square_data_gmm = []; r_square_label_gmm = [];
    for idx_trial = 1:ntrial
        start_trial = fixPOS(idx_trial);
        end_trial = cfPOS(idx_trial) + cfDUR(idx_trial)-1;
        c_power = log(signals(start_trial:end_trial,:));
        c_gmm_prob = gmm_prob(start_trial:end_trial, :);
        c_artifact = artifact(start_trial:end_trial);
        c_qda_prob = qda_prob(start_trial:end_trial,:);
        c_mask = mask(start_trial:end_trial);
        c_integrated = integrated_prob(start_trial:end_trial,:);
        trial_dur = size(c_power, 1);

        % for metrics r^2
        if ismember(cueTYP(idx_trial), [classes(1) classes(2)])
            r_square_data_all = [r_square_data_all; c_power];
            r_square_label_all = [r_square_label_all; repmat(cueTYP(idx_trial), trial_dur, 1)];

            tmp_c_power = c_power(c_gmm_prob(:,ic_index) > 0.5,:);
            r_square_data_gmm = [r_square_data_gmm; tmp_c_power];
            r_square_label_gmm = [r_square_label_gmm; repmat(cueTYP(idx_trial), sum(c_gmm_prob(:,ic_index) > 0.5), 1)];
        end

        if do_plot
            figure();
            subplot(411)
            imagesc(c_power')
            hold on;
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off;
            yticks(1:nchannels); yticklabels(channels_label)
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            title(['Log band power ' num2str(band(1)) '-' num2str(band(2))])

            subplot(412)
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


            subplot(413)
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

            subplot(414)
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
    [m] = computeMetrics(integratorCfg, artifact, gmm_prob, qda_prob, integrated_prob, events, event_start, cell2mat(gmmCfg.params.classes), classes);
    current_method = 'gmm';
    entry = struct();
    entry.File = filenames_my{idx_file};
    entry.Method = current_method; 
    
    if ~exist('sampleRate_ros', 'var'), sampleRate_ros = 16; end

    % ACTIVE METRICS
    entry.Act_Trial_Acc = m.act.acc.trial * 100;           % %
    entry.Act_Sample_QDA = m.act.acc.sample_qda * 100;     % %
    entry.Act_GMM_Conf = m.act.gmm.high_conf_ratio * 100;  % %
    entry.Act_Acc_NoTo  = m.act.acc.no_timeout * 100;
    entry.Act_Prog_Score = m.act.prog.timeout_score;       % Score (+1/-1)
    entry.Act_Time_Hit = m.act.time.hit / sampleRate_ros;  % Secondi
    entry.Act_Num_Timeout = m.act.num.timeout;
    
    % REST METRICS
    entry.Rest_Acc = m.rest.acc.trial * 100;               % %
    entry.Rest_Safe_Zone = m.rest.stab.safe_time_ratio * 100; % %
    entry.Rest_Time_Err = m.rest.time.err / sampleRate_ros;   % Secondi

    % OTHER METRICS
    entry.GMM_Avg_Active = m.gmm.avg_active; % Dovrebbe essere ALTO
    entry.GMM_Avg_Rest   = m.gmm.avg_rest;   % Dovrebbe essere BASSO
    entry.GMM_AUC        = m.gmm.auc;
    entry.Pareto_Curve_Rest = curve_rest; 
    entry.Pareto_Curve_Act  = curve_act;
    entry.Pareto_Curve_Act_NoTo = curve_act_noTimeout;
    entry.Pareto_Num_Timeouts = n_act_timeout;

    % CONTROL SNR
    entry.SNR_Active_Vel = m.snr.vel_active * 1000; % Scaling per leggibilità (es. 1e-3 -> 1)
    entry.SNR_Rest_Vel   = m.snr.vel_rest * 1000;
    entry.SNR_Ratio      = m.snr.ratio;

    if ~exist('Database','var'), Database = []; end
    Database = [Database; entry];

    [r2_values] = calc_r2_from_data(r_square_data_gmm, r_square_label_gmm, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | QDA data | ' num2str(size(r_square_data_gmm,1)) ' | band ' num2str(band(1)) '-' num2str(band(2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all, r_square_label_all, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all,1)) ' | band ' num2str(band(1)) '-' num2str(band(2))]);

end








%% -------------------------------------------------------------------------------
%                            TRADITIONAL METHOD
%  -------------------------------------------------------------------------------
for idx_file = 1:nFiles_trad
    fullpath_file_gdf = fullfile(pathname_trad, filenames_trad{idx_file});
    disp(['File ' num2str(idx_file) '/' num2str(nFiles_trad)]);
    disp(['   Loading gdf file : ', filenames_trad{idx_file}]);
    [c_signal,header] = sload(fullpath_file_gdf);
    c_signal = c_signal(:,1:nchannels);
    channels_label =  {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'Fp1', 'CP1', 'CPz', 'CP2', 'Fp2'};

    %% ----------------- load the file -----------------
    disp(['   Loading parameters file: ', filenames_trad{idx_file}(1:end-3) 'yaml'])
    fullpath_file_parameters = [pathname_trad(1:end-4) 'parameters/' filenames_trad{idx_file}(1:end-3) 'yaml'];
    [ringBufferCfg, artifactCfg, processingCfg, gmmCfg, qdaCfg, integratorCfg] = loadParameters(fullpath_file_parameters);

    disp(['   Loading QDA file: ', qdaCfg.file_name])
    qda_path = [pathname_trad(1:end-4) qdaCfg.file_name];
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
    band = processingCfg.bands;
    [signals, header_processed] = processing_onlineROS_CAR_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize, excl_chs);
    

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

    gmm_prob = zeros(size(signals,1),2);
    gmm_prob(:,1) = 1;

    gmmCfg.params.classes = [1, 0];

    %% ----------------- QDA -----------------
    disp('   applying QDA to all the signal') 
    X = [];
    if all(band == qdaCfg.model.bands)
        chs = qdaCfg.model.idchans{1};
        X = signals(:, chs);
    else
        disp('ERROR in the band!')
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
    r_square_data_all = []; r_square_data_all_2 = []; r_square_label_all = [];
    r_square_data_gmm = []; r_square_data_gmm_2 = []; r_square_label_gmm = [];
    for idx_trial = 1:ntrial
        start_trial = fixPOS(idx_trial);
        end_trial = cfPOS(idx_trial) + cfDUR(idx_trial)-1;
        c_power = log(signals(start_trial:end_trial,:));
        c_gmm_prob = gmm_prob(start_trial:end_trial, :);
        c_artifact = artifact(start_trial:end_trial);
        c_qda_prob = qda_prob(start_trial:end_trial,:);
        c_mask = mask(start_trial:end_trial);
        c_integrated = integrated_prob(start_trial:end_trial,:);
        trial_dur = size(c_power, 1);

        % for metrics r^2
        if ismember(cueTYP(idx_trial), [classes(1) classes(2)])
            r_square_data_all = [r_square_data_all; c_power];
            r_square_label_all = [r_square_label_all; repmat(cueTYP(idx_trial), trial_dur, 1)];
            tmp_c_power = c_power(c_gmm_prob(:,ic_index) > 0.5,:);
            r_square_data_gmm = [r_square_data_gmm; tmp_c_power];
            r_square_label_gmm = [r_square_label_gmm; repmat(cueTYP(idx_trial), sum(c_gmm_prob(:,ic_index) > 0.5), 1)];
        end

        if do_plot
            figure();
            subplot(411)
            imagesc(c_power')
            hold on;
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off;
            yticks(1:nchannels); yticklabels(channels_label)
            xlim([1 trial_dur])
            xticks(0:sampleRate_ros:trial_dur);
            xticklabels((0:sampleRate_ros:trial_dur)/sampleRate_ros)
            title(['Log band power ' num2str(band(1,1)) '-' num2str(band(1,2))])

            subplot(412)
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


            subplot(413)
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

            subplot(414)
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
    [m] = computeMetrics(integratorCfg, artifact, gmm_prob, qda_prob, integrated_prob, events, event_start, gmmCfg.params.classes, classes);

    current_method = 'traditional';
    entry = struct();
    entry.File = filenames_trad{idx_file};
    entry.Method = current_method; 
    
    if ~exist('sampleRate_ros', 'var'), sampleRate_ros = 16; end

    % ACTIVE METRICS
    entry.Act_Trial_Acc = m.act.acc.trial * 100;           % %
    entry.Act_Sample_QDA = m.act.acc.sample_qda * 100;     % %
    entry.Act_Acc_NoTo  = m.act.acc.no_timeout * 100;
    entry.Act_GMM_Conf = m.act.gmm.high_conf_ratio * 100;  % %
    entry.Act_Prog_Score = m.act.prog.timeout_score;       % Score (+1/-1)
    entry.Act_Time_Hit = m.act.time.hit / sampleRate_ros;  % Secondi
    entry.Act_Num_Timeout = m.act.num.timeout;
    
    % REST METRICS
    entry.Rest_Acc = m.rest.acc.trial * 100;               % %
    entry.Rest_Safe_Zone = m.rest.stab.safe_time_ratio * 100; % %
    entry.Rest_Time_Err = m.rest.time.err / sampleRate_ros;   % Secondi

    % OTHER METRICS
    entry.GMM_Avg_Active = m.gmm.avg_active; % Dovrebbe essere ALTO
    entry.GMM_Avg_Rest   = m.gmm.avg_rest;   % Dovrebbe essere BASSO
    entry.GMM_AUC        = m.gmm.auc;
    entry.Pareto_Curve_Rest = ''; 
    entry.Pareto_Curve_Act  = ''; 
    entry.Pareto_Curve_Act_NoTo = '';
    entry.Pareto_Num_Timeouts = '';

    % CONTROL SNR
    entry.SNR_Active_Vel = m.snr.vel_active * 1000; % Scaling per leggibilità (es. 1e-3 -> 1)
    entry.SNR_Rest_Vel   = m.snr.vel_rest * 1000;
    entry.SNR_Ratio      = m.snr.ratio;

    if ~exist('Database','var'), Database = []; end
    Database = [Database; entry];

    [r2_values] = calc_r2_from_data(r_square_data_gmm, r_square_label_gmm, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | QDA data | ' num2str(size(r_square_data_gmm,1)) ' | band ' num2str(band(1)) '-' num2str(band(2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all, r_square_label_all, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all,1)) ' | band ' num2str(band(1)) '-' num2str(band(2))]);

end





%% ------------ PLOTS ------------
if ~exist('Database', 'var') || isempty(Database)
    error('ERRORE: La variabile "Database" è vuota o non esiste. Esegui prima i cicli di analisi!');
end

methods = {Database.Method};
is_gmm = strcmpi(methods, 'gmm');
is_trad = strcmpi(methods, 'traditional') | strcmpi(methods, 'trad');
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
get_d = @(field) extract_data_nan_zero(Database, is_gmm, is_trad, field);


% FIGURA 1: ACTIVE TASK - PERFORMANCE
figure('Name', 'Fig 1: Active Task Performance', 'Color', 'w', 'Position', [50 50 1400 400]);

subplot(1, 3, 1);
[d, g] = get_d('Act_Trial_Acc');
custom_boxplot(d, g, colors, 'Accuracy (%)', 'Raw Trial Acc (inc. Timeouts)', [0 105]);

subplot(1, 3, 2);
[d, g] = get_d('Act_Sample_QDA');
custom_boxplot(d, g, colors, 'Accuracy (%)', 'QDA Sample Accuracy', [50 100]);

subplot(1, 3, 3);
try
    [d, g] = get_d('Act_Acc_NoTo');
    custom_boxplot(d, g, colors, 'Accuracy (%)', 'Trial Acc (Decisions Only)', [0 105]);
catch
    text(0.5,0.5,'Dato Mancante','HorizontalAlignment','center');
end


% FIGURA 2: ACTIVE TASK - TIME & QUALITY
figure('Name', 'Fig 2: Time vs Quality', 'Color', 'w', 'Position', [100 200 800 400]);

subplot(1, 2, 1);
[d, g] = get_d('Act_Time_Hit');
custom_boxplot(d, g, colors, 'Time (s)', 'Time to Hit', []);

subplot(1, 2, 2);
[d, g] = get_d('Act_Prog_Score');
custom_boxplot(d, g, colors, 'Score (+1/-1)', 'Trials Timeout Goodness', [-1.1 1.1]);
yline(0, 'k--', 'LineWidth', 1);
subtitle('+1 = Near Target, -1 = Wrong Side');


% FIGURA 3: REST TASK - SAFETY
figure('Name', 'Fig 3: Rest Safety', 'Color', 'w', 'Position', [150 350 1200 400]);

subplot(1, 3, 1);
[d, g] = get_d('Rest_Acc');
custom_boxplot(d, g, colors, 'Accuracy (%)', 'Rest Accuracy', [-5 105]);

subplot(1, 3, 2);
[d, g] = get_d('Rest_Safe_Zone');
custom_boxplot(d, g, colors, '% Time', 'Time in Dead Zone (0.4-0.6)', [-5 105]);

subplot(1, 3, 3);
[d, g] = get_d('Rest_Time_Err');
custom_boxplot(d, g, colors, 'Time (s)', 'Duration of False Positives', []);


% FIGURA 4: GMM INTENT DETECTION
db_gmm = Database(is_gmm);

if ~isempty(db_gmm) && isfield(db_gmm, 'GMM_AUC')
    figure('Name', 'Fig 4: GMM Intent Detection', 'Color', 'w', 'Position', [100 100 800 500]);
    
    data_act = [db_gmm.GMM_Avg_Active]';
    data_rst = [db_gmm.GMM_Avg_Rest]';
    data_auc = [db_gmm.GMM_AUC]';

    subplot(1, 2, 1);
    boxplot([data_act, data_rst], 'Labels', {'Active Tasks', 'Rest Task'}, 'Colors', 'k', 'Symbol', 'o');
    title('GMM Confidence Levels', 'FontSize', 12);
    ylabel('Avg Probability (0-1)', 'FontWeight', 'bold');
    subtitle('Active should be > Rest', 'FontSize', 10);
    grid on;
    h = findobj(gca, 'Tag', 'Box');
    if length(h) >= 2
        patch(get(h(2),'XData'), get(h(2),'YData'), [0 0.4470 0.7410], 'FaceAlpha', 0.5); % Blu
        patch(get(h(1),'XData'), get(h(1),'YData'), [0.5 0.5 0.5], 'FaceAlpha', 0.5);    % Grigio
    end
    
    subplot(1, 2, 2);
    boxplot(data_auc, 'Labels', {'Active vs Rest'}, 'Colors', 'k', 'Symbol', 'o');
    title('Discrimination Ability (AUC)', 'FontSize', 12);
    ylabel('AUC Score', 'FontWeight', 'bold');
    yline(0.5, 'r--', 'Chance Level');
    yline(0.7, 'g--', 'Good Separation');
    ylim([0.4 1.05]); grid on;
    h = findobj(gca, 'Tag', 'Box');
    if ~isempty(h)
        patch(get(h(1),'XData'), get(h(1),'YData'), [0.4660 0.6740 0.1880], 'FaceAlpha', 0.5); % Verde
    end
else
    warning('Dati GMM (AUC/Confidenza) non trovati o incompleti. Salto Figura 4.');
end


% FIGURA 5: CONTROL SNR
if isfield(Database, 'SNR_Ratio')
    figure('Name', 'Fig 5: Control Signal Quality (SNR)', 'Color', 'w', 'Position', [100 100 1000 400]);

    subplot(1, 3, 1);
    [d, g] = get_d('SNR_Active_Vel');
    custom_boxplot(d, g, colors, 'Velocity (a.u.)', 'Active Tasks', []);
    subtitle('Higher is Faster');

    subplot(1, 3, 2);
    [d, g] = get_d('SNR_Rest_Vel');
    custom_boxplot(d, g, colors, 'Velocity (a.u.)', 'Rest Task', []);
    subtitle('Lower is Better');

    subplot(1, 3, 3);
    [d, g] = get_d('SNR_Ratio');
    custom_boxplot(d, g, colors, 'SNR (Act/Rest)', 'Control Quality Ratio', []);
    yline(1, 'k--');
    subtitle('Higher is Better (>1 means Control)');
else
    warning('Dati SNR non trovati. Salto Figura 5.');
end


% FIGURE 6: VISUALIZZAZIONE "THRESHOLD TRENDS"
db_target = Database(is_gmm);
if isempty(db_target), error('Nessun file GMM trovato nel Database'); end

TOTAL_TRIALS = sum(ismember(cueTYP, [classes(1) classes(2)])); 
raw_rest_3d     = cat(3, db_target.Pareto_Curve_Rest);
raw_act_3d      = cat(3, db_target.Pareto_Curve_Act); 
raw_act_noto_3d = cat(3, db_target.Pareto_Curve_Act_NoTo);
raw_timeouts_3d = cat(3, db_target.Pareto_Num_Timeouts);
mean_rest_mat     = mean(raw_rest_3d, 3, 'omitnan') * 100;
mean_act_mat      = mean(raw_act_3d, 3, 'omitnan') * 100;      
mean_act_noto_mat = mean(raw_act_noto_3d, 3, 'omitnan') * 100; 
mean_timeouts_abs = mean(raw_timeouts_3d, 3, 'omitnan');
mean_timeouts_perc_mat = (mean_timeouts_abs / TOTAL_TRIALS) * 100;

gmm_act_val  = mean([db_target.Act_Trial_Acc], 'omitnan'); 
gmm_rest_val = mean([db_target.Rest_Acc], 'omitnan');
gmm_to_abs   = mean([db_target.Act_Num_Timeout], 'omitnan');
gmm_to_perc  = (gmm_to_abs / TOTAL_TRIALS) * 100;

figure('Name', 'Analysis vs GMM', 'Color', 'w', 'Position', [50 50 1600 550]);
x_axis = linspace(0.55, 1.0, size(mean_rest_mat, 2)); 
thresholds = x_axis; 

subplot(1, 3, 1);
imagesc(thresholds, thresholds, mean_act_mat); 
axis xy; colormap(gca, 'jet'); 
c = colorbar; c.Label.String = 'Accuracy [%]';
hold on;
plot([min(thresholds) max(thresholds)], [min(thresholds) max(thresholds)], 'k:', 'LineWidth', 1.5);
[~, ~] = contour(thresholds, thresholds, mean_act_mat, [gmm_act_val gmm_act_val], 'm-', 'LineWidth', 3); 
[r_good, c_good] = find(mean_act_mat >= gmm_act_val);
if ~isempty(r_good) && length(r_good) < 5 
    plot(thresholds(c_good), thresholds(r_good), 'm.', 'MarkerSize', 15);
end
xlabel(['Threshold ' num2str(classes(1))], 'FontSize', 11, 'FontWeight','bold');
ylabel(['Threshold ' num2str(classes(2))], 'FontSize', 11, 'FontWeight','bold');
title({'\bf ACTIVE Accuracy', sprintf('GMM Baseline: %.1f%%', gmm_act_val)}, 'FontSize', 14);

subplot(1, 3, 2);
imagesc(thresholds, thresholds, mean_rest_mat);
axis xy; colormap(gca, 'jet');
c = colorbar; c.Label.String = 'Accuracy [%]';
hold on;
plot([min(thresholds) max(thresholds)], [min(thresholds) max(thresholds)], 'k:', 'LineWidth', 1.5);
[~, ~] = contour(thresholds, thresholds, mean_rest_mat, [gmm_rest_val gmm_rest_val], 'm-', 'LineWidth', 3);
[r_good, c_good] = find(mean_rest_mat >= gmm_rest_val);
if ~isempty(r_good) && length(r_good) < 5
    plot(thresholds(c_good), thresholds(r_good), 'm.', 'MarkerSize', 15);
end
xlabel(['Threshold ' num2str(classes(1))], 'FontSize', 11, 'FontWeight','bold');
ylabel(['Threshold ' num2str(classes(2))], 'FontSize', 11, 'FontWeight','bold');
title({'\bf REST Accuracy', sprintf('GMM Baseline: %.1f%%', gmm_rest_val)}, 'FontSize', 14);

subplot(1, 3, 3);
imagesc(thresholds, thresholds, mean_timeouts_perc_mat);
axis xy; 
colormap(gca, flipud(parula)); 
c = colorbar; c.Label.String = 'Timeout Rate [%]';
caxis([0 50]); 
hold on;
plot([min(thresholds) max(thresholds)], [min(thresholds) max(thresholds)], 'k:', 'LineWidth', 1.5);
[~, ~] = contour(thresholds, thresholds, mean_timeouts_perc_mat, [gmm_to_perc gmm_to_perc], 'm-', 'LineWidth', 3); 
[r_good, c_good] = find(mean_timeouts_perc_mat <= gmm_to_perc);
if ~isempty(r_good) && length(r_good) < 5
    plot(thresholds(c_good), thresholds(r_good), 'm.', 'MarkerSize', 15);
end
xlabel(['Threshold ' num2str(classes(1))], 'FontSize', 11, 'FontWeight','bold');
ylabel(['Threshold ' num2str(classes(2))], 'FontSize', 11, 'FontWeight','bold');
title({'\bf TIMEOUT Rate', sprintf('GMM Value: %.1f%%', gmm_to_perc)}, 'FontSize', 14);







%% --- FUNZIONI DI SUPPORTO ---

function [data, groups] = extract_data_nan_zero(db, idx_gmm, idx_trad, field)
    % Verifica esistenza campo
    if ~isfield(db, field)
        warning(['Campo "' field '" mancante nel Database.']);
        data = []; groups = []; return;
    end
    
    d_gmm = [db(idx_gmm).(field)]';
    d_trad = [db(idx_trad).(field)]';
    
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
    
    % Gestione dinamica dei patch colorati
    if idx >= 1
        try
            % Cerca di colorare in base ai gruppi presenti
            group_names = unique(groups);
            has_gmm = any(strcmpi(group_names, 'GMM'));
            has_trad = any(strcmpi(group_names, 'Traditional'));
            
            % Nota: Boxplot disegna da destra a sinistra o viceversa a seconda della versione
            % Assegnazione colori semplice basata sugli handle
            if idx == 2
                 % Assumiamo ordine standard: ultimo handle = primo gruppo
                 patch(get(h_box(2),'XData'), get(h_box(2),'YData'), colors(1,:), 'FaceAlpha', 0.5); % GMM
                 patch(get(h_box(1),'XData'), get(h_box(1),'YData'), colors(2,:), 'FaceAlpha', 0.5); % Trad
            elseif idx == 1
                if has_gmm
                    patch(get(h_box(1),'XData'), get(h_box(1),'YData'), colors(1,:), 'FaceAlpha', 0.5);
                else
                    patch(get(h_box(1),'XData'), get(h_box(1),'YData'), colors(2,:), 'FaceAlpha', 0.5);
                end
            end
        catch
            % Fallback silenzioso
        end
    end
end

function [mat_rest, mat_act, mat_act_noTimeout, mat_act_timeout] = calc_pareto_matrix(integrated_prob, events, task_classes, start_task)
    % INPUT:
    % integrated_prob: vettore (o matrice) delle probabilità integrate. 
    %                  Assumiamo Colonna 1: 0=DX, 1=SX.
    % events: struct con .TYP, .POS, .DUR
    % task_classes: [CODE_SX, CODE_DX, CODE_REST]
    % start_task: codice evento inizio trial (es. 781)
    
    % Definiamo il range delle threshold (Confidenza)
    % Nota: Questa è la "forza" della soglia. 
    % Per SX (High) la soglia reale è val.
    % Per DX (Low)  la soglia reale è (1 - val).
    thresholds = 0.55 : 0.05 : 1.0; 
    
    n_th = length(thresholds);
    
    % OUTPUT: Ora sono MATRICI (n_th x n_th)
    % Riga (i) -> Soglia SX
    % Colonna (j) -> Soglia DX
    mat_rest            = zeros(n_th, n_th);
    mat_act             = zeros(n_th, n_th);
    mat_act_noTimeout   = zeros(n_th, n_th);
    mat_act_timeout     = zeros(n_th, n_th);
    
    CODE_SX   = task_classes(1); % Target: > Soglia Alta
    CODE_DX   = task_classes(2); % Target: < Soglia Bassa
    CODE_REST = task_classes(3); % Target: Nessun attraversamento
    
    % Estrazione eventi
    cfPOS = events.POS(events.TYP == start_task);
    cfDUR = events.DUR(events.TYP == start_task);
    cues = events.TYP(ismember(events.TYP, task_classes));
    
    ntrials = length(cfPOS);
    if length(cues) ~= ntrials
        warning('Disallineamento Cues/Start events. Controllo indici.');
        ntrials = min(length(cues), length(cfPOS));
    end
    
    % --- DOPPIO LOOP PER SOGLIE INDIPENDENTI ---
    
    % Loop 1: Varia la soglia per SX (Soglia Alta)
    for i_sx = 1:n_th
        th_val_sx = thresholds(i_sx); 
        th_high   = th_val_sx;       % Es. 0.70
        
        % Loop 2: Varia la soglia per DX (Soglia Bassa)
        for i_dx = 1:n_th
            th_val_dx = thresholds(i_dx);
            th_low    = 1.0 - th_val_dx; % Es. se val=0.60 -> th_low = 0.40
            
            % Contatori per questa specifica combinazione (SX=i, DX=j)
            cnt_act_hit = 0;
            cnt_act_timeout = 0;
            cnt_act_tot = 0;
            cnt_rest_ok = 0;
            cnt_rest_tot = 0;
            
            % Loop sui Trials
            for i_tr = 1:ntrials
                idx_s = cfPOS(i_tr);
                idx_e = idx_s + cfDUR(i_tr) - 1;
                
                % Protezione indici array
                if idx_e > length(integrated_prob)
                    idx_e = length(integrated_prob);
                end
                
                sig = integrated_prob(idx_s:idx_e, 1); 
                current_cue = cues(i_tr);
                
                % Cerca il primo attraversamento per SX (High)
                idx_cross_high = find(sig >= th_high, 1, 'first');
                
                % Cerca il primo attraversamento per DX (Low)
                idx_cross_low  = find(sig <= th_low, 1, 'first');
                
                res_high = ~isempty(idx_cross_high);
                res_low  = ~isempty(idx_cross_low);
                
                winner = 'none';
                
                % Logica per determinare chi vince
                if res_high && ~res_low
                    winner = 'high';
                elseif ~res_high && res_low
                    winner = 'low';
                elseif res_high && res_low
                    % Entrambe le soglie superate: vince chi accade prima
                    if idx_cross_high < idx_cross_low
                        winner = 'high';
                    else
                        winner = 'low';
                    end
                end
                
                % --- Valutazione Hit/Miss/Timeout ---
                if current_cue == CODE_SX || current_cue == CODE_DX
                    cnt_act_tot = cnt_act_tot + 1;
                    
                    if current_cue == CODE_SX
                        if strcmp(winner, 'high')
                            cnt_act_hit = cnt_act_hit + 1;
                        end
                    elseif current_cue == CODE_DX
                        if strcmp(winner, 'low')
                            cnt_act_hit = cnt_act_hit + 1;
                        end
                    end
                    
                    if strcmp(winner, 'none')
                        cnt_act_timeout = cnt_act_timeout + 1;
                    end
                    
                elseif current_cue == CODE_REST
                    cnt_rest_tot = cnt_rest_tot + 1;
                    
                    if strcmp(winner, 'none')
                        cnt_rest_ok = cnt_rest_ok + 1;
                    end
                end
            end
            
            % --- Assegnazione valori nella MATRICE (i_sx, i_dx) ---
            
            if cnt_act_tot > 0
                mat_act(i_sx, i_dx) = cnt_act_hit / cnt_act_tot;
                
                den_noTime = cnt_act_tot - cnt_act_timeout;
                if den_noTime > 0
                    mat_act_noTimeout(i_sx, i_dx) = cnt_act_hit / den_noTime;
                else
                    mat_act_noTimeout(i_sx, i_dx) = 0; % O NaN, a preferenza
                end
            end
            
            if cnt_rest_tot > 0
                mat_rest(i_sx, i_dx) = cnt_rest_ok / cnt_rest_tot;
            end
            
            mat_act_timeout(i_sx, i_dx) = cnt_act_timeout;
        end
    end
end