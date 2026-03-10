% load of gmm, qda and for each trial shows the results. Save the dataset
% for the anlisis
clear all; % close all;

addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/equal_ros')
addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/utils')
addpath(genpath('/home/paolo/Local/Matlab/yamlmatlab'));


% --- inizialization ---
nchannels = 16;
classes = [771 773 783];   

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
    for idx_file_t = 1:nFiles_trad
        fullpath_file_parameters = [pathname_trad(1:end-4) 'parameters/' filenames_trad{idx_file_t}(1:end-3) 'yaml'];
        [~, ~, ~, ~, qdaCfg_trad, ~] = loadParameters(fullpath_file_parameters);
    end

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

    % --- qda trad ---
    X_trad = [];
    for idx_band=1:nbands
        for idx_band2=1:nbands
            if all(bands(idx_band,:) == qdaCfg_trad.model.bands(idx_band2,:))
                chs = qdaCfg_trad.model.idchans{idx_band2};
                tmp = signals{idx_band}(:,chs); 
            end
        end
        X_trad = [X_trad, tmp];
    end
    qda_prob_trad = apply_qda_matrix(qdaCfg_trad.model, log(X_trad));

    %% ----------------- integrated prob -----------------
    event_start = 781;
    rejection_look_gmm = false;
    [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, gmm_prob, qda_prob, events, event_start, cell2mat(gmmCfg.params.classes), rejection_look_gmm);

    % --- integrated prob simulation of traditional ---
    gmm_prob_fake = zeros(size(signals{1},1),2);
    gmm_prob_fake(:,1) = 1;
    [integrated_prob_fake, ~] = applyIntegration(integratorCfg, artifact, gmm_prob_fake, qda_prob_trad, events, event_start, [1, 0], rejection_look_gmm);
    [curve_rest, curve_act, curve_act_noTimeout, n_act_timeout_matrix] = calc_pareto_matrix(integrated_prob_fake, events, classes, event_start);

    [acc_rest, acc_act, acc_act_noTimeout, n_act_timeout] = calc_pareto(integrated_prob_fake, events, classes, event_start, integratorCfg.feedbackThs);

    %% ----------------- plot prob integrated -----------------
    ic_index = find(cell2mat(gmmCfg.params.classes) == integratorCfg.ic_class_label);
    sampleRate_ros = bufferSize/chunkSize;
    do_plot = true;
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
        if ismember(cueTYP(idx_trial), [classes(1) classes(2)])
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
    m_sig = computeMetrics_signal(gmm_prob, qda_prob, qda_prob_trad, artifact, events, event_start, classes, cell2mat(gmmCfg.params.classes), integratorCfg);
    [m] = computeMetrics(integratorCfg, artifact, gmm_prob, qda_prob, integrated_prob, events, event_start, cell2mat(gmmCfg.params.classes), classes);
    current_method = 'gmm';
    entry = struct();
    entry.m_sig = m_sig;
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
    entry.Act_Time_Miss = m.act.time.miss / sampleRate_ros;
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
    entry.Pareto_Num_Timeouts = n_act_timeout_matrix;

    % SIMULATION
    entry.simulation.Rest_Acc = acc_rest * 100;
    entry.simulation.Act_Trial_Acc = acc_act * 100;
    entry.simulation.Act_Acc_NoTo = acc_act_noTimeout * 100;
    entry.simulation.Act_Num_Timeout = n_act_timeout;

    % CONTROL SNR
    entry.SNR_Active_Vel = m.snr.vel_active * 1000; % Scaling per leggibilità (es. 1e-3 -> 1)
    entry.SNR_Rest_Vel   = m.snr.vel_rest * 1000;
    entry.SNR_Ratio      = m.snr.ratio;

    % kappa
    entry.Act_Kappa_Sample = m.act.kappa.sample_qda;
    entry.Act_Kappa_Trial  = m.act.kappa.trial;

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
    [ringBufferCfg, artifactCfg, processingCfg, gmmCfg_trad, qdaCfg_trad, integratorCfg] = loadParameters(fullpath_file_parameters);

    disp(['   Loading QDA file: ', qdaCfg_trad.file_name])
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

    gmm_prob_trad = zeros(size(signals{1},1),2);
    gmm_prob_trad(:,1) = 1;

    gmmCfg_trad.params.classes = [1, 0];

    % my
    disp('   applying the GMM to all the signal') 
    mu_data = cellfun(@(x) x(1), gmmCfg.params.mu); 
    sigma_data = cellfun(@(x) x(1), gmmCfg.params.sigma);
    c_l = cell2mat(gmmCfg.params.central_left_idx);
    c_r = cell2mat(gmmCfg.params.central_right_idx);
    c_c = cell2mat(gmmCfg.params.central_idx);
    excl_chs = cell2mat(gmmCfg.params.excluded_idx);

    type = gmmCfg.params.type;
    nfeatures = gmmCfg.params.nfeatures;
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

    gmm_prob_my = posterior(gmmCfg.model, data_standardized);

    %% ----------------- QDA -----------------
    disp('   applying QDA to all the signal') 
    X = [];
    for idx_band=1:nbands
        for idx_band2=1:nbands
            if all(bands(idx_band,:) == qdaCfg_trad.model.bands(idx_band2,:))
                chs = qdaCfg_trad.model.idchans{idx_band2};
                tmp = signals{idx_band}(:,chs); 
            end
        end
        X = [X, tmp];
    end
    qda_prob_trad = apply_qda_matrix(qdaCfg_trad.model, log(X));

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
    qda_prob_my = apply_qda_matrix(qdaCfg.model, log(X)); 

    %% ----------------- integrated prob -----------------
    event_start = 781;
    rejection_look_gmm = false;
    [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, gmm_prob_trad, qda_prob_trad, events, event_start, gmmCfg_trad.params.classes, rejection_look_gmm);

    [integrated_prob_fake, ~] = applyIntegration(integratorCfg, artifact, gmm_prob_my, qda_prob_my, events, event_start, cell2mat(gmmCfg.params.classes), rejection_look_gmm);
    [curve_rest, curve_act, curve_act_noTimeout, n_act_timeout_matrix] = calc_pareto_matrix(integrated_prob_fake, events, classes, event_start);

    [acc_rest, acc_act, acc_act_noTimeout, n_act_timeout] = calc_pareto(integrated_prob_fake, events, classes, event_start, integratorCfg.feedbackThs);

    %% ----------------- plot prob integrated -----------------
    ic_index = find(gmmCfg_trad.params.classes == integratorCfg.ic_class_label);
    sampleRate_ros = bufferSize/chunkSize;
    do_plot = false;
    r_square_data_all_1 = []; r_square_data_all_2 = []; r_square_label_all = [];
    r_square_data_gmm_1 = []; r_square_data_gmm_2 = []; r_square_label_gmm = [];
    for idx_trial = 1:ntrial
        start_trial = fixPOS(idx_trial);
        end_trial = cfPOS(idx_trial) + cfDUR(idx_trial)-1;
        c_power_1 = log(signals{1}(start_trial:end_trial,:));
        c_power_2 = log(signals{2}(start_trial:end_trial,:));
        c_gmm_prob = gmm_prob_trad(start_trial:end_trial, :);
        c_artifact = artifact(start_trial:end_trial);
        c_qda_prob = qda_prob_trad(start_trial:end_trial,:);
        c_mask = mask(start_trial:end_trial);
        c_integrated = integrated_prob(start_trial:end_trial,:);
        trial_dur = size(c_power_1, 1);

        % for metrics r^2
        if ismember(cueTYP(idx_trial), [classes(1) classes(2)])
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
    m_sig = computeMetrics_signal(gmm_prob_my, qda_prob_my, qda_prob_trad, artifact, events, event_start, classes, cell2mat(gmmCfg.params.classes), integratorCfg);
    [m] = computeMetrics(integratorCfg, artifact, gmm_prob_trad, qda_prob_trad, integrated_prob, events, event_start, gmmCfg_trad.params.classes, classes);

    current_method = 'traditional';
    entry = struct();
    entry.m_sig = m_sig;
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
    entry.Act_Time_Miss = m.act.time.miss / sampleRate_ros;
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
    entry.Pareto_Num_Timeouts = n_act_timeout_matrix;

    % SIMULATION
    entry.simulation.Rest_Acc = acc_rest * 100;
    entry.simulation.Act_Trial_Acc = acc_act * 100;
    entry.simulation.Act_Acc_NoTo = acc_act_noTimeout * 100;
    entry.simulation.Act_Num_Timeout = n_act_timeout;

    % CONTROL SNR
    entry.SNR_Active_Vel = m.snr.vel_active * 1000; % Scaling per leggibilità (es. 1e-3 -> 1)
    entry.SNR_Rest_Vel   = m.snr.vel_rest * 1000;
    entry.SNR_Ratio      = m.snr.ratio;

    % Kappa
    entry.Act_Kappa_Sample = m.act.kappa.sample_qda;
    entry.Act_Kappa_Trial  = m.act.kappa.trial;

    if ~exist('Database','var'), Database = []; end
    Database = [Database; entry];

    [r2_values] = calc_r2_from_data(r_square_data_gmm_1, r_square_label_gmm, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | QDA data | ' num2str(size(r_square_data_gmm_1,1)) ' | band ' num2str(bands(1,1)) '-' num2str(bands(1,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all_1, r_square_label_all, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all_1,1)) ' | band ' num2str(bands(1,1)) '-' num2str(bands(1,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_gmm_2, r_square_label_gmm, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_gmm_2,1)) ' | band ' num2str(bands(2,1)) '-' num2str(bands(2,2))]);
    [r2_values] = calc_r2_from_data(r_square_data_all_2, r_square_label_all, 'Plot', false, 'ChanLabels', channels_label, 'title_data', ['traditional | all data | ' num2str(size(r_square_data_all_2,1)) ' | band ' num2str(bands(2,1)) '-' num2str(bands(2,2))]);

end

%% --- SAVE ---
subject_name = filenames_my{1}(1:2); % Prende i primi due caratteri, es. 'c7'
save_path = ['/home/paolo/cvsa/ic_cvsa_ws/record_mi/results_SMC/results_' subject_name '.mat'];

save(save_path, 'Database');
fprintf('   [INFO] Risultati salvati correttamente in: %s\n', save_path);



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


labels_database = unique({Database.Method});
for i = 1:length(labels_database)
    % FIGURE 6: VISUALIZZAZIONE "THRESHOLD TRENDS"
    db_idx = strcmpi({Database.Method}, labels_database{i});
    db_target = Database(db_idx);

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

    gmm_act_val  = mean([db_target.Act_Acc_NoTo], 'omitnan');
    gmm_rest_val = mean([db_target.Rest_Acc], 'omitnan');
    gmm_to_abs   = mean([db_target.Act_Num_Timeout], 'omitnan');
    gmm_to_perc  = (gmm_to_abs / TOTAL_TRIALS) * 100;

    figure('Name', 'Cross Validation', 'Color', 'w', 'Position', [50 50 1600 550]);
    x_axis = linspace(0.55, 1.0, size(mean_rest_mat, 2));
    thresholds = x_axis;

    subplot(1, 3, 1);
    imagesc(thresholds, thresholds, mean_act_noto_mat);
    axis xy; colormap(gca, 'jet');
    c = colorbar; c.Label.String = 'Accuracy [%]';
    hold on;
    plot([min(thresholds) max(thresholds)], [min(thresholds) max(thresholds)], 'k:', 'LineWidth', 1.5);
    [~, ~] = contour(thresholds, thresholds, mean_act_noto_mat, [gmm_act_val gmm_act_val], 'm-', 'LineWidth', 3);
    [r_good, c_good] = find(mean_act_noto_mat >= gmm_act_val);
    if ~isempty(r_good) && length(r_good) < 5
        plot(thresholds(c_good), thresholds(r_good), 'm.', 'MarkerSize', 15);
    end
    xlabel(['Threshold ' num2str(classes(1))], 'FontSize', 11, 'FontWeight','bold');
    ylabel(['Threshold ' num2str(classes(2))], 'FontSize', 11, 'FontWeight','bold');
    title({'\bf ACTIVE Accuracy no T', sprintf('%s Baseline: %.1f%%', labels_database{i}, gmm_act_val)}, 'FontSize', 14);

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
    title({'\bf REST Accuracy', sprintf('%s Baseline: %.1f%%', labels_database{i}, gmm_rest_val)}, 'FontSize', 14);

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
    title({'\bf TIMEOUT Rate', sprintf('%s Value: %.1f%%', labels_database{i}, gmm_to_perc)}, 'FontSize', 14);
    sgtitle(['data from evaluation where ' labels_database{i} ' methods was used'])


    % FIGURE 7
    all_acc_classic = [];
    all_acc_gmm_inf = [];
    for j = 1:length(db_target)
        all_acc_classic = [all_acc_classic; db_target(j).m_sig.acc_classic];
        all_acc_gmm_inf = [all_acc_gmm_inf; db_target(j).m_sig.acc_gmm_inf];
    end
    mean_classic = mean(all_acc_classic, 1, 'omitnan');
    mean_gmm_inf = mean(all_acc_gmm_inf, 1, 'omitnan');
    sem_classic = std(all_acc_classic, 0, 1, 'omitnan') ./ sqrt(size(all_acc_classic, 1));
    sem_gmm_inf = std(all_acc_gmm_inf, 0, 1, 'omitnan') ./ sqrt(size(all_acc_gmm_inf, 1));
    ths_labels = db_target(1).m_sig.ths_labels;
    figure('Color', 'w', 'Name', 'Informed QDA vs Classic QDA', 'Position', [100 100 900 600]);
    subplot(111)
    hold on;
    x = 1:length(ths_labels);
    width = 0.35;

    b1 = bar(x - width/2, mean_classic, width, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
    b2 = bar(x + width/2, mean_gmm_inf, width, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
    errorbar(x - width/2, mean_classic, sem_classic, 'k.', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    errorbar(x + width/2, mean_gmm_inf, sem_gmm_inf, 'k.', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    set(gca, 'XTick', x, 'XTickLabel', ths_labels, 'FontSize', 11);
    ylabel('Sample-wise Accuracy (%)', 'FontWeight', 'bold');
    xlabel('GMM Focalization Threshold (Online Gating)', 'FontWeight', 'bold');
    title({'\bf BCI Reliability Enhancement', ['data from evaluation where ' labels_database{i} ' methods was used']}, 'FontSize', 14);
    legend('Classic QDA (Trained on All Data)', 'Informed QDA (Trained on Focal Samples)', ...
        'Location', 'southoutside', 'Orientation', 'horizontal');
    grid on;
    ylim([45 100]);
    yline(50, 'r--', 'Chance Level', 'LineWidth', 1.5);
    box on;
    hold off;
end

%% FIGURA: KAPPA RELIABILITY (Informed vs Classic)
db_gmm = Database(strcmpi({Database.Method}, 'gmm'));

% Estrazione e media dei vettori Kappa tra tutti i file
all_k_cls = []; 
all_k_gmm = [];
for i = 1:length(db_gmm)
    all_k_cls = [all_k_cls; db_gmm(i).m_sig.kappa_classic];
    all_k_gmm = [all_k_gmm; db_gmm(i).m_sig.kappa_gmm_inf];
end

mean_k_cls = mean(all_k_cls, 1, 'omitnan');
mean_k_gmm = mean(all_k_gmm, 1, 'omitnan');
sem_k_cls  = std(all_k_cls, 0, 1, 'omitnan') ./ sqrt(size(all_k_cls,1));
sem_k_gmm  = std(all_k_gmm, 0, 1, 'omitnan') ./ sqrt(size(all_k_gmm,1));

figure('Name', 'Kappa Reliability Analysis', 'Color', 'w', 'Position', [100 100 800 500]);
x = 1:length(db_gmm(1).m_sig.ths_labels);
hold on;

% Plot Informed QDA (Blu)
errorbar(x, mean_k_gmm, sem_k_gmm, 's-', 'Color', [0.2 0.6 0.8], 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.6 0.8]);
% Plot Classic QDA (Grigio)
errorbar(x, mean_k_cls, sem_k_cls, 'o--', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [0.6 0.6 0.6]);

set(gca, 'XTick', x, 'XTickLabel', db_gmm(1).m_sig.ths_labels, 'FontSize', 11);
ylabel('Cohen’s Kappa (\kappa)', 'FontWeight', 'bold');
xlabel('GMM Focalization Threshold (Gate)', 'FontWeight', 'bold');
title({'\bf Performance Reliability Curve', '\rm Kappa as a function of Neural Focus'}, 'FontSize', 14);
legend('Informed QDA (Focus-Trained)', 'Classic QDA (All-Trained)', 'Location', 'best');
grid on; ylim([0 1]);
yline(0.6, 'k:', 'Substantial Control', 'LabelVerticalAlignment','bottom');
box on;

%% --- FUNZIONI DI SUPPORTO ---
[curve_rest, curve_act, curve_act_noTimeout, n_act_timeout_matrix] = calc_pareto_matrix(integrated_prob_fake, events, classes, event_start);

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

function [mat_acc_rest, mat_acc_act, mat_acc_act_noTimeout, mat_n_act_timeout] = calc_pareto_matrix(integrated_prob, events, task_classes, start_task)
    % INPUT:
    % integrated_prob: vettore (o matrice) delle probabilità integrate. 
    %                  Assumiamo Colonna 1: 0=DX, 1=SX.
    % events: struct con .TYP, .POS, .DUR
    % task_classes: [CODE_SX, CODE_DX, CODE_REST]
    % start_task: codice evento inizio trial (es. 781)
    
    % Def threshold to use
    thresholds = 0.55 : 0.05 : 1.0; 
    n_th = length(thresholds);
    
    % OUTPUT: Ora sono MATRICI (n_th x n_th)
    % Riga (i) -> Soglia SX
    % Colonna (j) -> Soglia DX
    mat_acc_rest            = zeros(n_th, n_th);
    mat_acc_act             = zeros(n_th, n_th);
    mat_acc_act_noTimeout   = zeros(n_th, n_th);
    mat_n_act_timeout       = zeros(n_th, n_th);
    
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
    
    % Loop 1: Varia la soglia per SX (Soglia Alta)
    for i_sx = 1:n_th
        th_val_sx = thresholds(i_sx); 
        th_high   = th_val_sx;       % Es. 0.70
        
        % Loop 2: Varia la soglia per DX (Soglia Bassa)
        for i_dx = 1:n_th
            th_val_dx = thresholds(i_dx);
            th_low    = 1.0 - th_val_dx;
            
            % Contatori per questa specifica combinazione (SX=i, DX=j)
            cnt_act_hit = 0;
            cnt_act_miss = 0;
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
                    if current_cue == CODE_SX
                        if strcmp(winner, 'high')
                            cnt_act_hit = cnt_act_hit + 1;
                        elseif strcmp(winner, 'low')
                            cnt_act_miss = cnt_act_miss + 1;
                        end
                    elseif current_cue == CODE_DX
                        if strcmp(winner, 'low')
                            cnt_act_hit = cnt_act_hit + 1;
                        elseif strcmp(winner, 'high')
                            cnt_act_miss = cnt_act_miss + 1;
                        end
                    end
                    
                    if strcmp(winner, 'none')
                        cnt_act_timeout = cnt_act_timeout + 1;
                    end
                    cnt_act_tot = cnt_act_tot + 1;
                    
                elseif current_cue == CODE_REST
                    cnt_rest_tot = cnt_rest_tot + 1;
                    
                    if strcmp(winner, 'none')
                        cnt_rest_ok = cnt_rest_ok + 1;
                    end
                end
            end
            
            % --- Assegnazione valori nella MATRICE (i_sx, i_dx) ---
            if cnt_act_tot > 0
                mat_acc_act(i_sx, i_dx) = cnt_act_hit / cnt_act_tot;
                
                den_noTime = cnt_act_tot + cnt_act_miss;
                if den_noTime > 0
                    mat_acc_act_noTimeout(i_sx, i_dx) = cnt_act_hit / den_noTime;
                else
                    mat_acc_act_noTimeout(i_sx, i_dx) = 0; % O NaN, a preferenza
                end
            end
            
            if cnt_rest_tot > 0
                mat_acc_rest(i_sx, i_dx) = cnt_rest_ok / cnt_rest_tot;
            end
            
            mat_n_act_timeout(i_sx, i_dx) = cnt_act_timeout;
        end
    end
end

function [acc_rest, acc_act, acc_act_noTimeout, n_act_timeout] = calc_pareto(integrated_prob, events, task_classes, start_task, ths)
    % INPUT:
    % integrated_prob: vettore (o matrice) delle probabilità integrate. 
    %                  Assumiamo Colonna 1: 0=DX, 1=SX.
    % events: struct con .TYP, .POS, .DUR
    % task_classes: [CODE_SX, CODE_DX, CODE_REST]
    % start_task: codice evento inizio trial (es. 781)
    % ths: threshold da dover usare

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
    
    th_high   = ths(1);
    th_low    = 1.0 - ths(2);
            
    % Contatori per questa specifica combinazione (SX=i, DX=j)
    cnt_act_hit = 0;
    cnt_act_miss = 0;
    cnt_act_timeout = 0;
    cnt_act_tot = 0;
    cnt_rest_ok = 0;
    cnt_rest_tot = 0;
            
    % Loop sui Trials
    for i_tr = 1:ntrials
        idx_s = cfPOS(i_tr);
        idx_e = idx_s + cfDUR(i_tr) - 1;
        if idx_e > length(integrated_prob)
            idx_e = length(integrated_prob);
        end

        sig = integrated_prob(idx_s:idx_e, 1);
        current_cue = cues(i_tr);

        % Cerca il primo attraversamento per SX (High)
        idx_cross_high = find(sig >= th_high, 1, 'first');
        idx_cross_low  = find(sig <= th_low, 1, 'first');

        res_high = ~isempty(idx_cross_high);
        res_low  = ~isempty(idx_cross_low);

        winner = 'none';
        if res_high && ~res_low
            winner = 'high';
        elseif ~res_high && res_low
            winner = 'low';
        elseif res_high && res_low
            if idx_cross_high < idx_cross_low
                winner = 'high';
            else
                winner = 'low';
            end
        end

        % --- Valutazione Hit/Miss/Timeout ---
        if current_cue == CODE_SX || current_cue == CODE_DX
            if current_cue == CODE_SX
                if strcmp(winner, 'high')
                    cnt_act_hit = cnt_act_hit + 1;
                elseif strcmp(winner, 'low')
                    cnt_act_miss = cnt_act_miss + 1;
                end
            elseif current_cue == CODE_DX
                if strcmp(winner, 'low')
                    cnt_act_hit = cnt_act_hit + 1;
                elseif strcmp(winner, 'high')
                    cnt_act_miss = cnt_act_miss + 1;
                end
            end

            if strcmp(winner, 'none')
                cnt_act_timeout = cnt_act_timeout + 1;
            end
            cnt_act_tot = cnt_act_tot + 1;

        elseif current_cue == CODE_REST
            cnt_rest_tot = cnt_rest_tot + 1;
            if strcmp(winner, 'none')
                cnt_rest_ok = cnt_rest_ok + 1;
            end
        end
    end

    acc_rest = cnt_rest_ok / cnt_rest_tot;
    acc_act = cnt_act_hit / cnt_act_tot;
    acc_act_noTimeout = cnt_act_hit / (cnt_act_hit + cnt_act_miss);
    n_act_timeout = cnt_act_timeout;
end

