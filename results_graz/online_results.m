% to save the metrics dataset online
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
subject = filenames_my{1}(1:2);

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
    CODE_SX = classes(1); CODE_DX = classes(2); CODE_REST = classes(3);
    CODE_HIT = 897; CODE_MISS = 898; CODE_TIMEOUT = 899;

    cfPOS = events.POS(events.TYP == event_start);
    cfDUR = events.DUR(events.TYP == event_start);
    ntrials = length(cfPOS);
    all_cues = events.TYP(ismember(events.TYP, classes));
    all_results = events.TYP(ismember(events.TYP, [CODE_HIT, CODE_MISS, CODE_TIMEOUT]));

    % counting variables
    cnt_act_hit = 0; cnt_act_miss = 0; cnt_act_timeout = 0; cnt_act_total = 0;
    cnt_rest_ok = 0; cnt_rest_err = 0; cnt_rest_total = 0;
    samp_corr_qda_my = 0; samp_tot_qda = 0;
    samp_corr_qda_trad = 0;

    % --- iterate over trials ---
    for i = 1:ntrials
        idx_start = cfPOS(i);
        idx_end = cfPOS(i) + cfDUR(i) - 1;
        dur = idx_end - idx_start + 1;

        cue = all_cues(i);
        res = all_results(i);

        s_art = artifact(idx_start:idx_end);
        s_qda_my = qda_prob(idx_start:idx_end, :);
        s_qda_trad = qda_prob_trad(idx_start:idx_end, :);

        % --- ACTIVE TASK ---
        if cue == CODE_SX || cue == CODE_DX
            cnt_act_total = cnt_act_total + 1;

            % Trial Result & Time
            if res == CODE_HIT
                cnt_act_hit = cnt_act_hit + 1;
            elseif res == CODE_MISS
                cnt_act_miss = cnt_act_miss + 1;
            elseif res == CODE_TIMEOUT
                cnt_act_timeout = cnt_act_timeout + 1;
            end

            % QDA Sample Accuracy
            if cue == CODE_SX, qda_col = 1; else, qda_col = 2; end
            if qda_col <= size(s_qda_my, 2)
                valid_mask = (s_art == 0);
                if any(valid_mask)
                    preds = s_qda_my(valid_mask, qda_col) >= 0.5;
                    samp_corr_qda_my = samp_corr_qda_my + sum(preds);
                    preds = s_qda_trad(valid_mask, qda_col) >= 0.5;
                    samp_corr_qda_trad = samp_corr_qda_trad + sum(preds);
                    samp_tot_qda = samp_tot_qda + length(preds);
                end
            end

            % --- REST TASK (783) ---
        elseif cue == CODE_REST
            cnt_rest_total = cnt_rest_total + 1;

            if res == CODE_HIT
                cnt_rest_ok = cnt_rest_ok + 1;
            elseif res == CODE_MISS
                cnt_rest_err = cnt_rest_err + 1;
            end
        end
    end

    % --- qda acc active ---
    entry = struct();
    entry.File = filenames_my{idx_file};
    entry.ActAccQDAMy = samp_corr_qda_my / samp_tot_qda * 100;
    entry.ActAccQDATrad = samp_corr_qda_trad / samp_tot_qda * 100;
    entry.ActAccTrial = cnt_act_hit / (cnt_act_total) *100; % Trial Accuracy
    entry.ActAccNo_timeout = cnt_act_hit / (cnt_act_miss + cnt_act_hit) *100;
    entry.RestAccTrial = cnt_rest_ok / (cnt_rest_total) * 100;
    entry.PercTimeout = cnt_act_timeout / cnt_act_total * 100;

    Database = [Database; entry];
end

%% save
res_dir = fullfile(pathname_my(1:38), 'results_graz/qda');
if ~exist(res_dir, 'dir'), mkdir(res_dir); end

save_name = fullfile(res_dir, ['online_' subject '.mat']);
save(save_name, 'Database', 'subject');

disp(['Dati salvati per soggetto ' subject ' in: ' save_name]);




