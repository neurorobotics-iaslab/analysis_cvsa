% OFFLINE SIMULATION per il caso Ibrido MI + CVSA
% script che emula il nodo ROS offline basandosi sul file GDF salvato.
clear all; close all;

addpath(genpath('/home/paolo/bci_vr_ws/src/analysis_bci/equal_ros'));
addpath(genpath('/home/paolo/bci_vr_ws/src/analysis_bci/utils'));

%% Initialization
nchannels = 32;

[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if isequal(filenames, 0)
    disp('No files selected.');
    return;
end
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
    channels_label = header.Label(1:nchannels);

    %% ----------------- load the file -----------------
    disp(['   Loading parameters file: ', filenames{idx_file}(1:end-3) 'yaml'])
    fullpath_file_parameters = [pathname(1:end-4) 'parameters/' filenames{idx_file}(1:end-3) 'yaml'];
    [ringBufferCfg, artifactCfg, processingCfg, qda_mi, qda_cvsa, integratorCfg, paradigm] = loadParameters(fullpath_file_parameters);

    if strcmp(paradigm, 'mi_lhrh') || strcmp(paradigm, 'hybrid')
        disp(['   Loading QDA mi file: ', qda_mi.file_name])
        qda_path = qda_mi.path_to_model;
        qda_mi.model = loadQDA(qda_path);
        classes = qda_mi.model.classes;
    end

    if strcmp(paradigm, 'cvsa_blbr') || strcmp(paradigm, 'hybrid')
        disp(['   Loading QDA cvsa file: ', qda_cvsa.file_name])
        qda_path = qda_cvsa.path_to_model;
        qda_cvsa.model = loadQDA(qda_path);
        classes = qda_cvsa.model.classes;
    end

    if strcmp(paradigm, 'hybrid')
        classes = [750 751];
    end

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
    [~, indici] = ismember(eog.label, channels_label);
    eog_channels = indici(indici > 0);
    filterOrder = processingCfg.filterOrder;
    bands = processingCfg.bands; % i knwo we are using one band
    nbands = size(bands, 1);
    do_hann = processingCfg.do_hann;
    signals = cell(nbands, 1);
    for idx_bands = 1:nbands
        c_band = bands(idx_bands,:);
        [signal_processed, header_processed] = processing_onlineROS_CAR_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, c_band, chunkSize, eog_channels, do_hann);
        signals{idx_bands} = signal_processed;
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

    %% ----------------- QDA CVSA -----------------
    if strcmp(paradigm, 'hybrid') || strcmp(paradigm, 'cvsa_blbr')
        disp('   applying QDA cvsa to all the signal') 
        X = [];
        bands = qda_cvsa.model.bands;
        nbands = qda_cvsa.model.nbands;
        for idx_band=1:nbands
            for idx_band2=1:nbands
                if all(bands(idx_band,:) == qda_cvsa.model.bands(idx_band2,:))
                    chs = qda_cvsa.model.idchans{idx_band2};
                    tmp = signals{idx_band}(:,chs);
                end
            end
            X = [X, tmp];
        end
        qda_prob_cvsa = apply_qda_matrix(qda_cvsa.model, log(X));
    end

    %% ----------------- QDA MI -----------------
    if strcmp(paradigm, 'hybrid') || strcmp(paradigm, 'mi_lhrh')
        disp('   applying QDA mi to all the signal')
        X = [];
        bands = qda_mi.model.bands;
        nbands = qda_mi.model.nbands;
        for idx_band=1:nbands
            for idx_band2=1:nbands
                if all(bands(idx_band,:) == qda_mi.model.bands(idx_band2,:))
                    chs = qda_mi.model.idchans{idx_band2};
                    tmp = signals{idx_band}(:,chs);
                end
            end
            X = [X, tmp];
        end
        qda_prob_mi = apply_qda_matrix(qda_mi.model, log(X));
    end

    %% ----------------- integrated prob -----------------
    event_start = 781;
    rejection = 0.5;
    fs_processed = sampleRate / chunkSize;
    if strcmp(paradigm, 'hybrid') 
        [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, qda_prob_cvsa, qda_prob_mi, header_processed.EVENT, event_start, rejection, paradigm, fs_processed);
    elseif strcmp(paradigm, 'mi_lhrh')
        [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, [], qda_prob_mi, header_processed.EVENT, event_start, rejection, paradigm, fs_processed);
    elseif strcmp(paradigm, 'cvsa_blbr')
        [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, qda_prob_cvsa, [], header_processed.EVENT, event_start, rejection, paradigm, fs_processed);
    end
    
    %% ----------------- plot prob integrated -----------------
    do_plot = true;
    r_square_data = []; r_square_label = [];
    cnt_hit = sum(boom == 897);
    cnt_miss = sum(boom == 898);
    cnt_tout = sum(boom == 899);
    time_hit = []; time_miss = []; time_tout = [];
    
    for idx_trial = 1:ntrial
        start_trial = fixPOS(idx_trial);
        end_trial = cfPOS(idx_trial) + cfDUR(idx_trial);
        trial_dur = end_trial - start_trial;
        c_power = log(signal_processed(start_trial:end_trial,:));
        c_artifact = artifact(start_trial:end_trial);
        c_mask = mask(start_trial:end_trial);
        c_integrated = integrated_prob(start_trial:end_trial,:);

        if strcmp(paradigm, 'mi_lhrh') || strcmp(paradigm, 'hybrid')
            c_qda_prob_mi = qda_prob_mi(start_trial:end_trial, :);
        end
        if strcmp(paradigm, 'cvsa_blbr') || strcmp(paradigm, 'hybrid')
            c_qda_prob_cvsa = qda_prob_cvsa(start_trial:end_trial,:);
        end

        % for metrics r^2
        r_square_data = [r_square_data; c_power];
        r_square_label = [r_square_label; repmat(cueTYP(idx_trial), trial_dur, 1)];

        if boom(idx_trial) == 897
            time_hit = [time_hit; size(c_power, 1) / fs_processed];
        elseif boom(idx_trial) == 898
            time_miss = [time_miss; size(c_power, 1) / fs_processed];
        elseif boom(idx_trial) == 899
            time_tout = [time_tout; size(c_power, 1) / fs_processed];
        end

        if do_plot
            figure();
            subplot(311)
            imagesc(c_power')
            hold on;
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off;
            yticks(1:nchannels); yticklabels(channels_label)
            title('Log band power')


            subplot(312)
            if strcmp(paradigm, 'cvsa_blbr')
                tmp_prob_cvsa = c_qda_prob_cvsa;
                tmp_prob_cvsa(c_mask == 0,1) = nan;
                plot(c_qda_prob_cvsa(:,1))
                hold on
            elseif strcmp(paradigm, 'mi_lhrh')
                tmp_prob_mi = c_qda_prob_mi;
                tmp_prob_mi(c_mask == 0,1) = nan;
                plot(c_qda_prob_mi(:,1))
                hold on
            elseif strcmp(paradigm, 'hybrid')
                tmp_prob_cvsa = c_qda_prob_cvsa;
                tmp_prob_cvsa(c_mask == 0,1) = nan;
                plot(c_qda_prob_cvsa(:,1))
                hold on
                tmp_prob_mi = c_qda_prob_mi;
                tmp_prob_mi(c_mask == 0,1) = nan;
                plot(c_qda_prob_mi(:,1))
            end
            plot(c_artifact);
            if strcmp(paradigm, 'mi_lhrh') || strcmp(paradigm, 'hybrid')
                scatter(1:size(c_qda_prob_mi, 1), c_qda_prob_mi(:,1), 15, 'black', 'filled')
                scatter(1:size(c_qda_prob_mi, 1), tmp_prob_mi(:,1), 15, 'green', 'filled')
            end
            if strcmp(paradigm, 'cvsa_blbr') || strcmp(paradigm, 'hybrid')
                scatter(1:size(c_qda_prob_cvsa, 1), c_qda_prob_cvsa(:,1), 15, 'black', 'filled')
                scatter(1:size(c_qda_prob_cvsa, 1), tmp_prob_cvsa(:,1), 15, 'green', 'filled')
            end
            yline(0.5, 'LineStyle','--');
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off
            ylim([0 1])
            xlim([1 trial_dur])
            if strcmp(paradigm, 'hybrid')
                legend('qda prob cvsa', 'qda prob mi','artifact','qda mi not used', 'qda mi prob used', 'qda cvsa not used', 'qda cvsa prob used')
            elseif strcmp(paradigm, 'cvsa_blbr')
                legend('qda prob cvsa','artifact', 'qda cvsa not used', 'qda cvsa prob used')
            elseif strcmp(paradigm, 'mi_lhrh')
                legend('qda prob mi', 'artifact','qda mi not used', 'qda mi prob used')
            end
            title('classifier probability')

            subplot(313)
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

    %% print accuracy
    disp('   metrics')
    disp(['      accuracy trial hit: ' num2str(cnt_hit/ntrial *100) '%'])
    disp(['      time mean hit: ' num2str(mean(time_hit)) 's'])
    disp(['      time mean miss: ' num2str(mean(time_miss)) 's'])
    disp(['      time mean tout: ' num2str(mean(time_tout)) 's'])
    [r2_values] = calc_r2_from_data(r_square_data, r_square_label, 'Plot', true, 'ChanLabels', channels_label, 'title_data', 'all data');

end
