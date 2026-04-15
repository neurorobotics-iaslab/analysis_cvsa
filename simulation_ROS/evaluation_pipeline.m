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
paradigm = filenames{1}(31:33);
if strcmp(paradigm, 'mi')
    classes = [769, 770];
    datapath = './src/qda_bci/';
    yaml_QDA_path_mi = [datapath 'cfg/mi/qda_test_mi.yaml'];
    disp('Loading QDA MI...');
    qda_mi = loadQDA(yaml_QDA_path_mi);
elseif strcmp(paradigm, 'cvsa')
    classes = [730, 731];
    datapath = './src/qda_bci/';
    yaml_QDA_path_cvsa = [datapath 'cfg/cvsa/qda_test_cvsa.yaml'];
    disp('Loading QDA CVSA...');
    qda_cvsa = loadQDA(yaml_QDA_path_cvsa);
elseif strcmp(paradigm, 'hybrid')
    classes = [750 751];
    datapath = './src/qda_bci/';
    yaml_QDA_path_mi = [datapath 'cfg/mi/qda_test_mi.yaml'];
    yaml_QDA_path_cvsa = [datapath 'cfg/cvsa/qda_test_cvsa.yaml'];
    disp('Loading QDA MI...');
    qda_mi = loadQDA(yaml_QDA_path_mi);
    disp('Loading QDA CVSA...');
    qda_cvsa = loadQDA(yaml_QDA_path_cvsa);
end
nclasses = length(classes);


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
    [ringBufferCfg, artifactCfg, processingCfg, qda_mi, qda_cvsa, integratorCfg] = loadParameters(fullpath_file_parameters);

    disp(['   Loading QDA mi file: ', qda_mi.file_name])
    gmm_path = [pathname(1:end-4) qda_mi.file_name];
    [qda_mi.model, qda_mi.params] = loadGMM(gmm_path);

    disp(['   Loading QDA cvsa file: ', qda_cvsa.file_name])
    qda_path = [pathname(1:end-4) qda_cvsa.file_name];
    qda_cvsa.model = loadQDA(qda_path);

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
    do_hann = processingCfg.do_hann;
    [signal_processed, header_processed] = processing_onlineROS_CSD_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize, do_hann);

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
    if strcmp(paradigm, 'hybrid') || strcmp(paradigm, 'cvsa')
        disp('   applying QDA cvsa to all the signal') 
        qda_prob_cvsa = apply_qda_matrix(qda_cvsa.model, log(signal_processed(:,qda_cvsa.model.idchans{1})));
    end

    %% ----------------- QDA MI -----------------
    if strcmp(paradigm, 'hybrid') || strcmp(paradigm, 'mi')
        disp('   applying QDA mi to all the signal')
        X = [];
        bands = qda_mi.bands;
        nbands = qda_mi.nbands;
        for idx_band=1:nbands
            for idx_band2=1:nbands
                if all(bands(idx_band,:) == qdaCfg.model.bands(idx_band2,:))
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
    if strcmp(paradigm, 'hybrid') 
        [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, qda_prob_cvsa, qda_prob_mi, event, event_start, rejection, paradigm, sampleRate);
    elseif strcmp(paradigm, 'mi')
        [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, [], qda_prob_mi, event, event_start, rejection, paradigm, sampleRate);
    elseif strcmp(paradigm, 'cvsa')
        [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, qda_prob_cvsa, [], event, event_start, rejection, paradigm, sampleRate);
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
        c_qda_prob_mi = qda_prob_mi(start_trial:end_trial, :);
        c_artifact = artifact(start_trial:end_trial);
        c_qda_prob_cvsa = qda_prob_cvsa(start_trial:end_trial,:);
        c_mask = mask(start_trial:end_trial);
        c_integrated = integrated_prob(start_trial:end_trial,:);

        % for metrics r^2
        r_square_data = [r_square_data; c_power];
        r_square_label = [r_square_label; repmat(cueTYP(idx_trial), trial_dur, 1)];

        if boom(idx_trial) == 897
            time_hit = [time_hit; size(c_power, 1) / sampleRate];
        elseif boom(idx_trial) == 898
            time_miss = [time_miss; size(c_power, 1) / sampleRate];
        elseif boom(idx_trial) == 899
            time_tout = [time_tout; size(c_power, 1) / sampleRate];
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
            tmp_prob_mi = c_qda_prob_mi;
            tmp_prob_mi(c_mask == 0,1) = nan;
            tmp_prob_cvsa = c_qda_prob_cvsa;
            tmp_prob_cvsa(c_mask == 0,1) = nan;
            plot(c_qda_prob_mi(:,1))
            hold on
            plot(c_qda_prob_cvsa(:,1))
            plot(c_artifact);
            scatter(1:size(c_qda_prob_mi, 1), c_qda_prob_mi(:,1), 15, 'black', 'filled')
            scatter(1:size(c_qda_prob_mi, 1), tmp_prob_mi(:,1), 15, 'green', 'filled')
            scatter(1:size(c_qda_prob_cvsa, 1), c_qda_prob_cvsa(:,1), 15, 'black', 'filled')
            scatter(1:size(c_qda_prob_cvsa, 1), tmp_prob_cvsa(:,1), 15, 'green', 'filled')
            yline(0.5, 'LineStyle','--');
            xline(cueDUR(idx_trial)+fixDUR(idx_trial), 'LineStyle','-');
            xline(fixDUR(idx_trial), 'LineStyle','-');
            hold off
            ylim([0 1])
            xlim([1 trial_dur])
            legend('qda prob mi', 'qda prob cvsa','artifact','qda mi not used', 'qda mi prob used', 'qda cvsa not used', 'qda cvsa prob used')
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
