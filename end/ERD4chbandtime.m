clc; clear all; close all;
%% Show the ERD/ERS for a file and for the channels selected
% compute and show the ERD/ERS as a single image for a selected channel

subject = 'h8';
day = '/20240926';

cal_eval = input('calibration (1) or evaluation (2): ');
if cal_eval == 1
    a = 'calibration';
    using_data_seelcted = input('all (1) or channels selected for the evaluation (2): ');
    if using_data_seelcted == 1
        channelSelected = [13 14 15 16 17 18 29 30 31 32 33 34 35 36 37 38 39];
        bands = [[6 8]; [8 10]; [10 12]; [12 14]; [14 16]; [16 18]; [6 9]; [9 12]; [12 15]; [15 18]; [8 14]];
        plot_narow = 2;
        plot_ncol = length(channelSelected);
    elseif using_data_seelcted == 2
        yaml_QDA_path = ['/home/paolo/cvsa_ws/record/' subject day '/qda_' subject '.yaml'];
        qda = loadQDA(yaml_QDA_path);
        channelSelected = unique(qda.idchans);
        plot_narow = length(channelSelected);
        plot_ncol = 2;
        bands = unique(qda.bands, 'rows');
    else
        disp('Error on the input, only 1 or 2 are allowd');
        return;
    end
elseif cal_eval == 2
    a = 'evaluation';
    yaml_QDA_path = ['/home/paolo/cvsa_ws/record/' subject day '/qda_' subject '.yaml'];
    qda = loadQDA(yaml_QDA_path);
    channelSelected = unique(qda.idchans);
    plot_narow = length(channelSelected);
    plot_ncol = 2;
    bands = unique(qda.bands, 'rows');
else
    disp('Error on the input, only 1 or 2 are allowd');
    return;
end
path = ['/home/paolo/cvsa_ws/record/' subject day '/gdf/' a];
files = dir(fullfile(path, '*.gdf'));

chanlog_path = '/home/paolo/new_chanlocs64.mat';

channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
        '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

% channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

eog_channels =  {'FP1', 'FP2', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'EOG', ...
        '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''};

classes = [730 731];
sampleRate = 512;
normalization = false; % if false compute the log power, if true the ERD/ERS
avg = 1;
filtOrder = 4;
th_eog = 3.0e4; %2.5e4;

% variable for concatenation
total.ERD = cell(1, size(bands, 1));
total.baseline = cell(1, size(bands, 1));
total.minDurTrial = Inf(1, size(bands, 1));
total.minDurBaseline = Inf(1, size(bands, 1));
total.cueTYP = [];
total.n_trial = 0;


% check all files
for idx_f = 1:length(files)
    file = fullfile(path, files(idx_f).name);
    disp(['loading file: ' file])
    [signal, header] = sload(file);
    nchannels = length(channels_label);

    % Extract infoo data
    signal = signal(:, 1:nchannels);
    [fixPOS, fixDUR, cuePOS, cueDUR, cueTYP, cfPOS, cfDUR, n_trial] = extract_info_data(header, 786, [730 731], 781);

    % Extract trial infoo
    [trialStart, trialStop, fixStart, fixStop, ck, tk] = extract_trial_infoo(signal, fixPOS, fixDUR, cfPOS, cfDUR, cueTYP, n_trial);

    % EOG check
    prev_n_trial = n_trial;
    for i=1:n_trial
        c_start = cfPOS(i);
        c_end = cfPOS(i) + cfDUR(i) - 1;
        data = signal(c_start:c_end,:);
        if eye_movement_check(data, eog_channels, th_eog, sampleRate)
            trialStart(i) = -1;
            trialStop(i) = -1;
            fixStop(i) = -1;
            fixStart(i) = -1;
            cueTYP(i) = -1;
        end
    end
    trialStart(trialStart == - 1) = [];
    trialStop(trialStop == - 1) = [];
    fixStart(fixStart == - 1) = [];
    fixStop(fixStop == -1) = [];
    cueTYP(cueTYP == -1) = [];
    n_trial = size(trialStart, 1);
    disp(['   [INFO] from ' num2str(prev_n_trial) ' to ' num2str(n_trial)]);
    if n_trial == 0
        disp('ERROR no trials');
        return;
    end
    total.n_trial = total.n_trial + n_trial;

    % concatenate cue
    total.cueTYP = cat(1, total.cueTYP, cueTYP);

    % processing over different bands
    for idx_band=1:size(bands, 1)
        band = bands(idx_band,:);
        disp(['   [INFO] band used: [' num2str(band(1)) ',' num2str(band(2)) ']']);

        % Processing
        signal_processed = processing_offline(signal, nchannels, sampleRate, band, filtOrder, avg);

        % ERD/ERS. Dimension of ERD are nsample x nchannels x ntrial
        [ERD, minDurTrial, minDurBaseline, baseline] = compute_ERDERS(signal_processed, trialStart, trialStop, fixStart, fixStop, nchannels, n_trial, ck, normalization);

        % concatenate ERD signal
        if idx_f == 1
            total.ERD{idx_band} = ERD;
            total.minDurTrial(idx_band) = minDurTrial;
            total.baseline{idx_band} = baseline;
            total.minDurBaseline(idx_band) = minDurBaseline;
        else
            if minDurTrial < total.minDurTrial(idx_band)
                tmp = total.ERD{idx_band};
                tmp = tmp(1:minDurTrial,:,:);
                total.ERD{idx_band} = cat(3, tmp, ERD);
                total.minDurTrial(idx_band) = minDurTrial;
            else
                tmp = ERD(1:total.minDurTrial(idx_band),:,:);
                total.ERD{idx_band} = cat(3, total.ERD{idx_band}, tmp);
            end

            if minDurBaseline < total.minDurBaseline(idx_band)
                tmp = total.baseline{idx_band};
                tmp = tmp(1:minDurBaseline,:,:);
                total.baseline{idx_band} = cat(3, tmp, baseline);
                total.minDurBaseline(idx_band) = minDurBaseline;
            else
                tmp = baseline(1:total.minDurBaseline(idx_band),:,:);
                total.baseline{idx_band} = cat(3, total.baseline{idx_band}, tmp);
            end
        end

        % Visualization
        if ~isempty(chanlog_path)
%             divisionSampleRate = 2;
            divisionSampleRate = 3;
%             showTopoplot_ERDERS(chanlog_path, ERD, sampleRate, divisionSampleRate, channels_label, band, minDurTrial, cueTYP, classes, nchannels, ['matlab' files(idx_f).name])
        end
    end
end

%% Calculate ERD/ERS 
trial_data = nan(total.minDurTrial(1), size(bands, 1), nchannels, size(total.ERD{1}, 3)); % samples x freqs x channels x trials
basline_data = nan(total.minDurBaseline(1), size(bands, 1), nchannels, size(total.baseline{1}, 3));
for idx_band=1:size(bands, 1)
    trial_data(:, idx_band, :, :) = total.ERD{idx_band}; 
    basline_data(:, idx_band, : ,:) = total.baseline{idx_band};
end
ERD = nan(size(trial_data));
for trial=1:total.n_trial 
    for ch=1:nchannels
        for freq=1:size(bands, 1)
            baseline = mean(basline_data(:,freq,ch,trial));
            ERD(:,freq,ch,trial) = trial_data(:,freq,ch,trial)/baseline;
        end
    end
end

%% Visualization divided for classes
band_strings = strings(size(bands, 1), 1);
for i = 1:size(bands, 1)
    band_strings(i) = sprintf('%d-%d', bands(i, 1), bands(i, 2));
end
figure;
t = linspace(0, minDurTrial/sampleRate, minDurTrial);
classes = [730 731];
selclassLb = {'Bottom left', 'Bottom right'};
chandles = [];
for cId = 1:length(classes)
    climits = nan(2, length(channelSelected));
    for chId = 1:length(channelSelected)
%         subplot(plot_narow, ceil(length(channelSelected)/3)*2, cId + (chId - 1)*2);
        subplot(plot_narow, plot_ncol, cId + (chId - 1)*2);
        cdata = mean(ERD(:, :, channelSelected(chId), total.cueTYP == classes(cId)), 4);
        imagesc(t, 1:size(bands, 1), cdata');
        set(gca,'YDir','normal');
        climits(:, chId) = get(gca, 'CLim');
        chandles = cat(1, chandles, gca);
        colormap(hot);
        colorbar;
        title(['Channel ' channels_label{channelSelected(chId)} ' | ' selclassLb{cId}]);
        xlabel('Time [s]');
        ylabel('Frequency [Hz]');
        yticks(1:size(bands, 1));
        yticklabels(band_strings)
        line([2 2],get(gca,'YLim'),'Color',[0 0 0])
        line([3 3],get(gca,'YLim'),'Color',[0 0 0])
    end
end

%% Visualization with the difference between classes
band_strings = strings(size(bands, 1), 1);
for i = 1:size(bands, 1)
    band_strings(i) = sprintf('%d-%d', bands(i, 1), bands(i, 2));
end
figure;
t = linspace(0, minDurTrial/sampleRate, minDurTrial);
title_diff = [num2str(classes(2)) ' - ' num2str(classes(1))];
nselclasses = length(classes);
chandles = [];
for cId = 1:nselclasses
    climits = nan(2, length(channelSelected));
    for chId = 1:length(channelSelected)
%         subplot(plot_narow, ceil(length(channelSelected)/3)*2, cId + (chId - 1)*2);
        subplot(plot_narow, plot_ncol, cId + (chId - 1)*2);
        cdata = mean(ERD(:, :, channelSelected(chId), total.cueTYP == classes(cId)), 4);
        imagesc(t, 1:size(bands, 1), cdata');
        set(gca,'YDir','normal');
        climits(:, chId) = get(gca, 'CLim');
        chandles = cat(1, chandles, gca);
        colormap(hot);
        colorbar;
        title(['Channel ' channels_label{channelSelected(chId)} ' | ' title_diff]);
        xlabel('Time [s]');
        ylabel('Frequency [Hz]');
        yticks(1:size(bands, 1));
        yticklabels(band_strings)
        line([2 2],get(gca,'YLim'),'Color',[0 0 0])
        line([3 3],get(gca,'YLim'),'Color',[0 0 0])
    end
end

