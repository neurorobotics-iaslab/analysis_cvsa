clc; clear all; close all;
% compute and show the log band of the file for 8-14 band in a topoplot
% works with the gdf

subject = 'c7';
cal_eval = input('calibration (1) or evaluation (2): ');
if cal_eval == 1
    a = 'calibration';
elseif cal_eval == 2
    a = 'evaluation';
else
    disp('Error on the input, only 1 or 2 are allowd');
    return;
end
path = ['/home/paolo/cvsa_ws/record/' subject '/gdf/' a];
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
band1 = (8:2:16)';
band2 = (10:2:18)';
bands = [band1, band2];
bands = cat(1, bands, [8,14]);
avg = 1;
filtOrder = 4;
th_eog = 2.5e4;

%% with all channel for the laplacian
total_signal = cell(1, size(bands, 1));
total_header = cell(1, size(bands, 1));
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
    for i=1:n_trial
        c_start = cfPOS(i);
        c_end = cfPOS(i) + cfDUR(i) - 1;
        data = signal(c_start:c_end,:);
        if eye_movement_check(data, eog_channels, th_eog, sampleRate)
            trialStart(i) = -1;
            trialStop(i) = -1;
            fixStop(i) = -1;
            fixStart(i) = -1;
        end
    end

    % processing over different bands
    for idx_band=1:size(bands, 1)
        band = bands(idx_band,:);
        disp(['   [INFO] band used: [' num2str(band(1)) ',' num2str(band(2)) ']']);

        % Processing
        signal_processed = processing_offline(signal, nchannels, sampleRate, band, filtOrder, avg);

        % ERD/ERS
        [ERD, minDurTrial, minDurFix] = compute_ERDERS(signal_processed, trialStart, trialStop, fixStart, fixStop, nchannels, n_trial, ck, normalization);

        % Visualization
        if ~isempty(chanlog_path)
            showTopoplot_ERDERS(chanlog_path, ERD, sampleRate, 2, channels_label, band, minDurTrial, cueTYP, classes, nchannels, files(idx_f).name)
        end

        if idx_f == 1
            tmp.EVENT.POS = header.EVENT.POS;
            tmp.EVENT.DUR = header.EVENT.DUR;
            tmp.EVENT.TYP = header.EVENT.TYP;
            total_header{idx_band} = tmp;
            total_signal{idx_band} = signal_processed;
        else
            tmp.EVENT.POS = cat(1, tmp.EVENT.POS, header.EVENT.POS + size(total_signal{idx_band}, 1));
            tmp.EVENT.DUR = cat(1, tmp.EVENT.DUR, header.EVENT.DUR);
            tmp.EVENT.TYP = cat(1, tmp.EVENT.TYP, header.EVENT.TYP);
            total_header{idx_band} = tmp;
            total_signal{idx_band} = cat(1, total_signal{idx_band}, signal_processed);
        end
    end
end

%% Same analysis for all the files
if input('Do you want to stop here? yes (1), no (): ') == 1
    return;
end
disp('All files concatenated' )

for idx_band=1:size(bands,1)
    band = bands(idx_band,:);
    disp(['   [INFO] band used: [' num2str(band(1)) ',' num2str(band(2)) ']']);

    header = total_header{idx_band};
    signal = total_signal{idx_band};

    % Extract infoo
    [fixPOS, fixDUR, cuePOS, cueDUR, cueTYP, cfPOS, cfDUR, n_trial] = extract_info_data(header, 786, [730 731], 781);
    [trialStart, trialStop, fixStart, fixStop, ck, tk] = extract_trial_infoo(signal, fixPOS, fixDUR, cfPOS, cfDUR, cueTYP, n_trial);

    % ERD/ERS
    [ERD, minDurTrial, minDurFix] = compute_ERDERS(signal, trialStart, trialStop, fixStart, fixStop, nchannels, n_trial, ck, normalization);

    % Visualization
    if ~isempty(chanlog_path)
        showTopoplot_ERDERS(chanlog_path, ERD, sampleRate, 2, channels_label, band, minDurTrial, cueTYP, classes, nchannels, 'all files concatenated')
    end
end


