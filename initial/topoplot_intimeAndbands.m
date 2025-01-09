% clc; clear all; close all;
%% compute and show the log band of the file for 8-14 band in a topoplot
% works with the gdf

subject = 'h8';
cal_eval = input('calibration (1) or evaluation (2): ');
if cal_eval == 1
    a = 'calibration';
elseif cal_eval == 2
    a = 'evaluation';
else
    disp('Error on the input, only 1 or 2 are allowd');
    return;
end
day = '/20241015';
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
% bands = [[6 8]; [8 10]; [10 12]; [12 14]; [14 16]; [8 14]];
bands = [[10 12]; [8 14]; [16 18]];
avg = 1;
filtOrder = 4;
th_eog = 3.0e4; %2.5e4;

% variable for concatenation
total.ERD = cell(1, size(bands, 1));
total.minDurTrial = Inf(1, size(bands, 1));
total.cueTYP = [];


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

    % concatenate cue
    total.cueTYP = cat(1, total.cueTYP, cueTYP);

    % processing over different bands
    for idx_band=1:size(bands, 1)
        band = bands(idx_band,:);
        disp(['   [INFO] band used: [' num2str(band(1)) ',' num2str(band(2)) ']']);

        % Processing
        signal_processed = processing_offline(signal, nchannels, sampleRate, band, filtOrder, avg);

        % ERD/ERS. Dimension of ERD are nsample x nchannels x ntrial
        [ERD, minDurTrial, minDurFix] = compute_ERDERS(signal_processed, trialStart, trialStop, fixStart, fixStop, nchannels, n_trial, ck, normalization);

        % concatenate ERD signal
        if idx_f == 1
            total.ERD{idx_band} = ERD;
            total.minDurTrial(idx_band) = minDurTrial;
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
        end

        % Visualization
        if ~isempty(chanlog_path)
            divisionSampleRate = 2;
%             divisionSampleRate = 3;
            showTopoplot_ERDERS(chanlog_path, ERD, sampleRate, divisionSampleRate, channels_label, band, minDurTrial, cueTYP, classes, nchannels, ['matlab' files(idx_f).name])
        end
    end
end

%% Same analysis for all the files
% if input('Do you want to stop here? yes (1), no (): ') == 1
%     return;
% end
disp('All files concatenated' )

for idx_band=1:size(bands,1)
    band = bands(idx_band,:);
    disp(['   [INFO] band used: [' num2str(band(1)) ',' num2str(band(2)) ']']);

    % Visualization
    if ~isempty(chanlog_path)
        divisionSampleRate = 2;
        showTopoplot_ERDERS(chanlog_path, total.ERD{idx_band}, sampleRate, divisionSampleRate, channels_label, band, total.minDurTrial(idx_band), total.cueTYP, classes, nchannels, 'all files concatenated')
    end
end


