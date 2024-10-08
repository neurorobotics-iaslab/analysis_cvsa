%% for each trial look the mean of emisphere r and l
clc; clearvars; close all;
% faccio right - left, se ho 730 vuol dire che ho chiesto di fare bottom
% left, quindi mi aspetto più segnale nell'emisfero destro, quindi nel plot
% the rimanga positivo

%% General infromation
subject = 'c7';
path = ['/home/paolo/cvsa_ws/record/' subject '/gdf/calibration'];
% path = ['/home/paolo/cvsa_ws/record/' subject '/gdf'];
files = dir(fullfile(path, '*.gdf'));

lap_path39 = '/home/paolo/laplacians/lap_39ch_CVSA.mat';
load(lap_path39);


% chs_l = {'P3', 'PZ', 'POZ', 'O1', 'P5', 'P1', 'PO5', 'PO3', 'PO7', 'OZ'};
% chs_r = {'PZ', 'P4', 'POZ', 'O2', 'P2', 'P6', 'PO4', 'PO6', 'PO8', 'OZ'};
chs_l = {'P3', 'O1', 'P5', 'P1', 'PO5', 'PO3', 'PO7'};
chs_r = {'P4', 'O2', 'P2', 'P6', 'PO4', 'PO6', 'PO8'};
channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
        'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
nchannels = size(channels_label, 2);
idx_ch_l = find(ismember(channels_label, chs_l));
idx_ch_r = find(ismember(channels_label, chs_r));

band = [8 14];
filterOrder = 4;

% Variables for data
TYP = [];
POS = [];
DUR = [];
TARGET = [];
classes = [730, 731];
sampleRate = 512;

for idx_f=1:length(files)
    %% Take data
    disp(['[INFO] Load file '  num2str(idx_f) '/' num2str(length(files))]);
    file = fullfile(path, files(idx_f).name);
    [signal, header] = sload(file);
%     load(file);
    signal = signal(:,1:nchannels);
    c_events = header.EVENT;

    %% Data processing
    sampleRate = 512;
    disp('   [PROC] Using log band power');
    % Apply car filter
%     disp('      [PROC] Apply car filter');
%     s_car = signal - mean(signal, 2);

    % Apply lap filter
%     disp('      [PROC] Apply lap filter');
%     s_lap = s_car * lap;

    % Apply filtering in band
    disp(['      [PROC] Apply filter ' num2str(band(1)) '-' num2str(band(2)) 'hz']);
    [b, a] = butter(filterOrder, band*2/sampleRate);
    s_band = filtfilt(b, a, signal);
%     s_band = filtfilt(b, a, s_lap);

    % Rect the signal
    disp('      [PROC] Rectifing signal');
    s_power = power(s_band, 2);

    % Apply average windows
    disp('      [PROC] Apply average windows');
    avg = 1;
    windowSize = avg * sampleRate;
    s_avg = zeros(size(s_power));
    for ch=1:nchannels
        s_avg(:,ch) = filter(ones(1, windowSize)/(windowSize), 1, s_power(:,ch));
    end

    % Apply log
    disp('      [PROC] Apply log to have log band power');
    s_log = log(s_avg);

    %% Extract informations
    nTrials = sum(c_events.TYP == 1);
    cfPOS = c_events.POS(c_events.TYP == 781);
    cfDUR = c_events.DUR(c_events.TYP == 781);

    cuePOS = c_events.POS(ismember(c_events.TYP, classes));
    cuePOS = cuePOS(length(cuePOS)-nTrials+1:end); % remove trial calibration eyes
    cueDUR = c_events.DUR(ismember(c_events.TYP, classes));
    cueDUR = cueDUR(length(cueDUR)-nTrials+1:end);
    cueTYP = c_events.TYP(ismember(c_events.TYP, classes));
    cueTYP = cueTYP(length(cueTYP)-nTrials+1:end);

    fixPOS = c_events.POS(c_events.TYP == 786);

    % extract trials
    trial_start = nan(nTrials, 1);
    trial_end   = nan(nTrials, 1);
    trial_cue   = nan(nTrials, 1);
    for idx_t=1:nTrials
        trial_start(idx_t) = fixPOS(idx_t);
        trial_end(idx_t)   = cfPOS(idx_t) + cfDUR(idx_t) - 1;
        trial_cue(idx_t)   = cueTYP(idx_t);
    end

    %% Compute the emisphere differences for each trial and show it
    figure();
    for idx_t=1:nTrials
        c_signal = s_log(trial_start(idx_t):trial_end(idx_t),:);

        c_signal_l = c_signal(:, idx_ch_l);
        c_signal_r = c_signal(:, idx_ch_r);

        dif = mean(c_signal_r, 2) - mean(c_signal_l,2);

        % show
        subplot(3, ceil(nTrials/3), idx_t);
        hold on;
        grid on;
        t_signal = linspace(0, size(dif, 1)/sampleRate, size(dif, 1));
        scatter(t_signal, dif, 2, "black", "filled");
        line([(cuePOS(idx_t) - fixPOS(idx_t)), (cuePOS(idx_t) - fixPOS(idx_t))]/sampleRate, [min(dif), max(dif)], 'Color', 'blue', 'LineWidth', 1);
        line([(cfPOS(idx_t) - fixPOS(idx_t)), (cfPOS(idx_t) - fixPOS(idx_t))]/sampleRate, [min(dif), max(dif)], 'Color', 'blue', 'LineWidth', 1);
        line([0, t_signal(end)], [0 0], 'Color', 'blue', 'LineWidth', 1);
        hold off;
        title([ 'Cue: ' num2str(trial_cue(idx_t))]);
    end

    sgtitle(['file: ' files(idx_f).name ' | log band power']);
    
    %% same as before but with psd
    disp('   [PROC] Using PSD');
    % data processing
    disp('      [PROC] Apply lap filter');
%     s_lap = s_car * lap;
    disp('      [PROC] Apply PSD');
    psd_wlength = 0.5;
    psd_wshift = 0.0625;
    psd_pshift = 0.25;
    psd_mlength =  1;
    [features, f] = proc_spectrogram(signal, psd_wlength, psd_wshift, psd_pshift, sampleRate, psd_mlength);
    disp('      [PROC] Apply log');
    features = log(features);
    c_events.POS = proc_pos2win(c_events.POS, psd_wshift*sampleRate, 'backward', psd_mlength*sampleRate);
    c_events.DUR = floor(c_events.DUR/(psd_wshift*sampleRate)) + 1;
    c_events.TYP = c_events.TYP;

    % Extract informations
    sampleRate = 16;
    nTrials = sum(c_events.TYP == 1);
    cfPOS = c_events.POS(c_events.TYP == 781);
    cfDUR = c_events.DUR(c_events.TYP == 781);

    cuePOS = c_events.POS(ismember(c_events.TYP, classes));
    cuePOS = cuePOS(length(cuePOS)-nTrials+1:end); % remove trial calibration eyes
    cueDUR = c_events.DUR(ismember(c_events.TYP, classes));
    cueDUR = cueDUR(length(cueDUR)-nTrials+1:end);
    cueTYP = c_events.TYP(ismember(c_events.TYP, classes));
    cueTYP = cueTYP(length(cueTYP)-nTrials+1:end);

    fixPOS = c_events.POS(c_events.TYP == 786);

    % extract trials
    trial_start = nan(nTrials, 1);
    trial_end   = nan(nTrials, 1);
    trial_cue   = nan(nTrials, 1);
    for idx_t=1:nTrials
        trial_start(idx_t) = fixPOS(idx_t);
        trial_end(idx_t)   = cfPOS(idx_t) + cfDUR(idx_t) - 1;
        trial_cue(idx_t)   = cueTYP(idx_t);
    end

    % Compute the emisphere differences for each trial and show it
    figure();
    for idx_t=1:nTrials
        c_signal = s_log(trial_start(idx_t):trial_end(idx_t),:);

        c_signal_l = c_signal(:, idx_ch_l);
        c_signal_r = c_signal(:, idx_ch_r);

        dif = mean(c_signal_r, 2) - mean(c_signal_l,2);

        % show
        subplot(3, ceil(nTrials/3), idx_t);
        hold on;
        grid on;
        t_signal = linspace(0, size(dif, 1)/sampleRate, size(dif, 1));
        scatter(t_signal, dif, 2, "black", "filled");
        line([(cuePOS(idx_t) - fixPOS(idx_t)), (cuePOS(idx_t) - fixPOS(idx_t))]/sampleRate, [min(dif), max(dif)], 'Color', 'blue', 'LineWidth', 1);
        line([(cfPOS(idx_t) - fixPOS(idx_t)), (cfPOS(idx_t) - fixPOS(idx_t))]/sampleRate, [min(dif), max(dif)], 'Color', 'blue', 'LineWidth', 1);
        text(t_signal(ceil(0.7*size(t_signal,2))), mean(dif(cfPOS(idx_t)-trial_start(idx_t):cfPOS(idx_t)-trial_start(idx_t)+cfDUR(idx_t)-1)), ['mean diff cf: ' num2str(mean(dif))])
        hold off;
        title([ 'Cue: ' num2str(trial_cue(idx_t))]);
    end

    sgtitle(['file: ' files(idx_f).name ' | PSD']);
end
