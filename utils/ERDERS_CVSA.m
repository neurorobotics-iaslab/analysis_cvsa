function [ERD, cueTYP, minDur] = ERDERS_CVSA(signal, header, channels_label_selected, lap_path, sampleRate, chanlog_path, classes, normalization, title)
nchannels = length(channels_label_selected);

%% Extract infoo data
signal = signal(:, 1:nchannels);
[fixPOS, fixDUR, cuePOS, cueDUR, cueTYP, cfPOS, cfDUR, n_trial] = extract_info_data(header, 786, [730 731], 781);

%% Extract trial infoo
[trialStart, trialStop, fixStart, fixStop, ck, tk] = extract_trial_infoo(signal, fixPOS, ...
    fixDUR, cfPOS, cfDUR, cueTYP, n_trial);

%% Processing
band = [8 14];
avg = 1;
filtOrder = 4;
signal = processing_offline(signal, lap_path, nchannels, sampleRate, band, filtOrder, avg);

%% ERD/ERS
[ERD, minDur, minDurFix] = compute_ERDERS(signal, trialStart, trialStop, fixStart, fixStop, nchannels, n_trial, ck, normalization);

%% Visualization
if ~isempty(chanlog_path)
    showTopoplot_ERDERS(chanlog_path, ERD, sampleRate, channels_label_selected, band, minDur, cueTYP, classes, nchannels, title)
end
end