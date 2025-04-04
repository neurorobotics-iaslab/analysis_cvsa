%% files for merging of two classifiers
%% compute the fischer to understand features using the log band
clc; clearvars;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/utils')

%% Initialization
a = 4:2:18;
b = a+2;
c = [a; b];
bands = cell(1, size(a, 2));
for i=1:length(a)
    bands{i} = [a(i), b(i)];
end
% bands = {[6 8], [8 10], [10 12], [12 14], [14 16]};
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
nbands = length(bands);
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
signals = cell(1, nbands);
headers = cell(1, nbands);
for idx_band = 1:nbands
    headers{idx_band}.TYP = [];
    headers{idx_band}.DUR = [];
    headers{idx_band}.POS = [];
    signals{idx_band} = [];
end
classes = [730 731];
cf_event = 781;
fix_event = 786;
nchannels = 39;
nclasses = length(classes);
filterOrder = 4;
avg = 1;

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);
%% concatenate the files
nFiles = length(filenames);
for idx_file= 1: nFiles
    fullpath_file_shift = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_shift);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;

     for idx_band = 1:nbands
        band = bands{idx_band};

        signal_processed = proc_512hz(c_signal, header.SampleRate, band, filterOrder, avg);
        
        c_header = headers{idx_band};
        c_header.sampleRate = header.SampleRate;
        c_header.channels_labels = header.Label;
        if isempty(find(header.EVENT.TYP == 2, 1)) % no eye calibration
            c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP);
            c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR);
            c_header.POS = cat(1, c_header.POS, header.EVENT.POS + size(signals{idx_band}, 1));
        else
            k = find(header.EVENT.TYP == 1, 1);
            c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP(k:end));
            c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR(k:end));
            c_header.POS = cat(1, c_header.POS, header.EVENT.POS(k:end) + size(signals{idx_band}, 1));
        end

        signals{idx_band} = cat(1, signals{idx_band}, signal_processed(:,:));
        headers{idx_band} = c_header;
    end
end


%% Labelling data 
[~, channelsSelected] = ismember(channels_select, channels_label);
nchannelsSelected = size(channelsSelected, 2);
events = headers{1};
sampleRate = events.sampleRate;
cuePOS = events.POS(ismember(events.TYP, classes));
cueDUR = events.DUR(ismember(events.TYP, classes));
cueTYP = events.TYP(ismember(events.TYP, classes));

fixPOS = events.POS(events.TYP == 786);
fixDUR = events.DUR(events.TYP == 786);

cfPOS = events.POS(events.TYP == 781);
cfDUR = events.DUR(events.TYP == 781);

%% Labeling data for the dataset
data = cell(1, nbands);
for idx_band=1:nbands
    c_signal = signals{idx_band};
    c_header = headers{idx_band};

    c_data = extract_all(c_signal, c_header, classes, fix_event, cf_event, 1);
    data{idx_band} = c_data;
end

%% using the qda perform the classification
[qda_filenames_shift, pathname_shift] = uigetfile('*.yaml', 'Select QDA SHIFT Files', 'MultiSelect', 'on');
if ischar(qda_filenames_shift)
    qda_filenames_shift = {qda_filenames_shift};
end
numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), qda_filenames_shift);
[~, sortedIdx] = sort(numbers);
qda_filenames_shift = qda_filenames_shift(sortedIdx);

[qda_filenames_sustained, pathname_manteined] = uigetfile('*.yaml', 'Select QDA SUSTAINED Files', 'MultiSelect', 'on');
if ischar(qda_filenames_sustained)
    qda_filenames_sustained = {qda_filenames_sustained};
end
numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), qda_filenames_sustained);
[~, sortedIdx] = sort(numbers);
qda_filenames_sustained = qda_filenames_sustained(sortedIdx);

%% compute the classificaitons
nqda = size(qda_filenames_shift, 2);
classified_data_shift = cell(1, nqda);
classified_data_sustained = cell(1, nqda);
classified_info_shift = cell(1, nqda);
classified_info_sustained = cell(1, nqda);

for idx_qda = 1:nqda
    % reasoning for shift classifier
    fullpath_file_shift = fullfile(pathname_shift, qda_filenames_shift{idx_qda});
    disp(['QDA shift: ' qda_filenames_shift{idx_qda}]);
    qda_shift = loadQDA(fullpath_file_shift);
    if qda_shift.nfeatures == 1
        c_idx_band = find(cellfun(@(x) isequal(x, qda_shift.bands), bands));
        selectedFeatures = [qda_shift.idchans, c_idx_band];
    else
        c_idx_band = [];
        for i = 1:length(qda_shift.bands) 
            c_idx_band = [c_idx_band find(ismember(cell2mat(bands'), qda_shift.bands(i,:), 'rows'))];
        end
        selectedFeatures = [qda_shift.idchans', c_idx_band'];
    end
    [X, ~, classified_info_shift{idx_qda}] = createDataset(filenames, data, selectedFeatures, bands, channels_label);

    classified_tmp = apply_qda_matrix(qda_shift, X);
    classified_data_shift{idx_qda} = classified_tmp;

    % reasoning for sustained classifier
    fullpath_file_sustained = fullfile(pathname_manteined, qda_filenames_sustained{idx_qda});
    disp(['QDA sustained ' qda_filenames_sustained{idx_qda}]);
    qda_sustained = loadQDA(fullpath_file_sustained);
    if qda_sustained.nfeatures == 1
        c_idx_band = find(cellfun(@(x) isequal(x, qda_sustained.bands), bands));
        selectedFeatures = [qda_sustained.idchans, c_idx_band];
    else
        c_idx_band = [];
        for i = 1:length(qda_sustained.bands) 
            c_idx_band = [c_idx_band find(ismember(cell2mat(bands'), qda_sustained.bands(i,:), 'rows'))];
        end
        selectedFeatures = [qda_sustained.idchans', c_idx_band'];
    end
    [X, ~, classified_info_sustained{idx_qda}] = createDataset(filenames, data, selectedFeatures, bands, channels_label);

    classified_tmp = apply_qda_matrix(qda_sustained, X);
    classified_data_sustained{idx_qda} = classified_tmp;
end

%% compute for the accuracy from start to x (only for the SHIFT)
time = 1; %s
samples = time * sampleRate;
acc_shift = nan(1, nqda);
start_trial_id = 41; % trial at which start the test set
train_trial = 1:start_trial_id-1;
ntrain_trial = length(train_trial);
for idx_qda = 1:nqda
    y_pred = nan(samples*ntrain_trial, 1);
    y_true = nan(samples*ntrain_trial, 1);
    qda_shift = nan(ntrain_trial, samples);
    for idx_test_trial = 1:ntrain_trial
        idx_trial = train_trial(idx_test_trial);
        c_start = classified_info_shift{idx_qda}.startTrial(idx_trial);
        c_end = c_start + samples - 1;
        c_data = classified_data_shift{idx_qda};
        qda_performance = c_data(c_start:c_end,1);
        qda_shift(idx_test_trial,:) = qda_performance;

        c_pred = ones(samples, 1) * classes(1);
        c_pred(qda_performance < 0.5) = classes(2);
        y_pred((idx_test_trial-1)*samples+1:idx_test_trial*samples) = c_pred;

        c_true = repmat(data{1}.typ(idx_trial), samples, 1);
        y_true((idx_test_trial-1)*samples+1:idx_test_trial*samples) = c_true;
    end

    acc_shift(idx_qda) = sum(y_pred == y_true) / length(y_true)*100;
end

%% compute for the accuracy form x to end (only for the sustained)
time = 1; %s
min_durCF = size(data{1}.cf, 1);
samples = min_durCF - time * sampleRate;
acc_sustained = nan(1, nqda);
for idx_qda = 1:nqda
    y_pred = nan(samples*ntrain_trial, 1);
    y_true = nan(samples*ntrain_trial, 1);
    for idx_test_trial = 1:ntrain_trial
        idx_trial = train_trial(idx_test_trial);
        c_start = classified_info_sustained{idx_qda}.startTrial(idx_trial) + time * sampleRate;
        c_end = c_start + samples - 1;
        c_data = classified_data_sustained{idx_qda};
        qda_performance = c_data(c_start:c_end,1);

        c_pred = ones(samples, 1) * classes(1);
        c_pred(qda_performance < 0.5) = classes(2);
        y_pred((idx_test_trial-1)*samples+1:idx_test_trial*samples) = c_pred;

        c_true = repmat(data{1}.typ(idx_trial), samples, 1);
        y_true((idx_test_trial-1)*samples+1:idx_test_trial*samples) = c_true;
    end

    acc_sustained(idx_qda) = sum(y_pred == y_true) / length(y_true)*100;
end

%% take the two best classifier (on train) and merge them
ntrial = size(classified_info_sustained{1}.startTrial, 1);
end_trial_id = ntrial;
test_trial = start_trial_id:end_trial_id;
ntest_trial = length(test_trial);
[val_qda_shift, idx_qda_shift] = max(acc_shift);
[val_qda_sustained, idx_qda_sustained] = max(acc_sustained);
disp('TRAIN: ')
disp(['acc shift: ' num2str(val_qda_shift) ' nfeatures: ' num2str(idx_qda_shift) ' | acc sustained: ' num2str(val_qda_sustained) ' nfeatures: ' num2str(idx_qda_sustained)]);

data_shift = classified_data_shift{idx_qda_shift};
data_sustained = classified_data_sustained{idx_qda_sustained};
qda_performance = nan(ntest_trial, min_durCF);
% first column is weight for shift, second for mantain
weight = [ones(time*sampleRate, 1), zeros(time*sampleRate, 1); 
          zeros(min_durCF - time*sampleRate, 1), ones(min_durCF - time*sampleRate, 1)];

y_label = cell(ntest_trial, 1);
y_true_total = nan(ntest_trial*min_durCF, 1);
y_pred_total = nan(ntest_trial*min_durCF, 1);
y_pred_shift = nan(ntest_trial*time*sampleRate, 1);
y_true_shift = nan(ntest_trial*time*sampleRate, 1);
y_pred_sustained = nan(ntest_trial*(min_durCF - time*sampleRate), 1);
y_true_sustained = nan(ntest_trial*(min_durCF - time*sampleRate), 1);
correctness = nan(ntest_trial, min_durCF);
for idx_test_trial = 1:ntest_trial
    idx_trial = test_trial(idx_test_trial);
    c_start = classified_info_shift{1}.startTrial(idx_trial);
    c_end = c_start + min_durCF - 1;

    qda_shift = weight(:,1) .* data_shift(c_start:c_end, 1);
    qda_sustained = weight(:,2) .* data_sustained(c_start:c_end, 1);

    qda_performance(idx_test_trial,:) = qda_shift + qda_sustained;

    % compute the new accuracy for the merge of the two classifiers
    y_pred = ones(min_durCF, 1) * classes(1);
    y_pred(qda_performance(idx_test_trial,:) < 0.5) = classes(2);
    y_true = repmat(data{1}.typ(idx_trial), min_durCF, 1);
    y_label{idx_test_trial} = [num2str(data{1}.typ(idx_trial)) ': ' sprintf('%.2f', sum(y_pred == y_true) / min_durCF*100)];

    % compute the total accuracy for the merge of the two classifiers
    y_true_total((idx_test_trial-1)*min_durCF+1:idx_test_trial*min_durCF) = y_true;
    y_pred_total((idx_test_trial-1)*min_durCF+1:idx_test_trial*min_durCF) = y_pred;

    % compute the accuracy for the shift and for the sustained part sparetly
    y_true_shift((idx_test_trial-1)*time*sampleRate+1:idx_test_trial*time*sampleRate) = repmat(data{1}.typ(idx_trial), time*sampleRate, 1);
    temp = qda_shift(qda_shift ~= 0);
    tmp_pred = ones(time*sampleRate, 1) * classes(1);
    tmp_pred(temp < 0.5) = classes(2);
    y_pred_shift((idx_test_trial-1)*time*sampleRate+1:idx_test_trial*time*sampleRate) = tmp_pred;
    y_true_sustained((idx_test_trial-1)*(min_durCF-time*sampleRate)+1:idx_test_trial*(min_durCF-time*sampleRate)) = repmat(data{1}.typ(idx_trial), min_durCF-time*sampleRate, 1);
    temp = qda_sustained(qda_sustained ~= 0);
    tmp_pred = ones(min_durCF-time*sampleRate, 1) * classes(1);
    tmp_pred(temp < 0.5) = classes(2);
    y_pred_sustained((idx_test_trial-1)*(min_durCF-time*sampleRate)+1:idx_test_trial*(min_durCF-time*sampleRate)) = tmp_pred;

    correctness(idx_test_trial,:) = y_pred == y_true;
end
acc_total = sum(y_pred_total == y_true_total) / length(y_pred_total) * 100;
acc_shift = sum(y_pred_shift == y_true_shift) / length(y_pred_shift) * 100;
acc_sustained = sum(y_pred_sustained == y_true_sustained) / length(y_pred_sustained) * 100;

disp('TEST:');
disp(['acc shift: ' num2str(acc_shift) ' | acc sustained: ' num2str(acc_sustained)]);

figure();
qda_performance(qda_performance < 0.3) = 0.3;
green_component = correctness .* qda_performance; % if prediction correct
red_component = (~correctness) .* qda_performance; %
color_map = cat(3, red_component, green_component, zeros(size(green_component))); % No blue component
imagesc(color_map)
yticks(1:ntrial)
yticklabels(y_label);
xticks(1:256:min_durCF)
xticklabels(((1:256:min_durCF) - 1)  / sampleRate)
title([subject ' | total: ' num2str(acc_total) ' | shift: ' num2str(acc_shift) ...
    ' | sustained: ' num2str(acc_sustained)])
