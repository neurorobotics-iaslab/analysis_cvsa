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

ntrial = size(cuePOS, 1);

%% Labeling data for the dataset
min_durCF = min(cfDUR);
min_durCUE = min(cueDUR);

trial_start = nan(ntrial, 1);
trial_end = nan(ntrial, 1);
trial_typ = nan(ntrial, 1);
for idx_trial = 1:ntrial
    trial_start(idx_trial) = cuePOS(idx_trial);
    trial_typ(idx_trial) = cueTYP(idx_trial);
    trial_end(idx_trial) = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;
end

min_trial_data = min(trial_end - trial_start+1);
trial_data = nan(min_trial_data, nbands, nchannels, ntrial);
for idx_band = 1:nbands
    c_signal = signals{idx_band};
    for trial = 1:ntrial
        c_start = trial_start(trial);
        c_end = trial_start(trial) + min_trial_data - 1;
        trial_data(:,idx_band,:,trial) = c_signal(c_start:c_end,:);
    end
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
classified_start_trial_shift = nan(ntrial, nqda);
classified_start_trial_sustained = nan(ntrial, nqda);

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

    % create the dataset and perform the classification
    X = nan(min_trial_data*ntrial, 1);
    for idx_trial = 1:ntrial
        for idx_f = 1:size(selectedFeatures, 1)
            c_channel = selectedFeatures(idx_f, 1);
            c_band = selectedFeatures(idx_f, 2);
            X((idx_trial-1)*min_trial_data+1:idx_trial*min_trial_data, idx_f) = trial_data(:,c_band,c_channel,idx_trial);
        end

        classified_start_trial_shift(idx_trial, idx_qda) = (idx_trial-1)*min_trial_data +1;
    end

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

    % create the dataset and perform the classification
    X = nan(min_trial_data*ntrial, 1);
    for idx_trial = 1:ntrial
        for idx_f = 1:size(selectedFeatures, 1)
            c_channel = selectedFeatures(idx_f, 1);
            c_band = selectedFeatures(idx_f, 2);
            X((idx_trial-1)*min_trial_data+1:idx_trial*min_trial_data, idx_f) = trial_data(:,c_band,c_channel,idx_trial);
        end

        classified_start_trial_sustained(idx_trial, idx_qda) = (idx_trial-1)*min_trial_data +1;
    end
    classified_tmp = apply_qda_matrix(qda_sustained, X);
    classified_data_sustained{idx_qda} = classified_tmp;
end

%% Use half of the third file as validation
time = 1; %s

start_trial_id_test_val = 41; % trial at which start the test set
ntrial_val = (ntrial - start_trial_id_test_val + 1) / 2;
ntrial_test = ntrial_val;

trial_val = start_trial_id_test_val:ntrial;
trial_val = trial_val(sort(cat(1, find(cueTYP(start_trial_id_test_val:ntrial) == classes(1), floor(ntrial_val/2), 'first'), find(cueTYP(start_trial_id_test_val:ntrial) == classes(2), floor(ntrial_val/2), 'first'))));
trial_test = start_trial_id_test_val:ntrial;
trial_test = trial_test(sort(cat(1, find(cueTYP(start_trial_id_test_val:ntrial) == classes(1), floor(ntrial_val/2), 'last'), find(cueTYP(start_trial_id_test_val:ntrial) == classes(2), floor(ntrial_val/2), 'last'))));

%% Select the best qda for the shift
samples = time * sampleRate;
accuracy_cv_shift = nan(nqda, 1);

for idx_qda = 1:nqda
    y_pred = nan(samples*length(trial_val), 1);
    y_true = nan(samples*length(trial_val), 1);
    for idx_trial_val = 1:length(trial_val)
        idx_trial = trial_val(idx_trial_val);
        c_start = classified_start_trial_shift(idx_trial, idx_qda) + min_durCUE;
        c_end = c_start + samples - 1;
        c_data = classified_data_shift{idx_qda};
        qda_performance = c_data(c_start:c_end,1);
        c_pred = ones(samples, 1) * classes(1);
        c_pred(qda_performance < 0.5) = classes(2);
        y_pred((idx_trial_val-1)*samples+1:idx_trial_val*samples) = c_pred;

        c_true = repmat(cueTYP(idx_trial), samples, 1);
        y_true((idx_trial_val-1)*samples+1:idx_trial_val*samples) = c_true;

    end
    accuracy_cv_shift(idx_qda) = sum(y_pred == y_true) / length(y_true)*100;
end

%% compute for the accuracy form x to end (only for the sustained)
samples = min_durCF - ceil(time * sampleRate);
accuracy_cv_sustained = nan(nqda, 1);
for idx_qda = 1:nqda
    y_pred = nan(samples*length(trial_val), 1);
    y_true = nan(samples*length(trial_val), 1);
    for idx_trial_val = 1:length(trial_val)
        idx_trial = trial_val(idx_trial_val);
        c_start = classified_start_trial_sustained(idx_trial, idx_qda) + ceil(time * sampleRate) + min_durCUE;
        c_end = c_start + samples - 1;
        c_data = classified_data_sustained{idx_qda};
        qda_performance = c_data(c_start:c_end,1);
        c_pred = ones(samples, 1) * classes(1);
        c_pred(qda_performance < 0.5) = classes(2);
        y_pred((idx_trial_val-1)*samples+1:idx_trial_val*samples) = c_pred;

        c_true = repmat(cueTYP(idx_trial), samples, 1);
        y_true((idx_trial_val-1)*samples+1:idx_trial_val*samples) = c_true;

    end
    accuracy_cv_sustained(idx_qda) = sum(y_pred == y_true) / length(y_true)*100;

end

%% take the two best classifier (on train) and merge them
[val_qda_shift, idx_qda_shift] = max(accuracy_cv_shift);
[val_qda_sustained, idx_qda_sustained] = max(accuracy_cv_sustained);
disp('With validation set: ')
disp(['acc shift: ' num2str(val_qda_shift) ' nfeatures: ' num2str(idx_qda_shift) ' | acc sustained: ' num2str(val_qda_sustained) ' nfeatures: ' num2str(idx_qda_sustained)]);

data_shift = classified_data_shift{idx_qda_shift};
data_sustained = classified_data_sustained{idx_qda_sustained};
qda_performance = nan(ntrial_test, min_durCF);
% first column is weight for shift, second for mantain
weight = [ones(time*sampleRate, 1), zeros(time*sampleRate, 1); 
          zeros(min_durCF - time*sampleRate, 1), ones(min_durCF - time*sampleRate, 1)];

y_label = cell(ntrial_test, 1);
y_true_total = nan(ntrial_test*min_durCF, 1);
y_pred_total = nan(ntrial_test*min_durCF, 1);
y_pred_shift = nan(ntrial_test*time*sampleRate, 1);
y_true_shift = nan(ntrial_test*time*sampleRate, 1);
y_pred_sustained = nan(ntrial_test*(min_durCF - time*sampleRate), 1);
y_true_sustained = nan(ntrial_test*(min_durCF - time*sampleRate), 1);
correctness = nan(ntrial_test, min_durCF);
for idx_test_trial = 1:ntrial_test
    idx_trial = trial_test(idx_test_trial);
    c_start = classified_start_trial_shift(idx_trial, 1) + min_durCUE;
    c_end = c_start + min_durCF - 1;

    qda_shift = weight(:,1) .* data_shift(c_start:c_end, 1);
    qda_sustained = weight(:,2) .* data_sustained(c_start:c_end, 1);

    qda_performance(idx_test_trial,:) = qda_shift + qda_sustained;

    % compute the new accuracy for the merge of the two classifiers
    y_pred = ones(min_durCF, 1) * classes(1);
    y_pred(qda_performance(idx_test_trial,:) < 0.5) = classes(2);
    y_true = repmat(cueTYP(idx_trial), min_durCF, 1);
    y_label{idx_test_trial} = [num2str(cueTYP(idx_trial)) ': ' sprintf('%.2f', sum(y_pred == y_true) / min_durCF*100)];

    % compute the total accuracy for the merge of the two classifiers
    y_true_total((idx_test_trial-1)*min_durCF+1:idx_test_trial*min_durCF) = y_true;
    y_pred_total((idx_test_trial-1)*min_durCF+1:idx_test_trial*min_durCF) = y_pred;

    % compute the accuracy for the shift and for the sustained part sparetly
    y_true_shift((idx_test_trial-1)*time*sampleRate+1:idx_test_trial*time*sampleRate) = repmat(cueTYP(idx_trial), time*sampleRate, 1);
    temp = qda_shift(qda_shift ~= 0);
    tmp_pred = ones(time*sampleRate, 1) * classes(1);
    tmp_pred(temp < 0.5) = classes(2);
    y_pred_shift((idx_test_trial-1)*time*sampleRate+1:idx_test_trial*time*sampleRate) = tmp_pred;
    y_true_sustained((idx_test_trial-1)*(min_durCF-time*sampleRate)+1:idx_test_trial*(min_durCF-time*sampleRate)) = repmat(cueTYP(idx_trial), min_durCF-time*sampleRate, 1);
    temp = qda_sustained(qda_sustained ~= 0);
    tmp_pred = ones(min_durCF-time*sampleRate, 1) * classes(1);
    tmp_pred(temp < 0.5) = classes(2);
    y_pred_sustained((idx_test_trial-1)*(min_durCF-time*sampleRate)+1:idx_test_trial*(min_durCF-time*sampleRate)) = tmp_pred;

    correctness(idx_test_trial,:) = y_pred == y_true;
end
acc_total = sum(y_pred_total == y_true_total) / length(y_pred_total) * 100;
acc_shift = sum(y_pred_shift == y_true_shift) / length(y_pred_shift) * 100;
accuracy_cv_sustained = sum(y_pred_sustained == y_true_sustained) / length(y_pred_sustained) * 100;

disp('TEST:');
disp(['acc shift: ' num2str(acc_shift) ' | acc sustained: ' num2str(accuracy_cv_sustained)]);

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
    ' | sustained: ' num2str(accuracy_cv_sustained)])
