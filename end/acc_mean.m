%% files for see the mean of the raw probabilities on the test set, in addition look the integrated prob
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
eog_threshold = 500;

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);
%% concatenate the files
nFiles = length(filenames);
trial_with_eog = [];
for idx_file= 1: nFiles
    fullpath_file_shift = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_shift);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;

    c_trial_with_eog = eog_detection(c_signal, header, eog_threshold, {'FP1', 'FP2', 'EOG'});
    trial_with_eog = [trial_with_eog; c_trial_with_eog];

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

minDurCue = min(cueDUR);
ntrial = length(cuePOS);

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

%% refactoring the data
% now the data are placed alternatevely, so odd trial id is for 730, even for 731
balanced_trial_idx = balanced_data_afterEOG(trial_with_eog, trial_typ, classes);
balanced_trial_data = trial_data(:,:,:,logical(balanced_trial_idx));
trial_typ = trial_typ(logical(balanced_trial_idx));
ntrial = sum(balanced_trial_idx);
idx_classes_trial = nan(ntrial/2, nclasses);
for idx_class = 1:nclasses
    idx_classes_trial(:,idx_class) = find(trial_typ == classes(idx_class));
end

tmp = nan(size(balanced_trial_data));
trial_typ = nan(size(trial_typ));
i = 1;
for idx_trial_class = 1:2:ntrial
    for idx_class = 1:nclasses
        tmp(:,:,:,idx_trial_class + idx_class - 1) = balanced_trial_data(:,:,:,idx_classes_trial(i, idx_class));
        trial_typ(idx_trial_class + idx_class - 1) = classes(idx_class);
    end
    i = i + 1;
end
trial_data = tmp;
percentual_test = 0.2;
ntrial_test = ceil(percentual_test*ntrial/2) * 2; % total trial for test

ntrial_train = ntrial - ntrial_test;
trial_typ_train  = trial_typ(1:ntrial_train);
trial_typ_test = trial_typ(ntrial_train:end);
trial_train = 1:ntrial_train;

%% using the qda perform the classification
[qda_filenames_shift, pathname_shift] = uigetfile('*.yaml', 'Select QDA SHIFT Files', 'MultiSelect', 'on');
if ischar(qda_filenames_shift)
    qda_filenames_shift = {qda_filenames_shift};
else
    numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), qda_filenames_shift);
    [~, sortedIdx] = sort(numbers);
    qda_filenames_shift = qda_filenames_shift(sortedIdx);
end

pathname_sustained = [pathname_shift(1:end-6) 'sustained/'];
qda_filenames_sustained = {'cl_best_sustained.yaml'};

pathname_cf = [pathname_shift(1:end-6) 'cf/'];
qda_filenames_cf = {'cl_best_cf.yaml'};

% [qda_filenames_sustained, pathname_sustained] = uigetfile('*.yaml', 'Select QDA sustained Files', 'MultiSelect', 'on');
% if ischar(qda_filenames_sustained)
%     qda_filenames_sustained = {qda_filenames_sustained};
% else
% 
%     numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), qda_filenames_sustained);
%     [~, sortedIdx] = sort(numbers);
%     qda_filenames_sustained = qda_filenames_sustained(sortedIdx);
% end

%% compute the classification for each qda (shift and sustained)
nqda = size(qda_filenames_shift, 2);
classified_data_shift = cell(1, nqda);
classified_data_sustained = cell(1, nqda);
classified_data_cf = cell(1, nqda);
classified_start_trial_shift = nan(ntrial, nqda);
classified_start_trial_sustained = nan(ntrial, nqda);
classified_start_trial_cf = nan(ntrial, nqda);

for idx_qda = 1:nqda
    % reasoning for SHIFT classifier
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

    % reasoning for SUSTAINED classifier
    fullpath_file_sustained = fullfile(pathname_sustained, qda_filenames_sustained{idx_qda});
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

     % reasoning for CF classifier
    fullpath_file_cf = fullfile(pathname_cf, qda_filenames_cf{idx_qda});
    disp(['QDA cf ' qda_filenames_cf{idx_qda}]);
    qda_cf = loadQDA(fullpath_file_cf);
    if qda_cf.nfeatures == 1
        c_idx_band = find(cellfun(@(x) isequal(x, qda_cf.bands), bands));
        selectedFeatures = [qda_cf.idchans, c_idx_band];
    else
        c_idx_band = [];
        for i = 1:length(qda_cf.bands) 
            c_idx_band = [c_idx_band find(ismember(cell2mat(bands'), qda_cf.bands(i,:), 'rows'))];
        end
        selectedFeatures = [qda_cf.idchans', c_idx_band'];
    end

    % create the dataset and perform the classification
    X = nan(min_trial_data*ntrial, 1);
    for idx_trial = 1:ntrial
        for idx_f = 1:size(selectedFeatures, 1)
            c_channel = selectedFeatures(idx_f, 1);
            c_band = selectedFeatures(idx_f, 2);
            X((idx_trial-1)*min_trial_data+1:idx_trial*min_trial_data, idx_f) = trial_data(:,c_band,c_channel,idx_trial);
        end

        classified_start_trial_cf(idx_trial, idx_qda) = (idx_trial-1)*min_trial_data +1;
    end
    classified_tmp = apply_qda_matrix(qda_cf, X);
    classified_data_cf{idx_qda} = classified_tmp;
end

%% check the number of classifier
if nqda == 1
    idx_qda_shift = 1;
    idx_qda_sus = 1;
    idx_qda_cf = 1;
end

%% show the best qda AUC for the shift in the interest part
% interval = [50 700]; % ms
% interval_sample = ceil(interval * sampleRate / 1000);
interval_sample = [qda_shift.interval_samples{1} qda_shift.interval_samples{2}];
samples = (interval_sample(2) - interval_sample(1)) + 1;
y_true_train = nan(samples * ntrial_train,1);
y_pred_train = nan(samples * ntrial_train,1);
for idx_trial=1:ntrial_train
    c_start = classified_start_trial_shift(idx_trial) + interval_sample(1);
    c_end = c_start + samples - 1;
    qda_performance = classified_data_shift{idx_qda_shift}(c_start:c_end,:);

    c_y_true = repmat(trial_typ(idx_trial), samples, 1);
    c_y_true(c_y_true == classes(1)) = 0;
    c_y_true(c_y_true == classes(2)) = 1;
    y_true_train((idx_trial-1)*samples+1:idx_trial*samples) = c_y_true;

    % qda_performance of 2 since the qda uses 0 for calss 730 and 1 for class 731
    y_pred_train((idx_trial-1)*samples+1:idx_trial*samples) = qda_performance(:,2);
end

y_true_test = nan(samples * ntrial_test, 1);
y_pred_test = nan(samples * ntrial_test, 1);
trial_test = ntrial_train + 1:ntrial_train+ntrial_test;
for idx_trial=1:ntrial_test
    c_trial = trial_test(idx_trial);
    c_start = classified_start_trial_shift(c_trial) + interval_sample(1);
    c_end = c_start + samples - 1;
    qda_performance = classified_data_shift{idx_qda_shift}(c_start:c_end,:);

    c_y_true = repmat(trial_typ(c_trial), samples, 1);
    c_y_true(c_y_true == classes(1)) = 0;
    c_y_true(c_y_true == classes(2)) = 1;
    y_true_test((idx_trial-1)*samples+1:idx_trial*samples) = c_y_true;

    % qda_performance of 2 since the qda uses 0 for calss 730 and 1 for class 731
    y_pred_test((idx_trial-1)*samples+1:idx_trial*samples) = qda_performance(:,2);
end

fig = figure('Units','normalized','OuterPosition',[0 0 1 1]);
subplot(1,2,1);
[xROC_shift_train, yROC_shift_train, ~, AUC_shift_train] = perfcurve(y_true_train, y_pred_train, 1);
plot(xROC_shift_train, yROC_shift_train, 'b-', 'LineWidth', 2);
hold on;
plot([0 1], [0 1], 'r--'); % Random classifier line
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(['ROC Curve (AUC = ' num2str(AUC_shift_train) ') | TRAIN']);
grid on;

subplot(1,2,2);
[xROC_shift_test, yROC_shift_test, ~, AUC_shift_test] = perfcurve(y_true_test, y_pred_test, 1);
plot(xROC_shift_test, yROC_shift_test, 'b-', 'LineWidth', 2);
hold on;
plot([0 1], [0 1], 'r--'); % Random classifier line
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(['ROC Curve (AUC = ' num2str(AUC_shift_test) ') | TEST']);
grid on;
sgtitle('SHIFT')
path_figure = ['/home/paolo/cvsa_ws/img_results/27_features/' subject '/calibFiles_acc/auc_shift.svg'];
print(fig, path_figure, '-dsvg');


%% show the best qda auc for the sustained in the interest part
% only the periodo taken into account
% interval = [1500 2700]; % ms --> keep into consideration the cue period (wxw)
% interval_sample = ceil(interval * sampleRate / 1000); % convert ms to sample
interval_sample = [qda_sustained.interval_samples{1} qda_sustained.interval_samples{2}];
samples = interval_sample(2) - interval_sample(1) + 1;
y_true_train = nan(samples * ntrial_train,1);
y_pred_train = nan(samples * ntrial_train,1);
for idx_trial=1:ntrial_train
    c_start = classified_start_trial_sustained(idx_trial) + interval_sample(1);
    c_end = c_start + samples - 1;
    qda_performance = classified_data_sustained{idx_qda_sus}(c_start:c_end,:);

    c_y_true = repmat(trial_typ(idx_trial), samples, 1);
    c_y_true(c_y_true == classes(1)) = 0;
    c_y_true(c_y_true == classes(2)) = 1;
    y_true_train((idx_trial-1)*samples+1:idx_trial*samples) = c_y_true;

    % qda_performance of 2 since the qda uses 0 for calss 730 and 1 for class 731
    y_pred_train((idx_trial-1)*samples+1:idx_trial*samples) = qda_performance(:,2);
end

y_true_test = nan(samples * ntrial_test, 1);
y_pred_test = nan(samples * ntrial_test, 1);
trial_test = ntrial_train + 1:ntrial_train+ntrial_test;
for idx_trial=1:ntrial_test
    c_trial = trial_test(idx_trial);
    c_start = classified_start_trial_sustained(c_trial) + interval_sample(1);
    c_end = c_start + samples - 1;
    qda_performance = classified_data_sustained{idx_qda_sus}(c_start:c_end,:);

    c_y_true = repmat(trial_typ(c_trial), samples, 1);
    c_y_true(c_y_true == classes(1)) = 0;
    c_y_true(c_y_true == classes(2)) = 1;
    y_true_test((idx_trial-1)*samples+1:idx_trial*samples) = c_y_true;

    % qda_performance of 2 since the qda uses 0 for calss 730 and 1 for class 731
    y_pred_test((idx_trial-1)*samples+1:idx_trial*samples) = qda_performance(:,2);
end

fig = figure('Units','normalized','OuterPosition',[0 0 1 1]);
subplot(1,2,1);
[xROC_shift_train, yROC_shift_train, ~, AUC_shift_train] = perfcurve(y_true_train, y_pred_train, 1);
plot(xROC_shift_train, yROC_shift_train, 'b-', 'LineWidth', 2);
hold on;
plot([0 1], [0 1], 'r--'); % Random classifier line
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(['ROC Curve (AUC = ' num2str(AUC_shift_train) ') | TRAIN']);
grid on;

subplot(1,2,2);
[xROC_shift_test, yROC_shift_test, ~, AUC_shift_test] = perfcurve(y_true_test, y_pred_test, 1);
plot(xROC_shift_test, yROC_shift_test, 'b-', 'LineWidth', 2);
hold on;
plot([0 1], [0 1], 'r--'); % Random classifier line
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(['ROC Curve (AUC = ' num2str(AUC_shift_test) ') | TEST']);
grid on;
sgtitle('SUSTAINED')
path_figure = ['/home/paolo/cvsa_ws/img_results/27_features/' subject '/calibFiles_acc/auc_sostained.svg'];
print(fig, path_figure, '-dsvg');

%% integrate the two best qda performances 
alpha = 0.98;
start_integrator = 0.5;
qda_best_shift_integrated_soft = ones(ntrial_train, min_trial_data + 1) * start_integrator;
qda_best_sustained_integrated_soft = ones(ntrial_train, min_trial_data + 1) * start_integrator;
qda_best_cf_integrated_soft = ones(ntrial_train, min_trial_data + 1) * start_integrator;
qda_best_shift_integrated_hard = ones(ntrial_train, min_trial_data + 1) * start_integrator;
qda_best_sustained_integrated_hard = ones(ntrial_train, min_trial_data + 1) * start_integrator;
qda_best_cf_integrated_hard = ones(ntrial_train, min_trial_data + 1) * start_integrator;
for idx_train_trial=1:ntrial_train
    idx_trial = trial_train(idx_train_trial);
    c_start = classified_start_trial_shift(idx_trial, idx_qda_shift); 
    c_end = c_start + min_trial_data - 1;
    idx_class = find(trial_typ(idx_trial) == classes);

    c_data_shift = classified_data_shift{idx_qda_shift};
    c_shift_soft = c_data_shift(c_start:c_end,idx_class);
    c_shift_hard = ones(min_trial_data, 1);
    c_shift_hard(c_shift_soft < 0.5) = 0;
    c_shift_hard(c_shift_soft == 0.5) = 0.5;

    c_data_sustained = classified_data_sustained{idx_qda_sus};
    c_sustained_soft = c_data_sustained(c_start:c_end,idx_class);
    c_sustained_hard = ones(min_trial_data, 1);
    c_sustained_hard(c_sustained_soft < 0.5) = 0;
    c_sustained_hard(c_sustained_soft == 0.5) = 0.5;

    c_data_cf = classified_data_cf{idx_qda_cf};
    c_cf_soft = c_data_cf(c_start:c_end, idx_class);
    c_cf_hard = ones(min_trial_data, 1);
    c_cf_hard(c_cf_soft < 0.5) = 0;
    c_cf_hard(c_cf_soft == 0.5) = 0.5;

    % compute the integration
    for idx_sample=1:min_trial_data
        % compute with soft predict
        qda_best_shift_integrated_soft(idx_train_trial, idx_sample + 1) = ...
            qda_best_shift_integrated_soft(idx_train_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_soft(idx_sample);
        qda_best_sustained_integrated_soft(idx_train_trial, idx_sample + 1) = ...
            qda_best_sustained_integrated_soft(idx_train_trial, idx_sample)* alpha + (1.0-alpha)*c_sustained_soft(idx_sample);
        qda_best_cf_integrated_soft(idx_train_trial, idx_sample + 1) = ...
            qda_best_cf_integrated_soft(idx_train_trial, idx_sample)*alpha + (1.0-alpha)*c_cf_soft(idx_sample);

        % compute for hard predict
        qda_best_shift_integrated_hard(idx_train_trial, idx_sample + 1) = ...
            qda_best_shift_integrated_hard(idx_train_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_hard(idx_sample);
        qda_best_sustained_integrated_hard(idx_train_trial, idx_sample + 1) = ...
            qda_best_sustained_integrated_hard(idx_train_trial, idx_sample)* alpha + (1.0-alpha)*c_sustained_hard(idx_sample);
        qda_best_cf_integrated_hard(idx_train_trial, idx_sample + 1) = ...
            qda_best_cf_integrated_hard(idx_train_trial, idx_sample)*alpha + (1.0-alpha)*c_cf_hard(idx_sample);
    end
end


%% Save the raws probabilities of both classifiers for the train set
qda_shift_raw_train = nan(ntrial_train, min_trial_data);
qda_cf_raw_train = nan(ntrial_train, min_trial_data);
qda_sus_raw_train = nan(ntrial_train, min_trial_data);
c_data_shift = classified_data_shift{idx_qda_shift};
c_data_cf = classified_data_cf{idx_qda_cf};
c_data_sus = classified_data_sustained{idx_qda_sus};
for idx_trial_train = 1:ntrial_train
    idx_trial = trial_train(idx_trial_train);
    c_start = classified_start_trial_shift(idx_trial, idx_qda_shift);
    c_end = c_start + min_trial_data - 1;
    idx_class = find(trial_typ(idx_trial) == classes);
    qda_shift_raw_train(idx_trial_train, :) = c_data_shift(c_start:c_end,idx_class);

    c_start = classified_start_trial_sustained(idx_trial, idx_qda_sus);
    c_end = c_start + min_trial_data - 1;
    qda_sus_raw_train(idx_trial_train, :) = c_data_sus(c_start:c_end,idx_class);

    c_start = classified_start_trial_cf(idx_trial, idx_qda_cf);
    c_end = c_start + min_trial_data - 1;
    qda_cf_raw_train(idx_trial_train, :) = c_data_cf(c_start:c_end,idx_class);
end

%% plot the best QDA for the two period
% mean and std of the raw probabilities
fig1 = figure('Units','normalized','OuterPosition',[0 0 1 1]);
subplot(1,3,1)
hold on;
plot(mean(qda_shift_raw_train, 1), 'Color','b');
plot(mean(qda_shift_raw_train, 1) + 0.2*std(qda_shift_raw_train, 0, 1), '--', 'Color','b');
plot(mean(qda_shift_raw_train, 1) - 0.2*std(qda_shift_raw_train, 0, 1), '--', 'Color','b');
plot(mean(qda_sus_raw_train, 1), 'Color','r');
plot(mean(qda_sus_raw_train, 1) + 0.2*std(qda_sus_raw_train, 0, 1), '--', 'Color','r');
plot(mean(qda_sus_raw_train, 1) - 0.2*std(qda_sus_raw_train, 0, 1), '--', 'Color','r');
plot(mean(qda_cf_raw_train, 1), 'Color','k');
plot(mean(qda_cf_raw_train, 1) + 0.2*std(qda_cf_raw_train, 0, 1), '--', 'Color','k');
plot(mean(qda_cf_raw_train, 1) - 0.2*std(qda_cf_raw_train, 0, 1), '--', 'Color','k');
hold off;
legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', 'QDA sustained mean', ...
    'QDA sustained + std', 'QDA sustained - std', 'QDA cf mean', 'QDA cf + std', 'Qda cf - std'});
xlim([0 min_trial_data])
ylim([0.0 1.0])
xticks(1:256:min_trial_data)
xticklabels(((1:256:min_trial_data) - 1)  / sampleRate)
title('raws probabilities')

% mean and std of integrated probabilities + the merged integration
subplot(1,3,2)
hold on;
plot(mean(qda_best_shift_integrated_soft, 1), 'Color','b');
plot(mean(qda_best_shift_integrated_soft, 1) + 0.2*std(qda_best_shift_integrated_soft,0,1), '--', 'Color','b');
plot(mean(qda_best_shift_integrated_soft, 1) - 0.2*std(qda_best_shift_integrated_soft,0,1), '--', 'Color','b');

plot(mean(qda_best_sustained_integrated_soft, 1), 'Color','r');
plot(mean(qda_best_sustained_integrated_soft, 1) + 0.2*std(qda_best_sustained_integrated_soft,0,1), '--', 'Color','r');
plot(mean(qda_best_sustained_integrated_soft, 1) - 0.2*std(qda_best_sustained_integrated_soft,0,1), '--', 'Color','r');

plot(mean(qda_best_cf_integrated_soft, 1), 'Color','k');
plot(mean(qda_best_cf_integrated_soft, 1) + 0.2*std(qda_best_cf_integrated_soft,0,1), '--', 'Color','k');
plot(mean(qda_best_cf_integrated_soft, 1) - 0.2*std(qda_best_cf_integrated_soft,0,1), '--', 'Color','k');
hold off;
legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', ...
    'QDA sustained mean', 'QDA sustained + std', 'QDA sustained - std', ...
    'QDA cf mean', 'QDA cf + std', 'QDA cf - std'});
xlim([0 min_trial_data])
ylim([0.0 1.0])
xticks(1:256:min_trial_data)
xticklabels(((1:256:min_trial_data) - 1)  / sampleRate)
title('integrated with soft probabilities')

% mean and std for the hard (before was only soft)
subplot(1,3,3)
hold on;
plot(mean(qda_best_shift_integrated_hard, 1), 'Color','b');
% plot(mean(qda_best_shift_integrated_hard, 1) + 0.3*std(qda_best_shift_integrated_hard,0,1), '--', 'Color','b');
% plot(mean(qda_best_shift_integrated_hard, 1) - 0.3*std(qda_best_shift_integrated_hard,0,1), '--', 'Color','b');

plot(mean(qda_best_sustained_integrated_hard, 1), 'Color','r');
% plot(mean(qda_best_sustained_integrated_hard, 1) + 0.3*std(qda_best_sustained_integrated_hard,0,1), '--', 'Color','r');
% plot(mean(qda_best_sustained_integrated_hard, 1) - 0.3*std(qda_best_sustained_integrated_hard,0,1), '--', 'Color','r');

plot(mean(qda_best_cf_integrated_hard, 1), 'Color','k');
% plot(mean(qda_merge_integrated_hard, 1) + 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','k');
% plot(mean(qda_merge_integrated_hard, 1) - 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','k');
hold off;
% legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', ...
%     'QDA sustained mean', 'QDA sustained + std', 'QDA sustained - std', ...
%     'QDA merge mean', 'QDA merge + std', 'QDA merge - std'});
legend({'QDA shift mean', ...
    'QDA sustained mean', ...
    'QDA cf mean'});
xlim([0 min_trial_data])
ylim([0.0 1.0])
xticks(1:256:min_trial_data)
xticklabels(((1:256:min_trial_data) - 1)  / sampleRate)
title('integrated with hard probabilities')

sgtitle([subject ' | raw probabilities and integration | TRAIN'])
path_figure = ['/home/paolo/cvsa_ws/img_results/27_features/' subject '/calibFiles_acc/results_train.svg'];
print(fig1, path_figure, '-dsvg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --> TEST set <--
% raw probabilities
c_data_shift = classified_data_shift{idx_qda_shift};
c_data_sus = classified_data_sustained{idx_qda_sus};
c_data_cf = classified_data_cf{idx_qda_cf};
shift_raw_best = nan(ntrial_test, min_trial_data);
sustained_raw_best = nan(ntrial_test, min_trial_data);
cf_raw_best = nan(ntrial_test, min_trial_data);
for idx_test_trial = 1:ntrial_test
    idx_trial = trial_test(idx_test_trial);
    idx_class = find(trial_typ(idx_trial) == classes);

    c_start = classified_start_trial_shift(idx_trial, idx_qda_shift);
    c_end = c_start + min_trial_data - 1;
    shift_raw_best(idx_test_trial,:) = c_data_shift(c_start:c_end, idx_class);

    c_start = classified_start_trial_sustained(idx_trial, idx_qda_sus);
    c_end = c_start + min_trial_data - 1;
    sustained_raw_best(idx_test_trial,:) = c_data_sus(c_start:c_end, idx_class);

    c_start = classified_start_trial_cf(idx_trial, idx_qda_cf);
    c_end = c_start + min_trial_data - 1;
    cf_raw_best(idx_test_trial,:) = c_data_cf(c_start:c_end, idx_class);
end

% integration
qda_best_shift_integrated_soft = ones(ntrial_test, min_trial_data + 1) * start_integrator;
qda_best_sustained_integrated_soft = ones(ntrial_test, min_trial_data + 1) * start_integrator;
qda_best_cf_integrated_soft = ones(ntrial_test, min_trial_data + 1) * start_integrator;
qda_best_shift_integrated_hard = ones(ntrial_test, min_trial_data + 1) * start_integrator;
qda_best_sustained_integrated_hard = ones(ntrial_test, min_trial_data + 1) * start_integrator;
qda_best_cf_integrated_hard = ones(ntrial_test, min_trial_data + 1) * start_integrator;
for idx_test_trial=1:ntrial_test
    idx_trial = trial_test(idx_test_trial);
    c_start = classified_start_trial_shift(idx_trial, idx_qda_shift); 
    c_end = c_start + min_trial_data - 1;
    idx_class = find(trial_typ(idx_trial) == classes);

    c_data_shift = classified_data_shift{idx_qda_shift};
    c_shift_soft = c_data_shift(c_start:c_end,idx_class);
    c_shift_hard = ones(min_trial_data, 1);
    c_shift_hard(c_shift_soft < 0.5) = 0;
    c_shift_hard(c_shift_soft == 0.5) = 0.5;

    c_data_sustained = classified_data_sustained{idx_qda_sus};
    c_sustained_soft = c_data_sustained(c_start:c_end,idx_class);
    c_sustained_hard = ones(min_trial_data, 1);
    c_sustained_hard(c_sustained_soft < 0.5) = 0;
    c_sustained_hard(c_sustained_soft == 0.5) = 0.5;

    c_data_cf = classified_data_cf{idx_qda_cf};
    c_cf_soft = c_data_cf(c_start:c_end,idx_class);
    c_cf_hard = ones(min_trial_data, 1);
    c_cf_hard(c_cf_soft < 0.5) = 0;
    c_cf_hard(c_cf_soft == 0.5) = 0.5;

    % compute the integration
    for idx_sample=1:min_trial_data
        % compute with soft predict
        qda_best_shift_integrated_soft(idx_test_trial, idx_sample + 1) = ...
            qda_best_shift_integrated_soft(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_soft(idx_sample);
        qda_best_sustained_integrated_soft(idx_test_trial, idx_sample + 1) = ...
            qda_best_sustained_integrated_soft(idx_test_trial, idx_sample)* alpha + (1.0-alpha)*c_sustained_soft(idx_sample);
        qda_best_cf_integrated_soft(idx_test_trial, idx_sample + 1) = ...
            qda_best_cf_integrated_soft(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_cf_soft(idx_sample);
        
        % compute for hard predict
        qda_best_shift_integrated_hard(idx_test_trial, idx_sample + 1) = ...
            qda_best_shift_integrated_hard(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_hard(idx_sample);
        qda_best_sustained_integrated_hard(idx_test_trial, idx_sample + 1) = ...
            qda_best_sustained_integrated_hard(idx_test_trial, idx_sample)* alpha + (1.0-alpha)*c_sustained_hard(idx_sample);
        qda_best_cf_integrated_hard(idx_test_trial, idx_sample + 1) = ...
            qda_best_cf_integrated_hard(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_cf_hard(idx_sample);
    end
end

% plot
fig2 = figure('Units','normalized','OuterPosition',[0 0 1 1]);
subplot(1,3,1)
hold on;
plot(mean(shift_raw_best, 1), 'Color','b');
plot(mean(shift_raw_best, 1) + 0.3*std(shift_raw_best, 0, 1), '--', 'Color','b');
plot(mean(shift_raw_best, 1) - 0.3*std(shift_raw_best, 0, 1), '--', 'Color','b');
plot(mean(sustained_raw_best, 1), 'Color','r');
plot(mean(sustained_raw_best, 1) + 0.3*std(sustained_raw_best, 0, 1), '--', 'Color','r');
plot(mean(sustained_raw_best, 1) - 0.3*std(sustained_raw_best, 0, 1), '--', 'Color','r');
plot(mean(cf_raw_best, 1), 'Color','k');
plot(mean(cf_raw_best, 1) + 0.3*std(cf_raw_best, 0, 1), '--', 'Color','k');
plot(mean(cf_raw_best, 1) - 0.3*std(cf_raw_best, 0, 1), '--', 'Color','k');
hold off;
legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', 'QDA sustained mean', ...
    'QDA sustained + std', 'QDA sustained - std', 'QDA cf mean', 'QDA cf + std', 'Qda cf - std'});
xlim([0 min_trial_data])
ylim([0.0 1.0])
xticks(1:256:min_trial_data)
xticklabels(((1:256:min_trial_data) - 1)  / sampleRate)
title('raws probabilities')

% mean and std of integrated probabilities + the merged integration
subplot(1,3,2)
hold on;
plot(mean(qda_best_shift_integrated_soft, 1), 'Color','b');
plot(mean(qda_best_shift_integrated_soft, 1) + 0.3*std(qda_best_shift_integrated_soft,0,1), '--', 'Color','b');
plot(mean(qda_best_shift_integrated_soft, 1) - 0.3*std(qda_best_shift_integrated_soft,0,1), '--', 'Color','b');

plot(mean(qda_best_sustained_integrated_soft, 1), 'Color','r');
plot(mean(qda_best_sustained_integrated_soft, 1) + 0.3*std(qda_best_sustained_integrated_soft,0,1), '--', 'Color','r');
plot(mean(qda_best_sustained_integrated_soft, 1) - 0.3*std(qda_best_sustained_integrated_soft,0,1), '--', 'Color','r');

plot(mean(qda_best_cf_integrated_soft, 1), 'Color','k');
plot(mean(qda_best_cf_integrated_soft, 1) + 0.3*std(qda_best_cf_integrated_soft,0,1), '--', 'Color','k');
plot(mean(qda_best_cf_integrated_soft, 1) - 0.3*std(qda_best_cf_integrated_soft,0,1), '--', 'Color','k');
hold off;
legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', ...
    'QDA sustained mean', 'QDA sustained + std', 'QDA sustained - std', ...
    'QDA merge mean', 'QDA merge + std', 'QDA merge - std'});
xlim([0 min_trial_data])
ylim([0.0 1.0])
xticks(1:256:min_trial_data)
xticklabels(((1:256:min_trial_data) - 1)  / sampleRate)
title('integrated with soft probabilities')

% mean and std for the hard (before was only soft)
subplot(1,3,3)
hold on;
plot(mean(qda_best_shift_integrated_hard, 1), 'Color','b');
% plot(mean(qda_best_shift_integrated_hard, 1) + 0.3*std(qda_best_shift_integrated_hard,0,1), '--', 'Color','b');
% plot(mean(qda_best_shift_integrated_hard, 1) - 0.3*std(qda_best_shift_integrated_hard,0,1), '--', 'Color','b');

plot(mean(qda_best_sustained_integrated_hard, 1), 'Color','r');
% plot(mean(qda_best_sustained_integrated_hard, 1) + 0.3*std(qda_best_sustained_integrated_hard,0,1), '--', 'Color','r');
% plot(mean(qda_best_sustained_integrated_hard, 1) - 0.3*std(qda_best_sustained_integrated_hard,0,1), '--', 'Color','r');

plot(mean(qda_best_cf_integrated_hard, 1), 'Color','k');
% plot(mean(qda_merge_integrated_hard, 1) + 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','k');
% plot(mean(qda_merge_integrated_hard, 1) - 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','k');
hold off;
% legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', ...
%     'QDA sustained mean', 'QDA sustained + std', 'QDA sustained - std', ...
%     'QDA merge mean', 'QDA merge + std', 'QDA merge - std'});
legend({'QDA shift mean', ...
    'QDA sustained mean', ...
    'QDA cf mean'});
xlim([0 min_trial_data])
ylim([0.0 1.0])
xticks(1:256:min_trial_data)
xticklabels(((1:256:min_trial_data) - 1)  / sampleRate)
title('integrated with hard probabilities')

sgtitle([subject ' | raw probabilities and integration | TEST'])
path_figure = ['/home/paolo/cvsa_ws/img_results/27_features/' subject '/calibFiles_acc/results_test.svg'];
print(fig2, path_figure, '-dsvg');

%% plot in general the merge vs the classical appraoch, compute also AUC for both
% understand the point into which perform the switch of the classifiers (using trian)
diff_raw_prob = mean(qda_shift_raw_train, 1) - mean(qda_sus_raw_train, 1);  
% Find zero-crossings
sign_change = diff(sign(diff_raw_prob));
idx_cross = find(sign_change ~= 0);
a = 0; b = min_trial_data;    
idx_cross = idx_cross(idx_cross >= 0.5*sampleRate);
x0 = idx_cross(1);     % Desired intersection point
epsilon = 1e-6;     
% Compute t(x) bounds
Tmin = -2 * (x0 - a) / (b - a);  
Tmax =  2 * (b - x0) / (b - a); 
% Compute required alpha for both ends
alpha_min = -log((1 - epsilon) / epsilon) / Tmin;
alpha_max = -log((1 - epsilon) / epsilon) / Tmax;
% Use the maximum to guarantee both ends are within epsilon
alpha_union = max(alpha_min, alpha_max);
% Normalized sigmoid input
t = 2 * ((1:min_trial_data) - x0) / (b - a);
% Define the sigmoid functions
g = (1 ./ (1 + exp(-alpha_union * t)))';  % Increasing
f = 1 - g;                       % Decreasing

merge_raw_best = nan(ntrial_test, min_trial_data);
qda_merge_integrated_soft = ones(ntrial_test, min_trial_data + 1) * start_integrator;
qda_merge_integrated_hard = ones(ntrial_test, min_trial_data + 1) * start_integrator;
for idx_test_trial=1:ntrial_test
    c_merge = shift_raw_best(idx_test_trial,:)'.*f + sustained_raw_best(idx_test_trial,:)'.*g;
    merge_raw_best(idx_test_trial, :) = c_merge;

    c_merge_hard = ones(min_trial_data, 1);
    c_merge_hard(c_merge < 0.5) = 0;
    c_merge_hard(c_merge == 0.5) = 0.5;

    for idx_sample=1:min_trial_data
        qda_merge_integrated_soft(idx_test_trial, idx_sample + 1) = ...
            qda_merge_integrated_soft(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_merge(idx_sample);
        qda_merge_integrated_hard(idx_test_trial, idx_sample + 1) = ...
            qda_merge_integrated_hard(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_merge_hard(idx_sample);
    end
end

% plot the merge in test wrt classical approach
fig3 = figure();
subplot(1,2,1)
hold on;
plot(mean(qda_merge_integrated_soft, 1), 'Color','b');
plot(mean(qda_merge_integrated_soft, 1) + 0.3*std(qda_merge_integrated_soft,0,1), '--', 'Color','b');
plot(mean(qda_merge_integrated_soft, 1) - 0.3*std(qda_merge_integrated_soft,0,1), '--', 'Color','b');

plot(mean(qda_best_cf_integrated_soft, 1), 'Color','k');
plot(mean(qda_best_cf_integrated_soft, 1) + 0.3*std(qda_best_cf_integrated_soft,0,1), '--', 'Color','k');
plot(mean(qda_best_cf_integrated_soft, 1) - 0.3*std(qda_best_cf_integrated_soft,0,1), '--', 'Color','k');
hold off;
legend({'QDA merge mean', 'QDA merge + std', 'Qda merge - std', ...
    'QDA cf mean', 'QDA cf + std', 'QDA cf - std'});
xlim([0 min_trial_data])
ylim([0.0 1.0])
xticks(1:256:min_trial_data)
xticklabels(((1:256:min_trial_data) - 1)  / sampleRate)
title('integrated with soft probabilities')

% mean and std for the hard (before was only soft)
subplot(1,2,2)
hold on;
plot(mean(qda_merge_integrated_hard, 1), 'Color','b');
% plot(mean(qda_merge_integrated_hard, 1) + 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','b');
% plot(mean(qda_merge_integrated_hard, 1) - 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','b');

plot(mean(qda_best_cf_integrated_hard, 1), 'Color','k');
% plot(mean(qda_merge_integrated_hard, 1) + 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','k');
% plot(mean(qda_merge_integrated_hard, 1) - 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','k');
hold off;
% legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', ...
%     'QDA sustained mean', 'QDA sustained + std', 'QDA sustained - std', ...
%     'QDA merge mean', 'QDA merge + std', 'QDA merge - std'});
legend({'QDA merge mean', ...
    'QDA cf mean'});
xlim([0 min_trial_data])
ylim([0.0 1.0])
xticks(1:256:min_trial_data)
xticklabels(((1:256:min_trial_data) - 1)  / sampleRate)
title('integrated with hard probabilities')

%% save the probabilities
% create the variables to save the values of interest
parameters.g = g; parameters.f=f; parameters.alpha= alpha;
parameters.alpha_union= alpha_union; parameters.typ = trial_typ(trial_test); parameters.intersect = x0;
parameters.features.shift.bands = qda_shift.bands;
parameters.features.shift.channels = qda_shift.chans;
parameters.features.sus.bands = qda_sustained.bands;
parameters.features.sus.channels = qda_sustained.chans;
parameters.features.cf.bands = qda_cf.bands;
parameters.features.cf.channels = qda_cf.chans;
raw.test_cf = nan(ntrial_test, min_trial_data, nclasses);
raw.test_shift = nan(ntrial_test, min_trial_data, nclasses);
raw.test_sus = nan(ntrial_test, min_trial_data, nclasses);
for idx_trial=1:ntrial_test
    c_trial = trial_test(idx_trial);
    c_start = classified_start_trial_shift(c_trial);
    c_end = c_start + min_trial_data - 1;

    raw.test_shift(idx_trial,:,:) = classified_data_shift{idx_qda_shift}(c_start:c_end,:);
    raw.test_sus(idx_trial,:,:) = classified_data_sustained{idx_qda_sus}(c_start:c_end,:);
    raw.test_cf(idx_trial,:,:) = classified_data_cf{idx_qda_cf}(c_start:c_end,:);
end


save(['/home/paolo/cvsa_ws/record/probabilities/' subject '.mat'], 'raw',"parameters")
