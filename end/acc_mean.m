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
        c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP);
        c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR);
        c_header.POS = cat(1, c_header.POS, header.EVENT.POS + size(signals{idx_band}, 1));

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

[qda_filenames_mantained, pathname_mantained] = uigetfile('*.yaml', 'Select QDA MANTAINED Files', 'MultiSelect', 'on');
if ischar(qda_filenames_mantained)
    qda_filenames_mantained = {qda_filenames_mantained};
end
numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), qda_filenames_mantained);
[~, sortedIdx] = sort(numbers);
qda_filenames_mantained = qda_filenames_mantained(sortedIdx);

%% compute the classification for each qda (shift and mantained)
nqda = size(qda_filenames_shift, 2);
ntrial = size(data{1}.cf,3);
classified_data_shift = cell(1, nqda);
classified_data_mantained = cell(1, nqda);
for i = 1:nqda
    classified_data_shift{i} = nan(size(data{1}.cf, 1)*ntrial, 2);
    classified_data_mantained{i} = nan(size(data{1}.cf, 1)*ntrial, 2);
end
classified_info_shift = cell(1, nqda);
classified_info_mantained = cell(1, nqda);

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
    classified_data_shift{idx_qda}(:,:) = classified_tmp;

    % reasoning for mantained classifier
    fullpath_file_mantained = fullfile(pathname_mantained, qda_filenames_mantained{idx_qda});
    disp(['QDA mantained ' qda_filenames_mantained{idx_qda}]);
    qda_mantained = loadQDA(fullpath_file_mantained);
    if qda_mantained.nfeatures == 1
        c_idx_band = find(cellfun(@(x) isequal(x, qda_mantained.bands), bands));
        selectedFeatures = [qda_mantained.idchans, c_idx_band];
    else
        c_idx_band = [];
        for i = 1:length(qda_mantained.bands) 
            c_idx_band = [c_idx_band find(ismember(cell2mat(bands'), qda_mantained.bands(i,:), 'rows'))];
        end
        selectedFeatures = [qda_mantained.idchans', c_idx_band'];
    end
    [X, ~, classified_info_mantained{idx_qda}] = createDataset(filenames, data, selectedFeatures, bands, channels_label);

    classified_tmp = apply_qda_matrix(qda_mantained, X);
    classified_data_mantained{idx_qda}(:,:) = classified_tmp;
end

%% compute for the accuracy from start to x (only for the SHIFT) and select the best one
time = 1; %s
samples = time * sampleRate;
min_durCF = size(data{1}.cf, 1);
acc_shift = nan(1, nqda);
start_trial_id = 41; % trial at which start the test set
train_trial = 1:start_trial_id-1;
ntrain_trial = length(train_trial);
qda_shift_mean_acc = nan(nqda, min_durCF);
qda_shift_std_acc = nan(nqda, min_durCF);
for idx_qda = 1:nqda
    y_pred = nan(samples*ntrain_trial, 1);
    y_true = nan(samples*ntrain_trial, 1);
    c_mean_acc = nan(ntrain_trial, min_durCF);
    for idx_train_trial = 1:ntrain_trial
        idx_trial = train_trial(idx_train_trial);
        c_start = classified_info_shift{idx_qda}.startTrial(idx_trial);
        c_end = c_start + min_durCF - 1;
        c_data = classified_data_shift{idx_qda};
        qda_performance = c_data(c_start:c_end,1);

        c_pred = ones(samples, 1) * classes(1);
        c_pred(qda_performance(1:samples,1) < 0.5) = classes(2);
        c_pred(qda_performance(1:samples, 1) == 0.5) = nan;
        y_pred((idx_train_trial-1)*samples+1:idx_train_trial*samples) = c_pred;

        c_true = repmat(data{1}.typ(idx_trial), samples, 1);
        y_true((idx_train_trial-1)*samples+1:idx_train_trial*samples) = c_true;

        % take only the correct probs in order to have 1 as correct an 0 as wrong classification
        idx_class = find(data{1}.typ(idx_trial) == classes);
        c_mean_acc(idx_train_trial,:) = c_data(c_start:c_end, idx_class);
    end

    % compute the mean of all the trials
    qda_shift_mean_acc(idx_qda,:) = mean(c_mean_acc, 1);
    qda_shift_std_acc(idx_qda,:) = std(c_mean_acc, 0, 1);

    acc_shift(idx_qda) = sum(y_pred == y_true) / length(y_true)*100;
end

% show the linear plot for each qda (in bigger the best one)
[val_qda_shift, idx_qda_shift] = max(acc_shift);
disp(['SHIFT | Max acc: ' num2str(val_qda_shift) ' | QDA with ' num2str(idx_qda_shift) ' features']);
figure();
hold on
for idx_qda=1:nqda
    if idx_qda == idx_qda_shift
        plot(qda_shift_mean_acc(idx_qda,:), 'LineWidth',5);
    else
        plot(qda_shift_mean_acc(idx_qda, :), 'LineWidth', 1);
    end
end
hold off;
xticks(1:256:min_durCF)
xticklabels(((1:256:min_durCF) - 1)  / sampleRate)
xlim([0 min_durCF]); 
ylim([0.3 1.0]);
title('All QDA shift | raw prob mean trials | TRAIN')


%% compute for the accuracy form x to end (only for the MANTAINED) and select the best one
time = 1; %s
min_durCF = size(data{1}.cf, 1);
samples = min_durCF - time * sampleRate;
acc_mantained = nan(1, nqda);
qda_mantained_mean_acc = nan(nqda, min_durCF);
qda_mantained_std_acc = nan(nqda, min_durCF);
for idx_qda = 1:nqda
    y_pred = nan(samples*ntrain_trial, 1);
    y_true = nan(samples*ntrain_trial, 1);
    c_mean_acc = nan(ntrain_trial, min_durCF);
    for idx_train_trial = 1:ntrain_trial
        idx_trial = train_trial(idx_train_trial);
        c_start = classified_info_shift{idx_qda}.startTrial(idx_trial);
        c_end = c_start + min_durCF - 1;
        c_data = classified_data_mantained{idx_qda};
        qda_performance = c_data(c_start:c_end,1);

        c_pred = ones(samples, 1) * classes(1);
        c_pred(qda_performance(samples:end,1) < 0.5) = classes(2);
        c_pred(qda_performance(samples:end,1) == 0.5) = nan;
        y_pred((idx_train_trial-1)*samples+1:idx_train_trial*samples) = c_pred;

        c_true = repmat(data{1}.typ(idx_trial), samples, 1);
        y_true((idx_train_trial-1)*samples+1:idx_train_trial*samples) = c_true;

        % take only the correct probs in order to have 1 as correct an 0 as wrong classification
        idx_class = find(data{1}.typ(idx_trial) == classes);
        c_mean_acc(idx_train_trial,:) = c_data(c_start:c_end, idx_class);
    end

    % compute the mean of all the trials
    qda_mantained_mean_acc(idx_qda,:) = mean(c_mean_acc, 1);
    qda_mantained_std_acc(idx_qda,:) = std(c_mean_acc, 0, 1);

    acc_mantained(idx_qda) = sum(y_pred == y_true) / length(y_true)*100;
end

% show the linear plot for each qda (in bigger the best one)
[val_qda_mantained, idx_qda_mantained] = max(acc_mantained);
disp(['MANTAINED | Max acc: ' num2str(val_qda_mantained) ' | QDA with ' num2str(idx_qda_mantained) ' features']);
figure();
hold on
for idx_qda=1:nqda
    if idx_qda == idx_qda_mantained
        plot(qda_mantained_mean_acc(idx_qda,:), 'LineWidth',5);
    else
        plot(qda_mantained_mean_acc(idx_qda, :), 'LineWidth', 1);
    end
end
hold off;
xticks(1:256:min_durCF)
xticklabels(((1:256:min_durCF) - 1)  / sampleRate)
xlim([0 min_durCF]); 
ylim([0.3 1.0]); 
title('All QDA mantained | raw prob mean trials | TRAIN')

%% integrate the two best qda performances
alpha = 0.9;
time_merge = 1.0; % time in sec at which perform the merge
start_itegrator = 0.5;
qda_best_shift_integrated_soft = ones(ntrain_trial, min_durCF + 1) * start_itegrator;
qda_best_mantained_integrated_soft = ones(ntrain_trial, min_durCF + 1) * start_itegrator;
qda_merge_integrated_soft = ones(ntrain_trial, min_durCF + 1) * start_itegrator;
qda_best_shift_integrated_hard = ones(ntrain_trial, min_durCF + 1) * start_itegrator;
qda_best_mantained_integrated_hard = ones(ntrain_trial, min_durCF + 1) * start_itegrator;
qda_merge_integrated_hard = ones(ntrain_trial, min_durCF + 1) * start_itegrator;
for idx_train_trial=1:ntrain_trial
    idx_trial = train_trial(idx_train_trial);
    c_start = classified_info_shift{idx_qda_shift}.startTrial(idx_trial); % taking shift or maintain in equal
    c_end = c_start + min_durCF - 1;
    idx_class = find(data{1}.typ(idx_trial) == classes);

    c_data_shift = classified_data_shift{idx_qda_shift};
    c_shift_soft = c_data_shift(c_start:c_end,idx_class);
    c_shift_hard = ones(min_durCF, 1);
    c_shift_hard(c_shift_soft < 0.5) = 0;
    c_shift_hard(c_shift_soft == 0.5) = 0.5;

    c_data_mantained = classified_data_mantained{idx_qda_mantained};
    c_mantained_soft = c_data_mantained(c_start:c_end,idx_class);
    c_mantained_hard = ones(min_durCF, 1);
    c_mantained_hard(c_mantained_soft < 0.5) = 0;
    c_mantained_hard(c_mantained_soft == 0.5) = 0.5;

    % compute the integration
    for idx_sample=1:min_durCF
        % compute with soft predict
        qda_best_shift_integrated_soft(idx_train_trial, idx_sample + 1) = ...
           qda_best_shift_integrated_soft(idx_train_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_soft(idx_sample);
        qda_best_mantained_integrated_soft(idx_train_trial, idx_sample + 1) = ...
           qda_best_mantained_integrated_soft(idx_train_trial, idx_sample)* alpha + (1.0-alpha)*c_mantained_soft(idx_sample);

        if idx_sample < time_merge*sampleRate
            qda_merge_integrated_soft(idx_train_trial, idx_sample + 1) = ...
               qda_merge_integrated_soft(idx_train_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_soft(idx_sample);
        else
            qda_merge_integrated_soft(idx_train_trial, idx_sample + 1) = ...
               qda_merge_integrated_soft(idx_train_trial, idx_sample)*alpha + (1.0-alpha)*c_mantained_soft(idx_sample);
        end

        % compute for hard predict
        qda_best_shift_integrated_hard(idx_train_trial, idx_sample + 1) = ...
           qda_best_shift_integrated_hard(idx_train_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_hard(idx_sample);
        qda_best_mantained_integrated_hard(idx_train_trial, idx_sample + 1) = ...
           qda_best_mantained_integrated_hard(idx_train_trial, idx_sample)* alpha + (1.0-alpha)*c_mantained_hard(idx_sample);

        if idx_sample < time_merge*sampleRate
            qda_merge_integrated_hard(idx_train_trial, idx_sample + 1) = ...
               qda_merge_integrated_hard(idx_train_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_hard(idx_sample);
        else
            qda_merge_integrated_hard(idx_train_trial, idx_sample + 1) = ...
               qda_merge_integrated_hard(idx_train_trial, idx_sample)*alpha + (1.0-alpha)*c_mantained_hard(idx_sample);
        end
    end
end

figure();
subplot(1,2,1)
hold on
plot(mean(qda_best_shift_integrated_soft, 1), 'LineWidth', 5);
plot(mean(qda_best_mantained_integrated_soft, 1), 'LineWidth', 5);
plot(mean(qda_merge_integrated_soft, 1), 'LineWidth', 5);
hold off
legend({'QDA shift', 'QDA manteined', 'QDA merge'})
title('Using soft probabilities')
xlim([0 min_durCF])
ylim([0.3 1.0])
xticks(1:256:min_durCF)
xticklabels(((1:256:min_durCF) - 1)  / sampleRate)

subplot(1,2,2)
hold on
plot(mean(qda_best_shift_integrated_hard, 1), 'LineWidth', 5);
plot(mean(qda_best_mantained_integrated_hard, 1), 'LineWidth', 5);
plot(mean(qda_merge_integrated_hard, 1), 'LineWidth', 5);
hold off
legend({'QDA shift', 'QDA manteined', 'QDA merge'})
title( 'Using hard probabilities')
xlim([0 min_durCF])
ylim([0.3 1.0])
xticks(1:256:min_durCF)
xticklabels(((1:256:min_durCF) - 1)  / sampleRate)
sgtitle([subject ' | integrated probabilities mean among trials | TRAIN'])

%% plot the best QDA for the two period
% mean and std of the raw probabilities
figure();
subplot(1,3,1)
hold on;
plot(qda_shift_mean_acc(idx_qda_shift,:), 'Color','b');
plot(qda_shift_mean_acc(idx_qda_shift,:) + 0.3*qda_shift_std_acc(idx_qda_shift,:), '--', 'Color','b');
plot(qda_shift_mean_acc(idx_qda_shift,:) - 0.3*qda_shift_std_acc(idx_qda_shift,:), '--', 'Color','b');
plot(qda_mantained_mean_acc(idx_qda_shift,:), 'Color','r');
plot(qda_mantained_mean_acc(idx_qda_shift,:) + 0.3*qda_mantained_std_acc(idx_qda_shift,:), '--', 'Color','r');
plot(qda_mantained_mean_acc(idx_qda_shift,:) - 0.3*qda_mantained_std_acc(idx_qda_shift,:), '--', 'Color','r');
hold off;
legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', 'QDA mantained mean', ...
    'QDA mantained + std', 'QDA mantained - std'});
xlim([0 min_durCF])
ylim([0.0 1.0])
xticks(1:256:min_durCF)
xticklabels(((1:256:min_durCF) - 1)  / sampleRate)
title('raws probabilities')

% mean and std of integrated probabilities + the merged integration
subplot(1,3,2)
hold on;
plot(mean(qda_best_shift_integrated_soft, 1), 'Color','b');
plot(mean(qda_best_shift_integrated_soft, 1) + 0.3*std(qda_best_shift_integrated_soft,0,1), '--', 'Color','b');
plot(mean(qda_best_shift_integrated_soft, 1) - 0.3*std(qda_best_shift_integrated_soft,0,1), '--', 'Color','b');

plot(mean(qda_best_mantained_integrated_soft, 1), 'Color','r');
plot(mean(qda_best_mantained_integrated_soft, 1) + 0.3*std(qda_best_mantained_integrated_soft,0,1), '--', 'Color','r');
plot(mean(qda_best_mantained_integrated_soft, 1) - 0.3*std(qda_best_mantained_integrated_soft,0,1), '--', 'Color','r');

plot(mean(qda_merge_integrated_soft, 1), 'Color','k');
plot(mean(qda_merge_integrated_soft, 1) + 0.3*std(qda_merge_integrated_soft,0,1), '--', 'Color','k');
plot(mean(qda_merge_integrated_soft, 1) - 0.3*std(qda_merge_integrated_soft,0,1), '--', 'Color','k');
hold off;
legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', ...
    'QDA mantained mean', 'QDA mantained + std', 'QDA mantained - std', ...
    'QDA merge mean', 'QDA merge + std', 'QDA merge - std'});
xlim([0 min_durCF])
ylim([0.0 1.0])
xticks(1:256:min_durCF)
xticklabels(((1:256:min_durCF) - 1)  / sampleRate)
title('integrated with soft probabilities')

% mean and std for the hard (before was only soft)
subplot(1,3,3)
hold on;
plot(mean(qda_best_shift_integrated_hard, 1), 'Color','b');
% plot(mean(qda_best_shift_integrated_hard, 1) + 0.3*std(qda_best_shift_integrated_hard,0,1), '--', 'Color','b');
% plot(mean(qda_best_shift_integrated_hard, 1) - 0.3*std(qda_best_shift_integrated_hard,0,1), '--', 'Color','b');

plot(mean(qda_best_mantained_integrated_hard, 1), 'Color','r');
% plot(mean(qda_best_mantained_integrated_hard, 1) + 0.3*std(qda_best_mantained_integrated_hard,0,1), '--', 'Color','r');
% plot(mean(qda_best_mantained_integrated_hard, 1) - 0.3*std(qda_best_mantained_integrated_hard,0,1), '--', 'Color','r');

plot(mean(qda_merge_integrated_hard, 1), 'Color','k');
% plot(mean(qda_merge_integrated_hard, 1) + 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','k');
% plot(mean(qda_merge_integrated_hard, 1) - 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','k');
hold off;
% legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', ...
%     'QDA mantained mean', 'QDA mantained + std', 'QDA mantained - std', ...
%     'QDA merge mean', 'QDA merge + std', 'QDA merge - std'});
legend({'QDA shift mean', ...
    'QDA mantained mean', ...
    'QDA merge mean'});
xlim([0 min_durCF])
ylim([0.0 1.0])
xticks(1:256:min_durCF)
xticklabels(((1:256:min_durCF) - 1)  / sampleRate)
title('integrated with hard probabilities')

sgtitle([subject ' | raw probabilities and integration | TRAIN'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --> TEST set <--
% test the best classifier in the test set
time = 1; %s
samples = time * sampleRate;
min_durCF = size(data{1}.cf, 1);
acc_shift = nan(1, nqda);
start_trial_id = 41; % trial at which start the test set
end_trial_id = ntrial;
test_trial = start_trial_id:end_trial_id;
ntest_trial = length(test_trial);

% shift
c_data = classified_data_shift{idx_qda_shift};
shift_raw_best = nan(ntest_trial, min_durCF);
for idx_test_trial = 1:ntest_trial
    idx_trial = test_trial(idx_test_trial);
    c_start = classified_info_shift{idx_qda_shift}.startTrial(idx_trial);
    c_end = c_start + min_durCF - 1;
 
    qda_performance = c_data(c_start:c_end,1);

    idx_class = find(data{1}.typ(idx_trial) == classes);
    shift_raw_best(idx_test_trial,:) = c_data(c_start:c_end, idx_class);
end
qda_shift_mean_acc = mean(shift_raw_best, 1);
qda_shift_std_acc = std(shift_raw_best, 0, 1);

% mantained
mantained_raw_best = nan(ntest_trial, min_durCF);
c_data = classified_data_mantained{idx_qda_mantained};
for idx_test_trial = 1:ntest_trial
    idx_trial = test_trial(idx_test_trial);
    c_start = classified_info_mantained{idx_qda_mantained}.startTrial(idx_trial);
    c_end = c_start + min_durCF - 1;
 
    qda_performance = c_data(c_start:c_end,1);

    idx_class = find(data{1}.typ(idx_trial) == classes);
    mantained_raw_best(idx_test_trial,:) = c_data(c_start:c_end, idx_class);
end
qda_mantained_mean_acc = mean(mantained_raw_best, 1);
qda_mantained_std_acc = std(mantained_raw_best, 0, 1);

% integration
qda_best_shift_integrated_soft = ones(ntest_trial, min_durCF + 1) * start_itegrator;
qda_best_mantained_integrated_soft = ones(ntest_trial, min_durCF + 1) * start_itegrator;
qda_merge_integrated_soft = ones(ntest_trial, min_durCF + 1) * start_itegrator;
qda_best_shift_integrated_hard = ones(ntest_trial, min_durCF + 1) * start_itegrator;
qda_best_mantained_integrated_hard = ones(ntest_trial, min_durCF + 1) * start_itegrator;
qda_merge_integrated_hard = ones(ntest_trial, min_durCF + 1) * start_itegrator;
for idx_test_trial=1:ntest_trial
    idx_trial = test_trial(idx_test_trial);
    c_start = classified_info_shift{idx_qda_shift}.startTrial(idx_trial); % taking shift or maintain in equal
    c_end = c_start + min_durCF - 1;
    idx_class = find(data{1}.typ(idx_trial) == classes);

    c_data_shift = classified_data_shift{idx_qda_shift};
    c_shift_soft = c_data_shift(c_start:c_end,idx_class);
    c_shift_hard = ones(min_durCF, 1);
    c_shift_hard(c_shift_soft < 0.5) = 0;
    c_shift_hard(c_shift_soft == 0.5) = 0.5;

    c_data_mantained = classified_data_mantained{idx_qda_mantained};
    c_mantained_soft = c_data_mantained(c_start:c_end,idx_class);
    c_mantained_hard = ones(min_durCF, 1);
    c_mantained_hard(c_mantained_soft < 0.5) = 0;
    c_mantained_hard(c_mantained_soft == 0.5) = 0.5;

    % compute the integration
    for idx_sample=1:min_durCF
        % compute with soft predict
        qda_best_shift_integrated_soft(idx_test_trial, idx_sample + 1) = ...
           qda_best_shift_integrated_soft(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_soft(idx_sample);
        qda_best_mantained_integrated_soft(idx_test_trial, idx_sample + 1) = ...
           qda_best_mantained_integrated_soft(idx_test_trial, idx_sample)* alpha + (1.0-alpha)*c_mantained_soft(idx_sample);

        if idx_sample < time_merge*sampleRate
            qda_merge_integrated_soft(idx_test_trial, idx_sample + 1) = ...
               qda_merge_integrated_soft(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_soft(idx_sample);
        else
            qda_merge_integrated_soft(idx_test_trial, idx_sample + 1) = ...
               qda_merge_integrated_soft(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_mantained_soft(idx_sample);
        end

        % compute for hard predict
        qda_best_shift_integrated_hard(idx_test_trial, idx_sample + 1) = ...
           qda_best_shift_integrated_hard(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_hard(idx_sample);
        qda_best_mantained_integrated_hard(idx_test_trial, idx_sample + 1) = ...
           qda_best_mantained_integrated_hard(idx_test_trial, idx_sample)* alpha + (1.0-alpha)*c_mantained_hard(idx_sample);

        if idx_sample < time_merge*sampleRate
            qda_merge_integrated_hard(idx_test_trial, idx_sample + 1) = ...
               qda_merge_integrated_hard(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_shift_hard(idx_sample);
        else
            qda_merge_integrated_hard(idx_test_trial, idx_sample + 1) = ...
               qda_merge_integrated_hard(idx_test_trial, idx_sample)*alpha + (1.0-alpha)*c_mantained_hard(idx_sample);
        end
    end
end

% plot
figure();
subplot(1,3,1)
hold on;
plot(qda_shift_mean_acc, 'Color','b');
plot(qda_shift_mean_acc + 0.3*qda_shift_std_acc, '--', 'Color','b');
plot(qda_shift_mean_acc - 0.3*qda_shift_std_acc, '--', 'Color','b');
plot(qda_mantained_mean_acc, 'Color','r');
plot(qda_mantained_mean_acc + 0.3*qda_mantained_std_acc, '--', 'Color','r');
plot(qda_mantained_mean_acc - 0.3*qda_mantained_std_acc, '--', 'Color','r');
hold off;
legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', 'QDA mantained mean', ...
    'QDA mantained + std', 'QDA mantained - std'});
xlim([0 min_durCF])
ylim([0.0 1.0])
xticks(1:256:min_durCF)
xticklabels(((1:256:min_durCF) - 1)  / sampleRate)
title('raws probabilities')

% mean and std of integrated probabilities + the merged integration
subplot(1,3,2)
hold on;
plot(mean(qda_best_shift_integrated_soft, 1), 'Color','b');
plot(mean(qda_best_shift_integrated_soft, 1) + 0.3*std(qda_best_shift_integrated_soft,0,1), '--', 'Color','b');
plot(mean(qda_best_shift_integrated_soft, 1) - 0.3*std(qda_best_shift_integrated_soft,0,1), '--', 'Color','b');

plot(mean(qda_best_mantained_integrated_soft, 1), 'Color','r');
plot(mean(qda_best_mantained_integrated_soft, 1) + 0.3*std(qda_best_mantained_integrated_soft,0,1), '--', 'Color','r');
plot(mean(qda_best_mantained_integrated_soft, 1) - 0.3*std(qda_best_mantained_integrated_soft,0,1), '--', 'Color','r');

plot(mean(qda_merge_integrated_soft, 1), 'Color','k');
plot(mean(qda_merge_integrated_soft, 1) + 0.3*std(qda_merge_integrated_soft,0,1), '--', 'Color','k');
plot(mean(qda_merge_integrated_soft, 1) - 0.3*std(qda_merge_integrated_soft,0,1), '--', 'Color','k');
hold off;
legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', ...
    'QDA mantained mean', 'QDA mantained + std', 'QDA mantained - std', ...
    'QDA merge mean', 'QDA merge + std', 'QDA merge - std'});
xlim([0 min_durCF])
ylim([0.0 1.0])
xticks(1:256:min_durCF)
xticklabels(((1:256:min_durCF) - 1)  / sampleRate)
title('integrated with soft probabilities')

% mean and std for the hard (before was only soft)
subplot(1,3,3)
hold on;
plot(mean(qda_best_shift_integrated_hard, 1), 'Color','b');
% plot(mean(qda_best_shift_integrated_hard, 1) + 0.3*std(qda_best_shift_integrated_hard,0,1), '--', 'Color','b');
% plot(mean(qda_best_shift_integrated_hard, 1) - 0.3*std(qda_best_shift_integrated_hard,0,1), '--', 'Color','b');

plot(mean(qda_best_mantained_integrated_hard, 1), 'Color','r');
% plot(mean(qda_best_mantained_integrated_hard, 1) + 0.3*std(qda_best_mantained_integrated_hard,0,1), '--', 'Color','r');
% plot(mean(qda_best_mantained_integrated_hard, 1) - 0.3*std(qda_best_mantained_integrated_hard,0,1), '--', 'Color','r');

plot(mean(qda_merge_integrated_hard, 1), 'Color','k');
% plot(mean(qda_merge_integrated_hard, 1) + 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','k');
% plot(mean(qda_merge_integrated_hard, 1) - 0.3*std(qda_merge_integrated_hard,0,1), '--', 'Color','k');
hold off;
% legend({'QDA shift mean', 'QDA shift + std', 'Qda shift - std', ...
%     'QDA mantained mean', 'QDA mantained + std', 'QDA mantained - std', ...
%     'QDA merge mean', 'QDA merge + std', 'QDA merge - std'});
legend({'QDA shift mean', ...
    'QDA mantained mean', ...
    'QDA merge mean'});
xlim([0 min_durCF])
ylim([0.0 1.0])
xticks(1:256:min_durCF)
xticklabels(((1:256:min_durCF) - 1)  / sampleRate)
title('integrated with hard probabilities')

sgtitle([subject ' | raw probabilities and integration | TEST'])