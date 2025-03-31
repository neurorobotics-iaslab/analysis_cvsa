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
    fullpath_file = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file);
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
[qda_filenames, pathname] = uigetfile('*.yaml', 'Select QDA Files', 'MultiSelect', 'on');
if ischar(qda_filenames)
    qda_filenames = {qda_filenames};
end
%%
subject = pathname(28:29);
numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), qda_filenames);
[~, sortedIdx] = sort(numbers);
qda_filenames = qda_filenames(sortedIdx);
nqda = size(qda_filenames, 2);
classified_data = cell(1, nqda);
ntrial = size(data{1}.cf,3);
for i = 1:nqda
    classified_data{i} = nan(size(data{1}.cf, 1)*ntrial, 2);
end
classified_info = cell(1, nqda);
acc_total = zeros(nqda, 1);

for idx_qda = 1:nqda
    fullpath_file = fullfile(pathname, qda_filenames{idx_qda});
    disp(['QDA: ' qda_filenames{idx_qda}]);
    qda = loadQDA(fullpath_file);

    if qda.nfeatures == 1
        c_idx_band = find(cellfun(@(x) isequal(x, qda.bands), bands));
        selectedFeatures = [qda.idchans, c_idx_band];
    else
        c_idx_band = [];
        for i = 1:length(qda.bands) 
            c_idx_band = [c_idx_band find(ismember(cell2mat(bands'), qda.bands(i,:), 'rows'))];
        end
        selectedFeatures = [qda.idchans', c_idx_band'];
    end
    [X, y, classified_info{idx_qda}] = createDataset(filenames, data, selectedFeatures, bands, channels_label);

    classified_tmp = apply_qda_matrix(qda, X);
    classified_data{idx_qda}(:,:) = classified_tmp;

    y_pred = ones(size(classified_tmp, 1), 1) * classes(1);
    y_pred(classified_tmp(:,1) < 0.5) = classes(2);
    acc_total(idx_qda) = sum(y_pred == y) / length(y);
    disp(['  accuracy all runs: ' num2str(acc_total(idx_qda))])
end

%% plot the last 20 trials
n_last = 20;
nrows = 2;
ntrial = size(classified_info{1}.startTrial,1);
handles = [];
figure();
for idx_trial = ntrial - n_last + 1 : ntrial
    c_start = classified_info{1}.startTrial(idx_trial);
    if idx_trial == ntrial
        c_end = size(classified_data{1}, 1)-1;
    else 
        c_end = classified_info{1}.startTrial(idx_trial + 1)-1;
    end

    qda_performance = nan(nqda, c_end-c_start + 1);
    for idx_qda=1:nqda
        c_data = classified_data{idx_qda};
        qda_performance(idx_qda,:) = c_data(c_start:c_end,1);
    end
    x_label = 0:256:size(qda_performance, 2);

    subplot(nrows, ceil(n_last/nrows), idx_trial - (ntrial - n_last))
    imagesc(qda_performance)
    yticks(1:nqda)
    yticklabels(1:nqda)
    xticks(x_label)
    xticklabels(x_label / sampleRate)
    handles = [handles, gca];
    title(['cue: ' num2str(data{1}.typ(idx_trial))])
    drawnow;
end
set(handles, 'clim', [0 1])
sgtitle([ subject ' | Only the cf | 1 means 730, 0 means 731'])

%% plot each qda and each trial ------> ALL TRIAL DURATION
start_trial_id = 41;
end_trial_id = ntrial;
nrows = 2;
nfigure = 4;
idx_figure = 1;
ntrial = size(classified_info{1}.startTrial,1);
handles = [];
a = classified_info{1}.startTrial;
b = a(2:end);
b = [b; size(classified_data{1}+1, 1)]; b = b + 1;
min_durCF = min(b-a);
idx_max_acc = nan;
max_acc = -inf;
figure();
for idx_qda = 1:nqda
    needed_trial = start_trial_id:end_trial_id;
    n_needed_trial = length(needed_trial);
    qda_performance = nan(n_needed_trial, min_durCF);
    y_label = cell(n_needed_trial, 1);
    c_data = classified_data{idx_qda};
    y_pred_selTrial = nan(n_needed_trial*min_durCF, 1);
    y_true_selTrial = nan(n_needed_trial*min_durCF, 1);
    for idx_needed_trial = 1:n_needed_trial
        idx_trial = needed_trial(idx_needed_trial);
        c_start = classified_info{1}.startTrial(idx_trial);
        if idx_trial == ntrial
            c_end = size(classified_data{1}, 1);
        else
            c_end = classified_info{1}.startTrial(idx_trial + 1)-1;
        end
        qda_performance(idx_needed_trial,:) = c_data(c_start:c_end,1);

        y_pred = ones(min_durCF, 1) * classes(1);
        y_pred(qda_performance(idx_needed_trial,:) < 0.5) = classes(2);
        y = repmat(data{1}.typ(idx_trial), min_durCF, 1);
        y_label{idx_needed_trial} = [num2str(data{1}.typ(idx_trial)) ': ' sprintf('%.2f', sum(y_pred == y) / min_durCF*100)];

        y_pred_selTrial((idx_needed_trial-1)*min_durCF+1:idx_needed_trial*min_durCF) = y_pred;
        y_true_selTrial((idx_needed_trial-1)*min_durCF+1:idx_needed_trial*min_durCF) = y;
    end

    acc_selTrial = sum(y_pred_selTrial == y_true_selTrial) / length(y_true_selTrial)*100;
    acc_selTrial_src = sprintf('%.2f', acc_selTrial);
    if max_acc < acc_selTrial
        idx_max_acc = idx_qda;
        max_acc = acc_selTrial;
    end

    if idx_qda == round(nqda/nfigure)*idx_figure + 1
        set(handles, 'clim', [0 1])
        sgtitle([ subject ' | Only the cf | 1 means 730, 0 means 731'])
        handles = [];
        idx_figure = idx_figure + 1;
        figure();
    end

    correctness = (y_pred_selTrial == y_true_selTrial);
    correctness_matrix = reshape(correctness, [min_durCF, n_needed_trial])';
    green_component = correctness_matrix .* qda_performance; % if prediction correct
    red_component = (~correctness_matrix) .* qda_performance; % 
    color_map = cat(3, red_component, green_component, zeros(size(green_component))); % No blue component
    subplot(nrows, ceil(nqda/nrows/nfigure), mod(idx_qda-1, round(nqda/nfigure)) + 1)
    imagesc(color_map)
    yticks(1:ntrial)
    yticklabels(y_label);
    xticks(x_label)
    xticklabels(x_label / sampleRate)
    handles = [handles, gca];
    title(['qda ' num2str(idx_qda) ' | acc: ' acc_selTrial_src])
    drawnow;
end
set(handles, 'clim', [0 1])
sgtitle([ subject ' | Only the cf | 1 means 730, 0 means 731'])

%%
idx_qda = idx_max_acc;
c_data = classified_data{idx_qda};
qda_performance = nan(n_needed_trial, min_durCF);
for idx_needed_trial = 1:n_needed_trial
    idx_trial = needed_trial(idx_needed_trial);
    c_start = classified_info{1}.startTrial(idx_trial);
    if idx_trial == ntrial
        c_end = size(classified_data{1}, 1);
    else
        c_end = classified_info{1}.startTrial(idx_trial + 1) -1;
    end
    qda_performance(idx_needed_trial,:) = c_data(c_start:c_end,1);

    y_pred = ones(min_durCF, 1) * classes(1);
    y_pred(qda_performance(idx_needed_trial,:) < 0.5) = classes(2);
    y = repmat(data{1}.typ(idx_trial), min_durCF, 1);
    y_label{idx_needed_trial} = [num2str(data{1}.typ(idx_trial)) ': ' sprintf('%.2f', sum(y_pred == y) / min_durCF*100)];

    y_pred_selTrial((idx_needed_trial-1)*min_durCF+1:idx_needed_trial*min_durCF) = y_pred;
    y_true_selTrial((idx_needed_trial-1)*min_durCF+1:idx_needed_trial*min_durCF) = y;
end
acc_selTrial = sprintf('%.2f', sum(y_pred_selTrial == y_true_selTrial) / length(y_true_selTrial)*100);

figure();
correctness = (y_pred_selTrial == y_true_selTrial);
correctness_matrix = reshape(correctness, [min_durCF, n_needed_trial])';
green_component = correctness_matrix .* qda_performance; % if prediction correct
red_component = (~correctness_matrix) .* qda_performance; %
color_map = cat(3, red_component, green_component, zeros(size(green_component))); % No blue component
imagesc(color_map)
yticks(1:ntrial)
yticklabels(y_label);
xticks(x_label)
xticklabels(x_label / sampleRate)
title(['qda ' num2str(idx_qda) ' | acc: ' acc_selTrial])

%% plot each qda and each trial  -----> CHOSEN TRIAL DURATION
start_trial_id = 41;
end_trial_id = ntrial;
nrows = 2;
nfigure = 4;
idx_figure = 1;
ntrial = size(classified_info{1}.startTrial,1);
handles = [];
a = classified_info{1}.startTrial;
b = a(2:end);
b = [b; size(classified_data{1}+1, 1)]; b = b + 1;
sample_dur = 512;
idx_max_acc = nan;
max_acc = -inf;
figure();
for idx_qda = 1:nqda
    needed_trial = start_trial_id:end_trial_id;
    n_needed_trial = length(needed_trial);
    qda_performance = nan(n_needed_trial, sample_dur);
    y_label = cell(n_needed_trial, 1);
    c_data = classified_data{idx_qda};
    y_pred_selTrial = nan(n_needed_trial*sample_dur, 1);
    y_true_selTrial = nan(n_needed_trial*sample_dur, 1);
    for idx_needed_trial = 1:n_needed_trial
        idx_trial = needed_trial(idx_needed_trial);
        c_start = classified_info{1}.startTrial(idx_trial);

        c_end = c_start + sample_dur - 1;

        qda_performance(idx_needed_trial,:) = c_data(c_start:c_end,1);

        y_pred = ones(sample_dur, 1) * classes(1);
        y_pred(qda_performance(idx_needed_trial,:) < 0.5) = classes(2);
        y = repmat(data{1}.typ(idx_trial), sample_dur, 1);
        y_label{idx_needed_trial} = [num2str(data{1}.typ(idx_trial)) ': ' sprintf('%.2f', sum(y_pred == y) / sample_dur*100)];

        y_pred_selTrial((idx_needed_trial-1)*sample_dur+1:idx_needed_trial*sample_dur) = y_pred;
        y_true_selTrial((idx_needed_trial-1)*sample_dur+1:idx_needed_trial*sample_dur) = y;
    end

    acc_selTrial = sum(y_pred_selTrial == y_true_selTrial) / length(y_true_selTrial)*100;
    acc_selTrial_src = sprintf('%.2f', acc_selTrial);
    if max_acc < acc_selTrial
        idx_max_acc = idx_qda;
        max_acc = acc_selTrial;
    end

    if idx_qda == round(nqda/nfigure)*idx_figure + 1
        set(handles, 'clim', [0 1])
        sgtitle([ subject ' | Only the cf | 1 means 730, 0 means 731'])
        handles = [];
        idx_figure = idx_figure + 1;
        figure();
    end

    correctness = (y_pred_selTrial == y_true_selTrial);
    correctness_matrix = reshape(correctness, [sample_dur, n_needed_trial])';
    green_component = correctness_matrix .* qda_performance; % if prediction correct
    red_component = (~correctness_matrix) .* qda_performance; % 
    color_map = cat(3, red_component, green_component, zeros(size(green_component))); % No blue component
    subplot(nrows, ceil(nqda/nrows/nfigure), mod(idx_qda-1, round(nqda/nfigure)) + 1)
    imagesc(color_map)
    yticks(1:ntrial)
    yticklabels(y_label);
    xticks(x_label)
    xticklabels(x_label / sampleRate)
    handles = [handles, gca];
    title(['qda ' num2str(idx_qda) ' | acc: ' acc_selTrial_src])
    drawnow;
end
set(handles, 'clim', [0 1])
sgtitle([ subject ' | Only the cf | 1 means 730, 0 means 731'])

%%
idx_qda = idx_max_acc;
c_data = classified_data{idx_qda};
qda_performance = nan(n_needed_trial, sample_dur);
for idx_needed_trial = 1:n_needed_trial
    idx_trial = needed_trial(idx_needed_trial);
    c_start = classified_info{1}.startTrial(idx_trial);
    c_end = c_start + sample_dur - 1;
    qda_performance(idx_needed_trial,:) = c_data(c_start:c_end,1);

    y_pred = ones(sample_dur, 1) * classes(1);
    y_pred(qda_performance(idx_needed_trial,:) < 0.5) = classes(2);
    y = repmat(data{1}.typ(idx_trial), sample_dur, 1);
    y_label{idx_needed_trial} = [num2str(data{1}.typ(idx_trial)) ': ' sprintf('%.2f', sum(y_pred == y) / sample_dur*100)];

    y_pred_selTrial((idx_needed_trial-1)*sample_dur+1:idx_needed_trial*sample_dur) = y_pred;
    y_true_selTrial((idx_needed_trial-1)*sample_dur+1:idx_needed_trial*sample_dur) = y;
end
acc_selTrial = sprintf('%.2f', sum(y_pred_selTrial == y_true_selTrial) / length(y_true_selTrial)*100);

figure();
correctness = (y_pred_selTrial == y_true_selTrial);
correctness_matrix = reshape(correctness, [sample_dur, n_needed_trial])';
green_component = correctness_matrix .* qda_performance; % if prediction correct
red_component = (~correctness_matrix) .* qda_performance; %
color_map = cat(3, red_component, green_component, zeros(size(green_component))); % No blue component
imagesc(color_map)
yticks(1:ntrial)
yticklabels(y_label);
xticks(x_label)
xticklabels(x_label / sampleRate)
title(['qda ' num2str(idx_qda) ' | acc: ' acc_selTrial])


