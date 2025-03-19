%% compute the fischer to understand features using the log band
clc; clearvars;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/utils')

%% Initialization
a = 4:2:18;
b = a+2;
c = [a; b];
bands = [];
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
subject = pathname(28:29);
numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), qda_filenames);
[~, sortedIdx] = sort(numbers);
qda_filenames = qda_filenames(sortedIdx);
nqda = size(qda_filenames, 2);
classified_data = cell(1, nqda);
classified_info = cell(1, nqda);

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
            c_idx_band = [c_idx_band find(cellfun(@(x) isequal(x, qda.bands(i,:)), bands))];
        end
        selectedFeatures = [qda.idchans', c_idx_band'];
    end
    disp(['n features: ' num2str(size(selectedFeatures, 1))])
    [X, y, classified_info{idx_qda}] = createDataset(filenames, data, selectedFeatures, bands, channels_label);

    classified_tmp = nan(size(X, 1), nclasses);
    for idx_sample = 1:size(X, 1)
        classified_tmp(idx_sample,:) = apply_qda(qda, X(idx_sample,:));
    end
    classified_data{idx_qda} = classified_tmp;
end
%% plot the last 20 trials
n_last = 20;
nrows = 2;
ntrial = size(classified_info{1}.startTrial,1);
cl = -inf;
handles = [];
figure();
for idx_trial = ntrial - n_last + 1 : ntrial
    c_start = classified_info{1}.startTrial(idx_trial) + 1;
    if idx_trial == ntrial
        c_end = size(classified_data{1}, 1);
    else 
        c_end = classified_info{1}.startTrial(idx_trial + 1) - 1;
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
    cl = max(cl, max(abs(qda_performance), [], 'all'));
    title(['cue: ' num2str(data{1}.typ(idx_trial))])
    drawnow;
end
set(handles, 'clim', [0 cl])
sgtitle([ subject ' | Only the cf | 1 means 730, 0 means 731'])

%% plot the mean for each class
a =  classified_info{1}.startTrial;
b = a(2:end); b = [b; size(X,1)];
minDurCf = min(b-a);
data_cf = nan(minDurCf, nqda, nclasses);
nrows = 4;
for idx_qda = 1:nqda
    c_data = classified_data{idx_qda};
    c_data_cf = nan(minDurCf, ntrial);
    for idx_trial = 1:ntrial
        c_start = classified_info{1}.startTrial(idx_trial) + 1;
        c_end = c_start + minDurCf - 1;
        c_data_cf(:,idx_trial) = c_data(c_start:c_end,1);
    end

    data_cf(:,idx_qda, 1) = mean(c_data_cf(:,data{1}.typ == classes(1)), 2);
    data_cf(:,idx_qda, 2) = mean(c_data_cf(:,data{1}.typ == classes(2)), 2);
end

handles = [];
cl = max(data_cf, [], 'all');
x_label = 0:256:size(data_cf, 1);
figure()
subplot(1,2,1)
imagesc(data_cf(:,:,1)')
yticks(1:nqda)
yticklabels(1:nqda)
xticks(x_label)
xticklabels(x_label / sampleRate)
title([num2str(classes(1)) ' | close to 1']);

handles = [handles, gca];
subplot(1,2,2)
imagesc(data_cf(:,:,2)')
yticks(1:nqda)
yticklabels(1:nqda)
xticks(x_label)
xticklabels(x_label / sampleRate)
title([num2str(classes(2)) ' | close to 0']);
handles = [handles, gca];
set(handles, 'clim', [0 cl]);
sgtitle([subject ' | Trial mean depending on class']);

