%% check for the ERP
clear all; % close all;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/EOG')

%% Initialization
band = [2 10];
headers.TYP = [];
headers.DUR = [];
headers.POS = [];
signals = [];
classes = [730 731];
cf_event = 781;
fix_event = 786;
nchannels = 39;
nclasses = length(classes);
filterOrder = 4;
avg = 1;
eog_threshold = 500;
load('/home/paolo/lap_39.mat')
roi = {{'F3', 'F1', 'Fz', 'F2', 'F4'}; {'FC3', 'FC1', 'FCz', 'FC2', 'FC4'}; {'C3', 'C1', 'Cz', 'C2', 'C4'}; {'CP3', 'CP1', 'CPz', 'CP2', 'CP4'};
    {'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6'}; {'PO7', 'PO5', 'PO3', 'POz', 'PO4', 'PO6', 'PO8'}; {'O1', 'Oz', 'O2'}};
roi_label = {{'F channels'}, {'FC channels'}, {'C channels'}, {'CP channels'}, {'P channels'}, {'PO channels'}, {'O channels'}};
roi_indices = cell(size(roi));
nroi = length(roi);

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
    for i = 1:nroi
        roi_indices{i} = find(ismember(channels_label, roi{i}));
    end

%     c_signal = c_signal * lap;
    c_trial_with_eog = eog_detection(c_signal, header, eog_threshold, {'FP1', 'FP2', 'EOG'});
    trial_with_eog = [trial_with_eog; c_trial_with_eog];

    signal_roi = nan(size(c_signal, 1), nroi); % for roi
    for i = 1:nroi
        signal_roi(:,i) = mean(c_signal(:, roi_indices{i}), 2);
    end

    % filter alpha band
    disp('      [proc] applying filtering')
    [b, a] = butter(filterOrder, band(2)*(2/header.SampleRate),'low');
    s_low = filter(b,a,signal_roi);
    [b, a] = butter(filterOrder, band(1)*(2/header.SampleRate),'high');
    s_filt = filter(b,a,s_low);

    signal_roi_processed = s_filt;

    headers.sampleRate = header.SampleRate;
    headers.channels_labels = header.Label;
    if isempty(find(header.EVENT.TYP == 2, 1)) % no eye calibration
        headers.TYP = cat(1, headers.TYP, header.EVENT.TYP);
        headers.DUR = cat(1, headers.DUR, header.EVENT.DUR);
        headers.POS = cat(1, headers.POS, header.EVENT.POS + size(signals, 1));
    else
        k = find(header.EVENT.TYP == 1, 1);
        headers.TYP = cat(1, headers.TYP, header.EVENT.TYP(k:end));
        headers.DUR = cat(1, headers.DUR, header.EVENT.DUR(k:end));
        headers.POS = cat(1, headers.POS, header.EVENT.POS(k:end) + size(signals, 1));
    end

    signals = cat(1, signals, signal_roi_processed(:,:));
end

%% Labelling data 
sampleRate = headers.sampleRate;
cuePOS = headers.POS(ismember(headers.TYP, classes));
cueDUR = headers.DUR(ismember(headers.TYP, classes));
cueTYP = headers.TYP(ismember(headers.TYP, classes));

fixPOS = headers.POS(headers.TYP == 786);
fixDUR = headers.DUR(headers.TYP == 786);

cfPOS = headers.POS(headers.TYP == 781);
cfDUR = headers.DUR(headers.TYP == 781);

minDurCue = min(cueDUR);
ntrial = length(cuePOS);

%% Labeling data for the dataset
min_durFIX = min(fixDUR);
min_durCF = min(cfDUR);
min_durCUE = min(cueDUR);

trial_start = nan(ntrial, 1);
trial_end = nan(ntrial, 1);
trial_typ = nan(ntrial, 1);
for idx_trial = 1:ntrial
    trial_start(idx_trial) = fixPOS(idx_trial);
    trial_typ(idx_trial) = cueTYP(idx_trial);
    trial_end(idx_trial) = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;
end

min_trial_data = min(trial_end - trial_start+1);
trial_data = nan(min_trial_data, nroi, ntrial);

for trial = 1:ntrial
    c_start = trial_start(trial);
    c_end = trial_start(trial) + min_trial_data - 1;
    trial_data(:,:,trial) = signals(c_start:c_end,:);
end


%% refactoring the data
% now the data are placed alternatevely, so odd trial id is for 730, even for 731
balanced_trial_idx = balanced_data_afterEOG(trial_with_eog, trial_typ, classes);
balanced_trial_data = trial_data(:,:,logical(balanced_trial_idx));
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
        tmp(:,:,idx_trial_class + idx_class - 1) = balanced_trial_data(:,:,idx_classes_trial(i, idx_class));
        trial_typ(idx_trial_class + idx_class - 1) = classes(idx_class);
    end
    i = i + 1;
end
trial_data = tmp;

% trial_data(:,:,:,3) = [];

%% see the errp
for idx_roi = 1:nroi
    c_data = squeeze(trial_data(:, idx_roi, :));

    figure();
    subplot(1,2,1)
    plot(c_data, 'Color', 'k');
    hold on;
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    xline(min_durFIX + sampleRate*0.7, '--r', 'stop', 'LabelOrientation', 'horizontal');
    plot(mean(c_data, 2), 'Color', 'red');
    hold off
    xticks(sampleRate:sampleRate:size(c_data, 1))
    xticklabels(string((sampleRate:sampleRate:size(c_data, 1)) / sampleRate));
    title('mean channel across trials')

    subplot(1,2,2)
    imagesc(c_data')
    colormap("jet")
    xticks(sampleRate:sampleRate:size(c_data, 1))
    xticklabels(string((sampleRate:sampleRate:size(c_data, 1)) / sampleRate));
    title('channel signal for all trials')
    ylabel('trial')
    xlabel('time')

    sgtitle(['channel: ' roi_label{idx_roi}{1}])
end
