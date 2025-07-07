%% file to show the lateralized index and the FEF and PFC components according to the cue onset
clear all; % close all;

addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/EOG')

%% Initialization
a = 4:2:18;
b = a+2;
c = [a; b];
bands = cell(1, size(a, 2) + 1);
for i=1:length(a)
    bands{i} = [a(i), b(i)];
end
bands{i + 1} = [8 14];
bands = {[8 13]};
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
load('/home/paolo/lap_39.mat')

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

        % laplacian
        c_signal = c_signal * lap;

        % filter alpha band
        disp('      [proc] applying filtering')
        [b, a] = butter(filterOrder, band(2)*(2/header.SampleRate),'low');
        s_low = filter(b,a,c_signal);
        [b, a] = butter(filterOrder, band(1)*(2/header.SampleRate),'high');
        s_filt = filter(b,a,s_low);

        signal_processed = s_filt; % filtered signal

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



%% alpha-band amplitude lateralization
% define left and right channels
ch_left = {'FC1', 'F1'};
ch_right = {'FC2', 'F2'};
[~, channelsSelected_left] = ismember(ch_left, channels_label);
nchannelsSelected_left = size(channelsSelected_left, 2);
[~, channelsSelected_right] = ismember(ch_right, channels_label);
nchannelsSelected_right = size(channelsSelected_right, 2);
figure()
idx_plot = 1;

for idx_band = 1:nbands
    amp_left = squeeze(mean(trial_data(:,idx_band,channelsSelected_left,:), 3)); % samples x trials
    amp_right = squeeze(mean(trial_data(:,idx_band,channelsSelected_right,:), 3)); % samples x trials

%     baseline_win = 1:min(fixDUR); 
%     baseline_left = mean(amp_left(baseline_win, :), 1);
%     baseline_right = mean(amp_right(baseline_win, :), 1);
% 
%     amp_left = (amp_left - baseline_left) ./ baseline_left;
%     amp_right = (amp_right - baseline_right) ./ baseline_right;


    ALI_left = (amp_right(:,trial_typ == classes(1)) - amp_left(:,trial_typ == classes(1)) ) ./ ( amp_right(:,trial_typ == classes(1)) + amp_left(:,trial_typ == classes(1)) );
    ALI_left_mean = mean(ALI_left, 2);
    ALI_right = (amp_left(:,trial_typ == classes(2)) - amp_right(:,trial_typ == classes(2)) ) ./ ( amp_right(:,trial_typ == classes(2)) + amp_left(:,trial_typ == classes(2)) );
    ALI_right_mean = mean(ALI_right, 2);

    time_vector = (0:min_trial_data-1) / sampleRate;  % seconds, adjust if you have event-related timing

    subplot(nbands, 2, idx_plot)
    hold on;
    grid on;
    plot(time_vector, ALI_left_mean, 'b', 'LineWidth', 2);
    plot(time_vector, ALI_right_mean, 'r', 'LineWidth', 2);
    xlabel('Time (s)');
    legend({'Attend Left', 'Attend Right'});
    title(['band: ' bands_str{idx_band}]);
    hold off;
    idx_plot = idx_plot + 1;

    subplot(nbands, 2, idx_plot)
    hold on;
    grid on;
    plot(time_vector, ALI_left_mean - ALI_right_mean, 'b', 'LineWidth', 2);
    xlabel('Time (s)');
    title(['diff | band: ' bands_str{idx_band}]);
    hold off;
    idx_plot = idx_plot + 1;
end
ch_right_str = strjoin(ch_right, ', ');
ch_left_str = strjoin(ch_left, ', ');

sgtitle(['Alpha Lateralization Index | channels right: ' ch_right_str ' | channels left: ' ch_left_str]);

%% FEF
% define left and right channels
ch_left = {'FC1', 'F1'};
ch_right = {'FC2', 'F2'};
[~, channelsSelected_left] = ismember(ch_left, channels_label);
nchannelsSelected_left = size(channelsSelected_left, 2);
[~, channelsSelected_right] = ismember(ch_right, channels_label);
nchannelsSelected_right = size(channelsSelected_right, 2);
figure()
idx_plot = 1;

for idx_band = 1:nbands
    amp_left = squeeze(mean(trial_data(:,idx_band,channelsSelected_left,:), 3)); % samples x trials
    amp_right = squeeze(mean(trial_data(:,idx_band,channelsSelected_right,:), 3)); % samples x trials

%     baseline_win = 1:min(fixDUR); 
%     baseline_left = mean(amp_left(baseline_win, :), 1);
%     baseline_right = mean(amp_right(baseline_win, :), 1);
% 
%     amp_left = (amp_left - baseline_left) ./ baseline_left;
%     amp_right = (amp_right - baseline_right) ./ baseline_right;

    ALI_left_mean = mean(amp_left, 2);
    ALI_right_mean = mean(amp_right, 2);

    time_vector = (0:min_trial_data-1) / sampleRate; 

    subplot(nbands, 2, idx_plot)
    hold on;
    grid on;
    plot(time_vector, ALI_left_mean, 'b', 'LineWidth', 2);
    plot(time_vector, ALI_right_mean, 'r', 'LineWidth', 2);
    xlabel('Time (s)');
    legend({'Attend Left', 'Attend Right'});
    title(['band: ' bands_str{idx_band}]);
    hold off;
    idx_plot = idx_plot + 1;

    subplot(nbands, 2, idx_plot)
    hold on;
    grid on;
    plot(time_vector, ALI_left_mean - ALI_right_mean, 'b', 'LineWidth', 2);
    xlabel('Time (s)');
    title(['diff | band: ' bands_str{idx_band}]);
    hold off;
    idx_plot = idx_plot + 1;
end
ch_right_str = strjoin(ch_right, ', ');
ch_left_str = strjoin(ch_left, ', ');

sgtitle(['Alpha Lateralization Index | channels right: ' ch_right_str ' | channels left: ' ch_left_str]);