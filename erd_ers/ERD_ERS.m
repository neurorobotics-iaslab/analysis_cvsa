%% ERD/ERS ---> with baseline correction
% the file given the gdf compute the erd and ers for different bands and classes
% clear all; close all;

addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/EOG')

%% Initialization
a = 4:2:18;
b = a+2;
c = [a; b];
bands = cell(1, size(a, 2) + 1);
for i=1:length(a)
    bands{i} = [a(i), b(i)];
end
bands{i+1} = [8 14];
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

        signal_processed = proc_ERD(c_signal, header.SampleRate, band, filterOrder);
        
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

%% baseline correction ERD/ERS
baseline = trial_data(1:min(fixDUR),:,:,:);

nsamples = size(trial_data, 1);
correction = nan(nsamples, nbands, nchannelsSelected, ntrial);

figure();
idx_plot = 1;
handle = []; c_max = -inf;
baseline_correction = true;

for idx_band = 1:nbands
    c_baseline = mean(mean(baseline(:,idx_band, channelsSelected,:), 1), 4);
    if baseline_correction
        correction(:,idx_band,:,:) = (trial_data(:,idx_band, channelsSelected,:) - c_baseline) ./ c_baseline; % baseline correction
    else
        correction(:,idx_band,:,:) = trial_data(:,idx_band, channelsSelected,:);
    end

    avg_correction_1 = mean(squeeze(correction(:,idx_band,:,trial_typ == classes(1))), 3);
    avg_correction_2 = mean(squeeze(correction(:,idx_band,:,trial_typ == classes(2))), 3);

    subplot(nbands,nclasses + 1,idx_plot)
    imagesc(avg_correction_1(:,:)')
    hold on
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off
    yticks(1:size(channelsSelected,2))
    yticklabels(channels_select)
    xticks(sampleRate:sampleRate:size(avg_correction_1, 1))
    xticklabels(string((sampleRate:sampleRate:size(avg_correction_1, 1)) / sampleRate));
    title(['class: ' num2str(classes(1))])
    handle = [handle, gca];
    c_max = max(c_max, max(avg_correction_1, [], 'all'));
    idx_plot = idx_plot + 1;
    title(['band: ' bands_str{idx_band} ' | class: ' num2str(classes(1))])

    subplot(nbands,nclasses + 1,idx_plot)
    imagesc(avg_correction_2(:,:)')
    hold on
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off
    yticks(1:size(channelsSelected,2))
    yticklabels(channels_select)
    xticks(sampleRate:sampleRate:size(avg_correction_2, 1))
    xticklabels(string((sampleRate:sampleRate:size(avg_correction_2, 1)) / sampleRate));
    title(['class: ' num2str(classes(2))])
    handle = [handle, gca];
    c_max = max(c_max, max(avg_correction_2, [], 'all'));
    idx_plot = idx_plot + 1;

    title(['band: ' bands_str{idx_band} ' | class: ' num2str(classes(2))])
    set(handle, 'CLim', [-c_max, c_max]);

    c_max = -inf; handle = [];

    subplot(nbands,nclasses + 1,idx_plot)
    imagesc(avg_correction_1(:,:)' - avg_correction_2(:,:)')
    hold on
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off
    yticks(1:size(channelsSelected,2))
    yticklabels(channels_select)
    xticks(sampleRate:sampleRate:size(avg_correction_2, 1))
    xticklabels(string((sampleRate:sampleRate:size(avg_correction_2, 1)) / sampleRate));
    title(['Diff: ' num2str(classes(1)) '-' num2str(classes(2))])
    handle = [handle, gca];
    c_max = max(c_max, max(abs(avg_correction_1(:,:)' - avg_correction_2(:,:)'), [], 'all'));
    idx_plot = idx_plot + 1;
    set(handle, 'CLim', [-c_max, c_max]);
    handle = []; c_max = -inf;
end

if baseline_correction
    sgtitle('ERD/ERS with baseline correction' )
else
    sgtitle('ERD/ERS without baseline correction' )
end





