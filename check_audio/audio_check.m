clear all;

classes = [730 731];
nclasses = length(classes);
audio.data = cell(1, nclasses);
audio.fs = cell(1,nclasses);
[audio.data{1}, audio.fs{1}] = audioread("/home/paolo/cvsa_ws/src/feedback_cvsa/cfg/left.mp3");
[audio.data{2}, audio.fs{2}] = audioread("/home/paolo/cvsa_ws/src/feedback_cvsa/cfg/right.mp3");
bag = bag_read("/home/paolo/cvsa_ws/src/analysis_cvsa/check_audio/files/f2.20250709.183939.calibration.cvsa_lbrb.bag");

events = struct2table(bag_process_events(bag));
my_audio = bag.audio;

cue_pos = events(ismember(events.TYP, classes),:).POS;
cue_typ = events(ismember(events.TYP, classes),:).TYP;


%% check for the ERP
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
bands = {[2 10]};
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
nbands = length(bands);
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

        % filter alpha band
        disp('      [proc] applying filtering')
        [b, a] = butter(filterOrder, band(2)*(2/header.SampleRate),'low');
        s_low = filtfilt(b,a,c_signal);
        [b, a] = butter(filterOrder, band(1)*(2/header.SampleRate),'high');
        s_filt = filtfilt(b,a,s_low);

%         disp('      [proc] applying laplacian');
%         s_lap = s_filt * lap;
%         s_lap = s_filt - mean(s_filt, 2);
        s_lap = s_filt;

        signal_processed = s_lap; % filtered signal --> all channels

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

% trial_data(:,:,:,3) = [];

%% see the errp
channels_select = {'PO4', 'PO6', 'PO8', 'O1', 'O2'};
[~, channelsSelected] = ismember(channels_select, channels_label);
nchannelsSelected = size(channelsSelected, 2);

for i = 1:size(cue_pos,1)
    figure();
    
    if cue_typ(i) == classes(1)
        window_size = (length(audio.data{1})-1)*1/audio.fs{1}; %s
        f = my_audio.data(cue_pos(i):ceil(cue_pos(i)+my_audio.fs*window_size));
        t = linspace(0,window_size,length(f));

        subplot(5,1,1)
        plot(t, f)
        title('Audio registered');

        subplot(5,1,2)
        plot((0:(length(audio.data{1})-1))*1/audio.fs{1}, audio.data{1})
        title('original audio')

        subplot(5,1,5)
        melSpectrogram(double(audio.data{1}), audio.fs{1});
        title('Spect original')

        subplot(5,1,3)
        melSpectrogram(double(f), my_audio.fs)
        title('spectrogram')

        c_data = squeeze(mean(trial_data(1:ceil((length(audio.data{1})-1)/audio.fs{1} * sampleRate), 1, channelsSelected, trial_typ == 730),4));
        subplot(5,1,4)
        plot(c_data, 'Color', 'k');
        hold on;
        plot(mean(c_data, 2), 'Color', 'red');
        hold off
        xticks(0:0.1*sampleRate:size(c_data, 1))
        xticklabels(string((0:0.1*sampleRate:size(c_data, 1))/sampleRate));
        xlim([0 0.8*sampleRate])
        title('mean channel across trials')
    else
        window_size = (length(audio.data{2})-1)*1/audio.fs{2}; %s
        f = my_audio.data(cue_pos(i):ceil(cue_pos(i)+my_audio.fs*window_size));
        t = linspace(0,window_size,length(f));

        subplot(5,1,1)
        plot(t, f)
        title('Audio registered');

        subplot(5,1,2)
        plot((0:(length(audio.data{2})-1))*1/audio.fs{2}, audio.data{2})
        title('original audio')

        subplot(5,1,3)
        melSpectrogram(double(f), my_audio.fs)
        title('spectrogram')

        c_data = squeeze(mean(trial_data(1:ceil((length(audio.data{2})-1)/audio.fs{1} * sampleRate), 1, channelsSelected, trial_typ == 731),4));
        subplot(5,1,4)
        plot(c_data, 'Color', 'k');
        hold on;
        plot(mean(c_data, 2), 'Color', 'red');
        hold off
        xticks(0:0.1*sampleRate:size(c_data, 1))
        xticklabels(string((0:0.1*sampleRate:size(c_data, 1))/sampleRate));
        xlim([0 0.7*sampleRate])
        title('mean channel across trials')

        subplot(5,1,5)
        melSpectrogram(double(audio.data{2}), audio.fs{2});
        title('Spect original')
    end

    sgtitle(num2str(cue_typ(i)))
end