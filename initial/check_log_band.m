%% compute the fischer to understand features using the log band
%% REASONING IN ALL THE TRIAL
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
% channels_select = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
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
trial_with_eog = [];
nFiles = length(filenames);
for idx_file= 1: nFiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file);
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

%% Trial extraction
ntrial = size(fixPOS, 1);
trial_start = nan(ntrial, 1);
trial_end = nan(ntrial, 1);
trial_typ = nan(ntrial, 1);
for idx_trial = 1:ntrial
    trial_start(idx_trial) = fixPOS(idx_trial);
    trial_typ(idx_trial) = cueTYP(idx_trial);
    trial_end(idx_trial) = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;
end

min_trial_data = min(trial_end - trial_start);
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

trial_data = trial_data(:,:,:,1:end-ntrial_test);
trial_typ  = trial_typ(1:ntrial - ntrial_test);

% %%
% chanlocs_path = '/home/paolo/chanlocs39.mat';
% load(chanlocs_path);
% trial_data_1 = mean(trial_data(:,:,:,trial_typ == classes(1)), 4);
% trial_data_2 = mean(trial_data(:,:,:,trial_typ == classes(2)), 4);
% data = trial_data_2 - trial_data_1;
% 
% nsignal = size(data, 1);
% 
% 
% fig = figure('Name', 'Fisher Score Topoplots', 'NumberTitle', 'off', ...
%     'Position', [100 100 1200 800]);
% 
% axHandles = gobjects(nbands, 1);
% for b = 1:nbands
%     axHandles(b) = subplot(ceil(nbands/2), 2, b);  % adjustable layout
% end
% 
% % Slider
% uicontrol('Style', 'slider', ...
%     'Min', 1, 'Max', nsignal, 'Value', 1, ...
%     'Units', 'normalized', ...
%     'Position', [0.1 0.02 0.8 0.04], ...
%     'SliderStep', [1/(nsignal-1) , 10/(nsignal-1)], ...
%     'Callback', @(src, event) updatePlots(round(get(src, 'Value')), ...
%         nbands, chanlocs, data, axHandles, bands_str, channelsSelected));
% 
% % Initial plot
% updatePlots(1, nbands, chanlocs, data, axHandles, bands_str, channelsSelected);

%%
chanlocs_path = '/home/paolo/chanlocs39.mat';
load(chanlocs_path);
trial_data_1 = mean(trial_data(:,:,:,trial_typ == classes(1)), 4);
trial_data_2 = mean(trial_data(:,:,:,trial_typ == classes(2)), 4);
data = trial_data_2 - trial_data_1;

nsignal = size(data, 1);
outputVideo = VideoWriter('diff_class2_class1.avi');  % or .mp4 if supported
outputVideo.FrameRate = 16;  % adjust playback speed
open(outputVideo);

% Create invisible figure for faster rendering
fig = figure('Visible', 'off', 'Position', [100 100 1200 800]);

axHandles = gobjects(nbands, 1);
for b = 1:nbands
    axHandles(b) = subplot(ceil(nbands/2), 2, b);
end

for winIdx = 1:32:nsignal
    updatePlots(winIdx, nbands, chanlocs, data, axHandles, bands_str, channelsSelected);
    frame = getframe(fig);
    writeVideo(outputVideo, frame);
end

close(outputVideo);
disp('✅ Video saved successfully!');


% --- Function to update plots ---
% --- Optimized updatePlots ---
function updatePlots(winIdx, nBands, chanlocs, data_all, axHandles, bands_str, channelsSelected)
%     data = zeros(39, 1);
    for b = 1:nBands
        cla(axHandles(b));
        data = data_all(winIdx, b, :);
        axes(axHandles(b));
        topoplot(data, chanlocs, ...
            'electrodes', 'off', ...
            'style', 'map', ...
            'numcontour', 0, ...
            'gridscale', 32);  % lower resolution
        title(bands_str{b});
        caxis([-0.5, 0.5]);
        colorbar;
    end
    if winIdx < 1024
        sgtitle(['Fixation', ' - time ', num2str(winIdx/512)])
    elseif winIdx < 1024+512
        sgtitle(['cue', ' - time ', num2str(winIdx/512)])
    else
        sgtitle(['cf', ' - time ', num2str(winIdx/512)])
    end
    drawnow;
end


