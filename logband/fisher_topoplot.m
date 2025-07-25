%% show the fisher score and the topoplot of a log band signal
clear all; % close all;

addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/EOG')
addpath('/home/paolo/Local/cnbi-smrtrain/toolboxes/cva')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')

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
% load('/home/paolo/lap_39.mat')

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

%         c_signal = c_signal * lap;

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
trial_data = tmp; % samples x bands x channels x trials

%% Computing fisher
normalize_std = false;
fisher = nan(nbands, nchannels);
mu = nan(nclasses, nchannels, nbands);
sigma = nan(nclasses, nchannels, nbands);

for idx_band = 1:nbands
    for idx_ch=1:nchannels
        c_data = trial_data(min_durFIX + min_durCUE:end,idx_band, idx_ch, trial_typ == classes(1));
        mu(1, idx_ch, idx_band) = mean(squeeze(c_data), 'all');
        sigma(1, idx_ch, idx_band) = std(squeeze(c_data(:)));

        c_data = trial_data(min_durFIX + min_durCUE:end,idx_band, idx_ch, trial_typ == classes(2));
        mu(2, idx_ch, idx_band) = mean(squeeze(c_data), 'all');
        sigma(2, idx_ch, idx_band) = std(squeeze(c_data(:)));
    end
end

% compute fisher
if normalize_std
    for idx_band = 1:nbands
        for idx_ch = 1:nchannels
            fisher(idx_band, idx_ch) = abs(mu(1, idx_ch, idx_band) - mu(2, idx_ch, idx_band)) ./ sqrt(sigma(1, idx_ch, idx_band).^2 + sigma(2, idx_ch, idx_band).^2);
        end
    end
else
    for idx_band = 1:nbands
        for idx_ch = 1:nchannels
            fisher(idx_band, idx_ch) = abs(mu(1, idx_ch, idx_band) - mu(2, idx_ch, idx_band));
        end
    end
end

% with cva
cva = nan(nbands, nchannels);
reshaped_dur = min_trial_data - (min_durFIX + min_durCUE - 1);
for idx_band= 1:nbands
    reshaped_signal = nan(reshaped_dur * ntrial, nchannels);
    ck = [];
    for idx_trial = 1:ntrial
        reshaped_signal((idx_trial - 1)*reshaped_dur + 1 : idx_trial*reshaped_dur, :) = trial_data(min_durFIX + min_durCUE:end,idx_band,:,idx_trial);
        ck = cat(1, ck, repmat(trial_typ(idx_trial), [reshaped_dur,1]));
    end

    cva(idx_band, :) = cva_tun_opt(reshaped_signal, ck);
end

%% Computing the log band diff as topoplot
chanlocs_path = '/home/paolo/chanlocs39.mat';
load(chanlocs_path);
diff = nan(nbands, nchannels);
for idx_band = 1:nbands
    c_data_1 = mean(mean(squeeze(trial_data(min_durFIX + min_durCUE:end,idx_band,:,trial_typ == classes(1))), 3), 1); % mean on trials and then on samples
    c_data_2 = mean(mean(squeeze(trial_data(min_durFIX + min_durCUE:end,idx_band,:,trial_typ == classes(2))), 3), 1);
 
    diff(idx_band,:) = c_data_1 - c_data_2;
end


%% show fisher, cva and log band diff for all the channels
figure();
subplot(1,2,1)
imagesc(fisher');
yticks(1:nchannels);
yticklabels(channels_label(1:nchannels));
xticks(1:nbands);
xticklabels(bands_str);
title('Fisher score');

subplot(1,2,2)
imagesc(cva')
yticks(1:nchannels);
yticklabels(channels_label(1:nchannels));
xticks(1:nbands);
xticklabels(bands_str);
title('CVA')
sgtitle('Only cf');


figure();
rows_plot = 2;
handles = []; cl = -inf;
for idx_band = 1:nbands
    chanlocs_data = diff(idx_band, :);
    subplot(rows_plot, ceil(nbands/rows_plot), idx_band)
    topoplot(squeeze(chanlocs_data), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(chanlocs_data)) max(abs(chanlocs_data))]);
    cl = max(cl, max(abs(chanlocs_data)));
    handles = [handles gca];
    axis image;
    colorbar;
    title(['band: ' bands_str{idx_band}])
end
set(handles, 'clim', [-cl, cl]);
sgtitle('Diff log band (bl - br) | only cf');


figure();
idx_plot = 1; ncols = 2; handles = []; cl = -inf;
for idx_band = 1:nbands
    data_1 = squeeze(mean(trial_data(:,idx_band,:,trial_typ == classes(1)), 4));
    data_2 = squeeze(mean(trial_data(:,idx_band,:,trial_typ == classes(2)), 4));
    c_diff = abs(data_1 - data_2);

    subplot(ceil(nbands/ncols), 2, idx_plot)
    imagesc(c_diff')
    hold on
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off
    yticks(1:nchannels)
    yticklabels(channels_label(1:end-1))
    xticks(sampleRate:sampleRate:size(c_diff, 1))
    xticklabels(string((sampleRate:sampleRate:size(c_diff, 1)) / sampleRate));
    idx_plot = idx_plot  + 1;
    handles = [handles; gca];
    cl = max(cl, max(abs(c_diff), [], 'all'));
    title(['band: ' bands_str{idx_band}]);
end
set(handles, 'clim', [0, cl])
sgtitle('Diff (bl - br) in time')


%% show fisher, cva and log band diff for the P, PO and O channels
figure();
subplot(1,2,1)
imagesc(fisher(:,channelsSelected)');
yticks(1:nchannelsSelected);
yticklabels(channels_select);
xticks(1:nbands);
xticklabels(bands_str);
title('Fisher score');

subplot(1,2,2)
imagesc(cva(:,channelsSelected)')
yticks(1:nchannelsSelected);
yticklabels(channels_select);
xticks(1:nbands);
xticklabels(bands_str);
title('CVA')
sgtitle('Only cf');


figure();
rows_plot = 2;
handles = []; cl = -inf;
for idx_band = 1:nbands
    chanlocs_data = zeros(1, nchannels);
    chanlocs_data(channelsSelected) = diff(idx_band,channelsSelected);
    subplot(rows_plot, ceil(nbands/rows_plot), idx_band)
    topoplot(squeeze(chanlocs_data), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(chanlocs_data)) max(abs(chanlocs_data))]);
    cl = max(cl, max(abs(chanlocs_data)));
    handles = [handles gca];
    axis image;
    colorbar;
    title(['band: ' bands_str{idx_band}])
end
set(handles, 'clim', [-cl, cl]);
sgtitle('Diff log band (bl - br) | only cf');


figure();
idx_plot = 1; ncols = 2; handles = []; cl = -inf;
for idx_band = 1:nbands
    data_1 = squeeze(mean(squeeze(trial_data(:,idx_band,channelsSelected,trial_typ == classes(1))), 3));
    data_2 = squeeze(mean(squeeze(trial_data(:,idx_band,channelsSelected,trial_typ == classes(2))), 3));
    c_diff = abs(data_1 - data_2);

    subplot(ceil(nbands/ncols), 2, idx_plot)
    imagesc(c_diff')
    hold on
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off
    yticks(1:nchannels)
    yticklabels(channels_select)
    xticks(sampleRate:sampleRate:size(c_diff, 1))
    xticklabels(string((sampleRate:sampleRate:size(c_diff, 1)) / sampleRate));
    idx_plot = idx_plot  + 1;
    handles = [handles; gca];
    cl = max(cl, max(abs(c_diff), [], 'all'));
    title(['band: ' bands_str{idx_band}]);
end
set(handles, 'clim', [0, cl])
sgtitle('Diff (bl - br) in time')
