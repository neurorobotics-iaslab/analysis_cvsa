%% compute the fischer to understand features using the log band
%% REASONING FOR THE SUSTAINED PART --> create the dataset
clc; clearvars;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/utils')

%% Initialization
create_dataset = true;
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

%% Labeling data for the dataset
data = cell(1, nbands);
for idx_band=1:nbands
    c_signal = signals{idx_band};
    c_header = headers{idx_band};

    c_data = extract_all(c_signal, c_header, classes, fix_event, cf_event, 1);
    data{idx_band} = c_data;
end

%% Trial extraction
ntrial = size(fixPOS, 1);
trial_start = nan(ntrial, 1);
trial_end = nan(ntrial, 1);
trial_typ = nan(ntrial, 1);
fix_start = nan(ntrial, 1);
fix_end = nan(ntrial, 1);
for idx_trial = 1:ntrial
    trial_start(idx_trial) = cuePOS(idx_trial);
    trial_typ(idx_trial) = cueTYP(idx_trial);
    trial_end(idx_trial) = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;
    fix_start(idx_trial) = fixPOS(idx_trial);
    fix_end(idx_trial) = fixPOS(idx_trial) + fixDUR(idx_trial) - 1;
end

min_trial_data = min(trial_end - trial_start);
min_fix_data = min(fix_end - fix_start);
trial_data = nan(min_trial_data, nbands, nchannels, ntrial);
basline_data = nan(min_fix_data, nbands, nchannels, ntrial);
for idx_band = 1:nbands
    c_signal = signals{idx_band};
    for trial = 1:ntrial
        c_start = trial_start(trial);
        c_end = trial_start(trial) + min_trial_data - 1;
        trial_data(:,idx_band,:,trial) = c_signal(c_start:c_end,:);
        c_fix_start = fix_start(trial);
        c_fix_end = fix_start(trial) + min_fix_data - 1;
        basline_data(:,idx_band,:,trial) = c_signal(c_fix_start:c_fix_end,:);
    end
end
trial_data_train = trial_data(:,:,:,1:end-ntrial/nFiles);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% last file is for the test
ntrial_train = ntrial - ntrial/nFiles;
trial_typ  = trial_typ(1:ntrial_train);
% trial_data_train = trial_data;

%% fisher for the susteined part
time = 0.3; %s
minDurCue = min(cueDUR);
c_start = minDurCue + ceil(sampleRate * time) + 1;
mean_c1 = nan(nchannelsSelected, nbands);
mean_c2 = nan(nchannelsSelected, nbands);
std_c1 = nan(nchannelsSelected, nbands);
std_c2 = nan(nchannelsSelected, nbands);
for idx_ch=1:nchannelsSelected
    idx_chSel = channelsSelected(idx_ch);

    for idx_band=1:nbands
        mean_c1(idx_ch, idx_band) = mean(trial_data_train(c_start:end,idx_band,idx_chSel,trial_typ == classes(1)), 'all');
        mean_c2(idx_ch, idx_band) = mean(trial_data_train(c_start:end,idx_band,idx_chSel,trial_typ == classes(2)), 'all');
        std_c1(idx_ch, idx_band) = std(trial_data_train(c_start:end,idx_band,idx_chSel,trial_typ == classes(1)), 0, 'all');
        std_c2(idx_ch, idx_band) = std(trial_data_train(c_start:end,idx_band,idx_chSel,trial_typ == classes(2)), 0, 'all');
    end
end
figure();
fisher_all_calib = abs(mean_c1 - mean_c2) ./ sqrt(std_c1.^2 + std_c2.^2);
imagesc(fisher_all_calib);
ylabel('channels');
ynew = 1:nchannelsSelected;
set(gca, 'YTick',ynew, 'YTickLabel', channels_select);
xlabel('frequencies');
xnew = 1:nbands;
x_label = bands_str;
set(gca, 'XTick',xnew, 'XTickLabel',x_label);
title([subject ' | fisher score for cf | CALIBRATION FILES']);


%% Fisher score in all the trial (cue + cf) and for each window 
logband_hz = sampleRate;
wlength = 0.10; %s
wlength = ceil(wlength * logband_hz);
overlap = 0.05; %s
overlap = ceil(overlap * logband_hz);
nwindows = round((min_trial_data - wlength)/overlap) + 1;
windows_center = zeros(1, nwindows);
fischers = zeros(nchannelsSelected, nbands, nwindows);
typ = nan(nwindows, 1);
for idx_w = 1:nwindows
    mean_c1 = nan(nchannelsSelected, nbands);
    mean_c2 = nan(nchannelsSelected, nbands);
    std_c1 = nan(nchannelsSelected, nbands);
    std_c2 = nan(nchannelsSelected, nbands);
    start_w = (idx_w-1)*overlap + 1;
    end_w = start_w + wlength - 1;
    if end_w > size(trial_data_train, 1)
        end_w = size(trial_data_train, 1);
    end

    for idx_ch=1:nchannelsSelected
        idx_chSel = channelsSelected(idx_ch);
        for idx_band=1:nbands
            mean_c1(idx_ch, idx_band) = mean(trial_data_train(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(1)), 'all');
            mean_c2(idx_ch, idx_band) = mean(trial_data_train(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(2)), 'all');
            std_c1(idx_ch, idx_band) = std(trial_data_train(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(1)), 0, 'all');
            std_c2(idx_ch, idx_band) = std(trial_data_train(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(2)), 0, 'all');
        end
    end

    windows_center(idx_w) = (start_w + end_w) / 2;
    fischers(:,:,idx_w) = abs(mean_c1 - mean_c2) ./ sqrt(std_c1.^2 + std_c2.^2);

    if start_w > min(cueDUR)
        typ(idx_w) = cf_event;
    end
end

%% metrics for the features selection all x w
max_nfeatures_allxw = 27;
allxw_img_calib = zeros(max_nfeatures_allxw, nwindows);
for c_nf=1:max_nfeatures_allxw
    [tot_sortedValues, tot_linearIndices_calib] = maxk(fisher_all_calib(:), c_nf);

    for idx_w=1:nwindows
        c_fischer = fischers(:,:,idx_w); % take only interested values
        [c_sortedValues, c_linearIndices] = maxk(c_fischer(:), c_nf);
        [c_rows, c_cols] = ind2sub(size(c_fischer), c_linearIndices);

        n_intersection = length(intersect(c_linearIndices, tot_linearIndices_calib));
        n_union = length(union(c_linearIndices, tot_linearIndices_calib));
        allxw_img_calib(c_nf, idx_w) = n_intersection / c_nf;
    end
end
max_allxw_calib = max(abs(allxw_img_calib), [], 'all');
allxw_img_calib = allxw_img_calib / max_allxw_calib;
max_allxw_calib = max(abs(allxw_img_calib), [], 'all');

time_diff = zeros(1, nwindows);
time_diff(2:end) = ((windows_center(1:end-1) + windows_center(2:end)) / 2) / sampleRate;

figure();
cf_start = find(typ == cf_event, 1);
imagesc(allxw_img_calib);
hold on;
plot([cf_start cf_start], [0 max_nfeatures_allxw], 'k');
hold off
xticks(1:5:size(time_diff,2))
xticklabels(time_diff(1:5:end))
yticks(1:max_nfeatures_allxw)
yticklabels(1:max_nfeatures_allxw)
set(gca, 'clim', [0 max_allxw_calib])
xlabel('time [s]')
title([subject ' | window features vs all cf features | CALIBRATION FILES'])


%% save the datasets %%%%%%%%%%%%%%%%%%
if create_dataset
    min_durCF = min(cfDUR);
    start_period = 1; % s ---> time to start the creation of dataset (from the end of cue)
    period = [start_period*sampleRate, min_durCF-1]; % period of interest in samples
    period_dur = period(2) - period(1);
    info.files = filenames;
    info.channelsLabel = channels_label;
    info.startTrial    = nan(ntrial,1);
    for c_nf=1:max_nfeatures_allxw
        [tot_sortedValues, tot_linearIndices_calib] = maxk(fisher_all_calib(:), c_nf);
        [f_idx_selChans, f_idx_band] = ind2sub(size(fisher_all_calib), tot_linearIndices_calib);

        info.bandSelected = bands(f_idx_band);
        info.bandSelected = vertcat(info.bandSelected{:});
        info.chSelected = channelsSelected(f_idx_selChans);

        X = nan(period_dur, c_nf);
        y = nan(period_dur, 1);

        for idx_trial = 1:ntrial
            for idx_f = 1:c_nf
                c_channel = channelsSelected(f_idx_selChans(idx_f));
                c_band = f_idx_band(idx_f);
                X((idx_trial-1)*period_dur+1:idx_trial*period_dur, idx_f) = ...
                    trial_data(minDurCue+period(1)+1:minDurCue+period(2),c_band,c_channel,idx_trial);
            end
            y((idx_trial-1)*period_dur+1:idx_trial*period_dur) = repmat(data{1}.typ(idx_trial), period_dur, 1);

            info.startTrial(idx_trial) = (idx_trial-1)*period_dur +1;
        end
        info.sampleRate = sampleRate;
        info.filterOrder = filterOrder;
        save_path = [pathname(1:end-4) 'two_classifier/dataset/sustained/sus_data' num2str(c_nf) '.mat'];
        save(save_path, 'X', 'y', 'info')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- reasoning for EVALUATION ---------------- %%
% Initialization
signals = cell(1, nbands);
headers = cell(1, nbands);
for idx_band = 1:nbands
    headers{idx_band}.TYP = [];
    headers{idx_band}.DUR = [];
    headers{idx_band}.POS = [];
    signals{idx_band} = [];
end

pathname = [pathname(1:end-16) 'evaluation/gdf/'];
files = dir(fullfile(pathname, '*.gdf'));
% concatenate the files
nFiles = length(files);
for idx_file= 1: nFiles
    fullpath_file = fullfile(pathname, files(1).name);
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', files(1).name]);
    [c_signal,header] = sload(fullpath_file);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;

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


% Labelling data 
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

% Trial extraction
ntrial = size(fixPOS, 1);
trial_start = nan(ntrial, 1);
trial_end = nan(ntrial, 1);
trial_typ = nan(ntrial, 1);
fix_start = nan(ntrial, 1);
fix_end = nan(ntrial, 1);
for idx_trial = 1:ntrial
    trial_start(idx_trial) = cuePOS(idx_trial);
    trial_typ(idx_trial) = cueTYP(idx_trial);
    trial_end(idx_trial) = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;
    fix_start(idx_trial) = fixPOS(idx_trial);
    fix_end(idx_trial) = fixPOS(idx_trial) + fixDUR(idx_trial) - 1;
end

min_trial_data = min(trial_end - trial_start);
min_fix_data = min(fix_end - fix_start);
trial_data = nan(min_trial_data, nbands, nchannels, ntrial);
basline_data = nan(min_fix_data, nbands, nchannels, ntrial);
for idx_band = 1:nbands
    c_signal = signals{idx_band};
    for trial = 1:ntrial
        c_start = trial_start(trial);
        c_end = trial_start(trial) + min_trial_data - 1;
        trial_data(:,idx_band,:,trial) = c_signal(c_start:c_end,:);
        c_fix_start = fix_start(trial);
        c_fix_end = fix_start(trial) + min_fix_data - 1;
        basline_data(:,idx_band,:,trial) = c_signal(c_fix_start:c_fix_end,:);
    end
end

%% Fisher score in all the trial (cue + cf) and for each window
nwindows = round((min_trial_data - wlength)/overlap) + 1;
windows_center = zeros(1, nwindows);
fischers = zeros(nchannelsSelected, nbands, nwindows);
typ = nan(nwindows, 1);
for idx_w = 1:nwindows
    mean_c1 = nan(nchannelsSelected, nbands);
    mean_c2 = nan(nchannelsSelected, nbands);
    std_c1 = nan(nchannelsSelected, nbands);
    std_c2 = nan(nchannelsSelected, nbands);
    start_w = (idx_w-1)*overlap + 1;
    end_w = start_w + wlength - 1;
    if end_w > size(trial_data, 1)
        end_w = size(trial_data, 1);
    end

    for idx_ch=1:nchannelsSelected
        idx_chSel = channelsSelected(idx_ch);
        for idx_band=1:nbands
            mean_c1(idx_ch, idx_band) = mean(trial_data(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(1)), 'all');
            mean_c2(idx_ch, idx_band) = mean(trial_data(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(2)), 'all');
            std_c1(idx_ch, idx_band) = std(trial_data(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(1)), 0, 'all');
            std_c2(idx_ch, idx_band) = std(trial_data(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(2)), 0, 'all');
        end
    end

    windows_center(idx_w) = (start_w + end_w) / 2;
    fischers(:,:,idx_w) = abs(mean_c1 - mean_c2) ./ sqrt(std_c1.^2 + std_c2.^2);

    if start_w > min(cueDUR)
        typ(idx_w) = cf_event;
    end
end

% fisher for all cf trial no cue
minDurCue = min(cueDUR);
c_start = ceil(time*sampleRate) + 1 + minDurCue;
for idx_ch=1:nchannelsSelected
    idx_chSel = channelsSelected(idx_ch);
    for idx_band=1:nbands
        mean_c1(idx_ch, idx_band) = mean(trial_data(c_start:end,idx_band,idx_chSel,trial_typ == classes(1)), 'all');
        mean_c2(idx_ch, idx_band) = mean(trial_data(c_start:end,idx_band,idx_chSel,trial_typ == classes(2)), 'all');
        std_c1(idx_ch, idx_band) = std(trial_data(c_start:end,idx_band,idx_chSel,trial_typ == classes(1)), 0, 'all');
        std_c2(idx_ch, idx_band) = std(trial_data(c_start:end,idx_band,idx_chSel,trial_typ == classes(2)), 0, 'all');
    end
end
figure();
fisher_all_eval = abs(mean_c1 - mean_c2) ./ sqrt(std_c1.^2 + std_c2.^2);
imagesc(fisher_all_eval);
ylabel('channels');
ynew = 1:nchannelsSelected;
set(gca, 'YTick',ynew, 'YTickLabel', channels_select);
xlabel('frequencies');
xnew = 1:nbands;
x_label = bands_str;
set(gca, 'XTick',xnew, 'XTickLabel',x_label);
title([subject ' | fisher score for cf | EVALUATION FILES ']);

%% metrics for the features selection all x w
allxw_img_calib = zeros(max_nfeatures_allxw, nwindows);
allxw_img_eval = zeros(max_nfeatures_allxw, nwindows);
for c_nf=1:max_nfeatures_allxw
    [~, tot_linearIndices_calib] = maxk(fisher_all_calib(:), c_nf);
%     [tot_rows, tot_cols] = ind2sub(size(fisher_all_calib), tot_linearIndices_calib);
    for idx_w=1:nwindows
        c_fischer = fischers(:,:,idx_w); % take only interested values
        [c_sortedValues, c_linearIndices] = maxk(c_fischer(:), c_nf);
        [c_rows, c_cols] = ind2sub(size(c_fischer), c_linearIndices);

        n_intersection = length(intersect(c_linearIndices, tot_linearIndices_calib));
        n_union = length(union(c_linearIndices, tot_linearIndices_calib));
        allxw_img_calib(c_nf, idx_w) = n_intersection / c_nf;
    end

    [~, tot_linearIndices_eval] = maxk(fisher_all_eval(:), c_nf);
%     [tot_rows, tot_cols] = ind2sub(size(fisher_all_calib), tot_linearIndices_calib);
    for idx_w=1:nwindows
        c_fischer = fischers(:,:,idx_w); % take only interested values
        [c_sortedValues, c_linearIndices] = maxk(c_fischer(:), c_nf);
%         [c_rows, c_cols] = ind2sub(size(c_fischer), c_linearIndices);

        n_intersection = length(intersect(c_linearIndices, tot_linearIndices_eval));
        n_union = length(union(c_linearIndices, tot_linearIndices_eval));
        allxw_img_eval(c_nf, idx_w) = n_intersection / c_nf;
    end
end
max_allxw_calib = max(abs(allxw_img_calib), [], 'all');
allxw_img_calib = allxw_img_calib / max_allxw_calib;
max_allxw_calib = max(abs(allxw_img_calib), [], 'all');

max_allxw_eval = max(abs(allxw_img_eval), [], 'all');
allxw_img_eval = allxw_img_eval / max_allxw_eval;
max_allxw_eval = max(abs(allxw_img_eval), [], 'all');

time_diff = zeros(1, nwindows);
time_diff(2:end) = ((windows_center(1:end-1) + windows_center(2:end)) / 2) / sampleRate;

figure();
subplot(1,2,1)
cf_start = find(typ == cf_event, 1);
imagesc(allxw_img_calib);
hold on;
plot([cf_start cf_start], [0 max_nfeatures_allxw], 'k');
hold off
xticks(1:5:size(time_diff,2))
xticklabels(time_diff(1:5:end))
yticks(1:max_nfeatures_allxw)
yticklabels(1:max_nfeatures_allxw)
set(gca, 'clim', [0 max_allxw_calib])
xlabel('time [s]')
title('window features (EVA) vs all cf features (CALIB)')

subplot(1,2,2)
cf_start = find(typ == cf_event, 1);
imagesc(allxw_img_eval);
hold on;
plot([cf_start cf_start], [0 max_nfeatures_allxw], 'k');
hold off
xticks(1:5:size(time_diff,2))
xticklabels(time_diff(1:5:end))
yticks(1:max_nfeatures_allxw)
yticklabels(1:max_nfeatures_allxw)
set(gca, 'clim', [0 max_allxw_eval])
xlabel('time [s]')
title('window features (EVAL) vs all cf features (EVAL)')

sgtitle(subject)
