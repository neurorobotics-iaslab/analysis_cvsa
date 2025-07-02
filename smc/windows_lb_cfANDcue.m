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
    trial_start(idx_trial) = cuePOS(idx_trial);
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

%% Fisher score in all the trial (cue + cf) and for each window
nrows = 4;
logband_hz = sampleRate;
wlength = 0.10; %s
wlength = ceil(wlength * logband_hz);
overlap = 0.05; %s
overlap = ceil(overlap * logband_hz);
nwindows_calib = round((min_trial_data - wlength)/overlap) + 1;
windows_center = zeros(1, nwindows_calib);plot_downsampling = 2;
handles = [];
cl = -inf;
fig1 = figure('Units','normalized','OuterPosition',[0 0 1 1], 'Visible','off');
fischers_calib = zeros(nchannelsSelected, nbands, nwindows_calib);
typ = nan(nwindows_calib, 1);
for idx_w = 1:nwindows_calib
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
    fischers_calib(:,:,idx_w) = abs(mean_c1 - mean_c2) ./ sqrt(std_c1.^2 + std_c2.^2);
    if mod(idx_w, plot_downsampling) == 0
        % calculate fisher
        subplot(nrows, ceil(nwindows_calib/nrows/plot_downsampling), idx_w/plot_downsampling)
        c_fisher = fischers_calib(:,:,idx_w);

        % show the fisher features
        imagesc(c_fisher);
        ylabel('channels');
        ynew = 1:nchannelsSelected;
        set(gca, 'YTick',ynew, 'YTickLabel', channels_select);
        xlabel('frequencies');
        xnew = 1:nbands;
        x_label = bands_str;
        set(gca, 'XTick',xnew, 'XTickLabel',x_label);
        handles = [handles gca];
        cl = max(cl, max(abs(c_fisher), [], 'all'));
        if start_w > min(cueDUR)
            title(['cf | ' num2str(idx_w)])
        else
            title(['cue | ' num2str(idx_w)])
        end
    end

    if start_w > min(cueDUR)
        typ(idx_w) = cf_event;
    end
end
set(handles, 'clim', [0 cl])
sgtitle([subject ' | fisher score for some windows | CALIBRATION FILES']);
path_figure = ['/home/paolo/cvsa_ws/img_results/27_features/' subject '/allFile_allCf/windows_fisher.svg'];
print(fig1, path_figure, '-dsvg');
close(fig1);

%% fisher for all cf trial no cue
minDurCue = min(cueDUR);
% trial_data_train = trial_data(:,:,:,1:end-20);
% trial_typ_train = trial_typ(1:end-20);
for idx_ch=1:nchannelsSelected
    idx_chSel = channelsSelected(idx_ch);
    for idx_band=1:nbands
        mean_c1(idx_ch, idx_band) = mean(trial_data(minDurCue:end,idx_band,idx_chSel,trial_typ == classes(1)), 'all');
        mean_c2(idx_ch, idx_band) = mean(trial_data(minDurCue:end,idx_band,idx_chSel,trial_typ == classes(2)), 'all');
        std_c1(idx_ch, idx_band) = std(trial_data(minDurCue:end,idx_band,idx_chSel,trial_typ == classes(1)), 0, 'all');
        std_c2(idx_ch, idx_band) = std(trial_data(minDurCue:end,idx_band,idx_chSel,trial_typ == classes(2)), 0, 'all');
    end
end
fig2 = figure('Visible','off');
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
path_figure = ['/home/paolo/cvsa_ws/img_results/27_features/' subject '/allFile_allCf/fisher_allCF.svg'];
print(fig2, path_figure, '-dsvg');
% close(fig2);

%% metrics for the features selection all x w
max_nfeatures_allxw = 27;
allxw_img_calib = zeros(max_nfeatures_allxw, nwindows_calib);
for c_nf=1:max_nfeatures_allxw
    [tot_sortedValues, tot_linearIndices_calib] = maxk(fisher_all_calib(:), c_nf);

    for idx_w=1:nwindows_calib
        c_fischer = fischers_calib(:,:,idx_w); % take only interested values
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

time_diff = zeros(1, nwindows_calib);
time_diff(2:end) = ((windows_center(1:end-1) + windows_center(2:end)) / 2) / sampleRate;

fig3 = figure('Units','normalized','OuterPosition',[0 0 1 1]);
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
path_figure = ['/home/paolo/cvsa_ws/img_results/27_features/' subject '/allFile_allCf/fisher_allxw.svg'];
print(fig3, path_figure, '-dsvg');


%% metrics for the features selection w x w
max_nfeatures_wxw = 27;
min_nfeatures_wxw = 1;
nfeatures_wxw = min_nfeatures_wxw:max_nfeatures_wxw;
nrows = 3;
fig4 = figure('Units','normalized','OuterPosition',[0 0 1 1]);
for idx_c_nf = 1:3:length(nfeatures_wxw)
    c_nf = nfeatures_wxw(idx_c_nf);
    wxw_img = zeros(nwindows_calib);
    for idx_w1=1:nwindows_calib
        c1_fischer = fischers_calib(:,:,idx_w1);
        [c1_sortedValues, c1_linearIndices] = maxk(c1_fischer(:), c_nf);
        for idx_w2=1:nwindows_calib
            c2_fischer = fischers_calib(:,:,idx_w2);
            [c2_sortedValues, c2_linearIndices] = maxk(c2_fischer(:), c_nf);

            n_wxw_intersection = length(intersect(c1_linearIndices, c2_linearIndices));
            n_wxw_union = length(union(c1_linearIndices, c2_linearIndices));

            wxw_img(idx_w1,idx_w2) = n_wxw_intersection / c_nf;
        end
    end
    max_wxw = max(abs(wxw_img), [], 'all');
    wxw_img = wxw_img / max_wxw;
    max_wxw = max(abs(wxw_img), [], 'all');

    subplot(nrows, ceil((max_nfeatures_wxw/3-min_nfeatures_wxw + 1)/nrows), (idx_c_nf-1)/3 + 1)
    hold on;
    imagesc(wxw_img);
    plot([cf_start, cf_start], [0 nwindows_calib], 'k');
    hold off
    set(gca, 'clim', [0 max_wxw])
    xticks(1:7:size(windows_center,2))
    xticklabels(windows_center(1:7:end)/sampleRate)
    yticks(1:7:size(windows_center,2))
    yticklabels(windows_center(1:7:end)/sampleRate)
    title(['n features: ' num2str(c_nf)])
end
sgtitle([subject ' | overlap ratio between windows | CALIBRATION FILES'])
path_figure = ['/home/paolo/cvsa_ws/img_results/27_features/' subject '/allFile_allCf/fisher_wxw.svg'];
print(fig4, path_figure, '-dsvg');
