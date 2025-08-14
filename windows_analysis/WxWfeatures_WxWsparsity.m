%% See the features selected in windows are they are sparse in the brain
clear all; % close all;

addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/EOG')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')

%% Initialization
bands = [{[4 8]} {[8 14]} {[14 30]}];
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
avg = 0.5;
eog_threshold = 500;

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);
day = filenames{1}(4:11);

%% concatenate the files
nFiles = length(filenames);
trial_with_eog = [];
for idx_file= 1: nFiles
    fullpath_file_shift = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_shift);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;
    sampleRate = header.SampleRate;

    c_trial_with_eog = eog_detection(c_signal, header, eog_threshold, {'FP1', 'FP2', 'EOG'});
    trial_with_eog = [trial_with_eog; c_trial_with_eog];

    for idx_band = 1:nbands
        band = bands{idx_band};

        % for power band
        disp('   [proc] power band');
        disp('      [proc] applying filtering')
        [b, a] = butter(filterOrder, band(2)*(2/sampleRate),'low');
        s_low = filter(b,a,c_signal);
        [b, a] = butter(filterOrder, band(1)*(2/sampleRate),'high');
        s_filt = filter(b,a,s_low);
        disp('      [proc] applying power')
%         s_rect = power(s_filt, 2);
        analytic = hilbert(s_filt);
        s_rect = abs(analytic).^2;
        disp('      [proc] applying average window')
        s_out = zeros(size(c_signal));
        nchannels = size(c_signal, 2);
        for idx_ch=1:nchannels
            s_out(:, idx_ch) = (filter(ones(1,avg*sampleRate)/avg/sampleRate, 1, s_rect(:, idx_ch)));
        end
       

        signal_processed = s_out;
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
minDurFix = min(fixDUR);
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

%% variables
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
[~, channelsSelected] = ismember(channels_select, channels_label);
nchannelsSelected = size(channelsSelected, 2);
LblMontage = cell(12,13);
LblMontage(1,:) = {'','','','','','FP1','FPZ','FP2','','','','',''};
LblMontage(2,:) = {'','AF9','AF7','','AF3','','AFZ','','AF4','','AF8','AF10',''};
LblMontage(3,:) = {'','F9','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F10',''};
LblMontage(4,:) = {'','FT9','FT7','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT8','FT10',''};
LblMontage(5,:) = {'','T9','T7','C5','C3','C1','CZ','C2','C4','C6','T8','T10',''};
LblMontage(6,:) = {'M1','TP9','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','TP10','M2'};
LblMontage(7,:) = {'','P9','P7','P5','P3','P1','PZ','P2','P4','P6','P8','P10',''};
LblMontage(8,:) = {'','PO9','PO7','PO5','PO3','PO1','POZ','PO2','PO4','PO6','PO8','PO10',''};
LblMontage(9,:) = {'','O9','','','','O1','OZ','O2','','','','O10',''};
LblMontage(10,:) = {'','','','','','','IZ','','','','','',''};
LblMontage(11,:) = {'','','','','','','','','','','','',''};
LblMontage(12,:) = {'','','','','','','EOG','','','','','',''};

trial_data(:,:,[1 2 19],:) = 0; % remove eog


%% windowing reasoning
win_size = 100; %ms
win_size = win_size * sampleRate /1000;
overlap = 100; %ms
overlap = overlap * sampleRate / 1000;
nwin = round((min_trial_data - win_size)/overlap) + 1;

% fisher score windows and the log band mean
fishers_windows = zeros(nchannels, nbands, nwin); % channels x bands x nwindows
logBand_win = zeros(nchannels, nbands, nwin, nclasses); % channels x bands x nwindows
windows_label = cell(1, nwin);
windows_center = zeros(1, nwin);

gini_series_all= zeros(nbands,nwin);
s_entropy_logband_all = nan(nwin, nbands);
s_m1_logband_all = nan(nwin, nbands);
s_m2_logband_all = nan(nwin, nbands);
s_m3_logband_all = nan(nwin, nbands);

s_entropy_logband_occipital = nan(nwin, nbands);
gini_series_occipital = zeros(nbands, nwin);

for idx_w = 1:nwin
    mean_c1 = nan(nchannels, nbands);
    mean_c2 = nan(nchannels, nbands);
    std_c1 = nan(nchannels, nbands);
    std_c2 = nan(nchannels, nbands);
    start_w = ceil((idx_w-1)*overlap + 1);
    end_w = ceil(start_w + win_size - 1);
    if end_w > min_trial_data
        end_w = min_trial_data;
    end
    windows_label{idx_w} = [num2str((start_w-1)/512) '-' num2str(end_w/512)];

    for idx_ch=1:nchannels
        for idx_band=1:nbands
            mean_c1(idx_ch, idx_band) = squeeze(mean(mean(trial_data(start_w:end_w,idx_band,idx_ch,trial_typ == classes(1)), 4), 1)); % mean trial then signal
            mean_c2(idx_ch, idx_band) = squeeze(mean(mean(trial_data(start_w:end_w,idx_band,idx_ch,trial_typ == classes(2)), 4), 1));
            std_c1(idx_ch, idx_band) = squeeze(mean(std(trial_data(start_w:end_w,idx_band,idx_ch,trial_typ == classes(1)), 0, 4), 1));
            std_c2(idx_ch, idx_band) = squeeze(mean(std(trial_data(start_w:end_w,idx_band,idx_ch,trial_typ == classes(2)), 0, 4), 1));
        end
    end

    windows_center(idx_w) = round((start_w-1 + end_w) / 2 / sampleRate,2);
    fishers_windows(:,:,idx_w) = abs(mean_c1 - mean_c2) ./ sqrt(std_c1.^2 + std_c2.^2);

    logBand_win(:,:,idx_w, 1) = mean_c1;
    logBand_win(:,:,idx_w, 2) = mean_c2;
    
    for idx_band = 1:nbands
        fishers_windows(isnan(fishers_windows(:,idx_band,idx_w)), idx_band, idx_w) = 0; % remove nan value for eog
        gini_series_all(idx_band, idx_w) = gini(squeeze(fishers_windows(:,idx_band, idx_w)));
        [s_entropy_logband_all(idx_w,idx_band), s_m1_logband_all(idx_w,idx_band), s_m2_logband_all(idx_w,idx_band), s_m3_logband_all(idx_w,idx_band)] = compute_sparsity(squeeze(fishers_windows(:,idx_band, idx_w)));

        gini_series_occipital(idx_band, idx_w) = gini(squeeze(fishers_windows(channelsSelected,idx_band, idx_w)));
        [s_entropy_logband_occipital(idx_w,idx_band), ~, ~, ~] = compute_sparsity(squeeze(fishers_windows(channelsSelected,idx_band, idx_w)));
    end
end
diff_logBand_win = abs(logBand_win(:,:,:,1) - logBand_win(:,:,:,2)); % diff bl - br


%% compute the wxw
nfeatures_occipital = ceil(nchannelsSelected*0.2); % number features of fisher to select --------------- we use 1 band only
nfeatures_all = ceil(nchannels * 0.2);

wxw_img_occipital = zeros(nwin, nwin, nbands);
wxw_img_all = zeros(nwin, nwin, nbands);

w_nearity_occipital = zeros(nwin, nbands);
w_nearity_all = zeros(nwin, nbands);
w_my_nearity_all = zeros(nwin, nbands);
for idx_band = 1:nbands
    for idx_w1 = 1:nwin
        c1_fisher_occipital = squeeze(fishers_windows(channelsSelected,idx_band, idx_w1)); % channels x 1
        [v, c1_linearIndices_occipital] = maxk(c1_fisher_occipital(:), nfeatures_occipital);

        c1_fisher_all = squeeze(fishers_windows(:,idx_band,idx_w1));
        [~, c1_linearIndices_all] = maxk(c1_fisher_all(:), nfeatures_all);

        % wxw
        for idx_w2 = 1:nwin
            c2_fisher_occipital = squeeze(fishers_windows(channelsSelected,idx_band,idx_w2));
            [~, c2_linearIndices_occipital] = maxk(c2_fisher_occipital(:), nfeatures_occipital);

            c2_fisher_all = squeeze(fishers_windows(:,idx_band,idx_w2));
            [~, c2_linearIndices_all] = maxk(c2_fisher_all(:), nfeatures_all);

            n_wxw_intersection = length(intersect(c1_linearIndices_occipital, c2_linearIndices_occipital));
            wxw_img_occipital(idx_w1,idx_w2, idx_band) = n_wxw_intersection / nfeatures_occipital;

            n_wxw_intersection = length(intersect(c1_linearIndices_all, c2_linearIndices_all));
            wxw_img_all(idx_w1,idx_w2, idx_band) = n_wxw_intersection / nfeatures_all;
        end

        % compute the nearity of the selected window
        [~, distance, ~] = compute_distance(channels_select(c1_linearIndices_occipital), LblMontage, channels_label);
        w_nearity_occipital(idx_w1, idx_band) = distance ;%/ sum(c1_fisher_occipital(c1_linearIndices_occipital));

        [~, distance, my_distance] = compute_distance(channels_label(c1_linearIndices_all), LblMontage, channels_label);
        w_nearity_all(idx_w1, idx_band) = distance;% / sum(c1_fisher_all(c1_linearIndices_all));
        w_my_nearity_all(idx_w1, idx_band) = my_distance;
    end

    % normalize the value of the metric
%     w_nearity_all(:,idx_band) = (w_nearity_all(:,idx_band) - min(w_nearity_all(:,idx_band))) / max(w_nearity_all(:,idx_band));
%     w_nearity_occipital(:,idx_band) = (w_nearity_occipital(:,idx_band) - min(w_nearity_occipital(:,idx_band))) / max(w_nearity_occipital(:,idx_band));
end

%% plot wxw
cue_startWIN = round((min_durFIX - win_size)/overlap) + 1;
cf_startWIN  = round((min_durCUE + min_durFIX - win_size)/overlap) + 1;
xticks_labels = arrayfun(@(x) sprintf('%.2f', x), windows_center, 'UniformOutput', false);

figure()
for idx_band = 1:nbands
    subplot(5,nbands,idx_band)
    imagesc(squeeze(wxw_img_all(:,:,idx_band)));
    hold on;
    xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    yticks(1:10:nwin); yticklabels(xticks_labels(1:10:end));
    xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
    title(['window x window features intersection | ' bands_str{idx_band}])

    subplot(5,nbands,idx_band+nbands)
    tmp = squeeze(w_nearity_all(:,idx_band)); 
%     tmp = tmp - min(tmp); tmp = tmp / max(tmp);
    plot(tmp) 
    tmp = squeeze(w_my_nearity_all(:,idx_band)); 
%     tmp = tmp - min(tmp); tmp = tmp / max(tmp);
    hold on;
    plot(tmp) 
    xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    legend([{'distance'}, {'sparsity roi'}])
    ylim([0 1])
    xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
    xlim([1, length(w_nearity_all(:,idx_band))])
    title(['window features nearity | ' bands_str{idx_band}])

    subplot(5,nbands,idx_band+2*nbands)
    plot(squeeze(gini_series_all(idx_band,:)))
    hold on;
    xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    ylim([0 1])
    xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
    xlim([0 length(squeeze(gini_series_all(idx_band,:  )))])
    title(['Gini values | ' bands_str{idx_band}])

    subplot(5,nbands,idx_band+3*nbands)
    plot(squeeze(s_entropy_logband_all(:,idx_band)))
    hold on;
    xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    ylim([0 1])
    xlim([1 length(squeeze(s_entropy_logband_all(:, idx_band)))])
    xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
    title(['entropy values | ' bands_str{idx_band}])

    subplot(5,nbands,idx_band+4*nbands)
    plot(squeeze(s_m1_logband_all(:,idx_band)))
    hold on;
    plot(squeeze(s_m2_logband_all(:,idx_band)))
    plot(squeeze(s_m3_logband_all(:,idx_band)))
    xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    legend([{'m1'}, {'m2'}, {'m3'}])
    ylim([0 1])
    xlim([1 length(squeeze(s_m1_logband_all(:, idx_band)))])
    xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
    title(['my metrics values | ' bands_str{idx_band}])
end
sgtitle([subject ' ' day ' | windowing | all channels']);

figure()
for idx_band = 1:nbands
    subplot(4,nbands,idx_band)
    imagesc(squeeze(wxw_img_occipital(:,:,idx_band)));
    hold on;
    xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    yticks(1:10:nwin); yticklabels(xticks_labels(1:10:end));
    xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
    title(['window x window features intersection | ' bands_str{idx_band}])

    subplot(4,nbands,idx_band+nbands)
    tmp = squeeze(w_nearity_occipital(:,idx_band));
    tmp = tmp - min(tmp); tmp = tmp / max(tmp);
    plot(tmp)
    hold on;
    xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    ylim([0 1])
    xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
    title(['window features nearity | ' bands_str{idx_band}])

    subplot(4,nbands,idx_band+2*nbands)
    plot(squeeze(gini_series_occipital(idx_band,:)))
    hold on;
    xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    ylim([0 1])
    xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
    title(['Gini values | ' bands_str{idx_band}])

    subplot(4,nbands,idx_band+3*nbands)
    plot(squeeze(s_entropy_logband_occipital(:,idx_band)))
    hold on;
    xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    ylim([0 1])
    xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
    title(['entropy values | ' bands_str{idx_band}])

end
sgtitle([subject ' ' day ' | windowing | occipital channels']);

%% take one of the metric to divide the data into shift-sustained-nothing
% use the nearity to divide the data into 2 or more clusters
K = 2; % number of states/clusters
chosen_band = 2;
all_data = squeeze(w_my_nearity_all(:,chosen_band));
[cluster_labels, C] = kmeans(all_data, K, 'Distance', 'sqeuclidean', 'Replicates', 10, 'Display', 'iter');

figure();
subplot(311)
imagesc(squeeze(diff_logBand_win(:,chosen_band,:)));
hold on;
xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
hold off;
yticks(1:nchannels); yticklabels(channels_label)
xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
title(['diff band power ' bands_str{chosen_band}])

subplot(312)
plot(cluster_labels)
hold on;
xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
hold off;
xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
title('clusters')

subplot(313)
plot(all_data)
hold on;
xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
hold off;
xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
title('nearity')

%% take the data according to the division and trian K qdas --> in this case 2
percentual = 0.75;
ntrial_train = ceil(ntrial*percentual); ntrial_test = ntrial - ntrial_train;
data_1_train = nan(sum(cluster_labels(cf_startWIN:end) == 1)*ceil(win_size)*ntrial_train, nchannels); i1_train = 1;
label_1_train = nan(sum(cluster_labels(cf_startWIN:end) == 1)*ceil(win_size)*ntrial_train, 1);
data_2_train = nan(sum(cluster_labels(cf_startWIN:end) == 2)*ceil(win_size)*ntrial_train, nchannels); i2_train = 1;
label_2_train = nan(sum(cluster_labels(cf_startWIN:end) == 2) * ceil(win_size)*ntrial_train, 1);

data_1_test = nan(sum(cluster_labels(cf_startWIN:end) == 1)*ceil(win_size)*ntrial_test, nchannels); i1_test = 1;
label_1_test = nan(sum(cluster_labels(cf_startWIN:end) == 1)*ceil(win_size)*ntrial_test, 1);
data_2_test = nan(sum(cluster_labels(cf_startWIN:end) == 2)*ceil(win_size)*ntrial_test, nchannels); i2_test = 1;
label_2_test = nan(sum(cluster_labels(cf_startWIN:end) == 2) * ceil(win_size)*ntrial_test, 1);
for idx_w = cf_startWIN:nwin
    start_w = ceil((idx_w-1)*overlap + 1);
    end_w = ceil(start_w + win_size - 1);
    if end_w > min_trial_data
        end_w = min_trial_data;
    end

    c_data = squeeze(trial_data(start_w:end_w, chosen_band,:,:));

    for t = 1:ntrial
        if t < ntrial_train
            if cluster_labels(idx_w) == 1
                data_1_train((i1_train-1)*ceil(win_size)+1:i1_train*ceil(win_size),:) = squeeze(c_data(:,:,t));
                label_1_train((i1_train-1)*ceil(win_size)+1:i1_train*ceil(win_size)) = repmat(trial_typ(t), 1, ceil(win_size));
                i1_train = i1_train + 1;
            elseif cluster_labels(idx_w) == 2
                data_2_train((i2_train-1)*ceil(win_size)+1:i2_train*ceil(win_size),:) = squeeze(c_data(:,:,t));
                label_2_train((i2_train-1)*ceil(win_size)+1:i2_train*ceil(win_size)) = repmat(trial_typ(t), 1, ceil(win_size));
                i2_train = i2_train + 1;
            end
        else
            if cluster_labels(idx_w) == 1
                data_1_test((i1_test-1)*ceil(win_size)+1:i1_test*ceil(win_size),:) = squeeze(c_data(:,:,t));
                label_1_test((i1_test-1)*ceil(win_size)+1:i1_test*ceil(win_size)) = repmat(trial_typ(t), 1, ceil(win_size));
                i1_test = i1_test + 1;
            elseif cluster_labels(idx_w) == 2
                data_2_test((i2_test-1)*ceil(win_size)+1:i2_test*ceil(win_size),:) = squeeze(c_data(:,:,t));
                label_2_test((i2_test-1)*ceil(win_size)+1:i2_test*ceil(win_size)) = repmat(trial_typ(t), 1, ceil(win_size));
                i2_test = i2_test + 1;
            end
        end
    end
end

%% compute the fisher score on th previous division
fisher_1 = nan(1, nchannels);
fisher_2 = nan(1, nchannels);
for idx_ch=1:nchannels
    % shift
    mu1 = mean(data_1_train(label_1_train == classes(1),idx_ch));
    sigma1 = std(data_1_train(label_1_train == classes(1),idx_ch));
    mu2 = mean(data_1_train(label_1_train == classes(2),idx_ch));
    sigma2 = std(data_1_train(label_1_train == classes(2),idx_ch));
    fisher_1(idx_ch) = abs(mu1 - mu2);

    % sus
    mu1 = mean(data_2_train(label_2_train == classes(1),idx_ch));
    sigma1 = std(data_2_train(label_2_train == classes(1),idx_ch));
    mu2 = mean(data_2_train(label_2_train == classes(2),idx_ch));
    sigma2 = std(data_2_train(label_2_train == classes(2),idx_ch));
    fisher_2(idx_ch) = abs(mu1 - mu2);
end

tmp = [fisher_1; fisher_2];
figure();
imagesc(tmp')
yticks(1:nchannels); yticklabels(channels_label)
xticks(1:3); xticklabels([{"shift"}, {"sustained"}])

% select the features
idx_features_1 = 1:nchannels; 
idx_features_1(ismember(idx_features_1, [1, 2, 19])) = [];

idx_features_2 = 1:nchannels; 
idx_features_2(ismember(idx_features_2, [1, 2, 19])) = [];

%% train the models
qdaModel_1 = fitcdiscr(data_1_train(:,idx_features_1), label_1_train, 'DiscrimType', 'quadratic');
CVModel = crossval(qdaModel_1, 'KFold', 5);
loss = kfoldLoss(CVModel);  % Average classification error
disp('--- class 1 ---')
disp('TRAIN')
fprintf("5-fold cross-validation loss: %.4f\n", loss);
p = predict(qdaModel_1, data_1_train(:,idx_features_1));
confusionmat(label_1_train, p)
disp(['Accuracy: ' num2str(sum(p == label_1_train)/size(p,1))])
disp('TEST')
p = predict(qdaModel_1, data_1_test(:,idx_features_1));
confusionmat(label_1_test, p)
disp(['Accuracy: ' num2str(sum(p == label_1_test)/size(p,1))])

disp('--- class 2 ---')
disp('TRAIN')
qdaModel_2 = fitcdiscr(data_2_train(:,idx_features_2), label_2_train, 'DiscrimType', 'quadratic');
CVModel = crossval(qdaModel_2, 'KFold', 5);
loss = kfoldLoss(CVModel);  % Average classification error
fprintf("5-fold cross-validation loss: %.4f\n", loss);
p = predict(qdaModel_2, data_2_train(:,idx_features_2));
confusionmat(label_2_train, p)
disp(['Accuracy: ' num2str(sum(p == label_2_train)/size(p,1))])
disp('TEST')
p = predict(qdaModel_2, data_2_test(:,idx_features_2));
confusionmat(label_2_test, p)
disp(['Accuracy: ' num2str(sum(p == label_2_test)/size(p,1))])

%% pseudo online test
for t = ntrial_train:ntrial
    c_data = squeeze(trial_data(:,chosen_band,:,t));

    [l,score1,~] = predict(qdaModel_1, c_data(:, idx_features_1));
    [~,score2,~] = predict(qdaModel_2, c_data(:, idx_features_2));

    figure();
    subplot(311)
    imagesc(c_data')
    hold on;
    xline(minDurCue, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(minDurCue + minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    yticks(1:nchannels); yticklabels(channels_label)
    xticks(sampleRate:sampleRate:size(c_data, 1))
    xticklabels(string((sampleRate:sampleRate:size(c_data, 1)) / sampleRate));
    title('band power')

    subplot(312)
    hold on;
    plot(score1(:,1))
    plot(score2(:,1))
    xline(minDurCue, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(minDurCue + minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    xlim([0 length(score1(:,1))])
    legend({'classe 1', 'classe 2'})
    title('probabilities')

    sgtitle(['class asked: ' num2str(trial_typ(t))])
end






%%
function [matrix_distance, sparsity_distance, sparsity_roi] = compute_distance(selectedChans, LblMontage, channels_label)
% --- compute sparsity distance ---
positions = zeros(length(selectedChans), 2);
for i = 1:length(selectedChans)
    [r, c] = find(strcmp(LblMontage, selectedChans{i}));
    if isempty(r)
        error(['Channel not found: ' selectedChans{i}]);
    end
    positions(i, :) = [r, c];
end

matrix_distance = pdist2(positions, positions, 'euclidean');
dist_vec = matrix_distance(triu(true(size(matrix_distance)), 1));
meanDist = mean(dist_vec);
medianDist = median(dist_vec);
[allRows, allCols] = find(~cellfun('isempty', LblMontage));
all_positions = [allRows, allCols];
maxDist = max(pdist2(all_positions, all_positions, 'euclidean'), [], 'all');
alpha = (meanDist/maxDist) / ((meanDist/maxDist) + (medianDist/maxDist));
normalizedDist = alpha*(meanDist/maxDist) + (1-alpha)*(medianDist/maxDist); % normalization
sparsity = 1 - normalizedDist; % 1 is more clustered, 0 more sparse
sparsity_distance = max(0, min(1, sparsity));

% --- sparsity roi ---
frontal = [1 2 3 4 5 20 21]; central = [6 7 8 9 10 11 12 22 23 24 25 26 27 28]; occipital = [13 14 15 16 17 18 29 30 31 32 33 34 35 36 37 38 39];
c_l = [6 22 25 8 11 27]; c_r = [7 24 26 10 12 28]; o_l = [29 13 30 37 33 34 17]; o_r = [31 15 32 35 36 38 18];
[~, idx_selected] = ismember(selectedChans, channels_label);

nFrontal = sum(ismember(idx_selected, frontal));
nC_l = sum(ismember(idx_selected, c_l));
nC_r = sum(ismember(idx_selected, c_r));
nO_l = sum(ismember(idx_selected, o_l));
nO_r = sum(ismember(idx_selected, o_r));
nCentral = sum(ismember(idx_selected, central));
nOccipital = sum(ismember(idx_selected, occipital));

tmp = [0 nFrontal/length(selectedChans) 0; 0 0 0; nC_l/length(selectedChans) 0 nC_r/length(selectedChans); 0 0 0; nO_l/length(selectedChans) 0 nO_r/length(selectedChans)];
tmp = tmp(:);
n = length(tmp);

% compute the first sparsity
l1 = norm(tmp, 1);      % Norma L1
l2 = norm(tmp, 2);      % Norma L2
if l2 == 0
    s_1 = 0;  % caso limite: vettore nullo
else
    s_1 = (sqrt(n) - (l1 / l2)) / (sqrt(n) - 1);
end

% compute the second sparsity
tmp = [nFrontal/length(selectedChans); nCentral/length(selectedChans); nOccipital/length(selectedChans)];
n = length(tmp);
l1 = norm(tmp, 1);      % Norma L1
l2 = norm(tmp, 2);      % Norma L2
if l2 == 0
    s_2 = 0;  % caso limite: vettore nullo
else
    s_2 = (sqrt(n) - (l1 / l2)) / (sqrt(n) - 1);
end

sparsity_roi = sqrt(s_1 * s_2);

end

%% sparsity and where
function [sparsityEntropy, m1, m2,m3] = compute_sparsity(signalVec)
    % Ensure column vector
    N = length(signalVec);
    % take the one which contribute at the 95% of the energy
    A2 = signalVec.^2;
    [A2_sorted, ~] = sort(A2, 'descend');
    cumulative = cumsum(A2_sorted); % Compute cumulative energy
    total = sum(A2_sorted);
    % Find how many components cover 95%
    threshold_index = find(cumulative >= 0.95 * total, 1, 'first');
    % Get the minimum value (in magnitude) that contributes to 95%
    value_cutoff = sqrt(A2_sorted(threshold_index));
    
    try
        activeIdx = find(signalVec >= value_cutoff);
    catch ME
        disp(ME.message);
    end
    % --- Entropy-Based Sparsity --- --> in general
    signalVec_entropy = signalVec(activeIdx).^2;
    total_entropy = sum(signalVec_entropy);
    if total_entropy == 0
        sparsityEntropy = NaN;
    else
        p = signalVec_entropy / total_entropy;
        H = -sum(p .* log2(p + eps));  % entropy
        Hmax = log2(N);
        sparsityEntropy = 1 - (H / Hmax);  % normalized: 0 to 1
    end
   
    % --- metodo 1 ---
    frontal = [1 2 3 4 5 20 21]; central = [6 7 8 9 10 11 12 22 23 24 25 26 27 28]; occipital = [13 14 15 16 17 18 29 30 31 32 33 34 35 36 37 38 39];
    c_l = [6 22 25 8 11 27]; c_r = [7 24 26 10 12 28]; o_l = [29 13 30 37 33 34 17]; o_r = [31 15 32 35 36 38 18];
    nFrontal = sum(ismember(activeIdx, frontal)); f_mean = mean(signalVec(ismember(activeIdx, frontal)));
    nCentral = sum(ismember(activeIdx, central)); c_mean = mean(signalVec(ismember(activeIdx, central)));
    nOccipital = sum(ismember(activeIdx, occipital)); o_mean = mean(signalVec(ismember(activeIdx, occipital)));
    total_mean = f_mean + c_mean + o_mean;
    n_total = nFrontal + nCentral + nOccipital;
    m1 = nFrontal/(nCentral + nOccipital) * f_mean/total_mean + nCentral/(nFrontal + nOccipital) * c_mean/total_mean + nOccipital/(nFrontal + nCentral) * o_mean/total_mean; 
    try
    m1 = m1 * (nchoosek(length(frontal), nFrontal)*nchoosek(length(central), nCentral) * nchoosek(length(occipital), nOccipital)) / nchoosek(length(signalVec), n_total);
    catch ME
        disp(ME.message)
    end

    % --- metodo 2 ---
    mean_m2 = [mean(signalVec(activeIdx(ismember(activeIdx,frontal)))), mean(signalVec(activeIdx(ismember(activeIdx,central)))), mean(signalVec(activeIdx(ismember(activeIdx,occipital))))];
    mean_m2(isnan(mean_m2)) = 0;
    n = length(mean_m2);
    l1 = norm(mean_m2, 1);      % Norma L1
    l2 = norm(mean_m2, 2);      % Norma L2
    if l2 == 0
        m2 = 0;  % caso limite: vettore nullo
    else
        m2 = (sqrt(n) - (l1 / l2)) / (sqrt(n) - 1);
    end

%     signalVec_entropy = mean_m2.^2;
%     total_entropy = sum(signalVec_entropy);
%     if total_entropy == 0
%         m2 = NaN;
%     else
%         p = signalVec_entropy / total_entropy;
%         H = -sum(p .* log2(p + eps));  % entropy
%         Hmax = log2(n);
%         m2 = 1 - (H / Hmax);  % normalized: 0 to 1
%     end

    % --- metodo 3 ---
    mean_m3 = [mean(signalVec(activeIdx(ismember(activeIdx,frontal)))), mean(signalVec(activeIdx(ismember(activeIdx,c_l)))), ...
        mean(signalVec(activeIdx(ismember(activeIdx,c_r)))), mean(signalVec(activeIdx(ismember(activeIdx,o_l)))), mean(signalVec(activeIdx(ismember(activeIdx,o_r))))];
    mean_m3(isnan(mean_m3)) = 0;
    n = length(mean_m3);
    l1 = norm(mean_m3, 1);      % Norma L1
    l2 = norm(mean_m3, 2);      % Norma L2
    if l2 == 0
        m3 = 0;  % caso limite: vettore nullo
    else
        m3 = (sqrt(n) - (l1 / l2)) / (sqrt(n) - 1);
    end

%     signalVec_entropy = mean_m3.^2;
%     total_entropy = sum(signalVec_entropy);
%     if total_entropy == 0
%         m3 = NaN;
%     else
%         p = signalVec_entropy / total_entropy;
%         H = -sum(p .* log2(p + eps));  % entropy
%         Hmax = log2(n);
%         m3 = 1 - (H / Hmax);  % normalized: 0 to 1
%     end
end

%%
function g = gini(x)
    x = abs(x(:));   % vectorize and take abs (in case)
    if sum(x) == 0
        g = 0;
        return;
    end
    x = sort(x);
    n = length(x);
    g = (2 * sum((1:n)' .* x)) / (n * sum(x)) - (n + 1) / n;
end

