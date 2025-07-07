%% show using the same windows: the wxw with 10 features, the log band, the phase
% use a windows of 100ms with overlap of 50ms, see this only for channels P, PO, O during all the trial.
% In the plot, report the start of the fixation and the cue.
% use the laplacian

clear all; % close all;

addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/EOG')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')

%% Initialization
a = 4:2:18;
b = a+2;
c = [a; b];
bands = cell(1, size(a, 2) + 1);
for i=1:length(a)
    bands{i} = [a(i), b(i)];
end
bands{i + 1} = [8 14];
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
nbands = length(bands);
signals = cell(2, nbands);
headers = cell(2, nbands);
for i = 1:2
    for idx_band = 1:nbands
        headers{i,idx_band}.TYP = [];
        headers{i,idx_band}.DUR = [];
        headers{i,idx_band}.POS = [];
        signals{i,idx_band} = [];
    end
end
erp_band = [2 10];
signals_erp = [];
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

    c_trial_with_eog = eog_detection(c_signal, header, eog_threshold, {'FP1', 'FP2', 'EOG'});
    trial_with_eog = [trial_with_eog; c_trial_with_eog];

%     % laplacian
%     disp('      [proc] applying laplacian')
%     c_signal = c_signal * lap;

    for idx_band = 1:nbands
        band = bands{idx_band};

        % for phase extraction
        disp('   [proc] Phase');
        disp('      [proc] applying filtering')
        [b, a] = butter(filterOrder, band(2)*(2/header.SampleRate),'low');
        s_low = filter(b,a,c_signal);
        [b, a] = butter(filterOrder, band(1)*(2/header.SampleRate),'high');
        s_filt = filter(b,a,s_low);
        signal_processed = s_filt;
        c_header = headers{1, idx_band};
        c_header.sampleRate = header.SampleRate;
        c_header.channels_labels = header.Label;
        if isempty(find(header.EVENT.TYP == 2, 1)) % no eye calibration
            c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP);
            c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR);
            c_header.POS = cat(1, c_header.POS, header.EVENT.POS + size(signals{1, idx_band}, 1));
        else
            k = find(header.EVENT.TYP == 1, 1);
            c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP(k:end));
            c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR(k:end));
            c_header.POS = cat(1, c_header.POS, header.EVENT.POS(k:end) + size(signals{1, idx_band}, 1));
        end
        signals{1, idx_band} = cat(1, signals{1, idx_band}, signal_processed(:,:));
        headers{1, idx_band} = c_header;

        % for log band extraction
        disp('   [proc] logband');

        signal_processed = proc_512hz(c_signal, header.SampleRate, band, filterOrder, avg);
        c_header = headers{2, idx_band};
        c_header.sampleRate = header.SampleRate;
        c_header.channels_labels = header.Label;
        if isempty(find(header.EVENT.TYP == 2, 1)) % no eye calibration
            c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP);
            c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR);
            c_header.POS = cat(1, c_header.POS, header.EVENT.POS + size(signals{2, idx_band}, 1));
        else
            k = find(header.EVENT.TYP == 1, 1);
            c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP(k:end));
            c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR(k:end));
            c_header.POS = cat(1, c_header.POS, header.EVENT.POS(k:end) + size(signals{2, idx_band}, 1));
        end
        signals{2, idx_band} = cat(1, signals{2, idx_band}, signal_processed(:,:));
        headers{2, idx_band} = c_header;
    end

    % for erp
    disp('   [proc] ERP');
    disp('      [proc] applying filtering')
    [b, a] = butter(filterOrder, erp_band(2)*(2/header.SampleRate),'low');
    s_low = filter(b,a,c_signal);
    [b, a] = butter(filterOrder, erp_band(1)*(2/header.SampleRate),'high');
    s_filt = filter(b,a,s_low);
    signal_processed = s_filt;
    
    signals_erp = cat(1, signals_erp, signal_processed(:,:));

end


%% Labelling data 
events = headers{1, 1};
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
trial_data4phase = nan(min_trial_data, nbands, nchannels, ntrial);
trial_data4logBand = nan(min_trial_data, nbands, nchannels, ntrial);
for idx_band = 1:nbands
    c_signal_phase = signals{1, idx_band};
    c_signal_logBand = signals{2, idx_band};
    for trial = 1:ntrial
        c_start = trial_start(trial);
        c_end = trial_start(trial) + min_trial_data - 1;
        trial_data4phase(:,idx_band,:,trial) = c_signal_phase(c_start:c_end,:);
        trial_data4logBand(:,idx_band,:,trial) = c_signal_logBand(c_start:c_end,:);
    end
end

trial_data4erp = nan(min_trial_data, nchannels, ntrial);
for trial = 1:ntrial
    c_start = trial_start(trial);
    c_end = trial_start(trial) + min_trial_data - 1;
    trial_data4erp(:,:,trial) = signals_erp(c_start:c_end,:);
end

%% refactoring the data
% now the data are placed alternatevely, so odd trial id is for 730, even for 731
balanced_trial_idx = balanced_data_afterEOG(trial_with_eog, trial_typ, classes);
balanced_trial_data_phase = trial_data4phase(:,:,:,logical(balanced_trial_idx));
balanced_trial_data_log = trial_data4logBand(:,:,:,logical(balanced_trial_idx));
balanced_trial_data_erp = trial_data4erp(:,:,logical(balanced_trial_idx));
trial_typ = trial_typ(logical(balanced_trial_idx));
ntrial = sum(balanced_trial_idx);
idx_classes_trial = nan(ntrial/2, nclasses);
for idx_class = 1:nclasses
    idx_classes_trial(:,idx_class) = find(trial_typ == classes(idx_class));
end

tmp_phase = nan(size(balanced_trial_data_phase));
tmp_log = nan(size(balanced_trial_data_log));
tmp_erp = nan(size(balanced_trial_data_erp));
trial_typ = nan(size(trial_typ));
i = 1;
for idx_trial_class = 1:2:ntrial
    for idx_class = 1:nclasses
        tmp_phase(:,:,:,idx_trial_class + idx_class - 1) = balanced_trial_data_phase(:,:,:,idx_classes_trial(i, idx_class));
        tmp_log(:,:,:,idx_trial_class + idx_class - 1) = balanced_trial_data_log(:,:,:,idx_classes_trial(i, idx_class));
        tmp_erp(:,:,idx_trial_class + idx_class - 1) = balanced_trial_data_erp(:,:,idx_classes_trial(i, idx_class));
        trial_typ(idx_trial_class + idx_class - 1) = classes(idx_class);
    end
    i = i + 1;
end
trial_data4phase = tmp_phase;
trial_data4logBand = tmp_log;
trial_data4erp = tmp_erp;

%% hilbert
hilbert_data = nan(min_trial_data, nbands, nchannels, ntrial);  % complex analytic signal

for b = 1:nbands
    for ch = 1:nchannels
        for tr = 1:ntrial
            % Extract signal vector for that band, channel, trial
            signal = trial_data4phase(:, b, ch, tr);
            
            % Compute analytic signal via Hilbert transform
            c_h = hilbert(signal);
            
            % Store result
            hilbert_data(:, b, ch, tr) = c_h;
        end
    end
end

phase_data = angle(hilbert_data); % samples x band x channels x trial
amplitude_data = abs(hilbert_data);
power_data = amplitude_data .^ 2;

%% windowing reasoning
win_size = 100; %ms
win_size = win_size * sampleRate /1000;
overlap = 50; %ms
overlap = overlap * sampleRate / 1000;
nwin = round((min_trial_data - win_size)/overlap) + 1;

% fisher score windows and the log band mean
fishers_windows = zeros(nchannels, nbands, nwin); % channels x bands x nwindows
logBand_win = zeros(nchannels, nbands, nwin, nclasses); % channels x bands x nwindows
windows_label = cell(1, nwin);
windows_center = zeros(1, nwin);
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
            mean_c1(idx_ch, idx_band) = squeeze(mean(mean(trial_data4logBand(start_w:end_w,idx_band,idx_ch,trial_typ == classes(1)), 4), 1));
            mean_c2(idx_ch, idx_band) = squeeze(mean(mean(trial_data4logBand(start_w:end_w,idx_band,idx_ch,trial_typ == classes(2)), 4), 1));
            std_c1(idx_ch, idx_band) = squeeze(mean(std(trial_data4logBand(start_w:end_w,idx_band,idx_ch,trial_typ == classes(1)), 0, 4), 1));
            std_c2(idx_ch, idx_band) = squeeze(mean(std(trial_data4logBand(start_w:end_w,idx_band,idx_ch,trial_typ == classes(2)), 0, 4), 1));
        end
    end

    windows_center(idx_w) = round((start_w-1 + end_w) / 2 / sampleRate,2);
    fishers_windows(:,:,idx_w) = abs(mean_c1 - mean_c2) ./ sqrt(std_c1.^2 + std_c2.^2);

    logBand_win(:,:,idx_w, 1) = mean_c1;
    logBand_win(:,:,idx_w, 2) = mean_c2;
end
diff_logBand_win = abs(logBand_win(:,:,:,1) - logBand_win(:,:,:,2)); % diff bl - br
diff_logBand_sample = mean(trial_data4logBand(:,:,:,trial_typ == classes(1)), 4) - mean(trial_data4logBand(:,:,:,trial_typ == classes(2)), 4);

% mean itpc
itpc_win = nan(nwin, nbands, nchannels);
itpc_sample = nan(min_trial_data, nbands, nchannels);

for b = 1:nbands
    for ch = 1:nchannels
        for idx_w = 1:nwin
            start_w = ceil((idx_w-1)*overlap + 1);
            end_w = ceil(start_w + win_size - 1);
            if end_w > min_trial_data
                end_w = min_trial_data;
            end

            tmp_phase = squeeze(phase_data(start_w:end_w,b, ch, :));
            tmp_itpc = zeros(1, size(tmp_phase, 1));

            for sample = 1:size(tmp_phase, 1)
                c_phase = squeeze(tmp_phase(sample, :)); % singal at each trial
                tmp_itpc(sample) = abs(mean(exp(1i * c_phase)));
            end
            itpc_win(idx_w,b,ch) = mean(tmp_itpc);
        end

        for sample = 1:min_trial_data
            c_phase = squeeze(phase_data(sample, b, ch,:)); % singal at each trial
            itpc_sample(sample, b, ch) = abs(mean(exp(1i * c_phase)));
        end
    end
end

% mean erp
erp_win = zeros(nwin, nchannels);
erp_sample = mean(trial_data4erp, 3);

for idx_w = 1:nwin
    start_w = ceil((idx_w-1)*overlap + 1);
    end_w = ceil(start_w + win_size - 1);
    if end_w > min_trial_data
        end_w = min_trial_data;
    end

    tmp_erp = squeeze(mean(mean(trial_data4erp(start_w:end_w,:,:), 3), 1));
    erp_win(idx_w, :) = tmp_erp;
end

%% compute the wxw
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
% channels_select = {'POZ', 'O1', 'O2', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
[~, channelsSelected] = ismember(channels_select, channels_label);
nchannelsSelected = size(channelsSelected, 2);
nfeatures = 27; % number features of fisher to select
wxw_img = zeros(nwin, nwin);
for idx_w1 = 1:nwin
    c1_fisher = squeeze(fishers_windows(channelsSelected,:,idx_w1));
    [c1_sortedValues, c1_linearIndices] = maxk(c1_fisher(:), nfeatures);
    for idx_w2 = 1:nwin
        c2_fisher = squeeze(fishers_windows(channelsSelected,:,idx_w2));
        [c2_sortedValues, c2_linearIndices] = maxk(c2_fisher(:), nfeatures);

        n_wxw_intersection = length(intersect(c1_linearIndices, c2_linearIndices));
        n_wxw_union = length(union(c1_linearIndices, c2_linearIndices));
        wxw_img(idx_w1,idx_w2) = n_wxw_intersection / nfeatures;
    end
end

%% plot wxw, log band diff and phase in windows
nfigure = 3;
idx_b = 1;
cue_startWIN = round((min_durFIX - win_size)/overlap) + 1;
cf_startWIN  = round((min_durCUE + min_durFIX - win_size)/overlap) + 1;
xticks_labels = arrayfun(@(x) sprintf('%.2f', x), windows_center, 'UniformOutput', false);
% handles_itpc = []; handles_logBand = [];
% cl_itpc = -inf; cl_logBand = -inf;
for idx_f = 1:nfigure
    figure()
    idx_plot = 1;
    subplot(nbands/nfigure+1,2, idx_plot)
    imagesc(wxw_img);
    hold on;
    xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    yticks(1:10:nwin); yticklabels(xticks_labels(1:10:end));
    xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
    title('window x window features intersection')
    idx_plot = idx_plot + 1;

    subplot(nbands/nfigure+1,2, idx_plot)
    plot(erp_win);
    hold on;
    xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    plot(mean(erp_win, 2), 'color', 'blue', 'LineWidth',3)
    hold off;
    yticks(1:10:nwin); yticklabels(xticks_labels(1:10:end));
    xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
    xlim([1, nwin]);
    title(['ERP | band: ' num2str(erp_band(1)) '-' num2str(erp_band(2))])
    idx_plot = idx_plot + 1;

    for i = 1:nbands/nfigure
        subplot(nbands/nfigure+1,2, idx_plot)
        imagesc(squeeze(itpc_win(:,idx_b,channelsSelected))');
        hold on;
        xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        yticks(1:nchannelsSelected); yticklabels(channels_label(channelsSelected));
        xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
        title(['itpc | band :' bands_str{idx_b}])
        idx_plot = idx_plot + 1;
%         cl_itpc = max(cl_itpc, max(abs(squeeze(itpc_win(cf_startWIN+ceil(win_size/overlap) + 1:end,idx_b,channelsSelected))), [], 'all'));
%         handles_itpc = [handles_itpc; gca];

        subplot(nbands/nfigure+1,2, idx_plot)
        imagesc(squeeze(diff_logBand_win(channelsSelected,idx_b,:)));
        hold on;
        xline(cue_startWIN, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(cf_startWIN, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        yticks(1:nchannelsSelected); yticklabels(channels_label(channelsSelected));
        xticks(1:10:nwin); xticklabels(xticks_labels(1:10:end));
        title(['bl - br log band | band: ' bands_str{idx_b}])
        idx_plot = idx_plot + 1;
%         cl_logBand = max(cl_logBand, max(abs(squeeze(diff_logBand_win(channelsSelected,idx_b,cf_startWIN:end))), [], 'all'));
%         handles_logBand = [handles_logBand; gca];

        idx_b = idx_b + 1;
    end

    sgtitle([subject ' ' day ' | windowing reasoning | figure #' num2str(idx_f)]);
end
% set(handles_itpc, 'clim', [0, cl_itpc])
% set(handles_logBand, 'clim', [0, cl_logBand])

%% plot the itpc and log band diff in time
idx_band = 1;
nfigure = 2;
% handles_itpc = []; handles_logBand = [];
% cl_itpc = -inf; cl_logBand = -inf;
for idx_f = 1:nfigure
    figure();
    idx_plot = 1;
    for i = 1:ceil(nbands)/nfigure
        subplot(ceil(nbands/nfigure), 2, idx_plot)
        imagesc(squeeze(itpc_sample(:,idx_band,channelsSelected))')
        hold on;
        xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        yticks(1:nchannelsSelected); yticklabels(channels_label(channelsSelected));
        xticks(sampleRate:sampleRate:size(itpc_sample, 1))
        xticklabels(string((sampleRate:sampleRate:size(itpc_sample, 1)) / sampleRate));
        title(['itpc | band :' bands_str{idx_band}])
        idx_plot = idx_plot + 1;
%         cl_itpc = max(cl_itpc, max(abs(squeeze(itpc_sample(cf_startWIN+ceil(win_size/overlap) + 1:end,idx_b,channelsSelected))), [], 'all'));
%         handles_itpc = [handles_itpc; gca];


        subplot(ceil(nbands/nfigure), 2, idx_plot)
        imagesc(squeeze(diff_logBand_sample(:,idx_band,channelsSelected))')
        hold on;
        xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        yticks(1:nchannelsSelected); yticklabels(channels_label(channelsSelected));
        xticks(sampleRate:sampleRate:size(itpc_sample, 1))
        xticklabels(string((sampleRate:sampleRate:size(itpc_sample, 1)) / sampleRate));
        title(['diff bl - br | band :' bands_str{idx_band}])
        idx_plot = idx_plot + 1;
%         cl_logBand = max(cl_logBand, max(abs(squeeze(diff_logBand_sample(channelsSelected,idx_b,cf_startWIN:end))), [], 'all'));
%         handles_logBand = [handles_logBand; gca];

        idx_band = idx_band + 1;
    end
    subplot(ceil(nbands/nfigure), 2, idx_plot)
    plot(erp_sample)
    hold on;
    plot(mean(erp_sample, 2), 'color', 'red', 'LineWidth',3)
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    xticks(sampleRate:sampleRate:size(itpc_sample, 1))
    xticklabels(string((sampleRate:sampleRate:size(itpc_sample, 1)) / sampleRate));
    xlim([1, min_trial_data]);
    title(['ERP | band: ' num2str(erp_band(1)) '-' num2str(erp_band(2))])
    idx_plot = idx_plot + 1;

    subplot(ceil(nbands/nfigure), 2, idx_plot)
    plot(erp_sample)
    hold on;
    plot(mean(erp_sample, 2), 'color', 'red', 'LineWidth',3)
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    xticks(sampleRate:sampleRate:size(itpc_sample, 1))
    xticklabels(string((sampleRate:sampleRate:size(itpc_sample, 1)) / sampleRate));
    xlim([1, min_trial_data]);
    title(['ERP | band: ' num2str(erp_band(1)) '-' num2str(erp_band(2))])

    sgtitle([subject ' ' day ' | time reasoning | figure #' num2str(idx_f)]);
end
% set(handles_itpc, 'clim', [0, cl_itpc])
% set(handles_logBand, 'clim', [0, cl_logBand])

