%% one class classifier using phase values and ROI in the channels
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
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
nbands = length(bands);
% channels_select = {'O1', 'O2', 'OZ'};
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
% roi = {{'F3', 'F1', 'Fz', 'F2', 'F4'}; {'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4'}; 
%     {'P5', 'P3', 'P1', 'PO7', 'PO5', 'PO3', 'O1'}; {'P6', 'P4', 'P2', 'PO8', 'PO6', 'PO4', 'O2'}};
% roi_label = {{'F channels'}, {'FC, C and CP channels'}, {'P, PO and O channels right'}, {'P, PO and O channels right'}};
% roi = {{'F3', 'F1', 'Fz', 'F2', 'F4'}; {'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4'}; 
%     {'P1', 'PO3', 'O1', 'P2', 'PO4', 'O2', 'Oz', 'POz', 'Pz'}};
% roi_label = {{'F channels'}, {'FC, C and CP channels'}, {'P, PO and O channels'}};
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

    c_trial_with_eog = eog_detection(c_signal, header, eog_threshold, {'FP1', 'FP2', 'EOG'});
    trial_with_eog = [trial_with_eog; c_trial_with_eog];

    for idx_band = 1:nbands
        band = bands{idx_band};

        % filter alpha band
        disp('      [proc] applying filtering')
        [b, a] = butter(filterOrder, band(2)*(2/header.SampleRate),'low');
        s_low = filter(b,a,c_signal);
        [b, a] = butter(filterOrder, band(1)*(2/header.SampleRate),'high');
        s_filt = filter(b,a,s_low);

%         disp('      [proc] applying laplacian');
%         s_filt = s_filt * lap;
        s_lap = s_filt;

        signal_processed = s_lap; % filtered signal --> all channels
        signal_roi_processed = nan(size(s_lap, 1), nroi); % for roi
        for i = 1:nroi
            signal_roi_processed(:,i) = mean(signal_processed(:, roi_indices{i}), 2);
        end


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

        signals{idx_band} = cat(1, signals{idx_band}, signal_roi_processed(:,:)); % for all channels
        headers{idx_band} = c_header;
    end
end
nchannels = size(signals{1}, 2);

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

%% hilbert
hilbert_data = nan(min_trial_data, nbands, nchannels, ntrial);  % complex analytic signal

for b = 1:nbands
    for ch = 1:nchannels
        for tr = 1:ntrial
            % Extract signal vector for that band, channel, trial
            signal = trial_data(:, b, ch, tr);
            
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

%% compute the ITPC Inter-Trial Phase Coherence
itcp_1 = nan(min_trial_data, nbands, nchannels);
itcp_2 = nan(min_trial_data, nbands, nchannels);
itcp = nan(min_trial_data, nbands, nchannels);

for b = 1:nbands
    for ch = 1:nchannels
        for sample = 1:min_trial_data
            c_phase_1 = squeeze(phase_data(sample, b, ch,trial_typ == classes(1))); % signal at each trial
            itcp_1(sample, b, ch) = abs(mean(exp(1i * c_phase_1)));

            c_phase_2 = squeeze(phase_data(sample, b, ch,trial_typ == classes(2))); % signal at each trial
            itcp_2(sample, b, ch) = abs(mean(exp(1i * c_phase_2))); % mean over trial of that signal value

            c_phase = squeeze(phase_data(sample,b,ch,:));
            itcp(sample, b, ch) = abs(mean(exp(1i * c_phase)));
        end
    end
end

%% divide the data for the two classes and see the ITPC value
figure();
idx_plot = 1;
handle = []; c_max = -inf;

for idx_band = 1:nbands
    subplot(nbands,2,idx_plot)
    hold on
    for idx_roi = 1:nroi
        plot(itcp(:,idx_band, idx_roi))
    end
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off
    legend(cellfun(@(c) c{1}, roi_label, 'UniformOutput', false))
    xticks(sampleRate:sampleRate:size(itcp_2, 1))
    xticklabels(string((sampleRate:sampleRate:size(itcp_2, 1)) / sampleRate));
    title(['band: ' bands_str{idx_band} ' | without trial separation'])
    idx_plot = idx_plot + 1;
    ylim([0 1])

    subplot(nbands,2,idx_plot)
    imagesc(squeeze(itcp(:,idx_band, :))')
    hold on
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durFIX + min_durCUE, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off
    yticks(1:nroi)
    yticklabels(cellfun(@(c) c{1}, roi_label, 'UniformOutput', false))
    xticks(sampleRate:sampleRate:size(itcp_2, 1))
    xticklabels(string((sampleRate:sampleRate:size(itcp_2, 1)) / sampleRate));
    handle = [handle, gca];
    c_max = max(c_max, max(abs(squeeze(itcp(:,idx_band, :))'), [], 'all'));
    title(['band: ' bands_str{idx_band} ' | without trial separation'])
    idx_plot = idx_plot + 1;
    set(handle, 'CLim', [-c_max, c_max]);
    handle = [];
    c_max = -inf;

end
sgtitle('Inter-Trial Phase Coherence (ITPC)')

%% window the phase and compute circular variance
window_size = 50; %ms
step_size = 25; %ms
wsize = window_size / 1000 * sampleRate;
ssize = step_size / 1000 * sampleRate;

num_windows = floor((min_trial_data - wsize) / ssize) + 1;

% circular variance
circ_var_windows = nan(num_windows, nbands, nchannels, ntrial); % compute the circular variance in that window
for b = 1:nbands
    for t = 1:ntrial
        for win = 1 : num_windows
            start_idx = round((win - 1) * ssize + 1);
            end_idx = round(start_idx + wsize - 1);
            c_phase = squeeze(phase_data(start_idx:end_idx,b,:,t));

            circ_var_windows(win,b,:,t) = 1-abs(mean(exp(1i * c_phase), 1)); % mean on the window
        end
    end
end

%% show for some channels how the circular variance is across the trials
for b = 1:nbands
    figure();
    idx_plot = 1; cl = -inf; handles = [];
    for ch = 1:nroi
        subplot(nroi, 1, idx_plot)
        ch_signal = nan(num_windows, ntrial);
        for t = 1:ntrial
            ch_signal(:,t) = circ_var_windows(:, b, ch, t);
        end
        imagesc(ch_signal');
        cl = max(cl, abs(max(ch_signal, [], 'all')));
        handles = [handles; gca];

        title(['roi: ' roi_label{ch}])
        idx_plot = idx_plot + 1;
    end
    cl = max(cl, max(abs(ch_signal), [], 'all'));
    set(handles, 'clim', [0 cl])
    sgtitle(['circular variance | band: ' bands_str{b}])
end

%% phase slope
phase_slope = nan(min_trial_data, nbands, nchannels, ntrial);

for b = 1:nbands
    for t = 1:ntrial
        for ch = 1:nchannels
            for sample = 1:min_trial_data
                if sample <= 5
                    prev_value = 0;
                else
                    prev_value = phase_data(sample - 5, b, ch,t);
                end
                current_value = phase_data(sample, b, ch, t);
                phase_slope(sample, b, ch, t) = current_value - prev_value;
            end
        end
    end
end

%% show
nrows = 4; ncols = 4;
for b = 1:nbands
    figure();
    idx_plot = 1; cl = -inf; handles = [];
    for ch = 1:nchannelsSelected - 1
        subplot(nrows, ncols, idx_plot)
        ch_signal = nan(min_trial_data, ntrial);
        for t = 1:ntrial
            ch_signal(:,t) = phase_slope(:, b, channelsSelected(ch), t);
        end
        imagesc(ch_signal');
        cl = max(cl, abs(max(ch_signal, [], 'all')));
        handles = [handles; gca];

        title(['channel: ' channels_select{ch}])
        idx_plot = idx_plot + 1;
    end
    cl = max(cl, max(abs(ch_signal), [], 'all'));
    set(handles, 'clim', [-cl cl])
    sgtitle(['circular variance | band: ' bands_str{b}])
end



%% SVM one class model  --> TODO select the features!
% compute the plv value for each signal
plv = nan(min_trial_data, nbands, nchannels, ntrial);
for b = 1:nbands
    for ch = 1:nchannels
        for t = 1:ntrial
            c_phase = squeeze(phase_data(:,b,ch,t));
            plv(:,b,ch,t) = real(exp(1i * c_phase));
        end
    end
end

%% take only the fixation
X = nan(min_durFIX * ntrial, nchannels);
for t = 1:ntrial
    X((t-1)*min_durFIX + 1: t*min_durFIX,:) = plv(1:min_durFIX,1,:,t);
end
Y = ones(size(X, 1), 1);

SVM_model = fitcsvm(X, Y, ...
    'KernelFunction', 'rbf', ...
    'KernelScale', 'auto', ...
    'Standardize', true, ...
    'OutlierFraction', 0.10);  % Assume 5% of training data may be outliers
%%
for t= 1:ntrial
    X_test = squeeze(plv(:,1,:,t));
    
    [label, score] = predict(SVM_model, X_test);
end
