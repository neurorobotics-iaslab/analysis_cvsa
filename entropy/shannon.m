%% for each trial compute the entropy entropy
clc; clearvars;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/utils')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/entropy/utils')

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

        signal_processed = proc_entropy_512hz(c_signal, header.SampleRate, band, filterOrder);
        
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

%% compute entropy single trial
window_size = floor(100 * sampleRate / 1000);
step_size = floor(50 * sampleRate / 1000); 

num_windows = floor((min_trial_data - window_size) / step_size) + 1;
entropy_matrix = zeros(num_windows, nchannels,  ntrial, nbands);  %% windows x channels x trials x bands

for idx_band = 1:nbands
    disp(['[INFO] computing entropy for band: ', bands_str{idx_band}])
    for t = 1:ntrial
        c_signal = trial_data(:,idx_band,:,t);

        trial_entropy = eeg_entropy_windowed(c_signal, 'shannon', window_size, step_size);
        entropy_matrix(:,:,t, idx_band) = trial_entropy;
    end
end

%% plotting
figure();
idx_plot = 1;
for idx_band = 1:nbands
    subplot(nbands,2,idx_plot)
    imagesc(squeeze(mean(entropy_matrix(:,:,trial_typ == classes(1), idx_band), 3))')
    hold on
    xline(floor((min(fixDUR) - window_size ) / step_size + 1), '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(floor((min(fixDUR) - window_size  + min_durCUE) / step_size + 1), '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off
    idx_plot = idx_plot + 1;
    yticks(1:nchannels)
    yticklabels(channels_label(1:nchannels))
    xlabel('window')
    title(['band: ' bands_str{idx_band} ' | class: ' num2str(classes(1))])

    subplot(nbands, 2, idx_plot)
    imagesc(squeeze(mean(entropy_matrix(:,:,trial_typ == classes(2), idx_band), 3))')
    hold on
    xline(floor((min(fixDUR) - window_size ) / step_size + 1), '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(floor((min(fixDUR) - window_size  + min_durCUE) / step_size + 1), '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off
    idx_plot = idx_plot + 1;
    yticks(1:nchannels)
    yticklabels(channels_label(1:nchannels))
    xlabel('window')
    title(['band: ' bands_str{idx_band} ' | class: ' num2str(classes(2))])
end
sgtitle('shannon entropy | mean over trial')



%% funtinos
function entropy_matrix = eeg_entropy_windowed(eeg_signal, method, window_size, step_size)
    % EEG_ENTROPY_WINDOWED Computes entropy in sliding windows
    %
    % entropy_matrix = eeg_entropy_windowed(eeg_signal, method, window_size, step_size)
    %
    % eeg_signal: matrix (time x channels) of EEG data
    % method: 'shannon' or 'approx'
    % window_size: number of samples per window
    % step_size: step between consecutive windows
    %
    % Returns: entropy_matrix (windows x channels)

    [num_samples, num_channels] = size(eeg_signal);
    num_windows = floor((num_samples - window_size) / step_size) + 1;
    entropy_matrix = zeros(num_windows, num_channels);

    for win = 1:num_windows
        start_idx = round((win - 1) * step_size + 1);
        end_idx = round(start_idx + window_size - 1);
        window_data = eeg_signal(start_idx:end_idx, :);

        a = eeg_entropy(window_data, method);
        entropy_matrix(win, :) = a;
    end
end
function entropy_values = eeg_entropy(eeg_signal, method)
    % EEG_ENTROPY Computes entropy for each channel of an EEG signal
    %
    % entropy_values = eeg_entropy(eeg_signal, method)
    %
    % eeg_signal: matrix (time x channels) of EEG data
    % method: 'shannon' or 'approx'
    %
    % Returns a vector of entropy values for each channel
    
    if nargin < 2
        method = 'shannon'; % Default method
    end
    
    [num_samples, num_channels] = size(eeg_signal);
    entropy_values = zeros(1, num_channels);
    
    for ch = 1:num_channels
        signal = eeg_signal(:, ch); % Extract individual channel
        
        switch lower(method)
            case 'shannon'
                % Normalize the signal between 0 and 1
                signal = signal - min(signal);
                signal = signal / max(signal);
                
                % Create probability histograms
%                 nbins = 32; % dividend the range 0-1 into 32 bins (same as paper stefano)
                nbins = ceil(log2(size(signal, 1))) + 1;  % --> try it
                p = histcounts(signal, nbins, 'Normalization', 'probability'); %probability to fall in a bin
                
                % Remove zeros to avoid log(0)
                p(p == 0) = [];
                
                % Compute Shannon entropy
                entropy_values(ch) = -sum(p .* log2(p));
            
            case 'approx'
                % Typical parameters for approximate entropy
                m = 2;
                r = 0.2 * std(signal);
                
                entropy_values(ch) = approx_entropy(signal, m, r);
            
            otherwise
                error('Unrecognized method. Use "shannon" or "approx".');
        end
    end
end

function ApEn = approx_entropy(signal, m, r)
    % APPROX_ENTROPY Computes approximate entropy of a signal
    N = length(signal);
    r = r * std(signal); % Tolerance factor
    
    C = zeros(1, N - m + 1);
    for i = 1:N - m + 1
        count = 0;
        for j = 1:N - m + 1
            if max(abs(signal(i:i+m-1) - signal(j:j+m-1))) <= r
                count = count + 1;
            end
        end
        C(i) = count / (N - m + 1);
    end
    
    phi_m = sum(log(C)) / (N - m + 1);
    
    % Repeat for m+1
    C = zeros(1, N - m);
    for i = 1:N - m
        count = 0;
        for j = 1:N - m
            if max(abs(signal(i:i+m) - signal(j:j+m))) <= r
                count = count + 1;
            end
        end
        C(i) = count / (N - m);
    end
    
    phi_m1 = sum(log(C)) / (N - m);
    
    % Compute approximate entropy
    ApEn = phi_m - phi_m1;
end