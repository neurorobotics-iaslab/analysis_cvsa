clc; clearvars;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/utils')
% addpath('/home/paolo/Local/cnbi-smrtrain')

%% Initialization
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
classes = [730 731];
cf_event = 781;
fix_event = 786;
nchannels = 39;
nclasses = length(classes);
filterOrder = 4;
avg = 1;

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files: select only one gdf', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);
%% concatenate the files
nFiles = length(filenames);
signals = [];
events.POS = [];
events.DUR = [];
events.TYP = [];
sampleRate = nan;
for idx_file= 1: nFiles
    fullpath_file_shift = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_shift);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;
    sampleRate = header.SampleRate;
    disp('      [proc] applying filtering')
    [b, a] = butter(4, 40*(2/sampleRate),'low');
    c_signal = filter(b,a,c_signal);
    [b, a] = butter(4, 1*(2/sampleRate),'high');
    c_signal = filter(b,a,c_signal);

    if isempty(find(header.EVENT.TYP == 2, 1)) % no eye calibration
        events.TYP = [events.TYP; header.EVENT.TYP];
        events.DUR = [events.DUR; header.EVENT.DUR];
        events.POS = [events.POS; header.EVENT.POS + size(signals, 1)];
        signals = [signals; c_signal];
    else
        k = find(header.EVENT.TYP == 1, 1);
        events.TYP = [events.TYP; header.EVENT.TYP(k:end)];
        events.DUR = [events.DUR; header.EVENT.DUR(k:end)];
        events.POS = [events.POS; header.EVENT.POS(k:end) + size(signals, 1)];
        signals = [signals; c_signal];
    end
end

%% Labelling data 
[~, channelsSelected] = ismember(channels_select, channels_label);
nchannelsSelected = size(channelsSelected, 2);
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
trial_data = nan(min_trial_data, nchannels, ntrial);
for trial = 1:ntrial
    c_start = trial_start(trial);
    c_end = trial_start(trial) + min_trial_data - 1;
    trial_data(:,:,trial) = signals(c_start:c_end,:);
end

%%
% Parameters
window = hamming(256);  % Or a suitable window length
nfft = 256;
f = 6:2:20;

% Allocate
spe_values = zeros(nchannels, ntrial);

for trial = 1:ntrial
    for ch = 1:nchannels
        signal = squeeze(trial_data(:, ch, trial));
        [Pxx, freqs] = pwelch(signal, window, 128, nfft, sampleRate);
        [~, ~, idx_freq] = intersect(f,freqs);


        % Normalize PSD
        Pxx_sub = Pxx(idx_freq) / sum(Pxx(idx_freq) + eps);
        
        % Compute Shannon entropy
        SpE = -sum(Pxx_sub .* log2(Pxx_sub + eps));  % use log2 or log depending on units
        spe_values(ch, trial) = SpE;
    end
end
imagesc(spe_values)