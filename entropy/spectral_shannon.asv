%% for each trial compute the entropy
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
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
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
%     disp('      [proc] applying filtering')
%     [b, a] = butter(4, 60*(2/sampleRate),'low');
%     c_signal = filter(b,a,c_signal);
%     [b, a] = butter(4, 1*(2/sampleRate),'high');
%     c_signal = filter(b,a,c_signal);

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
%% applying spectrogram
% disp('  [PROC] applaying spectrogram');
% win_length = 100; %ms
% win_overlap = 50; % ms
% window_length = ceil(win_length / 1000 * sampleRate);
% overlap = ceil(win_overlap / 1000 * sampleRate); % in samples
% num_windows = floor((size(signals, 1) - overlap) / (window_length - overlap));
% f = 6:2:20;
% nfreqs_selected = size(f,2);
% nfft      = 256; % default
% if(mod(nfft, 2) == 0)
%     nfreqs = (nfft/2) + 1;
% else
%     nfreqs = (nfft+1)/2;
% end
% 
% psd = nan(nfreqs, num_windows, nchannels); % dimension of psd
% for chId = 1:nchannels
%     [~,freqs,~,psd(:,:,chId)] = spectrogram(signals(:,chId), window_length, overlap, [], sampleRate);
% end
% [~, ~, idx_freq] = intersect(f,freqs);
% 
% psd_event.TYP = events.TYP;
% psd_event.POS = floor((events.POS - window_length) / (window_length - overlap)) + 1;
% psd_event.DUR = round(events.DUR/(window_length - overlap)) + 1; 

%%% equal in ros
wshift = 0.0625; % we have new signal after 0.0625 sec
pshift = 0.25;   % shifting inside the windows of 1 sec --> this should be 50 ms prev: 0.25
mlength = 1;     % moving average
wlength = 0.5;   % length of the internal windows --> this should be 100 ms prev: 0.5
wconv = 'backward';
f = 4:2:46; % in order to understand interesting features in a smaller
% number of frequencies
[PSD, freqs] = proc_spectrogram(signals, wlength, wshift, pshift, sampleRate, mlength); % windows x frequency x channels
[~, ~, idx_freq] = intersect(f,freqs);

disp('  [PROC] updating events from samples to windows');
psd_event.TYP = events.TYP;
psd_event.POS = proc_pos2win(events.POS, wshift*sampleRate, wconv, mlength*sampleRate); 
psd_event.conv = wconv;
psd_event.DUR = round(events.DUR/(wshift*sampleRate)) + 1; 



%% Labelling data 
[~, channelsSelected] = ismember(channels_select, channels_label);
nchannelsSelected = size(channelsSelected, 2);
cuePOS = psd_event.POS(ismember(psd_event.TYP, classes));
cueDUR = psd_event.DUR(ismember(psd_event.TYP, classes));
cueTYP = psd_event.TYP(ismember(psd_event.TYP, classes));

fixPOS = psd_event.POS(psd_event.TYP == 786);
fixDUR = psd_event.DUR(psd_event.TYP == 786);

cfPOS = psd_event.POS(psd_event.TYP == 781);
cfDUR = psd_event.DUR(psd_event.TYP == 781);

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
trial_data = nan(nfreqs_selected, min_trial_data, nchannels, ntrial);
for trial = 1:ntrial
    c_start = trial_start(trial);
    c_end = trial_start(trial) + min_trial_data - 1;
    trial_data(:,:,:,trial) = psd(idx_freq,c_start:c_end,:);
end

%%
spectral_entropy = nan(min_trial_data, nchannels, ntrial);

for t = 1:ntrial
    for ch = 1:nchannels
        for ti = 1:min_trial_data
            % Get power at all selected frequencies for this time point
            psd_slice = trial_data(:, ti, ch, t);
            
            % Normalize to sum to 1 (so it becomes a probability distribution)
            psd_norm = psd_slice / sum(psd_slice + eps);
            
            % Compute Shannon entropy
            spectral_entropy(ti, ch, t) = -sum(psd_norm .* log2(psd_norm + eps));
        end
    end
end

avg_entropy_1 = squeeze(mean(spectral_entropy(:,:,trial_typ == classes(1)), 3));  % [num_windows x channels]
avg_entropy_2 = squeeze(mean(spectral_entropy(:,:,trial_typ == classes(2)), 3));  % [num_windows x channels]

figure();
subplot(1,2,1)
imagesc(avg_entropy_1(:,channelsSelected)');
yticks(1:numel(channelsSelected));
yticklabels(channels_label(channelsSelected))
xlabel('windows')
hold on
xline(min(fixDUR), '--r', 'LineWidth', 1.5);  % Dashed red line
xline(min(fixDUR)+min(cueDUR), '--r', 'LineWidth', 1.5);  % Dashed blue line
hold off
title(['class: ' num2str(classes(1))]); colorbar;

subplot(1,2,2)
imagesc(avg_entropy_2(:,channelsSelected)');
yticks(1:numel(channelsSelected));
yticklabels(channels_label(channelsSelected))
xlabel('windows')
hold on
xline(min(fixDUR), '--r', 'LineWidth', 1.5);  % Dashed red line
xline(min(fixDUR)+min(cueDUR), '--r', 'LineWidth', 1.5);  % Dashed blue line
hold off
title(['class: ' num2str(classes(2))]); colorbar;
