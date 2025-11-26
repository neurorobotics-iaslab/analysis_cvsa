%% HMM to label the samples in time
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/EOG')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')
addpath(genpath('/home/paolo/Local/Matlab/HMM-MAR'))

%% === LOADING DATA ===
% Initialization
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
    trial_start(idx_trial) = cfPOS(idx_trial);
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
trial_data = tmp; % samples x bands x channels x trials --> only cf

%% show the log band diff in time for each channels and band. Then show the fisher score
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'}; 
[~, channelsSelected] = ismember(channels_select, channels_label);
nchannelsSelected = size(channelsSelected, 2);

fisher = nan(nbands, nchannelsSelected);
for b = 1:nbands
    c1 = squeeze(mean(mean(trial_data(:,b,channelsSelected,trial_typ == classes(1)), 4), 1));
    c2 = squeeze(mean(mean(trial_data(:,b,channelsSelected,trial_typ == classes(2)), 4), 1));

    std1 = squeeze(std(std(trial_data(:,b,channelsSelected,trial_typ==classes(1)), 0, 4), 0, 1));
    std2 = squeeze(std(std(trial_data(:,b,channelsSelected,trial_typ==classes(2)), 0, 4), 0, 1));

    fisher(b,:) = abs(c1 - c2) ./ sqrt(power(std1, 2) + power(std2, 2));
end

% plot fisher and diff log in mean
figure(); handles = []; cl=-inf;
for b = 1:nbands
    tmp1 = squeeze(mean(trial_data(:,b,channelsSelected, trial_typ==classes(1)), 4));
    tmp2 = squeeze(mean(trial_data(:,b,channelsSelected, trial_typ==classes(2)), 4));
    diff = abs(tmp1 - tmp2);

    subplot(ceil(nbands/2), 2, b)
    imagesc(diff')
    yticks(1:nchannelsSelected); yticklabels(channels_select);
    title(['band: ' bands_str{b}])
    handles = [handles; gca];
    cl = max(cl, max(diff, [], 'all'));
end
set(handles, 'clim', [0, cl]);
subplot(ceil(nbands/2), 2, nbands+1)
imagesc(fisher')
title('Fisher'); 
xticks(1:nbands); xticklabels(bands_str);
yticks(1:nchannelsSelected); yticklabels(channels_select);

%% ask the user the channels and the bands
select_band = {[8 10], [10 12], [12 14]}; 
channels_select = {{'PO5', 'PO3', 'PO7'}, {'P5', 'PO5', 'PO3', 'PO7'}, {'P5'}};
channels_select_flatten = [channels_select{:}];
total_electrodes = sum(cellfun(@(x) numel(x), channels_select));

idx_select_bands = find(cellfun(@(x) any(cellfun(@(y) isequal(x, y), select_band)), bands));
trial_data_select = nan(min_trial_data, total_electrodes, ntrial);
idx_ch = 1;
for i = 1:size(channels_select, 2)
    c_band = idx_select_bands(i);
    tmp_channels = channels_select{i};
    [~, c_channelsSelected] = ismember(tmp_channels, channels_label);

    c_data = squeeze(trial_data(:,c_band,:,:));
    trial_data_select(:,idx_ch:idx_ch+size(tmp_channels,2)-1,:) = c_data(:,c_channelsSelected,:);

    idx_ch = idx_ch + size(tmp_channels, 2); % update the channel for the trial_data_select
end

%% prepare the dataset
percentual = 0.7;
max_idx_trial_train = ceil(ntrial * percentual);
trial_data_train = trial_data_select(:,:,1:max_idx_trial_train);
trial_data_test = trial_data_select(:,:, max_idx_trial_train+1:end);

%% === HMM WITH SYMBOLS ===
% windowing parameters
winLength = round(0.1 * sampleRate);    % 100 ms window
stepSize  = round(0.05 * sampleRate);   % 50 ms step

% features extraction
nWins = floor((min_trial_data - winLength) / stepSize);
all_feats = nan(nWins * max_idx_trial_train, total_electrodes); % [total_windows x channels]

for t = 1:max_idx_trial_train
    data = squeeze(trial_data_train(:, :, t));  % [samples x chans]
    feat = zeros(nWins, total_electrodes);
    for w = 1:nWins
        idx = (w-1)*stepSize + (1:winLength);
        window = data(idx, :);
        feat(w,:) = mean(window, 1);  % mean on samples
    end

    all_feats((t-1)*nWins+1:t*nWins, :) = feat;
end

% applying the zscore
mu = mean(all_feats, 1);
sigma = std(all_feats, 0, 1);
all_feats = (all_feats - mu) ./ sigma; 

%% understand # symbols
maxK = 42;
distortions = zeros(1, maxK);
for k = 1:maxK
    [~, ~, sumd] = kmeans(all_feats, k, 'Replicates', 5, 'Distance', 'sqeuclidean');
    distortions(k) = sum(sumd);
end
figure();
plot(1:maxK, distortions); xlabel('k'); ylabel('Distortion'); title('Elbow method'); % look the graph to undersstan which is the best value for the numSymbols

%% quantize features into symbols and train the HMM model
numStates = 3;                           % HMM with 3 hidden states
numSymbols = 3;                         % Quantization bins

% computing the kmeans to have the symbols for each ssample
[idxs_cluster, C] = kmeans(all_feats, numSymbols, 'Distance', 'sqeuclidean', 'Replicates', 5); % try cosine or sqeuclidean

% initialize transitions and emissions randomly
TRANS = normalize(rand(numStates), 2);
EMIS = normalize(rand(numStates, numSymbols), 2);

% train the HMM model
[ESTTR, ESTEMIT] = hmmtrain(idxs_cluster', TRANS, EMIS, 'Maxiterations', 300, 'ALGORITHM', 'BaumWelch', 'VERBOSE',true); % try BaumWelch or Viterbi
[pred_train, ~] = hmmviterbi(idxs_cluster', ESTTR, ESTEMIT);

%% plotting the results with the discrete HMM model (only for the train)
% Posterior probabilities of hidden states at each time step
[posterior, logposterior] = hmmdecode(idxs_cluster', ESTTR, ESTEMIT);  % posterior is [numStates x T]

% Visualize for a specific trial
handles = []; cl = -inf;
for idx_trial = 1:10
    idx_start = (idx_trial-1)*nWins+1;
    idx_end = idx_trial*nWins;

    figure();
    subplot(3+size(posterior, 1),1,1)
    imagesc(all_feats((idx_trial-1)*nWins+1:idx_trial*nWins, :)')
    yticks(1:size(channels_select_flatten,2)); yticklabels(channels_select_flatten)
    title('Log band'); handles = [handles; gca]; cl = max(cl, abs(all_feats((idx_trial-1)*nWins+1:idx_trial*nWins, :)));

    subplot(3+size(posterior, 1),1,2)
    plot(pred_train((idx_trial-1)*nWins+1:idx_trial*nWins)); title('Decoded HMM States');
    ylim([0.5 numStates+0.5]); yticks(0:numStates); xlim([0 nWins])

    subplot(3+size(posterior, 1),1,3)
    plot(idxs_cluster((idx_trial-1)*nWins+1:idx_trial*nWins))
    title('Decoded kmean States');xlim([0 nWins])
    ylim([0.5 numSymbols+0.5])

    for i = 1:size(posterior, 1)
        subplot(3+size(posterior, 1),1,3+i)
        plot(posterior(i, idx_start:idx_end));
        title(['Posterior ' num2str(i)]);
        ylim([-0.1 1.1]); xlim([0 nWins])
    end
    sgtitle(['trial: ' num2str(idx_trial) ' | class: ' num2str(trial_typ(idx_trial))]);
end
set(handles, 'clim', [0 cl]);

