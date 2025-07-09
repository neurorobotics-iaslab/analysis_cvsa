%% for each trial compute the entropy
clc; clearvars;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/utils')
% addpath('/home/paolo/Local/cnbi-smrtrain') 
addpath(genpath('/home/paolo/Local/Matlab/MVGC1'))

%% Initialization
band = [1 40];
signals = [];
headers.TYP = [];
headers.DUR = [];
headers.POS = [];

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


    disp('      [proc] applying filtering')
    [b, a] = butter(filterOrder, band(2)*(2/header.SampleRate),'low');
    s_low = filter(b,a,c_signal);
    [b, a] = butter(filterOrder, band(1)*(2/header.SampleRate),'high');
    signal_processed = filter(b,a,s_low);

    headers.sampleRate = header.SampleRate;
    headers.channels_labels = header.Label;
    if isempty(find(header.EVENT.TYP == 2, 1)) % no eye calibration
        headers.TYP = cat(1, headers.TYP, header.EVENT.TYP);
        headers.DUR = cat(1, headers.DUR, header.EVENT.DUR);
        headers.POS = cat(1, headers.POS, header.EVENT.POS + size(signals, 1));
    else
        k = find(header.EVENT.TYP == 1, 1);
        headers.TYP = cat(1, headers.TYP, header.EVENT.TYP(k:end));
        headers.DUR = cat(1, headers.DUR, header.EVENT.DUR(k:end));
        headers.POS = cat(1, headers.POS, header.EVENT.POS(k:end) + size(signals, 1));
    end

    signals = cat(1, signals, signal_processed(:,:));
end


%% Labelling data
events = headers;
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
min_durFIX = min(fixDUR);

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

c_signal = signals;
for trial = 1:ntrial
    c_start = trial_start(trial);
    c_end = trial_start(trial) + min_trial_data - 1;
    trial_data(:,:,trial) = c_signal(c_start:c_end,:);
end


%% refactoring the data
% now the data are placed alternatevely, so odd trial id is for 730, even for 731
balanced_trial_idx = balanced_data_afterEOG(trial_with_eog, trial_typ, classes);
balanced_trial_data = trial_data(:,:,logical(balanced_trial_idx));
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
        tmp(:,:,idx_trial_class + idx_class - 1) = balanced_trial_data(:,:,idx_classes_trial(i, idx_class));
        trial_typ(idx_trial_class + idx_class - 1) = classes(idx_class);
    end
    i = i + 1;
end
trial_data = tmp;

%%
channels_select = {'FC1', 'FC2', 'F1', 'F2', 'O1', 'O2', 'PO8', 'PO7', 'PO6', 'PO5', 'PO4', 'PO3'};
[~, channelsSelected] = ismember(channels_select, channels_label);
nchannelsSelected = size(channelsSelected, 2);
trial_data_mvgc = permute(trial_data, [2 1 3]);
trial_data_mvgc = trial_data_mvgc(channelsSelected,:,:);

% Use one trial to estimate optimal model order
max_order = 15;
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(squeeze(trial_data_mvgc(:,:,1)), max_order, 'LWR');
morder = moAIC;  % Or use moBIC
fprintf('Selected model order: %d\n', morder);

% 3. Compute GC for baseline and task periods
regmode = 'OLS';
alpha   = 0.05;

% Baseline period
data_base = trial_data_mvgc(:,min_durFIX:min_durFIX + min_durCUE,:);
data_base = data_base(:,1:2:end,:);
[A_base, SIG_base] = tsdata_to_var(data_base, morder, regmode);
assert(~isbad(A_base), 'VAR estimation failed (baseline)');
F_base = var_to_autocov(A_base, SIG_base, morder);
assert(~isbad(F_base), 'Autocovariance failed (baseline)');
G_base = zeros(nchannelsSelected, nchannelsSelected);
for i = 1:nchannelsSelected
    for j = 1:nchannelsSelected
        if i == j
            G_base(i,j) = NaN;  
        else
            G_base(i,j) = autocov_to_mvgc(F_base, i, j);  % GC j→i
        end
    end
end
fprintf('Baseline GC computed.\n');

% Task period
data_task = trial_data_mvgc(:,min_durFIX + min_durCUE+1:min_trial_data,:);
[A_task, SIG_task] = tsdata_to_var(data_task, morder, regmode);
assert(~isbad(A_task), 'VAR estimation failed (task)');
F_task = var_to_autocov(A_task, SIG_task, morder);
assert(~isbad(F_task), 'Autocovariance failed (task)');
G_task = zeros(nchannelsSelected, nchannelsSelected);
for i = 1:nchannelsSelected
    for j = 1:nchannelsSelected
        if i == j
            G_task(i,j) = NaN;  
        else
            G_task(i,j) = autocov_to_mvgc(F_task, i, j);  % GC j→i
        end
    end
end
fprintf('Task GC computed.\n')

%% plot the matrix
figure;
subplot(1,2,1);
imagesc(G_base); colorbar;
title('Baseline Granger Causality'); xlabel('From'); ylabel('To');

subplot(1,2,2);
imagesc(G_task); colorbar;
title('Task Granger Causality'); xlabel('From'); ylabel('To');

G_diff = G_task - G_base;

figure;
imagesc(G_diff); colorbar;
title('Task - Baseline GC Difference');
xlabel('From'); ylabel('To');

%% time analysis
win_len = 100;   % number of timepoints in window
step = 50;
nwin = floor((min_trial_data - win_len)/step) + 1;

G_time = zeros(nchannelsSelected, nchannelsSelected, nwin);

for w = 1:nwin
    t1 = (w-1)*step + 1;
    t2 = t1 + win_len - 1;
    segment = trial_data_mvgc(:,t1:t2,:);
    
    [A_win, SIG_win] = tsdata_to_var(segment, morder, regmode);
    if isbad(A_win), continue; end
    
    F_win = var_to_autocov(A_win, SIG_win);
    if isbad(F_win), continue; end
    
    G_win = zeros(nchannelsSelected, nchannelsSelected);
    for i = 1:nchannelsSelected
        for j = 1:nchannelsSelected
            if i == j
                G_win(i,j) = NaN;
            else
                G_win(i,j) = autocov_to_mvgc(F_task, i, j);  % GC j→i
            end
        end
    end
    G_time(:,:,w) = G_win;
end

%% Plot GC over time between two channels
from_ch = 1; to_ch = 3;
figure;
plot(squeeze(G_time(to_ch, from_ch, :)));
xlabel('Time window'); ylabel('GC');
title(sprintf('GC from %d → %d over time', from_ch, to_ch));
