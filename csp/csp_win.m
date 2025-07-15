%% check the CSP as features
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
bands{i+1} = [8 14];
% bands = {[8 14]};
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
        disp('      [proc] applying filtering')
        [b, a] = butter(filterOrder, band(2)*(2/header.SampleRate),'low');
        s_low = filter(b,a,c_signal);
        [b, a] = butter(filterOrder, band(1)*(2/header.SampleRate),'high');
        s_filt = filter(b,a,s_low);
        signal_processed = s_filt;

        
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
trial_data = tmp;

channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
[~, channelsSelected] = ismember(channels_select, channels_label);
nchannelsSelected = size(channelsSelected, 2);
trial_data = trial_data(:,:,channelsSelected,:);

%% divide all into windows ---> work for only one band
win_size = 100; %ms
overlap = 50;
overlap = ceil(overlap * sampleRate / 1000);
win_size = ceil(win_size * sampleRate /1000);
nwin = floor((min_trial_data - win_size)/overlap) + 1;

trial_data_1 = trial_data(:,:,:,trial_typ == classes(1));
trial_data_2 = trial_data(:,:,:,trial_typ == classes(2));

percentual_train = 0.7;
trial_data_1_train = trial_data_1(:,:,:,1:size(trial_data_1, 4)*percentual_train);
trial_data_1_test = trial_data_1(:,:,:,size(trial_data_1, 4)*percentual_train + 1:end);
trial_data_2_train = trial_data_2(:,:,:,1:size(trial_data_1, 4)*percentual_train);
trial_data_2_test = trial_data_2(:,:,:,size(trial_data_1, 4)*percentual_train + 1:end);

W_sel = cell(1, nbands);
idxs_bands = [];
n_w = 3;
for idx_b = 1:nbands
    C1 = zeros(size(trial_data_1_train, 3), size(trial_data_1_train, 3));  % For class 1
    C2 = zeros(size(trial_data_1_train, 3), size(trial_data_1_train, 3));  % For class 2

    for t = 1:size(trial_data_1_train, 4)
        for idx_w = 1:nwin
            start_w = ceil((idx_w-1)*overlap + 1);
            end_w = ceil(start_w + win_size - 1);

            c_signal_1 = squeeze(trial_data_1_train(start_w:end_w,idx_b,:,t));
            c_signal_2 = squeeze(trial_data_2_train(start_w:end_w,idx_b,:,t));

            c_signal_1 = c_signal_1 - mean(c_signal_1, 1);
            c_signal_2 = c_signal_2 - mean(c_signal_2, 1);

            C1 = C1 + (c_signal_1' * c_signal_1) / trace(c_signal_1' * c_signal_1);
            C2 = C2 + (c_signal_2' * c_signal_2) / trace(c_signal_2' * c_signal_2);
        end
    end
    C1 = C1 / (nwin * size(trial_data_1_train, 4));
    C2 = C2 / (nwin * size(trial_data_2_train, 4));

    C = C1 + C2;
    [W, D] = eig(C1, C);
    % Sort eigenvalues (descending)
    [~, ind] = sort(diag(D), 'descend');
    W = W(:, ind);

    % Select first & last 2 filters
    c_W_sel = [W(:, 1:n_w), W(:, end-n_w+1:end)];
    W_sel{idx_b} = c_W_sel;

    idxs_bands = [idxs_bands; ones(1,size(c_W_sel, 2))*idx_b];
end


%% features extraction
X_feat = [];
Y_feat = [];

for idx_b = 1:nbands
    c_W_sel = W_sel{idx_b};
    c_X_feat = [];
    c_Y_feat = [];
    for t = 1:size(trial_data_1_train, 4)
        for idx_w = 1:nwin
            start_w = ceil((idx_w-1)*overlap + 1);
            end_w = ceil(start_w + win_size - 1);

            c_signal_1 = squeeze(trial_data_1_train(start_w:end_w,idx_b,:,t));
            c_signal_2 = squeeze(trial_data_2_train(start_w:end_w,idx_b,:,t));

            c_signal_1 = c_signal_1 - mean(c_signal_1, 1);
            c_signal_2 = c_signal_2 - mean(c_signal_2, 1);

            Z_1 = c_W_sel' * c_signal_1';
            Z_2 = c_W_sel' * c_signal_2'; 

            feat_1 = log(var(Z_1, 0, 2) / sum(var(Z_1, 0, 2))); 
            feat_2 = log(var(Z_2, 0, 2) / sum(var(Z_2, 0, 2))); 

            c_X_feat = [c_X_feat; feat_1'];
            c_Y_feat = [c_Y_feat; classes(1)];

            c_X_feat = [c_X_feat; feat_2'];
            c_Y_feat = [c_Y_feat; classes(2)];
        end
    end
    X_feat = [X_feat, c_X_feat];
    Y_feat = c_Y_feat;
end

%% features selection with fisher
fisher = nan(nbands,nclasses * n_w);
for idx_band = 1:nbands
    for idx_w = 1:n_w*nclasses
        c1 = X_feat(Y_feat == classes(1),(idx_band)*idx_w);
        c2 = X_feat(Y_feat == classes(2),(idx_band)*idx_w);

        mu1 = mean(c1, 'all');
        sigma1 = std(c1);

        mu2 = mean(c2, 'all');
        sigma2 = std(c2);

        fisher(idx_band, idx_w) = abs(mu1 - mu2) ./ sqrt(sigma1^2 + sigma2^2);
    end
end
figure();
imagesc(fisher')
xticks(1:nbands)
xticklabels(bands_str)
yticks(1:n_w*nclasses)

selected_rows = [5,5,6,6,6,4];
selected_cols = [5,6,5,6,7,9];

% X_feat = X_feat(:,selected_rows.*selected_cols);

%% tried a QDA
Mdl = fitcdiscr(X_feat, Y_feat, 'DiscrimType', 'linear'); % LDA
% Mdl = fitcsvm(X_feat, Y_feat, 'KernelFunction', 'linear');
cvMdl = crossval(Mdl, 'KFold', 5);
acc = 1 - kfoldLoss(cvMdl);
fprintf('Cross-validated accuracy: %.2f%%\n', acc*100);

%% test QDA
X_feat_test = [];
Y_feat_test = [];

for idx_b = 1:nbands
    c_W_sel = W_sel{idx_b};
    c_X_feat = [];
    c_Y_feat = [];
    for t = 1:size(trial_data_1_test, 4)
        for idx_w = 1:nwin
            start_w = ceil((idx_w-1)*overlap + 1);
            end_w = ceil(start_w + win_size - 1);

            c_signal_1 = squeeze(trial_data_1_test(start_w:end_w,idx_b,:,t));
            c_signal_2 = squeeze(trial_data_2_test(start_w:end_w,idx_b,:,t));

            c_signal_1 = c_signal_1 - mean(c_signal_1, 1);
            c_signal_2 = c_signal_2 - mean(c_signal_2, 1);

            Z_1 = c_W_sel' * c_signal_1';
            Z_2 = c_W_sel' * c_signal_2'; 

            feat_1 = log(var(Z_1, 0, 2) / sum(var(Z_1, 0, 2))); 
            feat_2 = log(var(Z_2, 0, 2) / sum(var(Z_2, 0, 2))); 

            c_X_feat = [c_X_feat; feat_1'];
            c_Y_feat = [c_Y_feat; classes(1)];

            c_X_feat = [c_X_feat; feat_2'];
            c_Y_feat = [c_Y_feat; classes(2)];
        end
    end
    X_feat_test = [X_feat_test, c_X_feat];
    Y_feat_test = c_Y_feat;
end

% X_feat_test = X_feat_test(:,selected_rows.*selected_cols);

pred = predict(Mdl, X_feat_test);
acc = sum(pred == Y_feat_test) / size(pred, 1);
fprintf('Cross-validated accuracy: %.2f%%\n', acc*100);
