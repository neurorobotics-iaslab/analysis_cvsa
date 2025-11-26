%% In this file i will work with cluster to labels the data
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

    for idx_band = 1:nbands
        band = bands{idx_band};

        % for log band extraction
        disp('   [proc] logband');

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
trial_data = tmp;

%% BASIC CASE (one band)
%% start clustering
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'}; 
[~, channelsSelected] = ismember(channels_select, channels_label);
nchannelsSelected = size(channelsSelected, 2);

c_trial_data = squeeze(trial_data(:,nbands,channelsSelected,:)); % take the 8-14 band and all 17 electrons
percentual_train = 0.7;
ntrial_train = ceil(ntrial*percentual_train);
ntrial_test = ntrial - ntrial_train;

all_data = nan(min_trial_data * ntrial_train, nchannelsSelected);
for i = 1:ntrial_train
    all_data(1+(i-1)*min_trial_data:i*min_trial_data,:) = squeeze(c_trial_data(:,:,i));
end
% all_data = zscore(all_data); % normalize data with zscore

K = 3; % number of states/clusters
[cluster_labels, C] = kmeans(all_data, K, 'Distance', 'sqeuclidean', 'Replicates', 10, 'Display', 'iter');

% Reshape cluster labels back into time x trial
cluster_per_trial =  nan(min_trial_data, ntrial_train);
for i = 1:ntrial_train
    cluster_per_trial(:,i) = cluster_labels(1+(i-1)*min_trial_data:i*min_trial_data);
end

% take centroid labels
[~, sortIdx] = sort(mean(C'));  % Sort clusters from low to high power
labels = ["Nothing", "Shift", "Sustained"];
mappedLabels = strings(size(mean(C')));
mappedLabels(sortIdx) = labels; % in this way we associate the correct label to the class

% cl = -inf; handles = [];
% for i = 1:ntrial
%     figure;
%     subplot(2,1,1)
%     plot(cluster_per_trial(:, i), 'LineWidth', 2);
%     xlabel('Time (samples)');
%     ylabel('Cluster ID');
%     title(sprintf('Temporal States - Trial %d', i));
%     ylim([0.5, K + 0.5]);
%     yticks(1:K);
%     yticklabels(mappedLabels)
%     xlim([0 min_trial_data]);
% 
%     subplot(2,1,2)
%     imagesc(squeeze(c_trial_data(:,:,i))')
%     handles = [handles; gca];
%     cl = max(cl, max(abs(squeeze(c_trial_data(:,:,i))), [], 'all'));
% end
% set(handles, 'clim', [0, cl])

figure(); hold on;
for i = 1:K
    plot(C(i,:))
end
title('centroids')
legend(mappedLabels)

%% take the data inside the shift and the sustained --> train
shfit_data = []; shift_label = [];
sustaiend_data = []; sustaiend_label = [];
noth_data = []; noth_label = [];
for i = 1:ntrial_train
    for sample = 1:min_trial_data
        if mappedLabels(cluster_per_trial(sample,i)) == "Shift"
            shfit_data = [shfit_data; squeeze(c_trial_data(sample,:,i))];
            shift_label = [shift_label; trial_typ(i)];
        elseif mappedLabels(cluster_per_trial(sample,i)) == "Sustained"
            sustaiend_data = [sustaiend_data; squeeze(c_trial_data(sample,:,i))];
            sustaiend_label = [sustaiend_label; trial_typ(i)];
        else
            noth_data = [noth_data; squeeze(c_trial_data(sample,:,i))];
            noth_label = [noth_label; trial_typ(i)];
        end
    end
end

%% compute the fisher score
fisher_shift = nan(1, nchannelsSelected);
fisher_sus = nan(1, nchannelsSelected);
fisher_noth = nan(1, nchannelsSelected);
for idx_ch=1:nchannelsSelected
    % shift
    mu1 = mean(shfit_data(shift_label == classes(1),idx_ch));
    sigma1 = std(shfit_data(shift_label == classes(1),idx_ch));
    mu2 = mean(shfit_data(shift_label == classes(2),idx_ch));
    sigma2 = std(shfit_data(shift_label == classes(2),idx_ch));
    fisher_shift(idx_ch) = abs(mu1 - mu2);

    % sus
    mu1 = mean(sustaiend_data(sustaiend_label == classes(1),idx_ch));
    sigma1 = std(sustaiend_data(sustaiend_label == classes(1),idx_ch));
    mu2 = mean(sustaiend_data(sustaiend_label == classes(2),idx_ch));
    sigma2 = std(sustaiend_data(sustaiend_label == classes(2),idx_ch));
    fisher_sus(idx_ch) = abs(mu1 - mu2);

    % noth
    mu1 = mean(noth_data(noth_label == classes(1),idx_ch));
    sigma1 = std(noth_data(noth_label == classes(1),idx_ch));
    mu2 = mean(noth_data(noth_label == classes(2),idx_ch));
    sigma2 = std(noth_data(noth_label == classes(2),idx_ch));
    fisher_noth(idx_ch) = abs(mu1 - mu2);
end

tmp = [fisher_shift; fisher_sus; fisher_noth];
figure();
imagesc(tmp')
xticks(1:3); xticklabels([{"shift"}, {"sustained"}, {"nothing"}])

%% train the models
idx_features_shift = 1:17; %[9 17];
idx_features_sus =   1:17; %[3 6 8 10 14 16 12 13];
idx_features_noth =  1:17; %[1 5 7 11 15];
% shift
qdaModel_shift = fitcdiscr(shfit_data(:,idx_features_shift), shift_label, 'DiscrimType', 'quadratic');
CVModel = crossval(qdaModel_shift, 'KFold', 5);
loss = kfoldLoss(CVModel);  % Average classification error
fprintf("5-fold cross-validation loss: %.4f\n", loss);
p = predict(qdaModel_shift, shfit_data(:,idx_features_shift));
confusionmat(shift_label, p)
disp(['Accuracy: ' num2str(sum(p == shift_label)/size(p,1))])

% sustained
qdaModel_sus = fitcdiscr(sustaiend_data(:,idx_features_sus), sustaiend_label, 'DiscrimType', 'quadratic');
CVModel = crossval(qdaModel_sus, 'KFold', 5);
loss = kfoldLoss(CVModel);  % Average classification error
fprintf("5-fold cross-validation loss: %.4f\n", loss);
p = predict(qdaModel_sus, sustaiend_data(:,idx_features_sus));
confusionmat(sustaiend_label, p)
disp(['Accuracy: ' num2str(sum(p == sustaiend_label)/size(p,1))])

% nothing
qdaModel_noth = fitcdiscr(noth_data(:,idx_features_noth), noth_label, 'DiscrimType','quadratic');
CVModel = crossval(qdaModel_noth, 'KFold', 5);
loss = kfoldLoss(CVModel);  % Average classification error
fprintf("5-fold cross-validation loss: %.4f\n", loss);
p = predict(qdaModel_noth, noth_data(:,idx_features_noth));
confusionmat(noth_label, p)
disp(['Accuracy: ' num2str(sum(p == noth_label)/size(p,1))])

%% See using all the thing what happens --> pseudo online
test_sus_pred = []; test_sus_label = [];
test_shift_pred = []; test_shift_label = [];
test_noth_pred = []; test_noth_label = [];
alpha = 0.02; threshold = 0.75;
for i = ntrial_train+1:ntrial
    c_data = squeeze(c_trial_data(:,:,i));

    [p_shift, p_shift_prob] = predict(qdaModel_shift, c_data(:,idx_features_shift));
    p_shift = p_shift == classes(2);
    [p_sus, p_sus_prob] = predict(qdaModel_sus, c_data(:,idx_features_sus));
    p_sus = p_sus == classes(2);
    [p_noth, p_noth_prob] = predict(qdaModel_noth, c_data(:,idx_features_noth));
    D = pdist2(c_data, C);
    [~, c_clusterLabels] = min(D, [], 2);

    current_state = ones(1, min_trial_data + 1) * 0.5;

    % pseudo online + accuracy all trial of all three classifiers
    for sample = 1:min_trial_data
        if mappedLabels(c_clusterLabels(sample)) == "Shift"
            % pseudo online part
            current_state(sample + 1) = current_state(sample)*(1-alpha) + alpha*p_shift(sample);

            % real time classification performances
            test_shift_pred = [test_shift_pred; p_shift(sample)];
            test_shift_label = [test_shift_label; trial_typ(i)];
        elseif mappedLabels(c_clusterLabels(sample)) == "Sustained"
            %pseudo online part
            current_state(sample + 1) = current_state(sample)*(1-alpha) + alpha*p_sus(sample);

            % real time classification performances
            test_sus_pred = [test_sus_pred; p_sus(sample)];
            test_sus_label = [test_sus_label; trial_typ(i)];
        else
            current_state(sample + 1) = current_state(sample);

            test_noth_pred = [test_noth_pred; p_noth(sample)];
            test_noth_label = [test_noth_label; trial_typ(i)];
        end
    end

    % plotting
    figure();
    subplot(6,1,1)
    imagesc(c_data')
    title('log band')

    subplot(6,1,2)
    plot(c_clusterLabels)
    ylim([0.5, K + 0.5]);
    yticks(1:K);
    yticklabels(mappedLabels)
    xlim([0 min_trial_data]);
    title('clustering')

    subplot(6,1,3)
    plot(p_shift)
    yticks(0:1); yticklabels([{'730'}, {'731'}])
    xlim([0 min_trial_data]);
    title('shift')

    subplot(6,1,4)
    plot(p_sus);
    yticks(0:1); yticklabels([{'730'}, {'731'}])
    xlim([0 min_trial_data]);
    title('sustained')

    subplot(6,1,5)
    plot(p_noth);
    yticks(0:1); yticklabels([{'730'}, {'731'}])
    xlim([0 min_trial_data]);
    title('nothing')

    subplot(6,1,6)
    plot(current_state);
    yticks(0:1); yticklabels([{'730'}, {'731'}])
    xlim([0 min_trial_data]);
    title('exponential')

    sgtitle(['asked: ' num2str(trial_typ(i))])
end
disp(['shift acc: ' num2str(sum(test_shift_label == test_shift_pred)/size(test_shift_pred, 1))])
disp(['sus acc: ' num2str(sum(test_sus_label == test_sus_pred)/size(test_sus_pred, 1))])
disp(['noth acc: ' num2str(sum(test_noth_label == test_noth_pred)/size(test_noth_pred, 1))])
disp(['In total: ' num2str((sum(test_shift_label == test_shift_pred) + sum(test_sus_label == test_sus_pred) + sum(test_noth_label == test_noth_pred) ) / (min_trial_data * ntrial_test))])




