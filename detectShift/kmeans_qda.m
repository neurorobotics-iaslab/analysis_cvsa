%% file to compute for each trial the gini value in a window of 100 ms
clear all; % close all;

addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/EOG')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/equal_ros')

%% Initialization
bands = [{[6 11]} {[8 14]} {[11 16]}];
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
avg = 1;% 0.75;
eog_threshold = 500;

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);

th_high = inf;

%% concatenate the files
nFiles = length(filenames);
trial_with_eog = [];
trial_with_high_voltage = [];
for idx_file= 1: nFiles
    fullpath_file_shift = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_shift);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;
    sampleRate = header.SampleRate;

    c_trial_with_eog = eog_detection(c_signal, header, eog_threshold, {'FP1', 'FP2', 'EOG'});
    trial_with_eog = [trial_with_eog; c_trial_with_eog];

    disp('   [proc] power band');
    for idx_band = 1:nbands
        band = bands{idx_band};

        % for power band using hilbert transformation
        bufferSize = floor(avg*sampleRate);
        chunkSize = 32;
        [signal_processed, header_processed] = processing_onlineROS_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize);

        if all(subject == 'h8')
            signal_processed(:,23) = 0;
        end

        if all(band == [8 14])
            c_trial_with_high_voltage = check_voltage_trial(signal_processed, header_processed, th_high);
            trial_with_high_voltage = [trial_with_high_voltage; c_trial_with_high_voltage];
        end

        c_header = headers{1, idx_band};
        c_header.sampleRate = header_processed.SampleRate/chunkSize;
        c_header.channels_labels = header_processed.Label;
        if isempty(find(header_processed.EVENT.TYP == 2, 1)) % no eye calibration
            c_header.TYP = cat(1, c_header.TYP, header_processed.EVENT.TYP);
            c_header.DUR = cat(1, c_header.DUR, header_processed.EVENT.DUR);
            c_header.POS = cat(1, c_header.POS, header_processed.EVENT.POS + size(signals{1, idx_band}, 1));
        else
            k = find(header_processed.EVENT.TYP == 1, 1);
            c_header.TYP = cat(1, c_header.TYP, header_processed.EVENT.TYP(k:end));
            c_header.DUR = cat(1, c_header.DUR, header_processed.EVENT.DUR(k:end));
            c_header.POS = cat(1, c_header.POS, header_processed.EVENT.POS(k:end) + size(signals{1, idx_band}, 1));
        end
        signals{1, idx_band} = cat(1, signals{1, idx_band}, signal_processed(:,:));
        headers{1, idx_band} = c_header;
    end
end


%% Labelling data 
events = headers{1,1};
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

% ncols = 5; nrows = 2; nfigures = ceil(ntrial / (ncols * nrows));
% t = 1;
% for idx_f = 1:nfigures
%     handles = [];
%     figure();
%     for idx_pos = 1:(ncols * nrows)
%         subplot(nrows, ncols, idx_pos)
%         imagesc(squeeze(trial_data(min_durCUE+min_durFIX+1:end,2,:,t))');
%         yticks(1:nchannels)
%         yticklabels(channels_label(1:end-1))
%         xticks(sampleRate:sampleRate:size(trial_data, 1))
%         xticklabels(string(((sampleRate:sampleRate:size(trial_data, 1)) / sampleRate) + ceil((min_durCUE+min_durFIX)/sampleRate)));
%         title(['trial only cf: ' num2str(t) ' class: ' num2str(trial_typ(t))])
%         colorbar;
%         handles = [handles; gca];
%         t = t + 1;
%     end
%     set(handles, 'clim', [0 500])
% end
%% refactoring the data
% now the data are placed alternatevely, so odd trial id is for 730, even for 731
trial_with_eog = trial_with_high_voltage | trial_with_eog;
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
trial_data = tmp; % samples x bands x channels x trials
trial_data(:,:,[1, 2, 19],:) = 0; % remove the power of the EOG channel, FP1 anf FP2 --> also in sparsity

%% compute sparsity
% define regions
occipital = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'}; [~, ch_occipital] = ismember(occipital, channels_label);
central = {'CP3', 'CP1', 'C3', 'FC3', 'C1', 'FC1', 'CP4', 'C4', 'FC4', 'CP2', 'C2', 'FC2', 'FCZ', 'CZ'}; [~, ch_central] = ismember(central, channels_label);
frontal = {'F3', 'F1', 'F2', 'F4', 'FP1', 'FP2', 'FZ'}; [~, ch_frontal] = ismember(frontal, channels_label);
nchannelsSelected = size(ch_occipital, 2);

nsparsity = 3;
sparsity = nan(min_trial_data, nbands, ntrial, nsparsity); % sample x band x trial x sparsity

for c = 1:ntrial
    c_data = squeeze(trial_data(:,:,:,c)); % samples x band x channels

    for sample = 1:min_trial_data
        c_sample = squeeze(c_data(sample,:,:)); % bands x channels

        for idx_band = 1:nbands
            tmp = squeeze(c_sample(idx_band,:)); % 1 x channels

            [sparsity(sample, idx_band, c,:), label_sparsity] = compute_features_kmeans(tmp);
        end
    end
end

%% show log band and index ---> all trial
handles = cell(1, nbands); cl = -inf(1, nbands);
for c = 1:ntrial
    figure();
    for idx_band = 1:nbands
        subplot(2,nbands,idx_band)
        imagesc(squeeze(trial_data(:,idx_band,:,c))')
        hold on;
        xline(minDurFix, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(minDurCue+minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        xticks(sampleRate:sampleRate:min_trial_data)
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        yticks(1:nchannels); yticklabels(channels_label)
        handles{idx_band} = [handles{idx_band}, gca];
        cl(idx_band) = max(cl(idx_band), ...
            max(abs(squeeze(trial_data(:, idx_band, :, c))), [], 'all'));
        title(['log band | ' bands_str{idx_band}])


        subplot(2,nbands, idx_band + nbands)
        plot(squeeze(sparsity(:,idx_band, c, 1)))
        hold on;
        for i = 2:nsparsity
            plot(squeeze(sparsity(:,idx_band, c, i)))
        end
        xline(minDurFix, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(minDurCue+minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        legend([label_sparsity,{'events'}])
        xticks(sampleRate:sampleRate:min_trial_data)
        xlim([0 min_trial_data])
        ylim([0, 1])
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        title(['sparsity (hoyer) | ' bands_str{idx_band}])
        
    end
    sgtitle(['task: ' num2str(trial_typ(c)) ' | trial ' num2str(c)])
end
for idx_band = 1:nbands
    set(handles{idx_band}, 'clim', [0, cl(idx_band)])
end

%% ----------------- KMEANS -----------------
K = 2;
choosen_band = 2;
percentual_train = 0.75; 
ntrial_train = floor((percentual_train * ntrial) / 2);
ntrial_train = ntrial_train * 2; % in this way even number for ntrial_train

sparsity_cf = squeeze(sparsity(minDurFix+minDurCue+1:end, choosen_band,:,:));

% z-score -> train and use the mu and var also for the test
train_data_3D = sparsity_cf(:, 1:ntrial_train, :);
train_data_2D = reshape(train_data_3D, size(train_data_3D, 1) * size(train_data_3D,2), size(train_data_3D,3));
mu_features = mean(train_data_2D, 1);
sigma_features = std(train_data_2D, 0, 1);
sigma_features(sigma_features == 0) = eps;
train_data_standardized_2D = (train_data_2D - mu_features) ./ sigma_features;
train_data_standardized_3D = reshape(train_data_standardized_2D, size(train_data_3D, 1), ntrial_train, size(train_data_3D,3));

test_data_3D = sparsity_cf(:, ntrial_train+1:end, :);
test_data_2D = reshape(test_data_3D, size(test_data_3D, 1) * (ntrial - ntrial_train), size(train_data_3D,3));
test_data_standardized_2D = (test_data_2D - mu_features) ./ sigma_features;
test_data_standardized_3D = reshape(test_data_standardized_2D, size(test_data_3D, 1), (ntrial - ntrial_train), size(train_data_3D,3));
sparsity_cf = cat(2, train_data_standardized_3D, test_data_standardized_3D);

disp('Esecuzione di K-means sui dati di training globali...');

[idx_train, C] = kmeans(train_data_standardized_2D, K, ...
    'Distance', 'sqeuclidean', ...
    'Replicates', 10, ...  % 'Replicates' è fondamentale per un risultato robusto
    'Display', 'final');

[~, sortIdx] = sort(C(:, 2), 'descend'); % Ordina in base alla feature 2 (OFI)
C = C(sortIdx, :);

% compute the labels for all the trials, with the selected centroids
train_kmeans = nan(ntrial_train * size(sparsity_cf, 1), nsparsity);
for c = 1:ntrial_train
    train_kmeans((c-1)*size(sparsity_cf, 1) + 1: c * size(sparsity_cf, 1),:) = sparsity_cf(:,c,:);
end
distances = pdist2(train_kmeans, C);  
[~, raw_labels] = min(distances, [], 2);

% save the cluster for each data
cluster_labels_train = nan(size(sparsity_cf, 1), ntrial_train);
for c = 1:ntrial_train
    cluster_labels_train(:,c) = raw_labels((c-1)*size(sparsity_cf, 1) + 1: c * size(sparsity_cf, 1));
end

% plot the C
disp('centroids: ')
disp(C)
disp('centroids diff features: ')
disp(C(1,:) - C(2,:))
%% ----------------- PLOT TRIALS -----------------
for c = 1:10
    figure();
    
    % --- log band ---
    subplot(3,1,1)
    imagesc(squeeze(trial_data(minDurFix+minDurCue+1:end,choosen_band,:,c))')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    yticks(1:nchannels); yticklabels(channels_label)
    title('log band')

    % --- sparsity indices ---
    subplot(3,1,2)
    plot(squeeze(sparsity_cf(:, c, 1)))
    hold on
    for i = 2:nsparsity
        plot(squeeze(sparsity_cf(:, c, i)))
    end
    hold off
    legend(label_sparsity)
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    xlim([1 size(sparsity_cf, 1)])
%     ylim([0, 1])
    title('sparsity (entropy)')

    % --- cluster labels k-means ---
    subplot(3,1,3)
    plot(squeeze(cluster_labels_train(:, c)), 'b')
    legend('kmeans')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    xlim([1 size(sparsity_cf, 1)])
    yticks(1:K);
    yticklabels([{'IC'}, {'NIC'}])   
    title('cluster comparison')
    
    sgtitle(['Trial: ' num2str(c) ' | Task: ' num2str(trial_typ(c)) ' | band: '  bands_str{choosen_band}])
end

%% test data kmeans
test_kmeans = nan((ntrial - ntrial_train) * size(sparsity_cf, 1), nsparsity);
for c = ntrial_train+1:ntrial
    test_kmeans((c-ntrial_train-1)*size(sparsity_cf, 1) + 1: (c-ntrial_train) * size(sparsity_cf, 1),:) = sparsity_cf(:,c,:);
end

distances = pdist2(test_kmeans, C);  
[~, raw_labels] = min(distances, [], 2);

cluster_labels_test = nan(size(sparsity_cf, 1), ntrial - ntrial_train);
for c = 1:ntrial-ntrial_train
    cluster_labels_test(:,c) = raw_labels((c - 1)*size(sparsity_cf, 1) + 1: c * size(sparsity_cf, 1));
end

%% extract train-test data
IC_train_data = []; IC_train_labels = [];
IC_test_data = []; IC_test_labels = [];

cl_train_data = []; cl_train_labels = [];
cl_test_data = []; cl_test_labels = [];

for c = 1:ntrial
    c_data = squeeze(trial_data(minDurCue+minDurFix+1:end,choosen_band,:,c)); % sample x channels

    if c <= ntrial_train
        for sample = 1:size(c_data,1)
            if cluster_labels_train(sample, c) == 1 % IC
                IC_train_data = [IC_train_data; c_data(sample,:)];
                IC_train_labels = [IC_train_labels; trial_typ(c)];
            end
        end
        cl_train_data = [cl_train_data; c_data];
        cl_train_labels = [cl_train_labels; repmat(trial_typ(c), size(c_data, 1), 1)];
    else
        for sample = 1:size(c_data,1)
            if cluster_labels_test(sample, c-ntrial_train) == 1 %"IC"
                IC_test_data = [IC_test_data; c_data(sample,:)];
                IC_test_labels = [IC_test_labels; trial_typ(c)];
            end
        end
        cl_test_data = [cl_test_data; c_data];
        cl_test_labels = [cl_test_labels; repmat(trial_typ(c), size(c_data, 1), 1)];
    end    
end

% balance the samples classes
if sum(IC_train_labels == classes(1)) < sum(IC_train_labels == classes(2))
    tmp = find(IC_train_labels == classes(2));
    IC_train_data(tmp(sum(IC_train_labels == classes(1)) + 1 :end),:) = [];
    IC_train_labels(tmp(sum(IC_train_labels == classes(1)) + 1 :end)) = [];
else
    tmp = find(IC_train_labels == classes(1));
    IC_train_data(tmp(sum(IC_train_labels == classes(2)) + 1 :end),:) = [];
    IC_train_labels(tmp(sum(IC_train_labels == classes(2)) + 1 :end)) = [];
end


% apply the log
IC_train_data = log(IC_train_data);
IC_test_data  = log(IC_test_data);
cl_train_data = log(cl_train_data);
cl_test_data  = log(cl_test_data);

%% fisher score on the train data
fisher_IC = nan(1, nchannels);
fisher_cl = nan(1, nchannels);

for idx_ch=1:nchannels-1
    % IC
    mu1 = mean(IC_train_data(IC_train_labels == classes(1),idx_ch));
    sigma1 = std(IC_train_data(IC_train_labels == classes(1),idx_ch));
    mu2 = mean(IC_train_data(IC_train_labels == classes(2),idx_ch));
    sigma2 = std(IC_train_data(IC_train_labels == classes(2),idx_ch));
    fisher_IC(idx_ch) = abs(mu1 - mu2);

    % cl
    mu1 = mean(cl_train_data(cl_train_labels == classes(1),idx_ch));
    sigma1 = std(cl_train_data(cl_train_labels == classes(1),idx_ch));
    mu2 = mean(cl_train_data(cl_train_labels == classes(2),idx_ch));
    sigma2 = std(cl_train_data(cl_train_labels == classes(2),idx_ch));
    fisher_cl(idx_ch) = abs(mu1 - mu2);
end

tmp = [fisher_IC; fisher_cl];
x_labels = ["IC", "classical"];
figure();
imagesc(tmp')
yticks(1:nchannels); yticklabels(channels_label)
xticks(1:4); xticklabels(x_labels)

%% train and test the qda
IC_select_channels = {}; %%%%% ----> features selection
cl_select_channels = {};

% features selection
if isempty(IC_select_channels) 
    IC_idx = ch_occipital;
else
    [~, IC_idx] = ismember(IC_select_channels, channels_label);
end
if isempty(cl_select_channels)
    cl_idx = ch_occipital;
else
    [~, cl_idx] = ismember(cl_select_channels, channels_label);
end

disp('------------------- QDA -------------------')
disp('--- IC ---')
qda_IC = fitcdiscr(IC_train_data(:,IC_idx), IC_train_labels, 'DiscrimType', 'quadratic');
CVModel = crossval(qda_IC, 'KFold', 5);
loss = kfoldLoss(CVModel);  % Average classification error
fprintf("5-fold cross-validation loss: %.4f\n", loss);
[p_IC_train, IC_train_score] = predict(qda_IC, IC_train_data(:,IC_idx));
confusionmat(IC_train_labels, p_IC_train)
disp(['Accuracy TRAIN: ' num2str(sum(p_IC_train == IC_train_labels)/size(p_IC_train,1))])
p = predict(qda_IC, IC_test_data(:,IC_idx));
confusionmat(IC_test_labels, p)
disp(['Accuracy TEST: ' num2str(sum(p == IC_test_labels)/size(p,1))])

disp('--- classical ---')
qda_cl = fitcdiscr(cl_train_data(:,cl_idx), cl_train_labels, 'DiscrimType', 'quadratic');
CVModel = crossval(qda_cl, 'KFold', 5);
loss = kfoldLoss(CVModel);  % Average classification error
fprintf("5-fold cross-validation loss: %.4f\n", loss);
[p_cl_train, cl_train_score] = predict(qda_cl, cl_train_data(:,cl_idx));
confusionmat(cl_train_labels, p_cl_train)
disp(['Accuracy TRAIN: ' num2str(sum(p_cl_train == cl_train_labels)/size(p_cl_train,1))])
p = predict(qda_cl, cl_test_data(:,cl_idx));
confusionmat(cl_test_labels, p)
disp(['Accuracy TEST: ' num2str(sum(p == cl_test_labels)/size(p,1))])

% plot the qda
IC_score1 = IC_train_score(IC_train_labels == classes(1),1);
IC_score2 = IC_train_score(IC_train_labels == classes(2),1);
[IC_x1,IC_f1] = ksdensity(IC_score1);
[IC_x2,IC_f2] = ksdensity(IC_score2);

cl_score1 = cl_train_score(cl_train_labels == classes(1),1);
cl_score2 = cl_train_score(cl_train_labels == classes(2),1);
[cl_x1,cl_f1] = ksdensity(cl_score1);
[cl_x2,cl_f2] = ksdensity(cl_score2);

figure();
subplot(121)
hold on;
grid on;
plot(IC_f1,IC_x1,'LineWidth',2);
plot(IC_f2,IC_x2,'LineWidth',2);
hold off;
title('IC qda in training');
legend(num2str(classes(1)),num2str(classes(2)));

subplot(122)
hold on;
grid on;
plot(cl_f1,cl_x1,'LineWidth',2);
plot(cl_f2,cl_x2,'LineWidth',2);
hold off;
title('cl qda in training');
legend(num2str(classes(1)),num2str(classes(2)));

%% show the trial test with the output of the classifiers and pseudo online
alpha = 0.96;
th = 0.8;
hit_trial_acc_onlineKmeans = nan(1, ntrial - ntrial_train);
hit_trial_acc_NOonlineKmeans = nan(1, ntrial - ntrial_train);
hit_trial_acc_classical = nan(1, ntrial - ntrial_train);

for c = ntrial_train+1:ntrial
    c_data = log(squeeze(trial_data(minDurFix+minDurCue+1:end, choosen_band,:,c))); % sample x channels

    [~, prob_IC_qda, ~] = predict(qda_IC, c_data(:,IC_idx));
    [~, prob_IC_lda, ~] = predict(lda_IC, c_data(:,IC_idx));
    [prob_IC_gmm, ~, ~] = predictGMM(IC_GMM_c1, IC_GMM_c2, c_data(:,IC_idx), classes);

    [~, prob_cl_qda, ~] = predict(qda_cl, c_data(:,cl_idx));
    [~, prob_cl_lda, ~] = predict(lda_cl, c_data(:,cl_idx));
    [prob_cl_gmm, ~, ~] = predictGMM(cl_GMM_c1, cl_GMM_c2, c_data(:,cl_idx), classes);
    

    figure();
    subplot(4,3,[1, 2, 3])
    imagesc(squeeze(trial_data(minDurFix+minDurCue+1:end,choosen_band,:,c))') % log bnd for the cf
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    yticks(1:nchannels); yticklabels(channels_label)
    title(['log band | ' bands_str{choosen_band}])

    subplot(4,3,4)
    if c <= ntrial_train
        c_cluster = squeeze(cluster_labels_train(:, c));
    else
        c_cluster = squeeze(cluster_labels_test(:, c-ntrial_train));
    end
    plot(c_cluster);
    legend('cluster')
    xlim([0 size(sparsity_cf, 1)])
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    title(['cluster | ' bands_str{choosen_band}])
    yticks(1:K);
    yticklabels({'IC', 'NIC'})

    % plot the prob of the classifiers
    subplot(4,3,7) % plot the prob for onlin kmeans
    plot(prob_IC_qda(:,1))
    hold on;
    yline(0.5, 'r--', 'LineWidth', 2);
    hold off;
    xlim([0 size(sparsity_cf, 1)])
    legend('qda')
    title('probabilities IC classifier')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    subplot(4,3,8) % plot the prob for NO online kmeans
    plot(prob_IC_qda(:,1))
    hold on;
    yline(0.5, 'r--', 'LineWidth', 2);
    hold off;
    xlim([0 size(sparsity_cf, 1)])
    legend('qda')
    title('probabilities IC classifier')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    subplot(4,3,9) % plot the prob for classcal approach
    plot(prob_cl_qda(:,1))
    hold on;
    yline(0.5, 'r--', 'LineWidth', 2);
    hold off;
    xlim([0 size(sparsity_cf, 1)])
    legend('qda')
    title('probabilities cl classifier')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));

    % plot the expo according to which classifier and the online kmeans
    subplot(4,3,10) 
    expo_qda = ones(size(c_cluster, 1)+1,1) * 0.5;
    for sample= 1:size(c_cluster,1)
        if c_cluster(sample) == 1 % IC
            expo_qda(sample+1) = expo_qda(sample) * alpha + (1-alpha)*prob_IC_qda(sample, 1);

            scores = [expo_qda(sample+1)];
            hit_trial_acc_onlineKmeans = check_threshold(trial_typ(c), c-ntrial_train, hit_trial_acc_onlineKmeans, scores, th);
        else
            expo_qda(sample+1) = expo_qda(sample);
        end
    end
    plot(expo_qda(2:end))
    hold on;
    yline(th, 'r--', 'LineWidth', 2);
    yline(1-th, 'r--', 'LineWidth', 2);
    hold off;
    xlim([0 size(sparsity_cf, 1)])
    title('exponential system | online kmeans')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    legend({'qda'})
    ylim([0 1])

    % plot the expo according to which classifier and the NO online kmeans
    subplot(4,3,11) 
    expo_qda = ones(size(c_cluster, 1)+1,1) * 0.5;
    for sample= 1:size(c_cluster,1)
        expo_qda(sample+1) = expo_qda(sample) * alpha + (1-alpha)*prob_IC_qda(sample, 1);

        scores = [expo_qda(sample+1)];
        hit_trial_acc_NOonlineKmeans = check_threshold(trial_typ(c), c-ntrial_train, hit_trial_acc_NOonlineKmeans, scores, th);
    end
    plot(expo_qda(2:end))
    hold on;
    yline(th, 'r--', 'LineWidth', 2);
    yline(1-th, 'r--', 'LineWidth', 2);
    hold off;
    xlim([0 size(sparsity_cf, 1)])
    title('exponential system | NO online kmeans')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    legend({'qda'})
    ylim([0 1])

    % plot the expo classical approach
    subplot(4,3,12) 
    expo_qda = ones(size(c_cluster, 1)+1,1) * 0.5;
    for sample= 1:size(c_cluster,1)
        expo_qda(sample+1) = expo_qda(sample) * alpha + (1-alpha)*prob_cl_qda(sample, 1);

        scores = [expo_qda(sample+1)];
        hit_trial_acc_classical = check_threshold(trial_typ(c), c-ntrial_train, hit_trial_acc_classical, scores, th);
    end
    plot(expo_qda(2:end))
    hold on;
    yline(th, 'r--', 'LineWidth', 2);
    yline(1-th, 'r--', 'LineWidth', 2);
    hold off;
    xlim([0 size(sparsity_cf, 1)])
    title('exponential system | classical')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    legend({'qda'})
    ylim([0 1])

    if c <= ntrial_train
        sgtitle(['trial ' num2str(c) ' | class asked: ' num2str(trial_typ(c)) ' | TRAIN'])
    else
        sgtitle(['trial ' num2str(c) ' | class asked: ' num2str(trial_typ(c)) ' | TEST'])
    end
end

% compute the accuracy
acc_hit = [nansum(hit_trial_acc_onlineKmeans, 2), nansum(hit_trial_acc_NOonlineKmeans, 2), nansum(hit_trial_acc_classical, 2)];
acc_miss = [nansum(1 - hit_trial_acc_onlineKmeans, 2), nansum(1 - hit_trial_acc_NOonlineKmeans, 2), nansum(1 - hit_trial_acc_classical, 2)];
acc_timeout = [sum(isnan(hit_trial_acc_onlineKmeans), 2), sum(isnan(hit_trial_acc_NOonlineKmeans), 2), sum(isnan(hit_trial_acc_classical), 2)];

classifiers_name = {'qda'};
res_order = {'online kmeans', 'NO online kmeans', 'classical'};
for j = 1:length(res_order)
    disp([res_order{j} ':']);
    for i = 1:length(classifiers_name)
        disp([classifiers_name{i} ' hit: ' num2str(acc_hit(i,j)) ' miss: ' num2str(acc_miss(i,j)) ' timeout: ' num2str(acc_timeout(i,j)) ' ntrial: ' num2str(ntrial- ntrial_train)])
        disp([classifiers_name{i} ' accuracy: ' num2str(acc_hit(i,j)/(ntrial-ntrial_train))])
    end
    fprintf('\n')
end


%%
function [prob_C1, prob_C2, predicted_class] = predictGMM(GMM_C1, GMM_C2, data, classes)
    % data deve essere un singolo campione [1 x N_canali] o 
    % una matrice [N_campioni x N_canali]
    likelihood_C1 = pdf(GMM_C1, data);
    likelihood_C2 = pdf(GMM_C2, data);

    total_likelihood = likelihood_C1 + likelihood_C2;
    total_likelihood(total_likelihood == 0) = 1; 

    prob_C1 = likelihood_C1 ./ total_likelihood;
    prob_C2 = likelihood_C2 ./ total_likelihood;

    predicted_class = zeros(size(prob_C1)); % Inizializza
    predicted_class(prob_C1 >= 0.5) = classes(1);
    predicted_class(prob_C1 < 0.5)  = classes(2);

end

function hit_trial_acc = check_threshold(trial_type, c_ntrial, hit_trial_acc, scores, th) % scores is qda, lda, svm, lr
    for idx_score = 1:length(scores)
        if isnan(hit_trial_acc(idx_score,c_ntrial))
            c_score = scores(idx_score);
            if trial_type == 730
                if c_score >= th
                    hit_trial_acc(idx_score, c_ntrial) = 1;
                elseif c_score <= 1-th
                    hit_trial_acc(idx_score, c_ntrial) = 0;
                end
            else
                if c_score >= th
                    hit_trial_acc(idx_score, c_ntrial) = 0;
                elseif c_score <= 1-th
                    hit_trial_acc(idx_score, c_ntrial) = 1;
                end
            end
        end
    end
end

% sparsity and where
function [sparsity, label_sparsity] = compute_features_kmeans(c_signal)
    % take the one which contribute at the 95% of the energy
    sparsity = nan(3,1);
    label_sparsity = [{'LI'},{'GI'},{'GB'}];
    o_l = sort([29 13 30 37 33 34 17]); o_r = sort([31 15 32 35 36 38 18]);
    frontal = sort([3 4 5 20 21]); c_l = sort([6 22 25 8 11 27]); c_r = sort([7 24 26 10 12 28]);

    % --- LAP --- Calcola il LAP per tutti i punti nella finestra passata (window_signal)
    % show the occipital lateralization that is strongand present during the CVSA
    P_left_window  = mean(c_signal(o_l));
    P_right_window = mean(c_signal(o_r));
    LAP_history = (P_right_window - P_left_window) ./ (P_right_window + P_left_window + eps);
    sparsity(1) = abs(LAP_history); % LAP_Mean

    % --- Gini Index + Occipital Power ---  -> when high there is a zone stronger, so IC
    % show the focusing is weighted in with the power in the occipital part, in this way ig CVSA then strong value 
    non_zeros_chs = setdiff(1:size(c_signal,2), [1, 2, 19]);
    global_mean = mean(c_signal(non_zeros_chs)); % car filter
    current_signal_normalized = c_signal - global_mean; % remove the global energy
    mean_roi = [mean(current_signal_normalized(frontal)), mean(current_signal_normalized(c_l)), ...
        mean(current_signal_normalized(c_r)), mean(current_signal_normalized(o_l)), ...
        mean(current_signal_normalized(o_r))];
    mean_roi = abs(mean_roi); % make sure the energy is positive--> we are using peak and valli with same significance
    mean_roi_ordered = sort(mean_roi);
    n = length(mean_roi_ordered);
    sum_roi_p = 0;
    for i = 1:n
        sum_roi_p = sum_roi_p + (n+1-i) * mean_roi_ordered(i);
    end
    total_sum = sum(mean_roi_ordered);
    if total_sum > 0
        gi = (1/n) * (n+1-2*sum_roi_p/total_sum);
    else
        gi = 0;
    end
    % compute the weight factor
    pot_occipital = mean_roi(4) + mean_roi(5);
    pot_total_roi = sum(mean_roi); % Somma di F, CL, CR, OL, OR
    if pot_total_roi > 0
        occipital_power = pot_occipital / pot_total_roi;
    else
        occipital_power = 0; 
    end
    sparsity(2) = occipital_power * gi;

    % --- GB ---
    % return the global power mean, 
    sparsity(3) = global_mean;
end