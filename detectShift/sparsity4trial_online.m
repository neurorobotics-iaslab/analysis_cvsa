%% file to compute for each trial the gini value in a window of 100 ms
clear all; % close all;

addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/EOG')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/equal_ros')

%% Initialization
bands = [{[14 22]} {[8 14]} {[22 30]}];
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
avg = 0.75;
eog_threshold = 500;

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);
if all(subject == 'c7')
    th_high = 500;
elseif all(subject == 'g2')
    th_high = 170;
elseif all(subject == 'g3')
    th_high = 3500;
elseif all(subject == 'f2')
    th_high = 150;
elseif all(subject == 'j2')
    th_high = 500;
elseif all(subject == 'd7')
    th_high = 350;
elseif all(subject == 'c8')
    th_high = 75;
elseif all(subject == 'h8')
    th_high = 75;
elseif all(subject == 'i7')
    th_high = 150;
else 
    th_high = 50;
end

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

sparsity = nan(min_trial_data, nbands, ntrial, 3); % sample x band x trial x sparsity

for c = 1:ntrial
    c_data = squeeze(trial_data(:,:,:,c)); % samples x band x channels

    for sample = 1:min_trial_data
        c_sample = squeeze(c_data(sample,:,:)); % bands x channels

        for idx_band = 1:nbands
            tmp = squeeze(c_sample(idx_band,:)); % 1 x channels
%             tmp = (tmp - min(tmp)); tmp = tmp / max(tmp); % normalize

            sparsity(sample, idx_band, c,:) = compute_sparsity(tmp);
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
        plot(squeeze(sparsity(:,idx_band, c, 2)))
        plot(squeeze(sparsity(:,idx_band, c, 3)))
        xline(minDurFix, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(minDurCue+minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        legend([{'bp mean'},{'bp max'},{'bp mean 2'},{'events'}])
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

% %% train data for the kmean (using only ban 8-14) --> all trial together
% --- rejected since there is the possibility to have trial with very high
% value so at the end the centroids are too far from them ---
% 
% train_kmeans = nan(ntrial_train * size(sparsity_cf, 1), 3);
% for t = 1:ntrial_train
%     train_kmeans((t-1)*size(sparsity_cf, 1) + 1: t * size(sparsity_cf, 1),:) = sparsity_cf(:,t,:);
% end
% K = 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [raw_labels, C] = kmeans(train_kmeans, K, 'Distance', 'sqeuclidean', 'Replicates', 10, 'Display', 'iter'); % K x p 
% 
% cluster_labels_train = nan(size(sparsity_cf, 1), ntrial_train);
% for t = 1:ntrial_train
%     cluster_labels_train(:,t) = raw_labels((t-1)*size(sparsity_cf, 1) + 1: t * size(sparsity_cf, 1));
% end
% 
% [~, sortIdx] = sort(vecnorm(C,2,2));  % Sort clusters from low to high power
% C = C(sortIdx,:);
% if K == 3
%     labels = ["Nothing", "Shift", "Sustained"];
% else
%     labels = ["Shift", "Sustained"];
% end
% 
% figure(); hold on;
% for i = 1:K
%     scatter(C(i,1), C(i,2), 'filled')
% end
% xlabel('mean'); ylabel('max');
% title('centroids with all the data together')
% legend(labels)

%% train data fro kmeans using 8-14 --> one kmeans for each trial and then keep the median/mean
K = 2;
choosen_band = 2;
percentual_train = 0.75; ntrial_train = ceil(percentual_train * ntrial);
sparsity_cf = squeeze(sparsity(minDurFix+minDurCue+1:end, choosen_band,:,:));
centroids = nan(ntrial_train,K,3);

for c = 1:ntrial_train
    c_data = squeeze(sparsity_cf(:,c,:));
    [~, c_C] = kmeans(c_data, K, 'Distance', 'sqeuclidean', 'Replicates', 10); %'Display', 'iter');

    [~, sortIdx] = sort(vecnorm(c_C,2,2), 'descend');  % sort centroids by norm 

    centroids(c,:,:) = c_C(sortIdx,:);
end

C = nan(K, size(sparsity_cf, 3));
for i = 1:K
    C(i,:) = median(centroids(:,i,:));
end

train_kmeans = nan(ntrial_train * size(sparsity_cf, 1), 3);
for c = 1:ntrial_train
    train_kmeans((c-1)*size(sparsity_cf, 1) + 1: c * size(sparsity_cf, 1),:) = sparsity_cf(:,c,:);
end

distances = pdist2(train_kmeans, C);  
[~, raw_labels] = min(distances, [], 2);
cluster_labels_train = nan(size(sparsity_cf, 1), ntrial_train);
for c = 1:ntrial_train
    cluster_labels_train(:,c) = raw_labels((c-1)*size(sparsity_cf, 1) + 1: c * size(sparsity_cf, 1));
end

if K == 3
    labels = ["Sustained", "Shift", "Nothing"];
else
    labels = ["Sustained", "Shift"];
end

figure(); hold on;
for i = 1:K
    scatter3(C(i,1), C(i,2), C(i,3), 'filled')
end
xlabel('mean'); ylabel('max'); zlabel('mean 2')
title('centroids with all the data together')
legend(labels)

% save the cluster label for the trian
% save('/home/paolo/cvsa_ws/src/analysis_cvsa/detectShift/kmeans_matlab.mat', 'raw_labels', 'sparsity_cf', 'K', 'C', 'ntrial_train')

%% show log band, index and cluster --> only cf
handles = []; cl = -inf;
for c = 1:ntrial_train % ntrial_train
    figure();
    subplot(4,1,1)
    imagesc(squeeze(trial_data(minDurFix+minDurCue+1:end,choosen_band,:,c))')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    yticks(1:nchannels); yticklabels(channels_label)
    handles = [handles, gca];
    cl = max(cl, ...
        max(abs(squeeze(trial_data(:, choosen_band, :, c))), [], 'all'));
    title(['log band | ' bands_str{choosen_band}])

    subplot(4,1,2)
    plot(squeeze(sparsity_cf(:, c, 1)))
    legend('band power')
    hold on
    plot(squeeze(sparsity_cf(:, c, 2)))
    plot(squeeze(sparsity_cf(:, c, 3)))
    hold off
    xticks(0:sampleRate:min_trial_data)
    xlim([0 size(sparsity_cf, 1)])
    ylim([0, 1])
    legend('bp mean', 'bp max', 'bp mean 2')
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    title(['sparsity (entropy) | ' bands_str{choosen_band}])

    subplot(4,1,3)
    plot(squeeze(cluster_labels_train(:, c)))
    legend('cluster')
    xticks(0:sampleRate:min_trial_data)
    xlim([0 size(sparsity_cf, 1)])
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    title(['cluster all trial together | ' bands_str{choosen_band}])
    yticks(1:K);
    yticklabels(labels)   

    subplot(414)
    scatter(sparsity_cf(:,c,1), sparsity_cf(:,c,3), 'filled')
    title('sparsity plot')

    sgtitle(['only CF | task: ' num2str(trial_typ(c)) ' | trial: ' num2str(c)])
end
% set(handles, 'clim', [0, cl])

%% test data kmeans
test_kmeans = nan((ntrial - ntrial_train) * size(sparsity_cf, 1), 3);
for c = ntrial_train+1:ntrial
    test_kmeans((c-ntrial_train-1)*size(sparsity_cf, 1) + 1: (c-ntrial_train) * size(sparsity_cf, 1),:) = sparsity_cf(:,c,:);
end

distances = pdist2(test_kmeans, C);  
[~, raw_labels] = min(distances, [], 2);

cluster_labels_test = nan(size(sparsity_cf, 1), ntrial - ntrial_train);
for c = 1:ntrial-ntrial_train
    cluster_labels_test(:,c) = raw_labels((c - 1)*size(sparsity_cf, 1) + 1: c * size(sparsity_cf, 1));
end

%% extract train-test qda
sus_train_data = []; sus_train_labels = [];
sus_test_data = []; sus_test_labels = [];

shi_train_data = []; shi_train_labels = [];
shi_test_data = []; shi_test_labels = [];

not_train_data = []; not_train_labels = [];
not_test_data = []; not_test_labels = [];

for c = 1:ntrial
    c_data = squeeze(trial_data(minDurCue+minDurFix+1:end,choosen_band,:,c)); % sample x channels

    if c <= ntrial_train
        for sample = 1:size(c_data,1)
            if cluster_labels_train(sample, c) == 2 % shift
                shi_train_data = [shi_train_data; c_data(sample,:)];
                shi_train_labels = [shi_train_labels; trial_typ(c)];
            elseif cluster_labels_train(sample, c) == 1 % sustained
                sus_train_data = [sus_train_data; c_data(sample,:)];
                sus_train_labels = [sus_train_labels; trial_typ(c)];
            elseif cluster_labels_train(sample, c) == 3 % nothing
                not_train_data = [not_train_data; c_data(sample,:)];
                not_train_labels = [not_train_labels; trial_typ(c)];
            end
        end
    else
        for sample = 1:size(c_data,1)
            if cluster_labels_test(sample, c-ntrial_train) == 2 %"Shift"
                shi_test_data = [shi_test_data; c_data(sample,:)];
                shi_test_labels = [shi_test_labels; trial_typ(c)];
            elseif cluster_labels_test(sample, c-ntrial_train) == 1 %"Sustained"
                sus_test_data = [sus_test_data; c_data(sample,:)];
                sus_test_labels = [sus_test_labels; trial_typ(c)];
            elseif cluster_labels_test(sample, c-ntrial_train) == 3 %"Nothing"
                not_test_data = [not_test_data; c_data(sample,:)];
                not_test_labels = [not_test_labels; trial_typ(c)];
            end
        end
    end    
end

% balance the samples classes
if sum(sus_train_labels == classes(1)) < sum(sus_train_labels == classes(2))
    tmp = find(sus_train_labels == classes(2));
    sus_train_data(tmp(sum(sus_train_labels == classes(1)) + 1 :end),:) = [];
    sus_train_labels(tmp(sum(sus_train_labels == classes(1)) + 1 :end)) = [];
else
    tmp = find(sus_train_labels == classes(1));
    sus_train_data(tmp(sum(sus_train_labels == classes(2)) + 1 :end),:) = [];
    sus_train_labels(tmp(sum(sus_train_labels == classes(2)) + 1 :end)) = [];
end
if sum(shi_train_labels == classes(1)) < sum(shi_train_labels == classes(2))
    tmp = find(shi_train_labels == classes(2));
    shi_train_data(tmp(sum(shi_train_labels == classes(1)) + 1 :end),:) = [];
    shi_train_labels(tmp(sum(shi_train_labels == classes(1)) + 1 :end)) = [];
else
    tmp = find(shi_train_labels == classes(1));
    shi_train_data(tmp(sum(shi_train_labels == classes(2)) + 1 :end),:) = [];
    shi_train_labels(tmp(sum(shi_train_labels == classes(2)) + 1 :end)) = [];
end

%% fisher score on the train data
fisher_shift = nan(1, nchannels);
fisher_sus = nan(1, nchannels);
if K == 3
    fisher_noth = nan(1, nchannels);
end
for idx_ch=1:nchannels-1
    % shift
    mu1 = mean(shi_train_data(shi_train_labels == classes(1),idx_ch));
    sigma1 = std(shi_train_data(shi_train_labels == classes(1),idx_ch));
    mu2 = mean(shi_train_data(shi_train_labels == classes(2),idx_ch));
    sigma2 = std(shi_train_data(shi_train_labels == classes(2),idx_ch));
    fisher_shift(idx_ch) = abs(mu1 - mu2);

    % sus
    mu1 = mean(sus_train_data(sus_train_labels == classes(1),idx_ch));
    sigma1 = std(sus_train_data(sus_train_labels == classes(1),idx_ch));
    mu2 = mean(sus_train_data(sus_train_labels == classes(2),idx_ch));
    sigma2 = std(sus_train_data(sus_train_labels == classes(2),idx_ch));
    fisher_sus(idx_ch) = abs(mu1 - mu2);

    % noth
    if K == 3
        mu1 = mean(not_train_data(not_train_labels == classes(1),idx_ch));
        sigma1 = std(not_train_data(not_train_labels == classes(1),idx_ch));
        mu2 = mean(not_train_data(not_train_labels == classes(2),idx_ch));
        sigma2 = std(not_train_data(not_train_labels == classes(2),idx_ch));
        fisher_noth(idx_ch) = abs(mu1 - mu2);
    end
end

if K == 3
    tmp = [fisher_shift; fisher_sus; fisher_noth];
    x_labels = ["shift", "sustained", "nothing"];
else
    tmp = [fisher_shift; fisher_sus];
    x_labels = ["shift", "sustained"];
end
figure();
imagesc(tmp')
yticks(1:nchannels); yticklabels(channels_label)
xticks(1:3); xticklabels(x_labels)

%% train and test the qda
sus_select_channels = {'PO7', 'PO3', 'PO5', 'O1'};
shi_select_channels = {};
not_select_channels = {};

% sustained
if isempty(sus_select_channels)
    % remove eog channel  --> all 0 values
    sus_test_data(:,all(sus_test_data==0,1)) = [];
    sus_train_data(:,all(sus_train_data==0,1)) = [];
    sus_idx = 1:size(sus_train_data, 2);
else
    [~, sus_idx] = ismember(sus_select_channels, channels_label);
end
disp('--- SUSTAINED ---')
qdaModel_sus = fitcdiscr(sus_train_data(:,sus_idx), sus_train_labels, 'DiscrimType', 'quadratic');
CVModel = crossval(qdaModel_sus, 'KFold', 5);
loss = kfoldLoss(CVModel);  % Average classification error
fprintf("5-fold cross-validation loss: %.4f\n", loss);
[p_sus_train, sus_train_score] = predict(qdaModel_sus, sus_train_data(:,sus_idx));
confusionmat(sus_train_labels, p_sus_train)
disp(['Accuracy TRAIN: ' num2str(sum(p_sus_train == sus_train_labels)/size(p_sus_train,1))])
p = predict(qdaModel_sus, sus_test_data(:,sus_idx));
confusionmat(sus_test_labels, p)
disp(['Accuracy TEST: ' num2str(sum(p == sus_test_labels)/size(p,1))])

% shift
if isempty(shi_select_channels)
    % remove eog channel  --> all 0 values
    shi_test_data(:,all(shi_test_data==0,1)) = [];
    shi_train_data(:,all(shi_train_data==0,1)) = [];
    shi_idx = 1:size(shi_train_data, 2);
else
    [~, shi_idx] = ismember(shi_select_channels, channels_label);
end
disp('--- SHIFT ---')
qdaModel_shi = fitcdiscr(shi_train_data(:,shi_idx), shi_train_labels, 'DiscrimType', 'quadratic');
CVModel = crossval(qdaModel_shi, 'KFold', 5);
loss = kfoldLoss(CVModel);  % Average classification error
fprintf("5-fold cross-validation loss: %.4f\n", loss);
[p_shi_train, shi_train_score] = predict(qdaModel_shi, shi_train_data(:,shi_idx));
confusionmat(shi_train_labels, p_shi_train)
disp(['Accuracy TRAIN: ' num2str(sum(p_shi_train == shi_train_labels)/size(p_shi_train,1))])
p = predict(qdaModel_shi, shi_test_data(:,shi_idx));
confusionmat(shi_test_labels, p)
disp(['Accuracy TEST: ' num2str(sum(p == shi_test_labels)/size(p,1))])

if K == 3
    if isempty(not_select_channels)
        % remove eog channel  --> all 0 values
        not_test_data(:,all(not_test_data==0,1)) = [];
        not_train_data(:,all(not_train_data==0,1)) = [];
        not_idx = 1:size(not_train_data, 2);
    else
        [~, not_idx] = ismember(not_select_channels, channels_label);
    end
    disp('--- NOTHING ---')
    qdaModel_not = fitcdiscr(not_train_data(:,not_idx), not_train_labels, 'DiscrimType', 'quadratic');
    CVModel = crossval(qdaModel_not, 'KFold', 5);
    loss = kfoldLoss(CVModel);  % Average classification error
    fprintf("5-fold cross-validation loss: %.4f\n", loss);
    p = predict(qdaModel_not, not_train_data(:,not_idx));
    confusionmat(not_train_labels, p)
    disp(['Accuracy TRAIN: ' num2str(sum(p == not_train_labels)/size(p,1))])
    p = predict(qdaModel_not, not_test_data(:,not_idx));
    confusionmat(not_test_labels, p)
    disp(['Accuracy TEST: ' num2str(sum(p == not_test_labels)/size(p,1))])
end

% plot the qda
shi_score1 = shi_train_score(shi_train_labels == classes(1));
shi_score2 = shi_train_score(shi_train_labels == classes(2));
[shi_x1,shi_f1] = ksdensity(shi_score1);
[shi_x2,shi_f2] = ksdensity(shi_score2);

sus_score1 = sus_train_score(sus_train_labels == classes(1));
sus_score2 = sus_train_score(sus_train_labels == classes(2));
[sus_x1,sus_f1] = ksdensity(sus_score1);
[sus_x2,sus_f2] = ksdensity(sus_score2);

figure();
subplot(121)
hold on;
grid on;
plot(shi_f1,shi_x1,'LineWidth',2);
plot(shi_f2,shi_x2,'LineWidth',2);
hold off;
title('Shift qda in training');
legend(num2str(classes(1)),num2str(classes(2)));

subplot(122)
hold on;
grid on;
plot(sus_f1,sus_x1,'LineWidth',2);
plot(sus_f2,sus_x2,'LineWidth',2);
hold off;
title('Sustained qda in training');
legend(num2str(classes(1)),num2str(classes(2)));

% save
% save('/home/paolo/cvsa_ws/src/analysis_cvsa/detectShift/qdaSUS_matlab.mat', 'sus_train_data', 'sus_train_labels', 'sus_idx', 'sus_train_score', 'qdaModel_sus')
% save('/home/paolo/cvsa_ws/src/analysis_cvsa/detectShift/qdaSHI_matlab.mat', 'shi_train_data', 'shi_train_labels', 'shi_idx', 'shi_train_score', 'qdaModel_shi')

%% train on overall the data --------------------------------------------------------
% sustained
sus_all_data = [sus_train_data; sus_test_data];
if isempty(sus_select_channels)
    % remove eog channel  --> all 0 values
    sus_all_data(:,all(sus_all_data==0,1)) = [];
    sus_idx = 1:size(sus_all_data, 2);
else
    [~, sus_idx] = ismember(sus_select_channels, channels_label);
end
disp('--- SUSTAINED ---')
qdaModel_sus = fitcdiscr(sus_all_data(:,sus_idx), sus_train_labels, 'DiscrimType', 'quadratic');
CVModel = crossval(qdaModel_sus, 'KFold', 5);
loss = kfoldLoss(CVModel);  % Average classification error
fprintf("5-fold cross-validation loss: %.4f\n", loss);
[p_sus_train, sus_train_score] = predict(qdaModel_sus, sus_all_data(:,sus_idx));
confusionmat(sus_train_labels, p_sus_train)
disp(['Accuracy TRAIN: ' num2str(sum(p_sus_train == sus_train_labels)/size(p_sus_train,1))])
p = predict(qdaModel_sus, sus_test_data(:,sus_idx));
confusionmat(sus_test_labels, p)
disp(['Accuracy TEST: ' num2str(sum(p == sus_test_labels)/size(p,1))])

% shift
shi_all = [shi_train_data; shi_test_data];
if isempty(shi_select_channels)
    % remove eog channel  --> all 0 values
    shi_all(:,all(shi_all==0,1)) = [];
    shi_idx = 1:size(shi_all, 2);
else
    [~, shi_idx] = ismember(shi_select_channels, channels_label);
end
disp('--- SHIFT ---')
qdaModel_shi = fitcdiscr(shi_all(:,shi_idx), shi_train_labels, 'DiscrimType', 'quadratic');
CVModel = crossval(qdaModel_shi, 'KFold', 5);
loss = kfoldLoss(CVModel);  % Average classification error
fprintf("5-fold cross-validation loss: %.4f\n", loss);
[p_shi_train, shi_train_score] = predict(qdaModel_shi, shi_all(:,shi_idx));
confusionmat(shi_train_labels, p_shi_train)
disp(['Accuracy TRAIN: ' num2str(sum(p_shi_train == shi_train_labels)/size(p_shi_train,1))])
p = predict(qdaModel_shi, shi_test_data(:,shi_idx));
confusionmat(shi_test_labels, p)
disp(['Accuracy TEST: ' num2str(sum(p == shi_test_labels)/size(p,1))])

if K == 3
    not_all = [not_train_data; not_test_data];
    if isempty(not_select_channels)
        % remove eog channel  --> all 0 values
        not_all(:,all(not_all==0,1)) = [];
        not_idx = 1:size(not_train_data, 2);
    else
        [~, not_idx] = ismember(not_select_channels, channels_label);
    end
    disp('--- NOTHING ---')
    qdaModel_not = fitcdiscr(not_all(:,not_idx), not_train_labels, 'DiscrimType', 'quadratic');
    CVModel = crossval(qdaModel_not, 'KFold', 5);
    loss = kfoldLoss(CVModel);  % Average classification error
    fprintf("5-fold cross-validation loss: %.4f\n", loss);
    p = predict(qdaModel_not, not_all(:,not_idx));
    confusionmat(not_train_labels, p)
    disp(['Accuracy TRAIN: ' num2str(sum(p == not_train_labels)/size(p,1))])
    p = predict(qdaModel_not, not_test_data(:,not_idx));
    confusionmat(not_test_labels, p)
    disp(['Accuracy TEST: ' num2str(sum(p == not_test_labels)/size(p,1))])
end

% plot the qda
shi_score1 = shi_train_score(shi_train_labels == classes(1));
shi_score2 = shi_train_score(shi_train_labels == classes(2));
[shi_x1,shi_f1] = ksdensity(shi_score1);
[shi_x2,shi_f2] = ksdensity(shi_score2);

sus_score1 = sus_train_score(sus_train_labels == classes(1));
sus_score2 = sus_train_score(sus_train_labels == classes(2));
[sus_x1,sus_f1] = ksdensity(sus_score1);
[sus_x2,sus_f2] = ksdensity(sus_score2);

figure();
subplot(121)
hold on;
grid on;
plot(shi_f1,shi_x1,'LineWidth',2);
plot(shi_f2,shi_x2,'LineWidth',2);
hold off;
title('Shift qda in training');
legend(num2str(classes(1)),num2str(classes(2)));

subplot(122)
hold on;
grid on;
plot(sus_f1,sus_x1,'LineWidth',2);
plot(sus_f2,sus_x2,'LineWidth',2);
hold off;
title('Sustained qda in training');
legend(num2str(classes(1)),num2str(classes(2)));

%% show the trial test with the output of the classifiers and pseudo online
alpha_sus = 0.85;
alpha_shi = 0.97;
for c = ntrial_train+1:ntrial
    c_data = squeeze(trial_data(minDurFix+minDurCue+1:end, choosen_band,:,c)); % sample x channels

    [~, prob_sus, ~] = predict(qdaModel_sus, c_data(:,sus_idx));
    [~, prob_shi, ~] = predict(qdaModel_shi, c_data(:,shi_idx));
    if K == 3
        [~, prob_not, ~] = predict(qdaModel_not, c_data(:,not_idx));
        legend_label = ["sustained", "shift", "nothing"];
    else
        legend_label = ["sustained", "shift"];
    end

    figure();
    subplot(4,1,1)
    imagesc(squeeze(trial_data(minDurFix+minDurCue+1:end,choosen_band,:,c))') % log bnd for the cf
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    yticks(1:nchannels); yticklabels(channels_label)
    title(['log band | ' bands_str{choosen_band}])

    subplot(4,1,2)
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
    yticklabels(labels)

    subplot(413)
    plot(prob_sus(:,1))
    hold on;
    plot(prob_shi(:,1))
    if K == 3
        plot(prob_not(:,1));
    end
    yline(0.5, 'r--', 'LineWidth', 2);
    hold off;
    xlim([0 size(sparsity_cf, 1)])
    legend(legend_label)
    title('probabilities')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));

    subplot(414)
    expo = ones(size(c_cluster, 1)+1,1) * 0.5;
    for sample= 1:size(c_cluster,1)
        if c_cluster(sample) == 1 % sustained
            expo(sample+1) = expo(sample) * alpha_sus + (1-alpha_sus)*prob_sus(sample, 1);
        elseif c_cluster(sample) == 10 % shift
            expo(sample+1) = expo(sample) * alpha_shi + (1-alpha_shi)*prob_shi(sample, 1);
        else % nothing
            expo(sample+1) = expo(sample);
        end
    end
    plot(expo(2:end))
    hold on;
    yline(0.9, 'r--', 'LineWidth', 2);
    yline(0.1, 'r--', 'LineWidth', 2);
    hold off;
    xlim([0 size(sparsity_cf, 1)])
    title('exponential system')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));

    if c <= ntrial_train
    sgtitle(['trial ' num2str(c) ' | class asked: ' num2str(trial_typ(c)) ' | TRAIN'])
    else
    sgtitle(['trial ' num2str(c) ' | class asked: ' num2str(trial_typ(c)) ' | TEST'])
    end

end



%%
c = 5;
prova = nan(min_trial_data, 3);
c_data = squeeze(trial_data(:,choosen_band,:,c)); % samples x channels

for sample = 1:min_trial_data
    c_sample = squeeze(c_data(sample,:)); % bands x channels

    %             tmp = (tmp - min(tmp)); tmp = tmp / max(tmp); % normalize

    prova(sample, :) = compute_sparsity(c_sample);
end
figure();
subplot(211)
imagesc(squeeze(trial_data(minDurFix+minDurCue+1:end,choosen_band,:,c))')

subplot(212)
plot(prova(minDurFix+minDurCue+1:end,:))



% sparsity and where
function sparsity = compute_sparsity(signalVec)
    % take the one which contribute at the 95% of the energy
    sparsity = nan(3,1);
    activeIdx = 1:length(signalVec);

    % --- compute the sparsity over 5 different ROIs ---
    frontal = [3 4 5 20 21]; c_l = [6 22 25 8 11 27]; c_r = [7 24 26 10 12 28]; o_l = [29 13 30 37 33 34 17]; o_r = [31 15 32 35 36 38 18];
    mean_roi = [mean(signalVec(activeIdx(ismember(activeIdx,frontal)))), mean(signalVec(activeIdx(ismember(activeIdx,c_l)))), ...
        mean(signalVec(activeIdx(ismember(activeIdx,c_r)))), mean(signalVec(activeIdx(ismember(activeIdx,o_l)))), ...
        mean(signalVec(activeIdx(ismember(activeIdx,o_r))))];
    mean_roi(isnan(mean_roi)) = 0;
    n = length(mean_roi);
    signalVec_entropy = mean_roi.^2;
    total_entropy = sum(signalVec_entropy);
    if total_entropy == 0
        sparsity(1) = 0;
    else
        p = signalVec_entropy / total_entropy;
        if (p(4) >= 1/n) || (p(5) >= 1/n)
            H = -sum(p .* log2(p + eps));  % entropy
            Hmax = log2(n);
            sparsity(1) = 1 - (H / Hmax);  % normalized: 0 to 1
        else 
            sparsity(1) = 0;
        end
    end

    max_roi = [max(signalVec(activeIdx(ismember(activeIdx,frontal)))), max(signalVec(activeIdx(ismember(activeIdx,c_l)))), ...
        max(signalVec(activeIdx(ismember(activeIdx,c_r)))), max(signalVec(activeIdx(ismember(activeIdx,o_l)))), ...
        max(signalVec(activeIdx(ismember(activeIdx,o_r))))];
    signalVec_entropy = max_roi.^2;
    total_entropy = sum(signalVec_entropy);
    if total_entropy == 0
        sparsity(2) = 0;
    else
        p = signalVec_entropy / total_entropy;
        if (p(4) >= 1/n) || (p(5) >= 1/n)
            H = -sum(p .* log2(p + eps));  % entropy
            Hmax = log2(n);
            sparsity(2) = 1 - (H / Hmax);  % normalized: 0 to 1
        else
            sparsity(2) = 0;
        end
    end

    frontal = [3 4 5 20 21]; central = [c_l, c_r, 9, 23]; occipital = [o_l, o_r, 14, 16, 39];
    mean_roi = [mean(signalVec(activeIdx(ismember(activeIdx,frontal)))), mean(signalVec(activeIdx(ismember(activeIdx,central)))), ...
        mean(signalVec(activeIdx(ismember(activeIdx,occipital))))];
    mean_roi(isnan(mean_roi)) = 0;
    n = length(mean_roi);
    signalVec_entropy = mean_roi.^2;
    total_entropy = sum(signalVec_entropy);
    if total_entropy == 0
        sparsity(3) = 0;
    else
        p = signalVec_entropy / total_entropy;
        if p(3) >= 1/n
            H = -sum(p .* log2(p + eps));  % entropy
            Hmax = log2(n);
            sparsity(3) = 1 - (H / Hmax);  % normalized: 0 to 1
        else
            sparsity(3) = 0;
        end
    end
        
%     l1 = norm(max_roi, 1);      % Norma L1
%     l2 = norm(max_roi, 2);      % Norma L2
%     if l2 == 0
%         sparsity(2) = 0;  % caso limite: vettore nullo
%     else
%         sparsity(2) = (sqrt(n) - (l1 / l2)) / (sqrt(n) - 1);
%     end
end