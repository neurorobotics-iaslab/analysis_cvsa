clear all; % close all;

addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_cvsa/equal_ros')

%% Initialization
threshold_gmm_ic = 0.7;
bands = [{[8 13]} {[18 24]}];
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
nbands = length(bands);
signals = cell(1, nbands);
artifacts = cell(1, nbands);
headers = cell(1, nbands);
for idx_band = 1:nbands
    headers{idx_band}.TYP = [];
    headers{idx_band}.DUR = [];
    headers{idx_band}.POS = [];
    signals{idx_band} = [];
    artifacts{idx_band} = [];
end
classes = [769 770];      
nchannels = 16;
nclasses = length(classes);
filterOrder = 4;
avg = 1;% 0.75;
channels_label = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'CP3', 'CP1', 'CP2', 'CP4', 'Pz'};

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);

th_high = inf;

%% concatenate the files
nFiles = length(filenames);
for idx_file= 1: nFiles
    fullpath_file_shift = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_shift);
    c_signal = c_signal(:,1:nchannels);
    sampleRate = header.SampleRate;

    [~, roi_idx_L] = ismember(upper({'C1', 'C3'}), upper(channels_label));
    [~, roi_idx_R] = ismember(upper({'C2', 'C4'}), upper(channels_label));

    excl_chs = [];

    for idx_band = 1:nbands
        band = bands{idx_band};

        % for power band using hilbert transformation
        bufferSize = floor(avg*sampleRate);
        chunkSize = 32;
        eog.label = [];
        muscle.filterOrder = 4;
        muscle.freq = 1; % remove antneuro problems
        muscle.threshold = 100;
        [signal_processed, header_processed] = processing_onlineROS_CAR_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize, excl_chs);
        artifact = artifact_rejection(c_signal, header, nchannels, bufferSize, chunkSize, eog, muscle);

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
        artifacts{1, idx_band} = cat(1, artifacts{1, idx_band}, artifact(:,:));
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
trial_start = nan(ntrial, 1);
trial_end = nan(ntrial, 1);
trial_typ = nan(ntrial, 1);
for idx_trial = 1:ntrial
    trial_start(idx_trial) = fixPOS(idx_trial);
    trial_typ(idx_trial) = cueTYP(idx_trial);
    trial_end(idx_trial) = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;
end

min_trial_data = min(trial_end - trial_start+1);
trial_data = nan(min_trial_data, nbands, nchannels, ntrial); % data x bands x channels x trial
artifacts_data = nan(min_trial_data, nbands, ntrial); % data x bands x trial
for idx_band = 1:nbands
    c_signal = signals{idx_band};
    c_artifact = artifacts{idx_band};
    for trial = 1:ntrial
        c_start = trial_start(trial);
        c_end = trial_start(trial) + min_trial_data - 1;
        trial_data(:,idx_band,:,trial) = c_signal(c_start:c_end,:);
        artifacts_data(:,idx_band,trial) = c_artifact(c_start:c_end,:);
    end
end

%% refactoring the data --> odd trial class 1 even class 2
idx_classes_trial = nan(ntrial/2, nclasses);
for idx_class = 1:nclasses
    idx_classes_trial(:,idx_class) = find(trial_typ == classes(idx_class));
end

tmp_data = nan(size(trial_data));
tmp_art = nan(size(artifacts_data));
trial_typ = nan(size(trial_typ));
i = 1;
for idx_trial_class = 1:2:ntrial
    for idx_class = 1:nclasses
        tmp_data(:,:,:,idx_trial_class + idx_class - 1) = trial_data(:,:,:,idx_classes_trial(i, idx_class));
        tmp_art(:,:,idx_trial_class + idx_class - 1) = artifacts_data(:,:,idx_classes_trial(i, idx_class));
        trial_typ(idx_trial_class + idx_class - 1) = classes(idx_class);
    end
    i = i + 1;
end
trial_data = tmp_data; % samples x bands x channels x trials
artifacts_data = tmp_art; % data x bands x trial
% trial_data(:,:,[1, 2, 19],:) = 0; % remove the power of the EOG channel, FP1 anf FP2 --> also in sparsity

%% compute sparsity
% define regions
nsparsity = 2;
sparsity = nan(min_trial_data, nbands, ntrial, nsparsity); % sample x band x trial x sparsity
% o_l_ch = {'P3', 'O1', 'P5', 'P1', 'PO5', 'PO3', 'PO7'};
% o_r_ch = {'P4', 'O2', 'P2', 'P6', 'PO4', 'PO6', 'PO8'};
% c_l_ch = {'FC1', 'C3', 'CP1', 'FC3', 'C1', 'CP3'};
% c_r_ch = {'FC2', 'C4', 'CP2', 'FC4', 'C2', 'CP4'};

% o_l_ch = {'O1', 'PO5', 'PO3', 'PO7'};
% o_r_ch = {'O2', 'PO4', 'PO6', 'PO8'};
o_l_ch = {'FC1', 'FC3'};
o_r_ch = {'FC2', 'FC4'};
c_l_ch = {'C3', 'C1'};
c_r_ch = {'C4', 'C2'};

[~, o_l] = ismember(o_l_ch, channels_label);
[~, o_r] = ismember(o_r_ch, channels_label);
[~, c_l] = ismember(c_l_ch, channels_label);
[~, c_r] = ismember(c_r_ch, channels_label);

for idx_trial = 1:ntrial
    c_data = squeeze(trial_data(:,:,:,idx_trial)); % samples x band x channels

    for sample = 1:min_trial_data
        c_sample = squeeze(c_data(sample,:,:)); % bands x channels

        for idx_band = 1:nbands
            tmp = squeeze(c_sample(idx_band,:)); % 1 x channels

            [feat_vals, label_sparsity] = compute_features_icnic(tmp, 'mi', o_l, o_r, c_l, c_r, nsparsity);
            
            % Salva tutte e 3 le feature
            sparsity(sample, idx_band, idx_trial, 1:nsparsity) = feat_vals;

        end
    end
end

%% show features
% sparsity_data = squeeze(sparsity(minDurFix+minDurCue+1:end, 1,:,:));
sparsity_data = squeeze(sparsity(:, 1,:,:));
% sparsity_data = squeeze(sparsity(1:minDurFix, 1,:,:));

data_3D = sparsity_data(:, 1:ntrial, :);
data_2D = reshape(data_3D, size(data_3D, 1) * size(data_3D,2), size(data_3D,3));

figure;
for i = 1:nsparsity
    for j = i+1:nsparsity
        subplot(nsparsity-1, nsparsity-1, (i-1)*(nsparsity-1)+j-i)
        scatter(data_2D(:,i), data_2D(:,j), 15, 'filled');
        xlabel(label_sparsity{i});
        ylabel(label_sparsity{j});
    end
end
sgtitle('Features comparison 2D');

figure
for i = 1:nsparsity
    subplot(nsparsity,1,i)
    histogram(data_2D(:, i));
    title(label_sparsity{i});
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
% for idx_band = 1:nbands
%     set(handles{idx_band}, 'clim', [0, cl(idx_band)])
% end

%% ----------------- GMM -----------------
choosen_band = 1;
percentual_train = 0.75; 
ntrial_train = floor((percentual_train * ntrial) / 2);
ntrial_train = ntrial_train * 2; % in this way even number for ntrial_train
K_range = 2:2; 
best_gmm = [];
min_bic = inf;
bics = [];

sparsity_cf = squeeze(sparsity(minDurFix+minDurCue+1:end, choosen_band,:,:));
artifacts_cf = squeeze(artifacts_data(minDurFix+minDurCue+1:end, choosen_band,:));
% sparsity_cf = squeeze(sparsity(:, choosen_band,:,:));
% artifacts_cf = squeeze(artifacts_data(:, choosen_band,:));

% z-score -> train and use the mu and var also for the test
train_data_3D = sparsity_cf(:, 1:ntrial_train, :);
train_data_2D = reshape(train_data_3D, size(train_data_3D, 1) * size(train_data_3D,2), size(train_data_3D,3));
artefact_train_2D = artifacts_cf(:, 1:ntrial_train);
artefact_train_1D = reshape(artefact_train_2D, size(artefact_train_2D, 1) * size(artefact_train_2D,2), 1);
train_data_2D_noArtif = train_data_2D(artefact_train_1D == 0,:);
mu_features = mean(train_data_2D_noArtif, 1);
sigma_features = std(train_data_2D_noArtif, 0, 1);
sigma_features(sigma_features == 0) = eps;
train_data_2D_noArtif = (train_data_2D_noArtif - mu_features) ./ sigma_features;
train_data_standardized_2D = (train_data_2D - mu_features) ./ sigma_features;
train_data_standardized_3D = reshape(train_data_standardized_2D, size(train_data_3D, 1), ntrial_train, size(train_data_3D,3));

test_data_3D = sparsity_cf(:, ntrial_train+1:end, :);
test_data_2D = reshape(test_data_3D, size(test_data_3D, 1) * (ntrial - ntrial_train), size(train_data_3D,3));
test_data_standardized_2D = (test_data_2D - mu_features) ./ sigma_features;
test_data_standardized_3D = reshape(test_data_standardized_2D, size(test_data_3D, 1), (ntrial - ntrial_train), size(train_data_3D,3));
sparsity_cf = cat(2, train_data_standardized_3D, test_data_standardized_3D);

disp('Esecuzione di GMM sui dati di training globali...');
options = statset('MaxIter', 1000, 'Display', 'off');

for k = K_range
    try
        % RegularizationValue = 1e-5 evita che le gaussiane collassino su un punto
        gmm_temp = fitgmdist(train_data_2D_noArtif, k, ...
                            'Options', options, ...
                            'CovarianceType', 'Full', ...
                            'SharedCovariance', false, ...
                            'RegularizationValue', 1e-5, ...
                            'Replicates', 25); 
        bics = [bics; gmm_temp.BIC];
        if gmm_temp.BIC < min_bic
            min_bic = gmm_temp.BIC;
            best_gmm = gmm_temp;
        end
    catch
        continue;
    end
end

gmm_model = best_gmm;
K = gmm_model.NumComponents;
disp(['GMM ottimizzato: K = ' num2str(K) ' (BIC = ' num2str(min_bic) ')']);

[~, sort_order] = sort(gmm_model.mu(:, 1), 'descend');
idx_ic = sort_order(1);  % Indice del cluster GMM per 'ic'
idx_nic = sort_order(2); % Indice del cluster GMM per 'nic'

% Crea le etichette finali
labels_gmm = {'NIC', 'IC'};

train_gmm = nan(ntrial_train * size(sparsity_cf, 1), nsparsity);
for c = 1:ntrial_train
    train_gmm((c-1)*size(sparsity_cf, 1) + 1: c * size(sparsity_cf, 1),:) = sparsity_cf(:,c,:);
end
P_soft = posterior(gmm_model, train_gmm);
cluster_labels_train = nan(size(sparsity_cf, 1), ntrial_train); % contains the prob to be ic
for c = 1:ntrial_train
    cluster_labels_train(:,c) = P_soft((c-1)*size(sparsity_cf, 1) + 1: c * size(sparsity_cf, 1),idx_ic);
end

fprintf('Mappatura: Cluster GMM %d -> "ic", Cluster GMM %d -> "nic"\n', idx_ic, idx_nic);

% plot the C
disp('centroids: ')
disp(gmm_model.mu)

%% --- VISUALIZZAZIONE GMM ---
figure('Color', 'w'); % Crea una figura con sfondo bianco
hold on;
scatter(train_data_2D_noArtif(:,1), train_data_2D_noArtif(:,2), 15, ...
        'MarkerFaceColor', [0.2 0.5 0.9], ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.4);

x_min = min(train_data_2D_noArtif(:,1)) - 1; x_max = max(train_data_2D_noArtif(:,1)) + 1;
y_min = min(train_data_2D_noArtif(:,2)) - 1; y_max = max(train_data_2D_noArtif(:,2)) + 1;
step = 0.05; 
[x1Grid, x2Grid] = meshgrid(x_min:step:x_max, y_min:step:y_max);
XGrid = [x1Grid(:), x2Grid(:)];

prob_GMM = pdf(gmm_model, XGrid);
prob_GMM = reshape(prob_GMM, size(x1Grid));
[C, h] = contour(x1Grid, x2Grid, prob_GMM, 10, 'LineWidth', 2, 'LineColor', [0.8 0.2 0.2]);
plot(gmm_model.mu(:,1), gmm_model.mu(:,2), 'k+', 'MarkerSize', 15, 'LineWidth', 3);

xlabel('Lateralization Index (Z-score)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Gini Index (Z-score)', 'FontSize', 12, 'FontWeight', 'bold');
title('GMM Fit: Cluster IC vs NIC', 'FontSize', 14);
legend({'Dati Reali', 'Ellissi GMM', 'Centroidi'}, 'Location', 'best');
grid on;
axis tight;
hold off;

% 1. Ottieni le etichette "Hard" dal GMM (assegna ogni punto al cluster più probabile)
cluster_idx = cluster(gmm_model, train_data_2D_noArtif);

% 2. Calcola la Silhouette
figure;
[s, h] = silhouette(train_data_2D_noArtif, cluster_idx);
mean_sil = mean(s);

disp(['Silhouette Score Medio: ' num2str(mean_sil)]);
% Salva il grafico per il paper
title(['Silhouette Plot (Score: ' num2str(mean_sil, '%.2f') ')']);

eva_ch = evalclusters(train_data_2D_noArtif, cluster_idx, 'CalinskiHarabasz');
disp(['Calinski-Harabasz Index: ' num2str(eva_ch.CriterionValues)]);

eva_db = evalclusters(train_data_2D_noArtif, cluster_idx, 'DaviesBouldin');
disp(['Davies-Bouldin Index: ' num2str(eva_db.CriterionValues)]);


%% ----------------- PLOT TRIALS -----------------
for c = 1:10
    figure();
    
    % --- log band ---
    subplot(4,1,1)
    imagesc(squeeze(trial_data(minDurFix+minDurCue+1:end,choosen_band,:,c))')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    yticks(1:nchannels); yticklabels(channels_label)
    title('log band')

    % --- sparsity indices ---
    subplot(4,1,2)
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
    title('sparsity')

    % --- cluster labels gmm ---
    subplot(4,1,3)
    plot(squeeze(cluster_labels_train(:, c)), 'r')
    hold on
    yline(threshold_gmm_ic, 'k--', 'LineWidth', 2);
    hold off
    legend('gmm', 'threshold gmm ic')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    xlim([1 size(sparsity_cf, 1)])
    yticks(0:K-1);
    yticklabels(labels_gmm)   
    title('cluster comparison')

    % --- plot if artefacts ---
    subplot(4,1,4)
    plot(squeeze(artifacts_data(minDurFix+minDurCue+1:end,choosen_band,c)), 'b')
    legend('artefact in the trial')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    xlim([1 size(sparsity_cf, 1)])
    yticks(0:1);
    ylim([0 1])
    yticklabels([{'no'}, {'yes'}])   
    title('artefacts')
    
    sgtitle(['Trial: ' num2str(c) ' | Task: ' num2str(trial_typ(c)) ' | band: '  bands_str{choosen_band}])
end

%% test data kmeans
test_gmm = nan((ntrial - ntrial_train) * size(sparsity_cf, 1), nsparsity);
for c = ntrial_train+1:ntrial
    test_gmm((c-ntrial_train-1)*size(sparsity_cf, 1) + 1: (c-ntrial_train) * size(sparsity_cf, 1),:) = sparsity_cf(:,c,:);
end

P_soft = posterior(gmm_model, test_gmm);

cluster_labels_test = nan(size(sparsity_cf, 1), ntrial - ntrial_train);
for c = 1:ntrial-ntrial_train
    cluster_labels_test(:,c) = P_soft((c - 1)*size(sparsity_cf, 1) + 1: c * size(sparsity_cf, 1), idx_ic);
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
            if artifacts_cf(sample,c) == 0
                if cluster_labels_train(sample, c) >= threshold_gmm_ic % IC
                    IC_train_data = [IC_train_data; c_data(sample,:)];
                    IC_train_labels = [IC_train_labels; trial_typ(c)];
                end
                cl_train_data = [cl_train_data; c_data(sample,:)];
                cl_train_labels = [cl_train_labels; trial_typ(c)];
            end
        end
    else
        for sample = 1:size(c_data,1)
            if artifacts_cf(sample,c) == 0
                if cluster_labels_test(sample, c-ntrial_train) >= threshold_gmm_ic %"IC"
                    IC_test_data = [IC_test_data; c_data(sample,:)];
                    IC_test_labels = [IC_test_labels; trial_typ(c)];
                end
                cl_test_data = [cl_test_data; c_data(sample,:)];
                cl_test_labels = [cl_test_labels; trial_typ(c)];
            end
        end

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
% occipital = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'}; 
occipital = channels_label;
[~, ch_occipital] = ismember(occipital, channels_label);
noccipital = size(ch_occipital, 2);

fisher_IC = nan(1, noccipital);
fisher_cl = nan(1, noccipital);

for idx_ch_occipital=1:noccipital
    idx_ch = ch_occipital(idx_ch_occipital);
    % IC
    mu1 = mean(IC_train_data(IC_train_labels == classes(1),idx_ch));
    sigma1 = std(IC_train_data(IC_train_labels == classes(1),idx_ch));
    mu2 = mean(IC_train_data(IC_train_labels == classes(2),idx_ch));
    sigma2 = std(IC_train_data(IC_train_labels == classes(2),idx_ch));
    fisher_IC(idx_ch_occipital) = abs(mu1 - mu2)^2 / (sigma1^2 + sigma2^2);

    % cl
    mu1 = mean(cl_train_data(cl_train_labels == classes(1),idx_ch));
    sigma1 = std(cl_train_data(cl_train_labels == classes(1),idx_ch));
    mu2 = mean(cl_train_data(cl_train_labels == classes(2),idx_ch));
    sigma2 = std(cl_train_data(cl_train_labels == classes(2),idx_ch));
    fisher_cl(idx_ch_occipital) = abs(mu1 - mu2)^2 / (sigma1^2 + sigma2^2);
end

tmp = [fisher_IC; fisher_cl];
x_labels = ["IC", "classical"];
figure();
imagesc(tmp')
title('gmm ic and classical fisher score')
colorbar;
yticks(1:noccipital); yticklabels(occipital)
xticks(1:4); xticklabels(x_labels)

% R^2
calc_r2_from_data(IC_train_data, IC_train_labels, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['QDA data | size data: ' num2str(size(IC_train_data,1))]);
calc_r2_from_data(cl_train_data, cl_train_labels, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['all data | size data: ' num2str(size(cl_train_data,1))]);

%% train and test the qda
IC_select_channels = {'C3', 'C4'}; %%%%% ----> features selection
cl_select_channels = {'C3', 'C4'};


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
hit_trial_acc_onlineGMM = nan(1, ntrial - ntrial_train);
hit_trial_acc_NOonlineGMM = nan(1, ntrial - ntrial_train);
hit_trial_acc_classical = nan(1, ntrial - ntrial_train);

for c = ntrial_train+1:ntrial
    c_data = log(squeeze(trial_data(minDurFix+minDurCue+1:end, choosen_band,:,c))); % sample x channels

    [~, prob_IC_qda, ~] = predict(qda_IC, c_data(:,IC_idx));

    [~, prob_cl_qda, ~] = predict(qda_cl, c_data(:,cl_idx));
    

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
    hold on
    yline(threshold_gmm_ic, 'k--', 'LineWidth', 2);
    hold off
    legend('cluster', 'threshold gmm ic')
    xlim([0 size(sparsity_cf, 1)])
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    title(['cluster | ' bands_str{choosen_band}])
    yticks(0:K-1);
    yticklabels(labels_gmm)

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
        if artifacts_cf(sample, c) == 0
            if c_cluster(sample) >= threshold_gmm_ic % IC
                expo_qda(sample+1) = expo_qda(sample) * alpha + (1-alpha)*prob_IC_qda(sample, 1);

                scores = [expo_qda(sample+1)];
                hit_trial_acc_onlineGMM = check_threshold(trial_typ(c), c-ntrial_train, hit_trial_acc_onlineGMM, scores, th);
            else
                expo_qda(sample+1) = expo_qda(sample);
            end
        else
            expo_qda(sample+1) = expo_qda(sample);
        end
    end
    plot(expo_qda(2:end))
    hold on;
    yline(th, 'r--', 'LineWidth', 2);
    yline(1-th, 'r--', 'LineWidth', 2);
    plot(artifacts_cf(:,c))
    hold off;
    xlim([0 size(sparsity_cf, 1)])
    title('exponential system | online kmeans')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    legend({'qda', 'artifact (0=no,1=yes)'})
    ylim([0 1])

    % plot the expo according to which classifier and the NO online kmeans
    subplot(4,3,11)
    expo_qda = ones(size(c_cluster, 1)+1,1) * 0.5;
    for sample= 1:size(c_cluster,1)
        if artifacts_cf(sample, c) == 0
            expo_qda(sample+1) = expo_qda(sample) * alpha + (1-alpha)*prob_IC_qda(sample, 1);

            scores = [expo_qda(sample+1)];
            hit_trial_acc_NOonlineGMM = check_threshold(trial_typ(c), c-ntrial_train, hit_trial_acc_NOonlineGMM, scores, th);
        else
            expo_qda(sample+1) = expo_qda(sample);
        end
    end
    plot(expo_qda(2:end))
    hold on;
    yline(th, 'r--', 'LineWidth', 2);
    yline(1-th, 'r--', 'LineWidth', 2);
    plot(artifacts_cf(:,c))
    hold off;
    xlim([0 size(sparsity_cf, 1)])
    title('exponential system | NO online kmeans')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    legend({'qda', 'artifact (0=no,1=yes)'})
    ylim([0 1])

    % plot the expo classical approach
    subplot(4,3,12)
    expo_qda = ones(size(c_cluster, 1)+1,1) * 0.5;
    for sample= 1:size(c_cluster,1)
        if artifacts_cf(sample, c) == 0
            expo_qda(sample+1) = expo_qda(sample) * alpha + (1-alpha)*prob_cl_qda(sample, 1);

            scores = [expo_qda(sample+1)];
            hit_trial_acc_classical = check_threshold(trial_typ(c), c-ntrial_train, hit_trial_acc_classical, scores, th);
        else
            expo_qda(sample+1) = expo_qda(sample);
        end
    end
    plot(expo_qda(2:end))
    hold on;
    yline(th, 'r--', 'LineWidth', 2);
    yline(1-th, 'r--', 'LineWidth', 2);
    plot(artifacts_cf(:,c))
    hold off;
    xlim([0 size(sparsity_cf, 1)])
    title('exponential system | classical')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue+1:sampleRate:min_trial_data) / sampleRate));
    legend({'qda', 'artifact (0=no,1=yes)'})
    ylim([0 1])

    if c <= ntrial_train
        sgtitle(['trial ' num2str(c) ' | class asked: ' num2str(trial_typ(c)) ' | TRAIN'])
    else
        sgtitle(['trial ' num2str(c) ' | class asked: ' num2str(trial_typ(c)) ' | TEST'])
    end
end

% compute the accuracy
acc_hit = [nansum(hit_trial_acc_onlineGMM, 2), nansum(hit_trial_acc_NOonlineGMM, 2), nansum(hit_trial_acc_classical, 2)];
acc_miss = [nansum(1 - hit_trial_acc_onlineGMM, 2), nansum(1 - hit_trial_acc_NOonlineGMM, 2), nansum(1 - hit_trial_acc_classical, 2)];
acc_timeout = [sum(isnan(hit_trial_acc_onlineGMM), 2), sum(isnan(hit_trial_acc_NOonlineGMM), 2), sum(isnan(hit_trial_acc_classical), 2)];

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