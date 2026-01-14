clear all; % close all;

addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/equal_ros')
addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/utils')

%% Initialization
DATAPAH = '/home/paolo/cvsa/ic_cvsa_ws/src/';
classes = [730 731];      
nchannels = 39;
nclasses = length(classes);
filterOrder = 4;
avg = 1;% 0.75;
threshold_gmm_ic = 0.7;

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);
time_str = datestr(now, 'ddmmyyyy_HHMMSS');
gmm_file = ['gmm_' subject '_' time_str '.yaml'];
save_path_gmm = [DATAPAH, 'gmm_cvsa/cfg/' gmm_file];
save_path_qda_dataset = [DATAPAH 'qda_cvsa/create_qda/datasets/gmm/data_' subject '_' time_str '.mat'];

%% understand the band
nFiles = length(filenames);
peaks = zeros(1, nFiles);
for idx_file = 1:nFiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    peaks(idx_file) = analyze_alpha_peak(fullpath_file, 'RestTrigger', 786, 'band', [8 14], ...
        'target_regions', {'O1', 'O2', 'OZ', 'PO7', 'PO8', 'PO3', 'PO4', 'PO5', 'PO6', 'POZ'});
end

%% start processing data
bands = [{[8 14]}];
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

for idx_file= 1: nFiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;
    sampleRate = header.SampleRate;

    excl_ch = {'FP1', 'FP2', 'EOG'};
    [~, excl_chs] = ismember(excl_ch, channels_label);

    disp('   [proc] power band');
    for idx_band = 1:nbands
        band = bands{idx_band};

        % for power band using hilbert transformation and artefact remotion -----------------------------------------------
        bufferSize = floor(avg*sampleRate);
        chunkSize = 32;
        eog.filterOrder = 4;
        eog.band = [1 7];
        eog.label = excl_ch;
        eog.h_threshold = 60;
        eog.v_threshold = 60;
        picks.filterOrder = 4;
        picks.freq = 1; % remove antneuro problems
        picks.threshold = 100;
        [signal_processed, header_processed] = processing_onlineROS_CAR_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize, excl_chs);
        artifact = artifact_rejection(c_signal, header, nchannels, bufferSize, chunkSize, eog, picks);

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
if sum(trial_typ == classes(1)) == sum(trial_typ == classes(2))
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
    artifacts_data = tmp_art;
end

%% compute sparsity
% define regions
nsparsity = 3;
sparsity = nan(min_trial_data, nbands, ntrial, nsparsity); % sample x band x trial x sparsity
% o_l_ch = {'P3', 'O1', 'P5', 'P1', 'PO5', 'PO3', 'PO7'};
% o_r_ch = {'P4', 'O2', 'P2', 'P6', 'PO4', 'PO6', 'PO8'};
% c_l_ch = {'FC1', 'C3', 'CP1', 'FC3', 'C1', 'CP3'};
% c_r_ch = {'FC2', 'C4', 'CP2', 'FC4', 'C2', 'CP4'};

o_l_ch = {'O1', 'PO5', 'PO3', 'PO7'};
o_r_ch = {'O2', 'PO4', 'PO6', 'PO8'};
c_l_ch = {'C3', 'CP1', 'C1', 'CP3'};
c_r_ch = {'C4', 'CP2', 'C2', 'CP4'};

[~, o_l] = ismember(o_l_ch, channels_label);
[~, o_r] = ismember(o_r_ch, channels_label);
[~, c_l] = ismember(c_l_ch, channels_label);
[~, c_r] = ismember(c_r_ch, channels_label);

type = 'cvsa';

for c = 1:ntrial
    c_data = squeeze(trial_data(:,:,:,c)); % samples x band x channels

    for sample = 1:min_trial_data
        c_sample = squeeze(c_data(sample,:,:)); % bands x channels

        for idx_band = 1:nbands
            tmp = squeeze(c_sample(idx_band,:)); % 1 x channels

            [sparsity(sample, idx_band, c,:), features_name] = compute_features_icnic(tmp, type, o_l, o_r, c_l, c_r, nsparsity);
        end
    end
end

%% ----------------- gmm -----------------
% update to work with subbands -> tesista
choosen_band = 1;
K_range = 2:2; 
best_gmm = [];
min_bic = inf;

sparsity_cf = squeeze(sparsity(minDurFix+minDurCue+1:end, choosen_band,:,:));
artifacts_cf = squeeze(artifacts_data(minDurFix+minDurCue+1:end, choosen_band,:));

% z-score -> train and use the mu and var also for the test
percentual_training = 0.75;
ntrial_train = 2 * round((ntrial * percentual_training) / 2);
data_3D = sparsity_cf(:, 1:ntrial_train, :);
data_2D = reshape(data_3D, size(data_3D, 1) * size(data_3D,2), size(data_3D,3));
artefact_2D = artifacts_cf(:, 1:ntrial_train);
artefact_1D = reshape(artefact_2D, size(artefact_2D, 1) * size(artefact_2D,2), 1);
data_2D_noArtif = data_2D(artefact_1D == 0,:);
mu_features = mean(data_2D_noArtif, 1);
sigma_features = std(data_2D_noArtif, 0, 1);
sigma_features(sigma_features == 0) = eps;
data_2D_noArtif = (data_2D_noArtif - mu_features) ./ sigma_features;

% extract data all trials
sparsity_cuecf = squeeze(sparsity(minDurFix+1:end, choosen_band,:,:));
data_2D = reshape(sparsity_cuecf, size(sparsity_cuecf, 1) * size(sparsity_cuecf,2), size(sparsity_cuecf,3));
data_standardized_2D = (data_2D - mu_features) ./ sigma_features;
sparsity_cuecf = reshape(data_standardized_2D, size(sparsity_cuecf, 1), ntrial, size(sparsity_cuecf,3));
sparsity_cf = sparsity_cuecf(minDurCue+1:end,:,:);


disp('Esecuzione di GMM sui dati di training globali...');
options = statset('MaxIter', 1000, 'Display', 'off');
for k = K_range
    try
        % RegularizationValue = 1e-5 evita che le gaussiane collassino su un punto
        gmm_temp = fitgmdist(data_2D_noArtif, k, ...
                            'Options', options, ...
                            'CovarianceType', 'Full', ...
                            'SharedCovariance', false, ...
                            'RegularizationValue', 1e-5, ...
                            'Replicates', 150); 
        
        if gmm_temp.BIC < min_bic
            min_bic = gmm_temp.BIC;
            best_gmm = gmm_temp;
        end
    catch ME
        fprintf(2, 'Errore durante il fit con k=%d:\n%s\n', k, ME.message);
        continue;
    end
end

gmm_model = best_gmm;
K = gmm_model.NumComponents;
disp(['GMM ottimizzato: K = ' num2str(K) ' (BIC = ' num2str(min_bic) ')']);

[~, sort_order] = sort(gmm_model.mu(:, 1), 'descend');
idx_nic = sort_order(2); 
idx_ic = sort_order(1); % strong lateralization = IC
classes_icnic = zeros(1,2);
classes_icnic(idx_ic) = 1;

% Crea le etichette finali
labels_gmm = {'NIC', 'IC'};

tmp_data_cuecf_gmm = nan(ntrial * size(sparsity_cuecf, 1), nsparsity);
for c = 1:ntrial
    tmp_data_cuecf_gmm((c-1)*size(sparsity_cuecf, 1) + 1: c * size(sparsity_cuecf, 1),:) = sparsity_cuecf(:,c,:);
end
P_soft = posterior(gmm_model, tmp_data_cuecf_gmm);
cluster_labels_cuecf = nan(size(sparsity_cuecf, 1), ntrial); % contains the prob to be ic
for c = 1:ntrial
    cluster_labels_cuecf(:,c) = P_soft((c-1)*size(sparsity_cuecf, 1) + 1: c * size(sparsity_cuecf, 1),idx_ic);
end
cluster_labels_cf = cluster_labels_cuecf(minDurCue+1:end,:);

fprintf('Mappatura: Cluster GMM %d -> "ic", Cluster GMM %d -> "nic"\n', idx_ic, idx_nic);

% plot the C
disp('centroids: ')
disp(gmm_model.mu)

%% --- VISUALIZZAZIONE GMM ---
cluster_idx = cluster(gmm_model, data_2D_noArtif);
colors = lines(K); 

% --- Scatter Plot 3D Interattivo ---
figure('Color', 'w', 'Name', 'GMM 3D Clustering'); 
hold on; grid on; rotate3d on;

view(3);
for k = 1:K
    idx_k = (cluster_idx == k);
    scatter3(data_2D_noArtif(idx_k, 1), data_2D_noArtif(idx_k, 2), data_2D_noArtif(idx_k, 3), ...
             20, colors(k,:), 'filled', 'MarkerFaceAlpha', 0.5);
end
plot3(gmm_model.mu(:,1), gmm_model.mu(:,2), gmm_model.mu(:,3), ...
      'k+', 'MarkerSize', 15, 'LineWidth', 3);

xlabel('Feature 1 (Z-score)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Feature 2 (Z-score)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Feature 3 (Z-score)', 'FontSize', 12, 'FontWeight', 'bold');
title(['GMM Fit 3D: K = ' num2str(K)], 'FontSize', 14);
legend({'Cluster 1', 'Cluster 2', 'Centroidi'}, 'Location', 'best');
hold off;

% --- Matrice di proiezioni 2D (Plotmatrix) ---
figure('Color', 'w', 'Name', 'GMM Feature Pairs');
[H,AX,BigAx,P,PAx] = plotmatrix(data_2D_noArtif);
for i = 1:size(AX,1)
    for j = 1:size(AX,2)
        if i ~= j
            cla(AX(i,j)); hold(AX(i,j), 'on');
            for k = 1:K
                idx_k = (cluster_idx == k);
                plot(AX(i,j), data_2D_noArtif(idx_k, j), data_2D_noArtif(idx_k, i), ...
                     '.', 'Color', colors(k,:), 'MarkerSize', 8);
            end
        end
    end
end
title(BigAx, 'Proiezioni 2D delle Feature (Pairwise Plot)');


% --- METRICS ---
cluster_idx = cluster(gmm_model, data_2D_noArtif);
figure;
[s, ~] = silhouette(data_2D_noArtif, cluster_idx);
mean_sil = mean(s);

disp(['Silhouette Score Medio: ' num2str(mean_sil)]);
title(['Silhouette Plot (Score: ' num2str(mean_sil, '%.2f') ')']);

eva_ch = evalclusters(data_2D_noArtif, cluster_idx, 'CalinskiHarabasz');
disp(['Calinski-Harabasz Index: ' num2str(eva_ch.CriterionValues)]);

eva_db = evalclusters(data_2D_noArtif, cluster_idx, 'DaviesBouldin');
disp(['Davies-Bouldin Index: ' num2str(eva_db.CriterionValues)]);

%% extract and save data for the QDA
test_trials  = ntrial_train+1:ntrial; % ---------------------------------------------------------------------------------------------------------------- chenage here for train/test
% test_trials = 1:ntrial_train;
ntrial_test = length(test_trials);
data = squeeze(trial_data(minDurCue+minDurFix+1:end,choosen_band,:,test_trials)); % take just the 8-14 band
data_artifact_cf = artifacts_cf(:,test_trials);
nsamples = size(data,1);
X_ic = []; X_traditional = []; X_nic = [];
y_ic = []; y_traditional = []; y_nic = [];
for idx_trial = 1:length(test_trials)
    for idx_sample = 1:nsamples
        if artifacts_cf(idx_sample,idx_trial) == 0 % no artifact
            X_traditional = [X_traditional; data(idx_sample,:,idx_trial)];
            y_traditional = [y_traditional; trial_typ(idx_trial)];
            if cluster_labels_cf(idx_sample, idx_trial) >= threshold_gmm_ic % IC state
                X_ic = [X_ic; data(idx_sample,:,idx_trial)];
                y_ic = [y_ic; trial_typ(idx_trial)];
            else
                X_nic = [X_nic; data(idx_sample,:, idx_trial)];
                y_nic = [y_nic; trial_typ(idx_trial)];
            end
        end
    end
end
X_ic = log(X_ic);
X_traditional = log(X_traditional);
X_nic = log(X_nic);

%% Comparison GMM: ALL vs IC vc NIC
% fisher score
occipital = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'}; 
[~, ch_occipital] = ismember(occipital, channels_label);
noccipital = size(ch_occipital, 2);

fisher = nan(3, noccipital);

for idx_ch_occipital=1:noccipital
    idx_ch = ch_occipital(idx_ch_occipital);
    % IC
    mu1 = mean(X_ic(y_ic == classes(1),idx_ch));
    sigma1 = std(X_ic(y_ic == classes(1),idx_ch));
    mu2 = mean(X_ic(y_ic == classes(2),idx_ch));
    sigma2 = std(X_ic(y_ic == classes(2),idx_ch));
    fisher(1, idx_ch_occipital) = abs(mu1 - mu2)^2 / (sigma1^2 + sigma2^2);

    % all
    mu1 = mean(X_traditional(y_traditional == classes(1),idx_ch));
    sigma1 = std(X_traditional(y_traditional == classes(1),idx_ch));
    mu2 = mean(X_traditional(y_traditional == classes(2),idx_ch));
    sigma2 = std(X_traditional(y_traditional == classes(2),idx_ch));
    fisher(2, idx_ch_occipital) = abs(mu1 - mu2)^2 / (sigma1^2 + sigma2^2);

    % nic
    mu1 = mean(X_nic(y_nic == classes(1),idx_ch));
    sigma1 = std(X_nic(y_nic == classes(1),idx_ch));
    mu2 = mean(X_nic(y_nic == classes(2),idx_ch));
    sigma2 = std(X_nic(y_nic == classes(2),idx_ch));
    fisher(3, idx_ch_occipital) = abs(mu1 - mu2)^2 / (sigma1^2 + sigma2^2);
end

figure();
imagesc(fisher')
colorbar;
yticks(1:noccipital); yticklabels(occipital)
xticks(1:3); xticklabels({'IC', 'traditional', 'NIC'})
sgtitle('gmm ic and classical fisher score')

figure('Color', 'w', 'Name', 'Paper Proof: Fisher Score Comparison');
b = bar(fisher');
b(1).FaceColor = [0.8 0.2 0.2]; 
b(2).FaceColor = [0.6 0.6 0.6]; 
b(3).FaceColor = [0.2 0.2 0.8]; 
legend({'IC (Selected)', 'Traditional', 'NIC (Rejected)'}, 'Location', 'best');
xticks(1:length(occipital));
xticklabels(occipital);
xtickangle(45);
ylabel('Fisher Score (Separability)');
title('Class Separability across GMM States');
grid on;

% R^2 
r2_ic_allch = calc_r2_from_data(X_ic, y_ic, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['QDA data | size data: ' num2str(size(X_ic,1))]);
r2_all_allch = calc_r2_from_data(X_traditional, y_traditional, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['all data | size data: ' num2str(size(X_traditional,1))]);
r2_nic_allch = calc_r2_from_data(X_nic, y_nic, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['nic data | size data: ' num2str(size(X_nic,1))]);

r2_ic_roi  = r2_ic_allch(ch_occipital);
r2_all_roi = r2_all_allch(ch_occipital);
r2_nic_roi = r2_nic_allch(ch_occipital);
figure('Color', 'w', 'Name', 'Paper Proof: R2 Separability Gain (Tri-state)', 'Position', [100, 100, 1000, 500]);
data_to_plot = [r2_ic_roi(:), r2_all_roi(:), r2_nic_roi(:)];
b = bar(data_to_plot);
b(1).FaceColor = [0.8 0.2 0.2]; % Rosso (IC)
b(2).FaceColor = [0.6 0.6 0.6]; % Grigio (All)
b(3).FaceColor = [0.2 0.2 0.8]; % Blu (NIC)
legend({'GMM IC (Signal)', 'Traditional (Mixed)', 'GMM NIC (Noise)'}, 'Location', 'best', 'FontSize', 10);
ylabel('Signed R^2 (Class Separability)', 'FontSize', 12, 'FontWeight', 'bold');
title('Separability Analysis: What does the GMM reject?', 'FontSize', 14);
grid on;
xticks(1:length(occipital));
xticklabels(occipital);
xtickangle(45);
yline(0, 'k-', 'LineWidth', 1);

% Gain Calculation (IC vs All) & Validation (NIC ~ 0)
gain = mean(abs(r2_ic_roi)) - mean(abs(r2_all_roi));
pct_gain = (gain / mean(abs(r2_all_roi))) * 100;
avg_nic_r2 = mean(abs(r2_nic_roi));

subtitle({['IC Gain over Traditional: +' num2str(pct_gain, '%.1f') '%'], ...
          ['Residual Separability in NIC: ' num2str(avg_nic_r2, '%.4f') ' (Should be near 0)']}, ...
          'FontSize', 10, 'Color', 'k');

%% --- TOPOPLOTS VISUALIZATION (Mean Difference) ---
path_locs = '/home/paolo/chanlocs39.mat'; 
if exist(path_locs, 'file')
    loc_data = load(path_locs);
    f_names = fieldnames(loc_data);
    chanlocs = loc_data.(f_names{1}); 
else
    error(['File chanlocs non trovato in: ' path_locs]);
end

% ALL data -- traditional approach
mu1_all = mean(X_traditional(y_traditional == classes(1), :), 1);
mu2_all = mean(X_traditional(y_traditional == classes(2), :), 1);
diff_all = mu1_all - mu2_all;

% IC data
mu1_ic = mean(X_ic(y_ic == classes(1), :), 1);
mu2_ic = mean(X_ic(y_ic == classes(2), :), 1);
diff_ic = mu1_ic - mu2_ic;

% NIC data
if ~isempty(X_nic)
    mu1_nic = mean(X_nic(y_nic == classes(1), :), 1);
    mu2_nic = mean(X_nic(y_nic == classes(2), :), 1);
    diff_nic = mu1_nic - mu2_nic;
else
    diff_nic = zeros(size(diff_ic));
    warning('Nessun dato NIC trovato.');
end

if exist('topoplot', 'file')
    figure('Color', 'w', 'Name', 'Paper Proof: Topoplot Contrast', 'Position', [100, 100, 1600, 500]);
    
    max_val = max([abs(diff_all), abs(diff_ic), abs(diff_nic)]);
    clim = [-max_val, max_val];
    
    % Plot 1: ALL
    subplot(1, 3, 1);
    topoplot(diff_all, chanlocs, 'maplimits', clim, 'electrodes', 'on', 'style', 'map', 'shading', 'interp');
    title({'ALL Data (Traditional)', ['samples: ' num2str(size(X_traditional,1))]}, 'FontSize', 12);
    colorbar;
    
    % Plot 2: IC (Il risultato buono)
    subplot(1, 3, 2);
    topoplot(diff_ic, chanlocs, 'maplimits', clim, 'electrodes', 'on', 'style', 'map', 'shading', 'interp');
    title({'IC Data (GMM Selected)',  ['samples: ' num2str(size(X_ic,1))]}, 'FontSize', 12, 'FontWeight', 'bold');
    colorbar;
    
    % Plot 3: NIC (Lo scarto)
    subplot(1, 3, 3);
    topoplot(diff_nic, chanlocs, 'maplimits', clim, 'electrodes', 'on', 'style', 'map', 'shading', 'interp');
    title({'NIC Data (GMM Rejected)',  ['samples: ' num2str(size(X_nic,1))]}, 'FontSize', 12);
    colorbar;
    
    colormap(jet); 
    sgtitle(['Mean Difference (Class ' num2str(classes(1)) ' - ' num2str(classes(2)) ')'], 'FontSize', 14);
else
    warning('EEGLAB topoplot non trovato.');
end

%% --- GMM GENERAL BEHAVIOUR ---
figure('Color', 'w', 'Name', 'Paper Proof: Attention Modulation');

prob_mean = mean(cluster_labels_cuecf(:,test_trials), 2);
prob_sem = std(cluster_labels_cuecf(:,test_trials), 0, 2) / sqrt(ntrial_test);

time_axis = (0:length(prob_mean)-1) / sampleRate; 
fill([time_axis fliplr(time_axis)], [prob_mean'+prob_sem' fliplr(prob_mean'-prob_sem')], ...
     [0.8 0.2 0.2], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
xline(minDurCue/sampleRate, 'k-', 'cue', 'LabelVerticalAlignment', 'bottom')
plot(time_axis, prob_mean, 'r-', 'LineWidth', 2);
yline(threshold_gmm_ic, 'k--', 'Threshold', 'LabelVerticalAlignment', 'bottom');
xlabel('Time from Task Start (s)');
ylabel('P(Attention | Data) by GMM');
title('Temporal Dynamics of Attention (Grand Average) | ondulatory behaviour');
grid on;
ylim([0 1]);
xlim([0 time_axis(end)]);

%% --- MULTI-FEATURE GMM VALIDATION ---
feat_LI_all = [];    
feat_LOG_all = [];   
feat_GINI_all = [];  
prob_gmm_all = [];

for tr = test_trials
    probs = cluster_labels_cf(:, tr);
    curr_LI   = squeeze(sparsity(minDurFix+minDurCue+1:end, choosen_band, tr, 1));
    curr_LOG  = squeeze(sparsity(minDurFix+minDurCue+1:end, choosen_band, tr, 2));
    curr_GINI = squeeze(sparsity(minDurFix+minDurCue+1:end, choosen_band, tr, 3));
    
    feat_LI_all   = [feat_LI_all; curr_LI];
    feat_LOG_all  = [feat_LOG_all; curr_LOG];
    feat_GINI_all = [feat_GINI_all; curr_GINI];
    prob_gmm_all  = [prob_gmm_all; probs];
end

[N,edges, bin_idx] = histcounts(prob_gmm_all,10);
feature_data = {feat_LI_all, feat_LOG_all, feat_GINI_all};

c_fig = figure('Color', 'w', 'Name', 'GMM Behavior on All Features', 'Position', [50, 50, 1400, 800]);
min_samples = 5; 
for f = 1:3
    set(0, 'CurrentFigure', c_fig);
    mean_bin = []; 
    sem_bin = []; 
    x_c = [];
    
    current_feat = feature_data{f};
    
    for b = 1:length(edges)-1
        idx = (bin_idx == b);
        
        if N(b) > min_samples
            vals = current_feat(idx);
            vals = vals(isfinite(vals)); 
            
            if ~isempty(vals)
                mean_bin = [mean_bin; mean(vals)];
                sem_bin = [sem_bin; std(vals) / sqrt(length(vals))];
                x_c = [x_c; (edges(b) + edges(b+1)) / 2];
            end
        end
    end
    
    % --- TREND PLOT ---
    ax_trend = subplot(2, 3, f);
    errorbar(x_c, mean_bin, sem_bin, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', 'r', 'Color', 'k');
    
    title(features_name{f}, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Feature Value');
    grid on; xlim([0 1]);
    
    % Calcolo Correlazione
    if length(mean_bin) > 2
        R = corr(x_c, mean_bin);
        subtitle(['Trend Corr: R = ' num2str(R, '%.2f')]);
    end
    
    % --- HISTOGRAM PLOT ---
    ax_hist = subplot(2, 3, f + 3); 
    bar(x_c, N, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
    
    xlabel('GMM Confidence');
    ylabel('# Samples');
    grid on; xlim([0 1]);
    title(['Samples per Bin (N=' num2str(sum(N)) ')']);
    
    linkaxes([ax_trend, ax_hist], 'x');
end
sgtitle('How does GMM Confidence relate to ALL Input Features?', 'FontSize', 16);

%% --- ONDULATORY PLOT ROI ---
% prepare the data
length_trial = size(trial_data,1) - minDurFix - minDurCue;
nroi = 2;
ond_logband_roi = nan(length_trial, floor(ntrial/nclasses), nclasses, nroi);
ond_cl = nan(length_trial * floor(ntrial/nclasses), nclasses);
for c = 1:nclasses
    ond_logband_roi(:,:,c, 1) = squeeze(mean(trial_data(minDurFix+minDurCue+1:end,choosen_band,o_l,trial_typ == classes(c)), 3));
    ond_logband_roi(:,:,c, 2) = squeeze(mean(trial_data(minDurFix+minDurCue+1:end,choosen_band,o_r,trial_typ == classes(c)), 3));
    
    tmp = cluster_labels_cuecf(minDurCue:end, trial_typ == classes(c));
    for i = 1:size(tmp, 2)
        start = (i-1)*size(tmp, 1) + 1;
        stop  = (i)*size(tmp,1);
        ond_cl(start:stop, c) = tmp(:,i);
    end
end

% plot
sat_factor = 1.5;
roi_names = {'ROI 1 (Left)', 'ROI 2 (Right)'};

for c = 1:nclasses
    figure('Name', ['Analisi Classe ' num2str(classes(c))], 'Color', 'w');
    handles = []; c_l = -inf;
    
    for r = 1:nroi
        subplot(1, 2, r);
        
        data_to_plot = ond_logband_roi(:, :, c, r)'; 
        imagesc(data_to_plot);
        
        mu_val = mean(data_to_plot(:), 'omitnan');
        std_val = std(data_to_plot(:), 'omitnan');
        
        limit_upper = mu_val + (sat_factor * std_val);
        
        handles = [handles; gca];
        c_l = max(limit_upper, c_l);
        colormap('jet'); 
        colorbar;
        axis tight;
        title([roi_names{r} ' - Class ' num2str(classes(c))]);
        xlabel('Time (samples)');
        ylabel('Trials');
    end

    set(handles, 'clim', [0 c_l]);
end

%% --- TRIAL PLOTS ---
for c = test_trials
    figure();
    
    % --- log band ---
    subplot(4,1,1)
    imagesc(squeeze(trial_data(minDurFix+minDurCue+1:end,choosen_band,:,c))')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue:sampleRate:min_trial_data) / sampleRate));
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
    legend(features_name)
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue:sampleRate:min_trial_data) / sampleRate));
    xlim([1 size(sparsity_cf, 1)])
%     ylim([0, 1])
    title('sparsity')

    % --- cluster labels gmm ---
    subplot(4,1,3)
    plot(squeeze(cluster_labels_cf(:, c)), 'r')
    hold on
    yline(threshold_gmm_ic, 'k--', 'LineWidth', 2);
    hold off
    legend('gmm', 'threshold gmm ic')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue:sampleRate:min_trial_data) / sampleRate));
    xlim([1 size(sparsity_cf, 1)])
    yticks(0:K-1);
    yticklabels(labels_gmm)   
    title('cluster comparison')

    % --- plot if artefacts ---
    subplot(4,1,4)
    plot(squeeze(artifacts_data(minDurFix+minDurCue+1:end,choosen_band,c)), 'b')
    legend('artefact in the trial')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue:sampleRate:min_trial_data) / sampleRate));
    xlim([1 size(sparsity_cf, 1)])
    yticks(0:1);
    ylim([0 1])
    yticklabels([{'no'}, {'yes'}])   
    title('artefacts')
    
    sgtitle(['Trial: ' num2str(c) ' | Task: ' num2str(trial_typ(c)) ' | band: '  bands_str{choosen_band}])
end
