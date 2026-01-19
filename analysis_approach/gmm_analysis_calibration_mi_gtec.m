clear all; % close all;

addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/equal_ros')
addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/utils')

%% Initialization
DATAPAH = '/home/paolo/cvsa/ic_cvsa_ws/src/';
classes = [769 770];
nchannels = 16;
nclasses = length(classes);
filterOrder = 4;
avg = 1;
threshold_gmm_ic = 0.5;
channels_label = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'Fp1', 'CP1', 'Pz', 'CP2', 'Fp2'};


%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);

%% understand the band
nFiles = length(filenames);
peaks = zeros(1, nFiles);
for idx_file = 1:nFiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    peaks(idx_file) = analyze_alpha_peak(fullpath_file, 'RestTrigger', 786, 'band', [8 14], ...
        'target_regions', {'C1', 'C3', 'C2', 'C4'});
end

%% start processing data
bands = [{[8 13]} {[18 24]}];
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
nbands = length(bands);
signals = cell(1, nbands);
artifacts = [];
headers = cell(1, nbands);
for idx_band = 1:nbands
    headers{idx_band}.TYP = [];
    headers{idx_band}.DUR = [];
    headers{idx_band}.POS = [];
    signals{idx_band} = [];
end

for idx_file= 1: nFiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file);
    c_signal = c_signal(:,1:nchannels);
    sampleRate = header.SampleRate;

    excl_ch = {'FP1', 'FP2', 'EOG'};
    [found, indices] = ismember(excl_ch, channels_label);
    excl_chs = indices(found);

    % for power band using hilbert transformation and artefact remotion -----------------------------------------------
    bufferSize = floor(avg*sampleRate);
    chunkSize = 32;
    eog.filterOrder = 4;
    eog.band = [1 10];
    eog.label = {'FP1', 'FP2'};
    eog.h_threshold = 75; %60;
    eog.v_threshold = 75; %60;
    picks.filterOrder = 4;
    picks.freq = 1; % remove antneuro problems
    picks.threshold = 120; %100;
    artifact = artifact_rejection(c_signal, header, nchannels, bufferSize, chunkSize, eog, picks);
    artifacts = cat(1, artifacts, artifact(:,:));

    disp('   [proc] power band');
    for idx_band = 1:nbands
        band = bands{idx_band};

        [signal_processed, header_processed] = processing_onlineROS_CAR_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize, excl_chs);
        
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
artifacts_data = nan(min_trial_data, ntrial); % data x trial
for idx_band = 1:nbands
    c_signal = signals{idx_band};
    c_artifact = artifacts;
    for trial = 1:ntrial
        c_start = trial_start(trial);
        c_end = trial_start(trial) + min_trial_data - 1;
        trial_data(:,idx_band,:,trial) = c_signal(c_start:c_end,:);
        artifacts_data(:,trial) = c_artifact(c_start:c_end,:);
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
        tmp_art(:,idx_trial_class + idx_class - 1) = artifacts_data(:,idx_classes_trial(i, idx_class));
        trial_typ(idx_trial_class + idx_class - 1) = classes(idx_class);
    end
    i = i + 1;
end
trial_data = tmp_data; % samples x bands x channels x trials
artifacts_data = tmp_art;

%% compute sparsity
% define regions
nsparsity = 1;
sparsity = nan(min_trial_data, ntrial, nsparsity*nbands); % sample x trial x sparsity*nbands
% o_l_ch = {'P3', 'O1', 'P5', 'P1', 'PO5', 'PO3', 'PO7'};
% o_r_ch = {'P4', 'O2', 'P2', 'P6', 'PO4', 'PO6', 'PO8'};
% c_l_ch = {'FC1', 'C3', 'CP1', 'FC3', 'C1', 'CP3'};
% c_r_ch = {'FC2', 'C4', 'CP2', 'FC4', 'C2', 'CP4'};

o_l_ch = {};
o_r_ch = {};
c_c_ch = {'Cz', 'Pz', 'FCz'};
c_l_ch = {'C3', 'CP1', 'C1'};
c_r_ch = {'C4', 'CP2', 'C2'};

[~, o_l] = ismember(o_l_ch, channels_label);
[~, o_r] = ismember(o_r_ch, channels_label);
[~, c_l] = ismember(c_l_ch, channels_label);
[~, c_r] = ismember(c_r_ch, channels_label);
[~, c_c] = ismember(c_c_ch, channels_label);

type = 'mi';

for c = 1:ntrial
    c_data = squeeze(trial_data(:,:,:,c)); % samples x band x channels

    for sample = 1:min_trial_data
        c_sample = squeeze(c_data(sample,:,:)); % bands x channels

        sparsity_vec = []; label_plot = [];

        for idx_band = 1:nbands
            tmp = squeeze(c_sample(idx_band,:)); % 1 x channels
            
            [tmp, label_plot_tmp] =  compute_features_icnic_mi(tmp, c_l, c_r, c_c);

            sparsity_vec = [sparsity_vec; tmp];
            label_plot = [label_plot, label_plot_tmp];
        end

        sparsity(sample, c, :) = sparsity_vec;
    end
end

for idx_band = 1:nbands
    for idx_s = 1:nsparsity
        idx = (idx_band-1)*nsparsity+idx_s;
        label_plot{idx} = [label_plot{idx}, ' ', bands_str{idx_band}];
    end
end

%% ----------------- gmm -----------------
K_range = 2:2; 
best_gmm = [];
min_bic = inf;

sparsity_cf = squeeze(sparsity(minDurFix+minDurCue+1:end, :,:));
artifacts_cf = squeeze(artifacts_data(minDurFix+minDurCue+1:end, :));

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
sparsity_cuecf = squeeze(sparsity(minDurFix+1:end, :,:));
data_2D = reshape(sparsity_cuecf, size(sparsity_cuecf, 1) * size(sparsity_cuecf,2), size(sparsity_cuecf,3));
data_standardized_2D = (data_2D - mu_features) ./ sigma_features;
sparsity_cuecf = reshape(data_standardized_2D, size(sparsity_cuecf, 1), ntrial, size(sparsity_cuecf,3));
sparsity_cf = sparsity_cuecf(minDurCue+1:end,:,:);

disp('Esecuzione di GMM sui dati di training globali...');
options = statset('MaxIter', 1000, 'Display', 'off');
BIC = [];
for k = K_range
    try
        % RegularizationValue = 1e-5 evita che le gaussiane collassino su un punto
        gmm_temp = fitgmdist(data_2D_noArtif, k, ...
                            'Options', options, ...
                            'CovarianceType', 'Full', ...
                            'SharedCovariance', false, ...
                            'RegularizationValue', 1e-5, ...
                            'Replicates', 50); 
        
        if gmm_temp.BIC < min_bic
            min_bic = gmm_temp.BIC;
            best_gmm = gmm_temp;
        end
        BIC =[BIC, gmm_temp.BIC];
    catch ME
        fprintf(2, 'Errore durante il fit con k=%d:\n%s\n', k, ME.message);
        continue;
    end
end

if length(BIC) > 1
    figure();
    plot(BIC)
end

gmm_model = best_gmm;
K = gmm_model.NumComponents;
disp(['GMM ottimizzato: K = ' num2str(K) ' (BIC = ' num2str(min_bic) ')']);

log_features = gmm_model.mu;
mean_score = mean(log_features, 2); 
[~, idx_ic] = max(mean_score); 
[~, idx_nic] = min(mean_score);
fprintf('Mappatura Automatica:\n');
fprintf('Cluster %d (Valore %.2f) -> IC (Task/High Focus)\n', idx_ic, mean_score(idx_ic));
fprintf('Cluster %d (Valore %.2f) -> NIC (Rest/Low Focus)\n', idx_nic, mean_score(idx_nic));

classes_icnic = zeros(1,2);
classes_icnic(idx_ic) = 1;

% Crea le etichette finali
train_gmm = nan(ntrial * size(sparsity_cuecf, 1), nsparsity*nbands);
for c = 1:ntrial
    train_gmm((c-1)*size(sparsity_cuecf, 1) + 1: c * size(sparsity_cuecf, 1),:) = sparsity_cuecf(:,c,:);
end
P_soft = posterior(gmm_model, train_gmm);
cluster_labels_cuecf = nan(size(sparsity_cuecf, 1), ntrial); % contains the prob to be ic
for c = 1:ntrial
    cluster_labels_cuecf(:,c) = P_soft((c-1)*size(sparsity_cuecf, 1) + 1: c * size(sparsity_cuecf, 1), idx_ic);
end
cluster_labels_cf = cluster_labels_cuecf(minDurCue+1:end,:);

% plot the C
disp('centroids: ')
disp(gmm_model.mu)

%% --- VISUALIZZAZIONE GMM ---
% Ottieni le etichette "Hard" dal GMM (assegna ogni punto al cluster più probabile)
cluster_idx = cluster(gmm_model, data_2D_noArtif);

% Calcola la Silhouette
figure;
[s, ~] = silhouette(data_2D_noArtif, cluster_idx);
mean_sil = mean(s);

disp(['Silhouette Score Medio: ' num2str(mean_sil)]);
% Salva il grafico per il paper
title(['Silhouette Plot (Score: ' num2str(mean_sil, '%.2f') ')']);

eva_ch = evalclusters(data_2D_noArtif, cluster_idx, 'CalinskiHarabasz');
disp(['Calinski-Harabasz Index: ' num2str(eva_ch.CriterionValues)]);

eva_db = evalclusters(data_2D_noArtif, cluster_idx, 'DaviesBouldin');
disp(['Davies-Bouldin Index: ' num2str(eva_db.CriterionValues)]);

%% extract data for the selection of the features
test_trials = 1:ntrial_train; % ------------------------------------------------------------------------------------------------------ here for test and train
ntrial_test = length(test_trials);
data = trial_data(minDurCue+minDurFix+1:end,:,:,test_trials); % data x bands x channels x trial
nsamples = size(data,1);
X_ic = []; X_all = []; X_nic = [];
count_artifact = 0; count_all = 0; count_rejected = 0;
for idx_band = 1:nbands
    tmp_X_ic = []; tmp_X_all = []; tmp_X_nic = [];
    y_ic = []; y_all = []; y_nic = [];
    trials = [];
    for idx_trial =  1:ntrial_test
        for idx_sample = 1:nsamples
            if artifacts_cf(idx_sample,idx_trial) == 0 % no artifact
                tmp_X_all = [tmp_X_all; data(idx_sample,idx_band,:,idx_trial)];
                y_all = [y_all; trial_typ(idx_trial)];
                if cluster_labels_cf(idx_sample, idx_trial) >= threshold_gmm_ic % IC state
                    tmp_X_ic = [tmp_X_ic; data(idx_sample,idx_band,:,idx_trial)];
                    y_ic = [y_ic; trial_typ(idx_trial)];
                    trials = [trials; idx_trial];
                else
                    tmp_X_nic = [tmp_X_nic; data(idx_sample,idx_band,:,idx_trial)];
                    y_nic = [y_nic; trial_typ(idx_trial)];
                    count_rejected = count_rejected + 1/nbands;
                end
            else
                count_artifact = count_artifact + 1/nbands;
            end
            count_all = count_all + 1/nbands;
        end
    end
    tmp_X_ic = log(tmp_X_ic);
    tmp_X_nic = log(tmp_X_nic);
    tmp_X_all = log(tmp_X_all);

    X_ic = [X_ic, tmp_X_ic];
    X_nic = [X_nic, tmp_X_nic];
    X_all = [X_all, tmp_X_all];
end

disp(['all sample for the trials: ' num2str(count_all) ', rejected for artifact: ' num2str(count_artifact) ', rejected gmm: ' num2str(count_rejected) ', acepted: ' num2str(size(X_ic, 1))])

%% Features selection QDA
% fisher score
fisher = nan(nbands*3, nchannels);
label_fisher = [];

for idx_ch=1:nchannels
    for idx_band = 1:nbands
        % IC
        mu1 = mean(X_ic(y_ic == classes(1),idx_band, idx_ch));
        sigma1 = std(X_ic(y_ic == classes(1),idx_band, idx_ch));
        mu2 = mean(X_ic(y_ic == classes(2),idx_band,idx_ch));
        sigma2 = std(X_ic(y_ic == classes(2),idx_band, idx_ch));
        fisher((idx_band-1)*3 + 1, idx_ch) = abs(mu1 - mu2)^2 / (sigma1^2 + sigma2^2);
        label_fisher = [label_fisher, {['IC', bands_str{idx_band}]}];

        % all
        mu1 = mean(X_all(y_all == classes(1), idx_band, idx_ch));
        sigma1 = std(X_all(y_all == classes(1),idx_band, idx_ch));
        mu2 = mean(X_all(y_all == classes(2),idx_band, idx_ch));
        sigma2 = std(X_all(y_all == classes(2),idx_band, idx_ch));
        fisher((idx_band-1)*3 +2, idx_ch) = abs(mu1 - mu2)^2 / (sigma1^2 + sigma2^2);
        label_fisher = [label_fisher, {['traditional', bands_str{idx_band}]}];

        % nic
        mu1 = mean(X_nic(y_nic == classes(1), idx_band, idx_ch));
        sigma1 = std(X_nic(y_nic == classes(1),idx_band, idx_ch));
        mu2 = mean(X_nic(y_nic == classes(2),idx_band, idx_ch));
        sigma2 = std(X_nic(y_nic == classes(2),idx_band, idx_ch));
        fisher((idx_band-1)*3 +3, idx_ch) = abs(mu1 - mu2)^2 / (sigma1^2 + sigma2^2);
        label_fisher = [label_fisher, {['nic', bands_str{idx_band}]}];
    end
end

figure();
imagesc(fisher')
colorbar;
yticks(1:nchannels); yticklabels(channels_label)
xticks(1:size(fisher, 1)); xticklabels(label_fisher)
sgtitle('gmm ic and classical fisher score')

figure('Color', 'w', 'Name', 'Fisher Score Comparison Bars', 'Position', [100, 100, 1200, 600]);
b = bar(fisher'); % Trasposto: X=Canali, Gruppi=Conditions
% Assegnazione colori ciclica (Rosso=IC, Grigio=All, Blu=NIC)
ic_color = [0.8 0.2 0.2];
all_color = [0.6 0.6 0.6];
nic_color = [0.2 0.2 0.8];
for k = 1:length(b)
    rem_k = mod(k-1, 3) + 1; % 1, 2, 3 ripetuto
    if rem_k == 1
        b(k).FaceColor = ic_color; % IC
    elseif rem_k == 2
        b(k).FaceColor = all_color; % All
    else
        b(k).FaceColor = nic_color; % NIC
    end
end
legend([b(1), b(2), b(3)], {'IC (Selected)', 'Traditional', 'NIC (Rejected)'}, 'Location', 'bestoutside');
xticks(1:nchannels); xticklabels(channels_label); xtickangle(45);
ylabel('Fisher Score');
title('Feature Separability: Does GMM select better data?');
grid on;


% R^2
for idx_band = 1:nbands
    r2_ic_allch = calc_r2_from_data(squeeze(X_ic(:,idx_band,:)), y_ic, 'Plot', false, 'ChanLabels', channels_label); 
    r2_all_allch = calc_r2_from_data(squeeze(X_all(:,idx_band,:)), y_all, 'Plot', false, 'ChanLabels', channels_label);
    r2_nic_allch = calc_r2_from_data(squeeze(X_nic(:,idx_band,:)), y_nic, 'Plot', false, 'ChanLabels', channels_label);


    figure('Color', 'w', 'Name', ['R2 Gain - Band ' bands_str{idx_band}], 'Position', [150, 150, 1000, 500]);
    data_to_plot = [r2_ic_allch(:), r2_all_allch(:), r2_nic_allch(:)];
    b_r2 = bar(data_to_plot);
    b_r2(1).FaceColor = ic_color;   % Rosso
    b_r2(2).FaceColor = all_color;  % Grigio
    b_r2(3).FaceColor = nic_color;  % Blu
    legend({'GMM IC (Signal)', 'Traditional (Mixed)', 'GMM NIC (Noise)'}, 'Location', 'best');
    ylabel('Signed R^2 (Class Separability)', 'FontSize', 12, 'FontWeight', 'bold');
    title(['Separability Analysis | Band: ' bands_str{idx_band}], 'FontSize', 14);
    grid on;
    xticks(1:nchannels); xticklabels(channels_label); xtickangle(45);
    yline(0, 'k-', 'LineWidth', 1);

    mean_r2_ic = mean(abs(r2_ic_allch));
    mean_r2_all = mean(abs(r2_all_allch));
    mean_r2_nic = mean(abs(r2_nic_allch));
    gain = mean_r2_ic - mean_r2_all;
    pct_gain = (gain / (mean_r2_all + eps)) * 100;
    subtitle({['IC vs All Improvement: +' num2str(pct_gain, '%.1f') '%'], ...
              ['Avg R^2: IC=' num2str(mean_r2_ic,'%.3f') ' | All=' num2str(mean_r2_all,'%.3f') ' | NIC=' num2str(mean_r2_nic,'%.3f')]}, ...
              'FontSize', 10, 'Color', 'k');
end

%% --- TOPOPLOTS VISUALIZATION (Mean Difference) ---
disp('Generazione Topoplot (Differenza Media)...');

path_locs = '/home/paolo/chanlocs39.mat'; 

if exist(path_locs, 'file')
    loc_data = load(path_locs);
    f_names = fieldnames(loc_data);
    chanlocs_full = loc_data.(f_names{1}); 
    disp(['[INFO] Chanlocs caricato.']);
else
    error(['File chanlocs non trovato in: ' path_locs]);
end

all_labels = {chanlocs_full.labels};
[is_present, idx_in_full] = ismember(channels_label, all_labels);
if ~all(is_present)
    missing_chans = all_labels(~is_present);
    error('I seguenti canali non sono stati trovati nel file chanlocs39: %s', strjoin(missing_chans, ', '));
end
chanlocs_subset = chanlocs_full(idx_in_full);
disp(['[INFO] Chanlocs ridotto a ' num2str(length(chanlocs_subset)) ' canali per il plotting.']);


figure('Color', 'w', 'Name', 'Topoplot Mean Difference', 'Position', [100, 100, 1200, 500]);
max_val = 0.1;
handles = [];
for idx_band = 1:nbands
    % --- X_ic ---
    X_c1 = squeeze(X_ic(y_ic == classes(1), idx_band, :)); % Classe 1 (es. 730)
    X_c2 = squeeze(X_ic(y_ic == classes(2), idx_band, :)); % Classe 2 (es. 731)

    mu1_ic = mean(X_c1, 1);
    mu2_ic = mean(X_c2, 1);
    diff_ic = mu1_ic - mu2_ic;

    % --- X_all ---
    Xall_c1 = squeeze(X_all(y_all == classes(1), idx_band, :));
    Xall_c2 = squeeze(X_all(y_all == classes(2), idx_band, :));

    mu1_all = mean(Xall_c1, 1);
    mu2_all = mean(Xall_c2, 1);
    diff_all = mu1_all - mu2_all;

    % --- X_nic ---
    Xnic_c1 = squeeze(X_nic(y_nic == classes(1), idx_band, :));
    Xnic_c2 = squeeze(X_nic(y_nic == classes(2), idx_band, :));

    mu1_nic = mean(Xnic_c1, 1);
    mu2_nic = mean(Xnic_c2, 1);
    diff_nic = mu1_nic - mu2_nic;

    % plotting
    max_val = max([abs(diff_ic), abs(diff_all), max_val, abs(diff_nic)]);

    subplot(nbands, 3, (idx_band-1)*3 + 1);
    topoplot(diff_ic, chanlocs_subset, 'electrodes', 'on', 'style', 'map', 'shading', 'interp');
    title(['GMM Selected (IC) | ', num2str(bands{idx_band}(1)) '-' num2str(bands{idx_band}(2)) ' Hz)'], 'FontSize', 12, 'FontWeight', 'bold');
    handles = [handles, gca];
    colorbar;

    subplot(nbands, 3, (idx_band-1)*3 + 2);
    topoplot(diff_all, chanlocs_subset, 'electrodes', 'on', 'style', 'map', 'shading', 'interp');
    title(['All Data | ', num2str(bands{idx_band}(1)) '-' num2str(bands{idx_band}(2)) ' Hz)'], 'FontSize', 12, 'FontWeight', 'bold');
    handles = [handles, gca];
    colorbar;

    subplot(nbands, 3, (idx_band-1)*3 + 3);
    topoplot(diff_nic, chanlocs_subset, 'electrodes', 'on', 'style', 'map', 'shading', 'interp');
    title(['NIC Data | ', num2str(bands{idx_band}(1)) '-' num2str(bands{idx_band}(2)) ' Hz)'], 'FontSize', 12, 'FontWeight', 'bold');
    handles = [handles, gca];
    colorbar;

    
    colormap(jet);
    sgtitle(['Mean Diff | ' num2str(classes(1)) ' - ' num2str(classes(2))], 'FontSize', 14);

end
set(handles, 'CLim', [-max_val, max_val])

%% --- MULTI-FEATURE GMM VALIDATION ---
feature_data = [];    
prob_gmm_all = [];

for tr = test_trials
    probs = cluster_labels_cf(:, tr);
    tmp_f = [];
    for f=1:size(sparsity_cf,3)
        tmp_f = [tmp_f, squeeze(sparsity_cf(:,tr,f))];
    end

    prob_gmm_all  = [prob_gmm_all; probs];
    feature_data = [feature_data; tmp_f];
end

[N,edges, bin_idx] = histcounts(prob_gmm_all,10);

c_fig = figure('Color', 'w', 'Name', 'GMM Behavior on All Features', 'Position', [50, 50, 1400, 800]);
min_samples = 5; 
for f = 1:size(feature_data, 2)
    set(0, 'CurrentFigure', c_fig);
    mean_bin = []; 
    sem_bin = []; 
    x_c = [];
    
    current_feat = feature_data(:,f);
    
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
    ax_trend = subplot(2, size(feature_data, 2), f);
    errorbar(x_c, mean_bin, sem_bin, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', 'r', 'Color', 'k');
    
    title(label_plot{f}, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Feature Value');
    grid on; xlim([0 1]);
    
    % Calcolo Correlazione
    if length(mean_bin) > 2
        R = corr(x_c, mean_bin);
        subtitle(['Trend Corr: R = ' num2str(R, '%.2f')]);
    end
    
    % --- HISTOGRAM PLOT ---
    ax_hist = subplot(2, size(feature_data, 2), f + size(feature_data, 2)); 
    bar(x_c, N, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
    
    xlabel('GMM Confidence');
    ylabel('# Samples');
    grid on; xlim([0 1]);
    title(['Samples per Bin (N=' num2str(sum(N)) ')']);
    
    linkaxes([ax_trend, ax_hist], 'x');
end
sgtitle('How does GMM Confidence relate to ALL Input Features?', 'FontSize', 16);


%% --- TRIAL PLOTS ---
for c = ntrial-10:ntrial
    figure();
    
    % --- log band band 1---
    subplot(4,1,1)
    imagesc(squeeze(trial_data(minDurFix+minDurCue+1:end,1,:,c))')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue:sampleRate:min_trial_data) / sampleRate));
    yticks(1:nchannels); yticklabels(channels_label)
    title(['log band ' bands_str{1}])

    % --- log band band 2 ---
    subplot(4,1,2)
    imagesc(squeeze(trial_data(minDurFix+minDurCue+1:end,2,:,c))')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue:sampleRate:min_trial_data) / sampleRate));
    yticks(1:nchannels); yticklabels(channels_label)
    title(['log band ' bands_str{2}])

    % --- sparsity indices ---
    subplot(4,1,3)
    plot(squeeze(sparsity_cf(:, c, 1)))
    hold on
    for i = 2:nsparsity*nbands
        plot(squeeze(sparsity_cf(:, c, i)))
    end
    hold off
    legend(label_plot)
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue:sampleRate:min_trial_data) / sampleRate));
    xlim([1 size(sparsity_cf, 1)])
    title('sparsity')

    % --- cluster labels gmm ---
    subplot(4,1,4)
    plot(squeeze(cluster_labels_cf(:, c)), 'r')
    hold on
    yline(threshold_gmm_ic, 'k--', 'LineWidth', 2);
    plot(squeeze(artifacts_data(minDurFix+minDurCue+1:end,c)), 'b')
    hold off
    legend('gmm', 'threshold gmm ic', 'artifact')
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue:sampleRate:min_trial_data) / sampleRate));
    xlim([1 size(sparsity_cf, 1)])
    yticks(0:K-1);
    yticklabels({'NIC/NO', 'IC/YES'})  
    ylim([0 1])
    title('cluster amd artifacts')
    
    sgtitle(['Trial: ' num2str(c) ' | Task: ' num2str(trial_typ(c))])
end

%% ERD ERS
choosen_band = 1;
data_for_plot = squeeze(trial_data(:, choosen_band, :, :)); % Diventa [Samples x Ch x Trials]
arts_for_plot = squeeze(artifacts_data(:, choosen_band, :)); % [Samples x Trials]

plot_header.SampleRate = sampleRate; % Assicurati che sia il rateo DOPO il chunking 
plot_header.Label = channels_label;  % O headers{1}.channels_labels

ch_list = {'C3', 'C1', 'C4', 'C2', 'Cz'};

disp('Plotting ERD/ERS...');
plot_erd_ers(log(data_for_plot), trial_typ, arts_for_plot, plot_header, minDurFix, ch_list, ['band: ' bands_str{choosen_band}]);


function plot_erd_ers(trial_data_band, trial_typ, artifacts_band, header, minDurFix, channels_to_plot, title_plot)
% PLOT_ERD_ERS Visualizza l'andamento temporale della potenza (ERD/ERS).
%
% INPUT:
%   trial_data_band: Matrice 3D [Samples x Channels x Trials] (già selezionata la banda!)
%   trial_typ:       Vettore [Trials x 1] con le classi (es. 730, 731)
%   artifacts_band:  Matrice 2D [Samples x Trials] maschera artefatti (0=clean, 1=artifact)
%   header:          Header struct (per SampleRate e Label)
%   minDurFix:       Numero di campioni della fissazione (per allineare lo 0 al Cue)
%   channels_to_plot: Cell array con i nomi dei canali da plottare (es. {'PO7', 'PO8'})

    fs = header.SampleRate;
    labels = header.Label;
    
    % Identifica le classi
    classes = unique(trial_typ);
    if length(classes) ~= 2
        warning('Trovate %d classi. Il plot è ottimizzato per 2 classi.', length(classes));
    end
    
    % Asse Temporale: 0 = Inizio Cue
    [n_samples, ~, n_trials] = size(trial_data_band);
    t_axis = ((1:n_samples) - minDurFix) / fs; % Secondi relativi al Cue
    
    % Baseline Period: Fixation
    base_start = 1; 
    base_end = minDurFix;
    
    figure('Color', 'w', 'Name', 'ERD/ERS Time Course', 'Position', [100, 100, 1200, 600]);
    
    for i = 1:length(channels_to_plot)
        ch_name = channels_to_plot{i};
        
        % Trova indice canale
        % Gestione label (rimuove 'EEG ' se presente per il confronto)
        clean_labels = strrep(labels, 'EEG ', '');
        ch_idx = find(strcmpi(clean_labels, ch_name), 1);
        
        if isempty(ch_idx)
            warning('Canale %s non trovato.', ch_name);
            continue;
        end
        
        % Estrai dati canale: [Samples x Trials]
        raw_ch_data = squeeze(trial_data_band(:, ch_idx, :));
        % artifacts_band è [Samples x Trials]
        raw_ch_data(artifacts_band == 1) = NaN;
        
        % --- CALCOLO ERD% ---
        % Formula: (Potenza(t) - Baseline) / Baseline * 100
        baseline = mean(mean(raw_ch_data(base_start:base_end, :), 1, 'omitnan'), 2, 'omitnan');
        
        erd_data = zeros(size(raw_ch_data));
        for tr = 1:n_trials
            erd_data(:, tr) = raw_ch_data(:, tr) - baseline; % dB change
        end
        
        % Media per Classe
        mean_c1 = mean(erd_data(:, trial_typ == classes(1)), 2, 'omitnan');
        mean_c2 = mean(erd_data(:, trial_typ == classes(2)), 2, 'omitnan');
        
        % Standard Error (per l'ombra)
        se_c1 = std(erd_data(:, trial_typ == classes(1)), 0, 2, 'omitnan') ./ sqrt(sum(trial_typ == classes(1)));
        se_c2 = std(erd_data(:, trial_typ == classes(2)), 0, 2, 'omitnan') ./ sqrt(sum(trial_typ == classes(2)));
        
        % Smoothing grafico (Moving average per pulire le linee)
        smooth_win = round(fs * 0.2); % 200ms
        mean_c1 = movmean(mean_c1, smooth_win);
        mean_c2 = movmean(mean_c2, smooth_win);
        
        % --- SUBPLOT ---
        subplot(ceil(length(channels_to_plot)/2), 2, i);
        hold on;
        
        % Area Fissazione (Grigio chiaro)
        xline(0, 'k-', 'LineWidth', 1.5); % Cue Onset
        xline(1, 'k', 'LineWidth',1.5);
        
        % Plot con ombra (funzione interna semplice o fill)
        plot_shaded(t_axis, mean_c1, se_c1, [0.2 0.6 1]); % Blu (Classe 1)
        plot_shaded(t_axis, mean_c2, se_c2, [1 0.4 0.4]); % Rosso (Classe 2)
        
        title(['Canale ' ch_name], 'FontSize', 12, 'FontWeight', 'bold');
        if i == 1
            legend({'cue start', 'cf start', ['std class ' num2str(classes(1))], ['mean class ' num2str(classes(1))], ...
                ['std class ' num2str(classes(2))], ['mean class ' num2str(classes(2))]}, 'Location', 'best');
        end
        xlabel('Time (s)');
        ylabel('Log Power Change (dB-like)');
        grid on;
        xlim([t_axis(1), t_axis(end)]);
    end
    sgtitle(title_plot)
end

function plot_shaded(x, y, err, color)
    % Helper per disegnare linea + errore
    x = x(:)'; y = y(:)'; err = err(:)';
    % Rimuovi NaN per il patch
    valid = ~isnan(y) & ~isnan(err);
    x = x(valid); y = y(valid); err = err(valid);
    
    fill([x fliplr(x)], [y+err fliplr(y-err)], color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(x, y, 'Color', color, 'LineWidth', 2);
end