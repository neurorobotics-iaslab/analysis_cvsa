clear all; % close all;

addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_cvsa/equal_ros')

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
artifacts_data = tmp_art;

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

            [sparsity(sample, idx_band, c,:), label_plot] = compute_features_icnic(tmp, type, o_l, o_r, c_l, c_r, nsparsity);
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
data_3D = sparsity_cf(:, :, :);
data_2D = reshape(data_3D, size(data_3D, 1) * size(data_3D,2), size(data_3D,3));
artefact_2D = artifacts_cf(:, :);
artefact_1D = reshape(artefact_2D, size(artefact_2D, 1) * size(artefact_2D,2), 1);
data_2D_noArtif = data_2D(artefact_1D == 0,:);
mu_features = mean(data_2D_noArtif, 1);
sigma_features = std(data_2D_noArtif, 0, 1);
sigma_features(sigma_features == 0) = eps;
data_2D_noArtif = (data_2D_noArtif - mu_features) ./ sigma_features;
data_standardized_2D = (data_2D - mu_features) ./ sigma_features;
sparsity_cf = reshape(data_standardized_2D, size(data_3D, 1), ntrial, size(data_3D,3));

disp('Esecuzione di GMM sui dati di training globali...');
options = statset('MaxIter', 1000, 'Display', 'off');
for k = K_range
    try
        % RegularizationValue = 1e-5 evita che le gaussiane collassino su un punto
        gm_temp = fitgmdist(data_2D_noArtif, k, ...
                            'Options', options, ...
                            'CovarianceType', 'Full', ...
                            'SharedCovariance', false, ...
                            'RegularizationValue', 1e-5, ...
                            'Replicates', 150); 
        
        if gm_temp.BIC < min_bic
            min_bic = gm_temp.BIC;
            best_gmm = gm_temp;
        end
    catch ME
        fprintf(2, 'Errore durante il fit con k=%d:\n%s\n', k, ME.message);
        continue;
    end
end

gmm_model = best_gmm;
K = gmm_model.NumComponents;
disp(['GMM ottimizzato: K = ' num2str(K) ' (BIC = ' num2str(min_bic) ')']);

[~, sort_order] = sort(gmm_model.mu(:, 2), 'descend');
idx_nic = sort_order(1); 
idx_ic = sort_order(2); % Potenza BASSA = Il cervello sta lavorando (ERD su uno o entrambi i lati) -> IC
classes_icnic = zeros(1,2);
classes_icnic(idx_ic) = 1;

% Crea le etichette finali
labels_gmm = {'IC', 'NIC'};

train_gmm = nan(ntrial * size(sparsity_cf, 1), nsparsity);
for c = 1:ntrial
    train_gmm((c-1)*size(sparsity_cf, 1) + 1: c * size(sparsity_cf, 1),:) = sparsity_cf(:,c,:);
end
P_soft = posterior(gmm_model, train_gmm);
cluster_labels = nan(size(sparsity_cf, 1), ntrial); % contains the prob to be ic
for c = 1:ntrial
    cluster_labels(:,c) = P_soft((c-1)*size(sparsity_cf, 1) + 1: c * size(sparsity_cf, 1),idx_ic);
end

fprintf('Mappatura: Cluster GMM %d -> "ic", Cluster GMM %d -> "nic"\n', idx_ic, idx_nic);

% plot the C
disp('centroids: ')
disp(gmm_model.mu)

%% --- VISUALIZZAZIONE GMM ---
figure('Color', 'w'); % Crea una figura con sfondo bianco
hold on;
scatter(data_2D_noArtif(:,1), data_2D_noArtif(:,2), 15, ...
        'MarkerFaceColor', [0.2 0.5 0.9], ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.4);

x_min = min(data_2D_noArtif(:,1)) - 1; x_max = max(data_2D_noArtif(:,1)) + 1;
y_min = min(data_2D_noArtif(:,2)) - 1; y_max = max(data_2D_noArtif(:,2)) + 1;
step = 0.05; 
[x1Grid, x2Grid] = meshgrid(x_min:step:x_max, y_min:step:y_max);
XGrid = [x1Grid(:), x2Grid(:)];

prob_GMM = pdf(gmm_model, XGrid);
prob_GMM = reshape(prob_GMM, size(x1Grid));
[C, ~] = contour(x1Grid, x2Grid, prob_GMM, 10, 'LineWidth', 2, 'LineColor', [0.8 0.2 0.2]);
plot(gmm_model.mu(:,1), gmm_model.mu(:,2), 'k+', 'MarkerSize', 15, 'LineWidth', 3);

xlabel(label_plot{1}, 'FontSize', 12, 'FontWeight', 'bold');
ylabel(label_plot{2}, 'FontSize', 12, 'FontWeight', 'bold');
title('GMM Fit: Cluster IC vs NIC', 'FontSize', 14);
legend({'Dati Reali', 'Ellissi GMM', 'Centroidi'}, 'Location', 'best');
grid on;
axis tight;
hold off;

% 1. Ottieni le etichette "Hard" dal GMM (assegna ogni punto al cluster più probabile)
cluster_idx = cluster(gmm_model, data_2D_noArtif);

% 2. Calcola la Silhouette
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

%% extract and save data for the QDA
data = squeeze(trial_data(minDurCue+minDurFix+1:end,choosen_band,:,:)); % take just the 8-14 band
nsamples = size(data,1);
X = []; X_all = [];
y = []; y_all = [];
trials = [];
for idx_trial =  1:ntrial
    for idx_sample = 1:nsamples
        if artifacts_cf(idx_sample,idx_trial) == 0 % no artifact
            X_all = [X_all; data(idx_sample,:,idx_trial)];
            y_all = [y_all; trial_typ(idx_trial)];
            if cluster_labels(idx_sample, idx_trial) >= threshold_gmm_ic % IC state
                X = [X; data(idx_sample,:,idx_trial)];
                y = [y; trial_typ(idx_trial)];
                trials = [trials; idx_trial];
            end
        end
    end
end
X = log(X);
X_all = log(X_all);


%% Check the data used for the QDA
% fisher score
occipital = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'}; 
[~, ch_occipital] = ismember(occipital, channels_label);
noccipital = size(ch_occipital, 2);

fisher = nan(2, noccipital);

for idx_ch_occipital=1:noccipital
    idx_ch = ch_occipital(idx_ch_occipital);
    % IC
    mu1 = mean(X(y == classes(1),idx_ch));
    sigma1 = std(X(y == classes(1),idx_ch));
    mu2 = mean(X(y == classes(2),idx_ch));
    sigma2 = std(X(y == classes(2),idx_ch));
    fisher(1, idx_ch_occipital) = abs(mu1 - mu2)^2 / (sigma1^2 + sigma2^2);

    % all
    mu1 = mean(X_all(y_all == classes(1),idx_ch));
    sigma1 = std(X_all(y_all == classes(1),idx_ch));
    mu2 = mean(X_all(y_all == classes(2),idx_ch));
    sigma2 = std(X_all(y_all == classes(2),idx_ch));
    fisher(2, idx_ch_occipital) = abs(mu1 - mu2)^2 / (sigma1^2 + sigma2^2);
end

figure();
imagesc(fisher')
colorbar;
yticks(1:noccipital); yticklabels(occipital)
xticks(1:2); xticklabels({'IC', 'traditional'})
sgtitle('gmm ic and classical fisher score')

% R^2
calc_r2_from_data(X, y, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['QDA data | size data: ' num2str(size(X,1))]);
calc_r2_from_data(X_all, y_all, 'Plot', true, 'ChanLabels', channels_label, 'title_data', ['all data | size data: ' num2str(size(X_all,1))]);

%% --- TOPOPLOTS VISUALIZATION (Mean Difference) ---
disp('Generazione Topoplot (Differenza Media)...');

path_locs = '/home/paolo/chanlocs39.mat'; 

if exist(path_locs, 'file')
    loc_data = load(path_locs);
    % Estrae la variabile in modo dinamico
    f_names = fieldnames(loc_data);
    chanlocs = loc_data.(f_names{1}); 
    disp(['[INFO] Chanlocs caricato.']);
else
    error(['File chanlocs non trovato in: ' path_locs]);
end

% --- A. Per X (GMM Selected) ---
X_c1 = X(y == classes(1), :); % Classe 1 (es. 730)
X_c2 = X(y == classes(2), :); % Classe 2 (es. 731)

% Media sulle colonne (canali)
mu1_ic = mean(X_c1, 1);
mu2_ic = mean(X_c2, 1);
diff_ic = mu1_ic - mu2_ic;

% --- B. Per X_all (All Data) ---
Xall_c1 = X_all(y_all == classes(1), :);
Xall_c2 = X_all(y_all == classes(2), :);

mu1_all = mean(Xall_c1, 1);
mu2_all = mean(Xall_c2, 1);
diff_all = mu1_all - mu2_all;

% 3. Plotting
if exist('topoplot', 'file')
    figure('Color', 'w', 'Name', 'Topoplot Mean Difference', 'Position', [100, 100, 1200, 500]);
    
    % Calcolo limiti colore comuni per confronto diretto
    % Troviamo il valore massimo di differenza in assoluto per centrare la colormap
    max_val = max([abs(diff_ic), abs(diff_all)]);
    if max_val == 0, max_val = 0.1; end
    clim = [-max_val, max_val];
    
    % Subplot 1: GMM (IC)
    subplot(1, 2, 1);
    topoplot(diff_ic, chanlocs, 'maplimits', clim, 'electrodes', 'on', 'style', 'map', 'shading', 'interp');
    title({'GMM Selected (IC)', ['Mean Diff (' num2str(classes(1)) ' - ' num2str(classes(2)) ')']}, 'FontSize', 12, 'FontWeight', 'bold');
    colorbar;
    
    % Subplot 2: All Data
    subplot(1, 2, 2);
    topoplot(diff_all, chanlocs, 'maplimits', clim, 'electrodes', 'on', 'style', 'map', 'shading', 'interp');
    title({'All Data ', ['Mean Diff (' num2str(classes(1)) ' - ' num2str(classes(2)) ')']}, 'FontSize', 12, 'FontWeight', 'bold');
    colorbar;
    
    % Colormap: Rosso = Classe 1 Maggiore, Blu = Classe 2 Maggiore
    colormap(jet); 
    sgtitle(['Difference power (' num2str(bands{choosen_band}(1)) '-' num2str(bands{choosen_band}(2)) ' Hz)'], 'FontSize', 14);
    
else
    warning('EEGLAB topoplot non trovato.');
end

%% --- TRIAL PLOTS ---
for c = 1:10
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
    legend(label_plot)
    xticks(0:sampleRate:min_trial_data)
    xticklabels(string((minDurFix+minDurCue:sampleRate:min_trial_data) / sampleRate));
    xlim([1 size(sparsity_cf, 1)])
%     ylim([0, 1])
    title('sparsity')

    % --- cluster labels gmm ---
    subplot(4,1,3)
    plot(squeeze(cluster_labels(:, c)), 'r')
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

%% ERD ERS
% 1. Prepara i dati per la banda scelta (es. choosen_band = 1)
data_for_plot = squeeze(trial_data(:, choosen_band, :, :)); % Diventa [Samples x Ch x Trials]
arts_for_plot = squeeze(artifacts_data(:, choosen_band, :)); % [Samples x Trials]

plot_header.SampleRate = sampleRate; % Assicurati che sia il rateo DOPO il chunking 
plot_header.Label = channels_label;  % O headers{1}.channels_labels

ch_list = {'PO7', 'PO8', 'O1', 'P5', 'PO5', 'O2'};

disp('Plotting ERD/ERS...');
plot_erd_ers(log(data_for_plot), trial_typ, arts_for_plot, plot_header, minDurFix, ch_list);


function plot_erd_ers(trial_data_band, trial_typ, artifacts_band, header, minDurFix, channels_to_plot)
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