clear; clc; close all;

%% load data
data_path = '/home/paolo/cvsa/ic_cvsa_ws/record_mi/results_graz/gmm';
files = dir(fullfile(data_path, 'GMM_Validation_*.mat'));
n_subjects = length(files);
if n_subjects == 0
    error('Nessun file trovato! Esegui prima lo script di analisi per salvare i dati.');
end
disp(['Trovati ' num2str(n_subjects) ' soggetti. Generazione plot...']);

%% inizialization
feat_idx_mu = 1;   % Colonna 1: 8-13 Hz
feat_idx_beta = 2; % Colonna 2: 18-24 Hz
min_samples_bin = 10;
n_bins = 10;
edges = linspace(0, 1, n_bins+1);
centers = (edges(1:end-1) + edges(2:end)) / 2;
colors = lines(n_subjects);
figure('Color', 'w', 'Name', 'GMM Final Validation', 'Position', [100, 50, 900, 800]);
% Variabili di supporto
hist_matrix = zeros(n_bins, n_subjects);
subj_names = cell(1, n_subjects);
legend_labels = cell(1, n_subjects);
legend_handles = []; 

%% --- RENDS ---
subplot(2, 1, 1); hold on;

for i = 1:n_subjects
    load(fullfile(files(i).folder, files(i).name));
    subj_names{i} = subject;

    %% (8-13 Hz)
    curr_feat = feature_data(:, feat_idx_mu);
    [N, ~, bin_idx] = histcounts(prob_gmm_all, edges);
    mean_bin = []; sem_bin = []; x_c = [];
    
    for b = 1:length(edges)-1
        if N(b) > min_samples_bin
            vals = curr_feat(bin_idx == b);
            vals = vals(isfinite(vals));
            if ~isempty(vals)
                mean_bin = [mean_bin; mean(vals)];
                sem_bin = [sem_bin; std(vals)/sqrt(length(vals))];
                x_c = [x_c; centers(b)];
            end
        end
    end
    
    % --- corr ---
    if length(mean_bin) > 2
        R_mu = corr(x_c, mean_bin);
    else
        R_mu = NaN;
    end
    
    h = errorbar(x_c, mean_bin, sem_bin, '-o', 'Color', colors(i,:), ...
        'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'CapSize', 0);
    legend_handles(end+1) = h;

    %% BETA (18-24 Hz)
    R_beta = NaN;
    if size(feature_data, 2) >= feat_idx_beta
        curr_feat = feature_data(:, feat_idx_beta);
        mean_bin = []; sem_bin = []; x_c = [];
        
        for b = 1:length(edges)-1
            if N(b) > min_samples_bin
                vals = curr_feat(bin_idx == b);
                vals = vals(isfinite(vals));
                if ~isempty(vals)
                    mean_bin = [mean_bin; mean(vals)];
                    sem_bin = [sem_bin; std(vals)/sqrt(length(vals))];
                    x_c = [x_c; centers(b)];
                end
            end
        end
        
        % --- Corr ---
        if length(mean_bin) > 2
            R_beta = corr(x_c, mean_bin);
        end
        
        errorbar(x_c, mean_bin, sem_bin, '--s', 'Color', colors(i,:), ...
            'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'w', ...
            'CapSize', 0, 'HandleVisibility', 'off'); % HandleVisibility off per non sporcare la legenda
    end
    
    legend_labels{i} = sprintf('%s', subject);
    hist_matrix(:, i) = histcounts(prob_gmm_all, edges, 'Normalization', 'probability')';
end

grid on;
title('A) Feature-Confidence Correlation Trends', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Focality Feature Value');
xlabel('GMM Confidence P(IC)');
xlim([0 1]);
h_solid = plot(nan, nan, '-k', 'LineWidth', 2);
h_dash = plot(nan, nan, '--k', 'LineWidth', 2);
legend([legend_handles, h_solid, h_dash], [legend_labels, {'\mu-band (8-13 Hz)', '\beta-band (18-24 Hz)'}], ...
    'Location', 'northwest', 'FontSize', 9, 'NumColumns', 1, 'Interpreter', 'tex');


%% --- ISTOG ---
subplot(2, 1, 2);
b_handle = bar(centers, hist_matrix, 'grouped');
for i = 1:n_subjects
    b_handle(i).FaceColor = colors(i,:);
    b_handle(i).EdgeColor = 'none';
end

grid on;
title('B) GMM Confidence Distribution (Grouped by Subject)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('GMM Confidence P(IC)');
ylabel('Probability Density');
xlim([0 1]);
ylim([0 max(hist_matrix(:))*1.15]);
legend(subj_names, 'Location', 'north', 'Orientation', 'horizontal', 'Box', 'off');
sgtitle('GMM Unsupervised Gating Validation', 'FontSize', 14, 'FontWeight', 'bold');