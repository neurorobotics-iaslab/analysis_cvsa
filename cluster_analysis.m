% Assumes trial_data is [samples × channels × trials]
[n_samples, n_channels, n_trials] = size(trial_data);

% STEP 1: Concatenate all trials into one big [samples*trials × channels] matrix
all_data = reshape(trial_data, [n_samples * n_trials, n_channels]);

% Optional: z-score each channel (across all time)
all_data = zscore(all_data);

% STEP 2: Dimensionality reduction with PCA
[coeff, score, ~, ~, explained] = pca(all_data);

% Keep enough components to explain ~95% variance
cum_expl = cumsum(explained);
n_components = find(cum_expl > 95, 1);
reduced_data = score(:, 1:n_components);

% STEP 3: (Optional) t-SNE for visualization
% tsne_data = tsne(reduced_data, 'NumDimensions', 2, 'Perplexity', 30);

% STEP 4: Temporal clustering using k-means
K = 3; % number of states/clusters
[cluster_labels, C] = kmeans(reduced_data, K, 'Replicates', 10);

% STEP 5: Reshape cluster labels back into time x trial
cluster_per_trial = reshape(cluster_labels, [n_samples, n_trials]);

%% STEP 6: Plot temporal evolution of states for one trial
trial_id = 17;
figure;
plot(cluster_per_trial(:, trial_id), 'LineWidth', 2);
xlabel('Time (samples)');
ylabel('Cluster ID');
title(sprintf('Temporal States - Trial %d', trial_id));
ylim([0.5, K + 0.5]);
yticks(1:K);

%%
% Assumes trial_data is [samples × channels × trials]
[n_samples, n_channels, n_trials] = size(trial_data);

% STEP 1: Concatenate all trials into one big [samples*trials × channels] matrix
all_data = reshape(trial_data, [n_samples * n_trials, n_channels]);

% Optional: Normalize each channel (z-score)
all_data = zscore(all_data);  % Fast, can also do sliding-window in real-time

% RMS over sliding windows
win_size = 50;  % in samples (e.g., 200 ms @ 250 Hz)
features = movmean(trial_data.^2, win_size, 1);  % samples x channels x trials
features = sqrt(features);

% Then reshape and cluster same way
all_data = reshape(features, [n_samples * n_trials, n_channels]);


% STEP 2: Clustering (K-means directly on raw or z-scored EEG)
K = 3;  % number of states/phases
[cluster_labels, C] = kmeans(all_data, K, 'Replicates', 10);

% STEP 3: Reshape cluster labels back into [samples × trials]
cluster_per_trial = reshape(cluster_labels, [n_samples, n_trials]);

%% STEP 4: Visualize cluster timecourse for one trial
trial_id = 17;
figure;
plot(cluster_per_trial(:, trial_id), 'LineWidth', 2);
xlabel('Time (samples)');
ylabel('Cluster ID');
title(sprintf('Temporal States - Trial %d (Raw Features)', trial_id));
ylim([0.5, K + 0.5]);
yticks(1:K);
