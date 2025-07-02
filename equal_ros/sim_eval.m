%% in this file we use the evaluation file to check if the cf will hit or miss. We use the same processing used in ros
clc; clearvars;% close all;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/utils');
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils');

%% variables
bands = {[8 10], [10 12], [12 14], [14 16], [8 14]};
avg = 1;
bufferSize_integrator = 48;
chunkSize = 32;
nchannels = 39;
bufferSize_processing = 512;
filterOrder = 4;
classes = [730 731];
fix_event = 786;
cf_event = 781;
hit_event = 897;
miss_event = 898;
timeout_event = 899;
startTrial_event = 786; %% we use only data inside fix to boom
nclasses = length(classes);
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
nbands = length(bands);
signals = cell(1, nbands);
headers = cell(1, nbands);
for idx_band = 1:nbands
    headers{idx_band}.TYP = [];
    headers{idx_band}.DUR = [];
    headers{idx_band}.POS = [];
    headers{idx_band}.startNewFile = [];
    headers{idx_band}.ths_rejection = [];
    headers{idx_band}.ths = [];
    headers{idx_band}.bufferSize_integrator = [];
    headers{idx_band}.alpha = [];
    signals{idx_band} = [];
end

%% load the gdfs
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF evaluation Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end

%% Concatenate all + processing + extract selected channels
nFiles = length(filenames);
hit = zeros(nFiles, 1); miss = zeros(nFiles, 1); timeout = zeros(nFiles, 1);
for idx_file= 1: nFiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [signal,header] = sload(fullpath_file);
    signal = signal(:,1:nchannels);
    channels_label = header.Label;
    [~, idx_channels_select] = ismember(channels_select, channels_label);
    % load the parameters used in ros
    parameters_file = [pathname(1:end-4) 'parameters/' filenames{idx_file}(1:end-3) 'yaml'];
    a = ReadYaml(parameters_file);
    ths_rejection = a.integrator.thresholds_rejection;
    ths = a.trainingCVSA_node.thresholds;
    if contains(filenames{idx_file}, 'expo')
        alpha = a.integrator.alpha;
    else
        bufferSize = a.integrator.buffer_size*512/16;
    end

    for idx_band = 1:nbands
        band = bands{idx_band};

        signal_processed = proc_512hz(signal, header.SampleRate, band, filterOrder, avg);
        
        c_header = headers{idx_band};
        c_header.sampleRate = header.SampleRate;
        c_header.channels_labels = header.Label;
        c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP);
        c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR);
        c_header.POS = cat(1, c_header.POS, header.EVENT.POS + size(signals{idx_band}, 1));
        c_header.startNewFile = cat(1, c_header.startNewFile, size(signals{idx_band}, 1) + 1);

        signals{idx_band} = cat(1, signals{idx_band}, signal_processed(:,:));
        c_header.ths = cat(1, c_header.ths, repmat(ths,size(signal_processed, 1), 1));
        c_header.ths_rejection = cat(1, c_header.ths_rejection, repmat(ths_rejection, size(signal_processed, 1), 1));
        if contains(filenames{idx_file}, 'expo')
            c_header.alpha = cat(1, c_header.alpha, repmat(alpha, size(signal_processed, 1), 1));
            c_header.bufferSize_integrator = cat(1, c_header.bufferSize_integrator, nan(size(signal_processed, 1), 1));
        else
            c_header.bufferSize_integrator = cat(1, c_header.bufferSize_integrator, repmat(bufferSize, size(signal_processed, 1), 1));
            c_header.alpha = cat(1, c_header.alpha, nan(size(signal_processed, 1), 1));
        end
        headers{idx_band} = c_header;
        headers{idx_band} = c_header;
    end

    % for the accuracy
    hit(idx_file) = sum(header.EVENT.TYP == hit_event);
    miss(idx_file) = sum(header.EVENT.TYP == miss_event);
    timeout(idx_file) = sum(header.EVENT.TYP == timeout_event);
end
subject = filenames{1}(1:2);

%% extract all data of the cf
% now data contains all the cf data and the info give the positions and
% durations
cf_data = cell(1, nbands); 
for idx_band = 1:nbands
    [cf_data{idx_band}, cf_info] = extract_event(signals{idx_band}, headers{idx_band}, classes, [cf_event]);
end

%% load QDA and extract features
[filename_qda, path_qda] = uigetfile('*.yaml;*.yml', 'Select a YAML of the QDA File');
path_qda = fullfile(path_qda, filename_qda);
qda = loadQDA(path_qda);

%% features extraction
% bands_interest =[[8 14]; [8 14];[8 14];[8 14];[8 14];[8 14];[8 14];[8 14];[8 14];[8 14];[8 14];[8 14];[8 14];[8 14];[8 14];[8 14];[8 14];];
% [~,idx_select] = ismember(channels_select, headers{1}.channels_labels);
% features = features_extraction(cf_data, bands, bands_interest, idx_select);
features = features_extraction(cf_data, bands, qda.bands, qda.idchans);
nsamples = size(features, 1);

%% compute the classification and the integrator
idx_trial = 1;
idx_integrator = 1;
integrator_buffer = [];
ntrial = size(cf_info.startEvent, 1);
prob_integrated_B = inf(nsamples, nclasses);
probs = inf(nsamples, nclasses);
for idx_feature=1:nsamples
    % reset
    if idx_feature == cf_info.startEvent(idx_trial)
        if idx_trial == ntrial
            idx_trial = ntrial;
        else
            idx_trial = idx_trial + 1;
        end
        if ~isnan(cf_info.bufferSize_integrator(idx_trial))
            bufferSize_integrator = cf_info.bufferSize_integrator(idx_trial);
        else
            bufferSize_integrator = 48;
        end
        idx_integrator = 1;
        integrator_buffer = repmat(classes, 1, bufferSize_integrator/nclasses);
        ths_rejection = cf_info.ths_rejection(idx_trial,:);
    end
    if idx_integrator > bufferSize_integrator
        idx_integrator = 1;
    end

    % classify
    feature = features(idx_feature,:);
    prob = apply_qda(qda, feature);
    [~, idxMax] = max(prob);
    if prob(idxMax) > ths_rejection(idxMax)
        integrator_buffer(idx_integrator) = classes(idxMax);
        idx_integrator = idx_integrator + 1;
    end
    prob_integ = nan(1, nclasses);
    for idx_class=1:nclasses
        count = sum(integrator_buffer == classes(idx_class));
        prob_integ(idx_class) = count/bufferSize_integrator;
    end
    probs(idx_feature,:) = prob;
    prob_integrated_B(idx_feature,:) = prob_integ;
end

%% Compute the integrator exponential
prob_integrated_E = inf(nsamples, nclasses);
idx_trial = 1;
for idx_feature=1:nsamples
    % reset
    if idx_feature == cf_info.startEvent(idx_trial)
        if idx_trial == ntrial
            idx_trial = ntrial;
        else
            idx_trial = idx_trial + 1;
        end
        if ~isnan(cf_info.alpha(idx_trial))
            alpha = cf_info.alpha(idx_trial);
        else
            alpha = 0.97;
        end
        c_prob_integ = [0.5 0.5];
        ths_rejection = cf_info.ths_rejection(idx_trial,:);
    end

    % classify
    feature = features(idx_feature,:);
    prob = apply_qda(qda, feature);
    [~, idxMax] = max(prob);

    if prob(idxMax) > ths_rejection(idxMax)
        prob = zeros(1, nclasses);
        prob(idxMax) = 1.0;
        c_prob_integ = (1.0 - alpha)*prob + alpha*c_prob_integ;
    end
    prob_integrated_E(idx_feature,:) = c_prob_integ;
end

%% plot the cf in time (prob integrated)
last_trials = 20;
time_plot = 0.5;
nrows = 3;
figure();
for idx_trial=ntrial-last_trials+1:ntrial
    ths = cf_info.ths(idx_trial,:);
    start_cf = cf_info.startEvent(idx_trial);
    end_cf = cf_info.endEvent(idx_trial);
    ths_rejection = cf_info.ths_rejection(idx_trial,:);
    
    c_prob = probs(start_cf:end_cf, :);
    c_prob = cat(1, nan(1, nclasses), c_prob);

    c_prob_integrated_B = prob_integrated_B(start_cf:end_cf,:);
    c_prob_integrated_B = cat(1, repmat(1/nclasses, 1, nclasses), c_prob_integrated_B);

    c_prob_integrated_E = prob_integrated_E(start_cf:end_cf,:);
    c_prob_integrated_E = cat(1, repmat(1/nclasses, 1, nclasses), c_prob_integrated_E);

    subplot(nrows, ceil(last_trials/nrows), idx_trial - (ntrial - last_trials));
    hold on;
    x = 1:end_cf-start_cf + 2;
    plot(c_prob_integrated_B(:,1));
    plot(c_prob_integrated_E(:,1));
    yline(ths_rejection(1), '--r')
    yline(1.0 - ths_rejection(2), '--r')
    yline(ths(1), '-g');
    yline(1.0 - ths(2), '-g');
    scatter(x, c_prob(:,1), 5, 'k', 'filled');
    title({['cue: ' num2str(cf_data{1}.typ(idx_trial)) ' | trial ' num2str(idx_trial)] [' | ' num2str(cf_info.hit(idx_trial))]});
    ylim([0,1])
    xticks_ = [(1:cf_info.sampleRate*time_plot:max(x))-1 max(x)];
    xticks(xticks_)
    xticklabels(string((xticks_)/cf_info.sampleRate))
    drawnow;
    hold off;
end
legend({'integrated prob of 730 with buffer', 'integrated prob of 730 with exponential', ...
        'rejection 730', 'rejection 731', 'threhsold 730', 'threhsold 731', 'qda probs'}, 'Location', 'best');
sgtitle(['subject: ' subject ' | cf probabilities'])

%% plot the features for each last trial
% plot with imagesc
figure();
% nfeatures = qda.nfeatures;
nfeatures = size(features, 2);
for idx_trial=ntrial-last_trials+1:ntrial
    start_trial = cf_info.startEvent(idx_trial);
    end_trial = cf_info.endEvent(idx_trial);
    c_features = features(start_trial:end_trial,:);

    x = 1:end_trial-start_trial + 1;
    leg = cell(nfeatures, 1);
    subplot(nrows, ceil(last_trials/nrows), idx_trial- (ntrial-last_trials))
    hold on
    imagesc(c_features')
    hold off
    xlim([0 size(c_features, 1)])
    xticks_ = [(1:cf_info.sampleRate*time_plot:max(x))-1 max(x)];
    xticks(xticks_)
    xticklabels(string((xticks_)/cf_info.sampleRate))
    ylim([0.5, nfeatures+0.5])
    yticks(1:nfeatures)
    yticklabels([headers{1}.channels_labels(qda.idchans)])
%     yticklabels([headers{1}.channels_labels(idx_select)])
    title({['cue: ' num2str(cf_data{1}.typ(idx_trial)) ' | trial ' num2str(idx_trial)] [' | ' num2str(cf_info.hit(idx_trial))]});
    drawnow;
end
sgtitle(['subject: ' subject ' | features received for the QDA']);

%% Compute accuracy raw prob
count_allProb = zeros(nclasses, 4); % general, hit, miss, timeout
count_correctProb = zeros(nclasses, 4);

for idx_trial=1:ntrial
    start_cf = cf_info.startEvent(idx_trial);
    end_cf = cf_info.endEvent(idx_trial);

    c_cueTYP = cf_data{1}.typ(idx_trial);
    c_prob = probs(start_cf:end_cf, :);
    idx_cue = find(classes == c_cueTYP);
    c_prob = c_prob(:,idx_cue);
    ths_rejection = cf_info.ths_rejection(idx_trial,idx_cue);

    count_allProb(idx_cue, 1) = count_allProb(idx_cue, 1) + size(c_prob, 1);
    count_correctProb(idx_cue, 1) = count_correctProb(idx_cue, 1) + sum(c_prob > ths_rejection);

    if cf_info.hit(idx_trial) == hit_event
        count_allProb(idx_cue, 2) = count_allProb(idx_cue, 2) + size(c_prob, 1);
        count_correctProb(idx_cue, 2) = count_correctProb(idx_cue, 2) + sum(c_prob > ths_rejection);
    elseif cf_info.hit(idx_trial) == miss_event
        count_allProb(idx_cue, 3) = count_allProb(idx_cue, 3) + size(c_prob, 1);
        count_correctProb(idx_cue, 3) = count_correctProb(idx_cue, 3) + sum(c_prob > ths_rejection);
    elseif cf_info.hit(idx_trial) == timeout_event
        count_allProb(idx_cue, 4) = count_allProb(idx_cue, 4) + size(c_prob, 1);
        count_correctProb(idx_cue, 4) = count_correctProb(idx_cue, 4) + sum(c_prob > ths_rejection);
    end
end
disp('Raw accuracy');
for idx_class=1:nclasses
    disp(['   class: ' num2str(classes(idx_class))]);
    disp(['      Accuracy: ' num2str(count_correctProb(idx_class, 1)/count_allProb(idx_class, 1)) ...
        ' | correct: ' num2str(count_correctProb(idx_class, 1)) ' all: ' num2str(count_allProb(idx_class, 1))]);
    disp(['      Accuracy while hit: ' num2str(count_correctProb(idx_class, 2) / count_allProb(idx_class, 2)) ' | ' ...
        'correct: ' num2str(count_correctProb(idx_class, 2)) ' all: ' num2str(count_allProb(idx_class, 2))]);
    disp(['      Accuracy while miss: ' num2str(count_correctProb(idx_class, 3) / count_allProb(idx_class, 3)) ' | ' ...
        'correct: ' num2str(count_correctProb(idx_class, 3)) ' all: ' num2str(count_allProb(idx_class, 3))]);
    disp(['      Accuracy while timeout: ' num2str(count_correctProb(idx_class, 4) / count_allProb(idx_class, 4)) ' | ' ...
        'correct: ' num2str(count_correctProb(idx_class, 4)) ' all: ' num2str(count_allProb(idx_class, 4))]);
end
disp(['   General ' num2str(sum(count_correctProb(:,1))/sum(count_allProb(:,1))) ]);

%% Compute accuracy trial
disp('Trial acuracy')
disp('   For file:')
for idx_file=1:nFiles
    c_hit = hit(idx_file); c_miss = miss(idx_file); c_timeout = timeout(idx_file);
    disp(['      Accuracy: ' num2str(c_hit/(c_miss + c_timeout + c_hit)) ' Accuracy rejection: ' num2str(c_hit/(c_hit + c_miss)) ...
        ' | hit: ' num2str(c_hit) ' miss: ' num2str(c_miss) ' timeout:' num2str(c_timeout)])
end
count_allProb = [sum(hit) sum(miss) sum(timeout)];
disp(['   Mean all file: Accuracy: ' num2str(count_allProb(1)/sum(count_allProb)) ...
    ' Accuracy rejection: ' num2str(count_allProb(1)/(count_allProb(1) + count_allProb(2))) ...
    ' | hit: ' num2str(count_allProb(1)) ' miss: ' num2str(count_allProb(2)) ' timeout: ' num2str(count_allProb(3))])

%% QDA analysis
y_true_label = nan(nsamples, 1);
for idx_trial=1:ntrial
    c_start = cf_info.startEvent(idx_trial);
    c_end = cf_info.endEvent(idx_trial);
    if cf_data{1}.typ(idx_trial) == classes(1)
        y_true_label(c_start:c_end) = classes(1);
    else
        y_true_label(c_start:c_end) = classes(2);
    end
end

y_predict_label = nan(nsamples, 1);
positive_class = 1;
negative_class = 2;
for idx_sample=1:nsamples
    feature = features(idx_sample,:);
    prob = apply_qda(qda,feature);
    if prob(positive_class) > 1/nclasses
        y_predict_label(idx_sample) = classes(positive_class);
    else
        y_predict_label(idx_sample) = classes(negative_class);
    end
end

conf_matrix = confusionmat(y_true_label, y_predict_label); % [TN FP; FN TP]
TN = conf_matrix(1,1); FP = conf_matrix(1,2); FN = conf_matrix(2,1); TP = conf_matrix(2,2);

% some metrics
disp(['Accuracy: ', num2str((TP+TN)/(TP + TN + FP + FN))]);
precision = TP/(TP + FP);
disp(['Precision: ', num2str(precision)]);
recall = TP/(TP + FN);
disp(['Recall: ', num2str(recall)]);
disp(['F1 score: ', num2str(2*((precision*recall)/(precision+recall)))])
disp(['True negative rate: ', num2str(TN/(TN + FP))])
disp(['False positive rate; ', num2str(FP/(FP + TN))])
disp(['False negative rate: ', num2str(FN/(FN + TP))])
disp(['Balanced Accuracy: ', num2str((recall + TN/(TN + FP))/2)]);
disp(['MCC: ', num2str((TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))]);

% F1-score
[precision, recall, ~] = perfcurve(y_true_label, probs(:, positive_class), classes(positive_class), 'XCrit', 'TPR', 'YCrit', 'PPV');
f1_score = 2 * (precision .* recall) ./ (precision + recall);
disp(['F1-Score: ', num2str(max(f1_score))]);


% ROC Curve
[roc_x, roc_y, ~, auc, optrocpt] = perfcurve(y_true_label, probs(:, positive_class), classes(positive_class));
figure;
subplot(1, 2, 1);
hold on;
plot(roc_x, roc_y, 'b-', 'LineWidth', 2);
plot([0 1], [0, 1], 'r-', 'LineWidth', 2)
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(['ROC Curve (AUC = ' num2str(auc) ')']);
hold off;
grid on;

% precision-recall curve
[prec, recall, ~, auc_pr] = perfcurve(y_true_label, probs(:, positive_class), classes(positive_class), 'XCrit', 'reca', 'YCrit', 'prec');
subplot(1, 2, 2);
plot(recall, prec, 'r-', 'LineWidth', 2);
xlabel('Recall');
ylabel('Precision');
title(['Precision-Recall Curve (AUC = ' num2str(auc_pr) ')']);
grid on;