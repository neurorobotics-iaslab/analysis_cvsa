%% in this file we use the evaluation file to check if the cf will hit or miss. We use the same processing used in ros
clc; clearvars; close all;

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
    signals{idx_band} = [];
end

%% load the gdfs
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF evaluation Files', 'MultiSelect', 'on');

% Concatenate all + processing + extract selected channels
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
    bufferSize = a.integrator.buffer_size;

    for idx_band = 1:nbands
        band = bands{idx_band};
        [signal_processed, header_processed] = processing_offline(signal, header, nchannels, bufferSize_processing, filterOrder, band, chunkSize);
        c_header = headers{idx_band};
        c_header.sampleRate = header_processed.SampleRate;
        c_header.channels_labels = header.Label;
        c_header.TYP = cat(1, c_header.TYP, header_processed.EVENT.TYP);
        c_header.DUR = cat(1, c_header.DUR, header_processed.EVENT.DUR);
        c_header.POS = cat(1, c_header.POS, header_processed.EVENT.POS + size(signals{idx_band}, 1));
        c_header.startNewFile = cat(1, c_header.startNewFile, size(signals{idx_band}, 1) + 1);

        signals{idx_band} = cat(1, signals{idx_band}, signal_processed(:,:));
        c_header.ths = cat(1, c_header.ths, repmat(ths,size(signal_processed, 1), 1));
        c_header.ths_rejection = cat(1, c_header.ths_rejection, repmat(ths_rejection, size(signal_processed, 1), 1));
        c_header.bufferSize_integrator = cat(1, c_header.bufferSize_integrator, repmat(bufferSize, size(signal_processed, 1), 1));
        headers{idx_band} = c_header;
    end

    % for the accuracy
    hit(idx_file) = sum(header.EVENT.TYP == hit_event);
    miss(idx_file) = sum(header.EVENT.TYP == miss_event);
    timeout(idx_file) = sum(header.EVENT.TYP == timeout_event);
end

%% extract all data
% for each band the signal are: samples x channels x trial
%   in data the signal keep this division for data.trial, data.cue data.cf data.fix
data = cell(1, nbands); 
for idx_band = 1:nbands
    [data{idx_band}, info] = extract_cf(signals{idx_band}, headers{idx_band}, classes, cf_event);
end

%% config the integrator
[filename_qda, path_qda] = uigetfile('*.yaml;*.yml', 'Select a YAML of the QDA File');
path_qda = fullfile(path_qda, filename_qda);
% path_qda = '/home/paolo/cvsa_ws/src/qda_cvsa/cfg/qda_c7_lbp_20250131.yaml';
qda = loadQDA(path_qda);

%% compute the classification and the integrator
features = features_extraction(data, bands, qda.bands, qda.idchans);
nfeature = size(features, 1);
idx_trial = 1;
idx_integrator = 1;
sintegrator_buffer = [];
ntrial = size(info.startCf, 1);
prob_integrated_B = inf(nfeature, nclasses);
probs = inf(nfeature, nclasses);
for idx_feature=1:nfeature
    % reset
    if idx_feature == info.startCf(idx_trial)
        if idx_trial == ntrial
            idx_trial = ntrial;
        else
            idx_trial = idx_trial + 1;
        end
        bufferSize_integrator = info.bufferSize_integrator(idx_trial);
        idx_integrator = 1;
        integrator_buffer = repmat(classes, 1, bufferSize_integrator/nclasses);
        ths_rejection = info.ths_rejection(idx_trial,:);
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
prob_integrated_E = inf(nfeature, nclasses);
alpha = 0.97;
idx_trial = 1;
for idx_feature=1:nfeature
    % reset
    if idx_feature == info.startCf(idx_trial)
        if idx_trial == ntrial
            idx_trial = ntrial;
        else
            idx_trial = idx_trial + 1;
        end
        c_prob_integ = [0.5 0.5];
        ths_rejection = info.ths_rejection(idx_trial,:);
    end
    if idx_integrator > bufferSize_integrator
        idx_integrator = 1;
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

%% plot the trials
last_trials = 20;
time_plot = 0.5;
for idx_trial=ntrial-last_trials+1:ntrial
    ths = info.ths(idx_trial,:);
    start_cf = info.startCf(idx_trial);
    end_cf = info.endCf(idx_trial);
    
    c_prob = probs(start_cf:end_cf, :);
    c_prob = cat(1, nan(1, nclasses), c_prob);

    c_prob_integrated_B = prob_integrated_B(start_cf:end_cf,:);
    c_prob_integrated_B = cat(1, repmat(1/nclasses, 1, nclasses), c_prob_integrated_B);

    c_prob_integrated_E = prob_integrated_E(start_cf:end_cf,:);
    c_prob_integrated_E = cat(1, repmat(1/nclasses, 1, nclasses), c_prob_integrated_E);

    figure();
    hold on;
    x = 1:end_cf-start_cf + 2;
    plot(c_prob_integrated_B(:,1));
    plot(c_prob_integrated_E(:,1));
    yline(ths_rejection(1), '--r')
    yline(1.0 - ths_rejection(2), '--r')
    yline(ths(1), '-g');
    yline(1.0 - ths(2), '-g');
    scatter(x, c_prob(:,1), 5, 'k', 'filled');
    legend({'integrated prob of 730 with buffer', 'integrated prob of 730 with exponential', ...
        'rejection 730', 'rejection 731', 'threhsold 730', 'threhsold 731', 'qda probs'}, 'Location', 'best');
    
    title(['class asked: ' num2str(data{1}.typ(idx_trial)) ' | trial ' num2str(idx_trial) ' | ' num2str(info.hit(idx_trial))]);
    ylim([0,1])

    xticks_ = [(1:info.sampleRate*time_plot:max(x))-1 max(x)];
    xticks(xticks_)
    xticklabels(string((xticks_)/info.sampleRate))
    hold off;
end

%% Compute accuracy
for idx_file=1:nFiles
    c_hit = hit(idx_file); c_miss = miss(idx_file); c_timeout = timeout(idx_file);
    disp(['Accuracy overall: ' num2str(c_hit/(c_miss + c_timeout + c_hit)) ' Accuracy rejection: ' num2str(c_hit/(c_hit + c_miss)) ...
        ' | hit: ' num2str(c_hit) ' miss: ' num2str(c_miss) ' timeout:' num2str(c_timeout)])
end
disp('ALL FILE');
all_hit = sum(hit); all_miss = sum(miss); all_timeout = sum(timeout);
disp(['Accuracy overall: ' num2str(all_hit/(all_miss + all_timeout + all_hit)) ' Accuracy rejection: ' num2str(all_hit/(all_hit + all_miss)) ...
    ' | hit: ' num2str(all_hit) ' miss: ' num2str(all_miss) ' timeout:' num2str(all_timeout)])