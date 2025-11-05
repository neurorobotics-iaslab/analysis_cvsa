%% file to check the power in a grand average fro the two classes
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
    fullpath_file = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file);
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

for t = 1:ntrial
    c_data = squeeze(trial_data(:,:,:,t)); % samples x band x channels

    for sample = 1:min_trial_data
        c_sample = squeeze(c_data(sample,:,:)); % bands x channels

        for idx_band = 1:nbands
            tmp = squeeze(c_sample(idx_band,:)); % 1 x channels

            sparsity(sample, idx_band, t,:) = compute_sparsity(tmp);
        end
    end
end

%% compute the sparsity for the mean trial divided by classes
sparsity_overall = nan(min_trial_data, nbands, nclasses, 3); 
for c = 1:nclasses
    t = trial_typ == classes(c);

    c_data = squeeze(mean(trial_data(:,:,:,t), 4)); % sample x band x channels

    for sample = 1:min_trial_data
        c_sample = squeeze(c_data(sample,:,:)); % bands x channels

        for idx_band = 1:nbands
            tmp = squeeze(c_sample(idx_band,:)); % 1 x channels
%             tmp = (tmp - min(tmp)); tmp = tmp / max(tmp); % normalize

            sparsity_overall(sample, idx_band, c,:) = compute_sparsity(tmp);
        end
    end
end

%% compute the sparsity for all the trial merged
sparsity_all = nan(min_trial_data, nbands, 3);
c_data = squeeze(mean(trial_data(:,:,:,:), 4));
for sample = 1:min_trial_data
    c_sample = squeeze(c_data(sample,:,:));

    for idx_band = 1:nbands
        tmp = squeeze(c_sample(idx_band,:)); % 1 x channels

        sparsity_all(sample, idx_band,:) = compute_sparsity(tmp);
    end
end


%% show log band and index ---> all trial
handles = cell(1, nbands); cl = -inf(1, nbands);
for c = 1:nclasses
    t = trial_typ == classes(c);
    figure();
    for idx_band = 1:nbands
        subplot(3,nbands,idx_band)
        imagesc(squeeze(mean(trial_data(:,idx_band,:,t), 4))')
        hold on;
        xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(min_durCUE+min_durFIX, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        xticks(sampleRate:sampleRate:min_trial_data)
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        yticks(1:nchannels); yticklabels(channels_label)
        handles{idx_band} = [handles{idx_band}, gca];
        cl(idx_band) = max(cl(idx_band), ...
            max(abs(squeeze(trial_data(:, idx_band, :, t))), [], 'all'));
        title(['log band | ' bands_str{idx_band}])


        subplot(3,nbands, idx_band + nbands)
        plot(squeeze(mean(sparsity(:,idx_band, t, 1),3)))
        hold on;
        plot(squeeze(mean(sparsity(:,idx_band, t, 2),3)))
        plot(squeeze(mean(sparsity(:,idx_band, t, 3),3)))
        xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(min_durCUE+min_durFIX, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        legend([{'bp mean'},{'bp max'},{'bp mean 2'},{'events'}])
        xticks(sampleRate:sampleRate:min_trial_data)
        xlim([0 min_trial_data])
        ylim([0, 1])
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        title(['mean of sparsity (entropy) over trials | ' bands_str{idx_band}])

        subplot(3,nbands, idx_band + 2*nbands)
        plot(squeeze(sparsity_overall(:,idx_band, c, 1)))
        hold on;
        plot(squeeze(sparsity_overall(:,idx_band, c, 2)))
        plot(squeeze(sparsity_overall(:,idx_band, c, 3)))
        xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(min_durCUE+min_durFIX, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        legend([{'bp mean'},{'bp max'},{'bp mean 2'},{'events'}])
        xticks(sampleRate:sampleRate:min_trial_data)
        xlim([0 min_trial_data])
        ylim([0, 1])
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        title(['sparsity (entropy) | ' bands_str{idx_band}])
        
    end
    sgtitle(['class: ' num2str(classes(c))])
end

figure();
for idx_band = 1:nbands
    subplot(3, nbands, idx_band)
    imagesc(squeeze(mean(trial_data(:,idx_band,:,:), 4))')
    hold on;
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durCUE+min_durFIX, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    xticks(sampleRate:sampleRate:min_trial_data)
    xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
    yticks(1:nchannels); yticklabels(channels_label)
    handles{idx_band} = [handles{idx_band}, gca];
    cl(idx_band) = max(cl(idx_band), ...
        max(abs(squeeze(trial_data(:, idx_band, :, :))), [], 'all'));
    title(['log band | ' bands_str{idx_band}])

    subplot(3,nbands, idx_band + nbands)
    plot(squeeze(mean(sparsity(:, idx_band, :, 1),3)))
    hold on;
    plot(squeeze(mean(sparsity(:, idx_band, :, 2),3)))
    plot(squeeze(mean(sparsity(:, idx_band, :, 3),3)))
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durCUE+min_durFIX, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    legend([{'bp mean'},{'bp max'},{'bp mean 2'},{'events'}])
    xticks(sampleRate:sampleRate:min_trial_data)
    xlim([0 min_trial_data])
    ylim([0, 1])
    xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
    title(['mean of sparsity (entropy) over trials | ' bands_str{idx_band}])

    subplot(3,nbands, idx_band + 2*nbands)
    plot(squeeze(sparsity_all(:, idx_band, 1)))
    hold on;
    plot(squeeze(sparsity_all(:, idx_band, 2)))
    plot(squeeze(sparsity_all(:, idx_band, 3)))
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durCUE+min_durFIX, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    legend([{'bp mean'},{'bp max'},{'bp mean 2'},{'events'}])
    xticks(sampleRate:sampleRate:min_trial_data)
    xlim([0 min_trial_data])
    ylim([0, 1])
    xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
    title(['sparsity (entropy) | ' bands_str{idx_band}])
end
sgtitle('Mean without class division')

figure();
for idx_band = 1:nbands
    subplot(1,nbands,idx_band)
    diff = (squeeze(mean(trial_data(:,idx_band,:,trial_typ == classes(1)), 4)) - squeeze(mean(trial_data(:,idx_band,:,trial_typ == classes(2)), 4)))';
    imagesc(diff)
    hold on;
    xline(min_durFIX, '--r', 'Cue', 'LabelOrientation', 'horizontal');
    xline(min_durCUE+min_durFIX, '--r', 'Cf', 'LabelOrientation', 'horizontal');
    hold off;
    xticks(sampleRate:sampleRate:min_trial_data)
    xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
    yticks(1:nchannels); yticklabels(channels_label)
    cl = max(abs(diff), [], 'all');
    set(gca, 'clim', [-cl cl])
    title(['log band | DIFF: ' num2str(classes(1)) '-' num2str(classes(2)) ' | ' bands_str{idx_band}])
end


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