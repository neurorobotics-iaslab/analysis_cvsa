%% file to compute for each trial the gini value in a window of 100 ms
clear all; % close all;

addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/EOG')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')

%% Initialization
bands = [{[8 14]} {[14 22]} {[22 30]}];
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
nbands = length(bands);
signals = cell(2, nbands);
headers = cell(2, nbands);
for i = 1:2
for idx_band = 1:nbands
    headers{i,idx_band}.TYP = [];
    headers{i,idx_band}.DUR = [];
    headers{i,idx_band}.POS = [];
    signals{i,idx_band} = [];
end
end
classes = [730 731];
cf_event = 781;
fix_event = 786;
nchannels = 39;
nclasses = length(classes);
filterOrder = 4;
avg = 1;
eog_threshold = 500;

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);
day = filenames{1}(4:11);

%% concatenate the files
nFiles = length(filenames);
trial_with_eog = [];
for idx_file= 1: nFiles
    fullpath_file_shift = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_shift);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;

    c_trial_with_eog = eog_detection(c_signal, header, eog_threshold, {'FP1', 'FP2', 'EOG'});
    trial_with_eog = [trial_with_eog; c_trial_with_eog];

    for idx_band = 1:nbands
        band = bands{idx_band};
        sampleRate = header.SampleRate;

        % for power band
        disp('   [proc] power band');
        disp('      [proc] applying filtering')
        [b, a] = butter(filterOrder, band(2)*(2/sampleRate),'low');
        s_low = filter(b,a,c_signal);
        [b, a] = butter(filterOrder, band(1)*(2/sampleRate),'high');
        s_filt = filter(b,a,s_low);
        disp('      [proc] applying power')
%         s_rect = power(s_filt, 2); ---> is better tot ake the power from hilbert
        analytic = hilbert(s_filt);
        s_rect = abs(analytic).^2;
        disp('      [proc] applying average window')
        s_out = zeros(size(c_signal));
        nchannels = size(c_signal, 2);
        for idx_ch=1:nchannels
            s_out(:, idx_ch) = (filter(ones(1,avg*sampleRate)/avg/sampleRate, 1, s_rect(:, idx_ch)));
        end

        signal_processed = s_out;

        c_header = headers{1, idx_band};
        c_header.sampleRate = header.SampleRate;
        c_header.channels_labels = header.Label;
        if isempty(find(header.EVENT.TYP == 2, 1)) % no eye calibration
            c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP);
            c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR);
            c_header.POS = cat(1, c_header.POS, header.EVENT.POS + size(signals{1, idx_band}, 1));
        else
            k = find(header.EVENT.TYP == 1, 1);
            c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP(k:end));
            c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR(k:end));
            c_header.POS = cat(1, c_header.POS, header.EVENT.POS(k:end) + size(signals{1, idx_band}, 1));
        end
        signals{1, idx_band} = cat(1, signals{1, idx_band}, signal_processed(:,:));
        headers{1, idx_band} = c_header;

        % for hilbert
        disp('      [proc] processing for hilbert')
        [b, a] = butter(filterOrder, band(2)*(2/header.SampleRate),'low');
        s_low = filter(b,a,c_signal);
        [b, a] = butter(filterOrder, band(1)*(2/header.SampleRate),'high');
        s_filt = filter(b,a,s_low);

        signal_processed = s_filt; % filtered signal

        c_header = headers{2,idx_band};
        c_header.sampleRate = header.SampleRate;
        c_header.channels_labels = header.Label;
        if isempty(find(header.EVENT.TYP == 2, 1)) % no eye calibration
            c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP);
            c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR);
            c_header.POS = cat(1, c_header.POS, header.EVENT.POS + size(signals{2,idx_band}, 1));
        else
            k = find(header.EVENT.TYP == 1, 1);
            c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP(k:end));
            c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR(k:end));
            c_header.POS = cat(1, c_header.POS, header.EVENT.POS(k:end) + size(signals{2,idx_band}, 1));
        end

        signals{2,idx_band} = cat(1, signals{2,idx_band}, signal_processed(:,:));
        headers{2,idx_band} = c_header;
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

%% computing the position of the electrons
LblMontage = cell(12,13);
LblMontage(1,:) = {'','','','','','FP1','FPZ','FP2','','','','',''};
LblMontage(2,:) = {'','AF9','AF7','','AF3','','AFZ','','AF4','','AF8','AF10',''};
LblMontage(3,:) = {'','F9','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F10',''};
LblMontage(4,:) = {'','FT9','FT7','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT8','FT10',''};
LblMontage(5,:) = {'','T9','T7','C5','C3','C1','CZ','C2','C4','C6','T8','T10',''};
LblMontage(6,:) = {'M1','TP9','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','TP10','M2'};
LblMontage(7,:) = {'','P9','P7','P5','P3','P1','PZ','P2','P4','P6','P8','P10',''};
LblMontage(8,:) = {'','PO9','PO7','PO5','PO3','PO1','POZ','PO2','PO4','PO6','PO8','PO10',''};
LblMontage(9,:) = {'','O9','','','','O1','OZ','O2','','','','O10',''};
LblMontage(10,:) = {'','','','','','','IZ','','','','','',''};
LblMontage(11,:) = {'','','','','','','','','','','','',''};
LblMontage(12,:) = {'','','','','','','EOG','','','','','',''};
[nRows, nCols] = size(LblMontage);
labelPos = containers.Map();
for r = 1:nRows
    for c = 1:nCols
        label = strtrim(LblMontage{r, c});
        if ~isempty(label)
            labelPos(label) = [c, r];  % (col, row)
        end
    end
end

% Get coordinates for each channel
coords_electrons = nan(nchannels, 2);
for i = 1:nchannels
    label = strtrim(channels_label{i});
    if isKey(labelPos, label)
        coords_electrons(i, :) = labelPos(label); % (y x)
    else
        warning('Label %s not found in LblMontage.', label);
    end
end

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
trial_data_logband = nan(min_trial_data, nbands, nchannels, ntrial);
trial_data_hilbert = nan(min_trial_data, nbands, nchannels, ntrial);
for idx_band = 1:nbands
    c_signal_logband = signals{1,idx_band};
    c_signal_hilbert = signals{2,idx_band};
    for trial = 1:ntrial
        c_start = trial_start(trial);
        c_end = trial_start(trial) + min_trial_data - 1;
        trial_data_logband(:,idx_band,:,trial) = c_signal_logband(c_start:c_end,:);
        trial_data_hilbert(:,idx_band,:,trial) = c_signal_hilbert(c_start:c_end,:);
    end
end

%% refactoring the data
% now the data are placed alternatevely, so odd trial id is for 730, even for 731
balanced_trial_idx = balanced_data_afterEOG(trial_with_eog, trial_typ, classes);
balanced_trial_data_logband = trial_data_logband(:,:,:,logical(balanced_trial_idx));
balanced_trial_data_hilbert = trial_data_hilbert(:,:,:,logical(balanced_trial_idx));
trial_typ = trial_typ(logical(balanced_trial_idx));
ntrial = sum(balanced_trial_idx);
idx_classes_trial = nan(ntrial/2, nclasses);
for idx_class = 1:nclasses
    idx_classes_trial(:,idx_class) = find(trial_typ == classes(idx_class));
end

tmp_logband = nan(size(balanced_trial_data_logband));
tmp_hilbert = nan(size(balanced_trial_data_hilbert));
trial_typ = nan(size(trial_typ));
i = 1;
for idx_trial_class = 1:2:ntrial
    for idx_class = 1:nclasses
        tmp_logband(:,:,:,idx_trial_class + idx_class - 1) = balanced_trial_data_logband(:,:,:,idx_classes_trial(i, idx_class));
        tmp_hilbert(:,:,:,idx_trial_class + idx_class - 1) = balanced_trial_data_hilbert(:,:,:,idx_classes_trial(i, idx_class));
        trial_typ(idx_trial_class + idx_class - 1) = classes(idx_class);
    end
    i = i + 1;
end
trial_data_logband = tmp_logband; % samples x bands x channels x trials
trial_data_hilbert = tmp_hilbert;
trial_data_logband(:,:,19,:) = 0; % remove the power of the EOG channel
trial_data_hilbert(:,:,19,:) = 0;

%% Applying hilbert
hilbert_data = nan(min_trial_data, nbands, nchannels, ntrial);  % complex analytic signal

for b = 1:nbands
    for ch = 1:nchannels
        for tr = 1:ntrial
            % Extract signal vector for that band, channel, trial
            signal = trial_data_hilbert(:, b, ch, tr);
            
            % Compute analytic signal via Hilbert transform
            c_h = hilbert(signal);
            
            % Store result
            hilbert_data(:, b, ch, tr) = c_h;
        end
    end
end

phase_data = angle(hilbert_data); % samples x band x channels x trial
amplitude_data = abs(hilbert_data);
power_data = amplitude_data .^ 2;

%% compute sparsity
% define regions
occipital = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'}; [~, ch_occipital] = ismember(occipital, channels_label);
central = {'CP3', 'CP1', 'C3', 'FC3', 'C1', 'FC1', 'CP4', 'C4', 'FC4', 'CP2', 'C2', 'FC2', 'FCZ', 'CZ'}; [~, ch_central] = ismember(central, channels_label);
frontal = {'F3', 'F1', 'F2', 'F4', 'FP1', 'FP2', 'FZ'}; [~, ch_frontal] = ismember(frontal, channels_label);
nchannelsSelected = size(ch_occipital, 2);

s_entropy_hilbert = nan(min_trial_data, nbands, ntrial); % samples x band x trial (for all metrics: 1 is focus in a point, 0 is sparse among the channels)
s_entropy_logband = nan(min_trial_data, nbands, ntrial);
s_m1_hilbert = nan(min_trial_data, nbands, ntrial); 
s_m1_logband = nan(min_trial_data, nbands, ntrial);
s_m2_logband = nan(min_trial_data, nbands, ntrial);
s_m2_hilbert = nan(min_trial_data, nbands, ntrial);
s_m3_logband = nan(min_trial_data, nbands, ntrial);
s_m3_hilbert = nan(min_trial_data, nbands, ntrial);

li_index_logband = nan(min_trial_data, nbands, ntrial, 3);
li_index_hilbert = nan(min_trial_data, nbands, ntrial, 3);
for t = 1:ntrial
    c_data_logband = squeeze(trial_data_logband(:,:,:,t)); % samples x band x channels
    c_data_hilbert = squeeze(power_data(:,:,:,t)); % samples x band x channels

    for sample = 1:min_trial_data
        c_sample_logband = squeeze(c_data_logband(sample,:,:)); % bands x channels
        c_sample_hilbert = squeeze(c_data_hilbert(sample,:,:)); % bands x channels

        for idx_band = 1:nbands
            tmp_logband = squeeze(c_sample_logband(idx_band,:)); % 1 x channels
%             tmp_logband = (tmp_logband - min(tmp_logband)); tmp_logband = tmp_logband / max(tmp_logband); % normalize

            tmp_hilbert = squeeze(c_sample_hilbert(idx_band,:));
%             tmp_hilbert = (tmp_hilbert - min(tmp_hilbert)); tmp_hilbert = tmp_hilbert / max(tmp_hilbert); % normalize
                
            [s_entropy_hilbert(sample,idx_band,t), s_m1_hilbert(sample,idx_band,t), s_m2_hilbert(sample, idx_band, t), s_m3_hilbert(sample, idx_band, t)] = compute_sparsity(tmp_hilbert);
            [s_entropy_logband(sample,idx_band,t), s_m1_logband(sample,idx_band,t), s_m2_logband(sample, idx_band, t), s_m3_logband(sample, idx_band, t)] = compute_sparsity(tmp_logband);

            data_li_logband = [mean(tmp_logband(ch_occipital)); mean(tmp_logband(ch_central)); mean(tmp_logband(ch_frontal));];
            data_li_hilbert = [mean(tmp_hilbert(ch_occipital)); mean(tmp_hilbert(ch_central)); mean(tmp_hilbert(ch_frontal));];
            i = 1;
            for idx_li1 = 1:length(data_li_logband)
                for idx_li2 = idx_li1+1:length(data_li_logband)
                    li_index_logband(sample, idx_band, t, i) = (data_li_logband(idx_li1) - data_li_logband(idx_li2)) / (data_li_logband(idx_li1) + data_li_logband(idx_li2));
                    li_index_hilbert(sample, idx_band, t, i) = (data_li_hilbert(idx_li1) - data_li_hilbert(idx_li2)) / (data_li_hilbert(idx_li1) + data_li_hilbert(idx_li2));
                    i = i + 1;
                end 
            end
        end
    end
end

%% show the log band and the gini val for both cases
handles = []; cl = -inf;
for t = 7:ntrial
    figure();
    for idx_band = 1:nbands
        subplot(5,nbands,idx_band)
        imagesc(squeeze(trial_data_logband(:,idx_band,:,t))')
        hold on;
        xline(minDurFix, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(minDurCue+minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        xticks(sampleRate:sampleRate:min_trial_data)
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        yticks(1:nchannels); yticklabels(channels_label)
        handles = [handles; gca];
        cl = max(cl, max(abs(squeeze(trial_data_logband(:,idx_band,:,t))), [], 'all'));
        title(['log band | ' bands_str{idx_band}])

        subplot(5,nbands,idx_band + nbands)
        imagesc(squeeze(power_data(:,idx_band,:,t))')
        hold on;
        xline(minDurFix, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(minDurCue+minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        xticks(sampleRate:sampleRate:min_trial_data)
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        yticks(1:nchannels); yticklabels(channels_label)
        handles = [handles; gca];
        cl = max(cl, max(abs(squeeze(trial_data_logband(:,idx_band,:,t))), [], 'all'));
        title(['hilbert power | ' bands_str{idx_band}])

        subplot(5,nbands,idx_band+2*nbands)
        plot(s_entropy_hilbert(:,idx_band,t))
        hold on;
        plot(s_entropy_logband(:,idx_band,t))
        hold off;
        legend([{'hilbert'},{'logband'}])
        hold on;
        xline(minDurFix, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(minDurCue+minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        xticks(sampleRate:sampleRate:min_trial_data)
        xlim([0 min_trial_data])
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        title(['sparsity entropy | ' bands_str{idx_band}])

        subplot(5,nbands,idx_band+3*nbands)
        tmp_h = s_m1_hilbert(:,idx_band,t);
%         tmp_h = tmp_h - min(tmp_h(minDurFix+ minDurCue + 1:end)); tmp_h = tmp_h / max(tmp_h(minDurFix+ minDurCue+1:end));
        tmp_l = s_m1_logband(:,idx_band,t);
%         tmp_l = tmp_l - min(tmp_l(minDurFix+ 1+minDurCue:end)); tmp_l = tmp_l / max(tmp_l(minDurFix+ minDurCue+1:end));
        plot(tmp_h)
        hold on;
        plot(tmp_l)
        hold off;
        legend([{'hilbert'},{'logband'}])
        hold on;
        xline(minDurFix, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(minDurCue+minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        xticks(sampleRate:sampleRate:min_trial_data)
        ylim([0 1])
        xlim([0 min_trial_data])
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        title(['sparsity m1 (count zone) | ' bands_str{idx_band}])

        subplot(5,nbands, idx_band + 4*nbands)
        plot(squeeze(s_m2_hilbert(:,idx_band, t)))
        hold on;
        plot(squeeze(s_m2_logband(:,idx_band, t)))
        xline(minDurFix, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(minDurCue+minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        legend([{'hilbert'},{'logband'}])
        xticks(sampleRate:sampleRate:min_trial_data)
        xlim([0 min_trial_data])
        ylim([0 1])
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        title(['sparsity m2 (hoyer) | ' bands_str{idx_band}])

        
    end
    sgtitle(['task: ' num2str(trial_typ(t))])

    figure();
    for idx_band = 1:nbands
        subplot(2,nbands, idx_band)
        plot(squeeze(li_index_logband(:,idx_band,t,:)))
        hold on;
        plot(squeeze(li_index_hilbert(:,idx_band,t,:)))
        hold off;
        legend([{'lb: occ-cent'}, {'lb: occ-fornt'}, {'lb cent-front'}, {'h: occ-cent'}, {'h: occ-fornt'}, {'h cent-front'}])
        hold on;
        xline(minDurFix, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(minDurCue+minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        xticks(sampleRate:sampleRate:min_trial_data)
        xlim([0 min_trial_data])
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        title(['li index | ' bands_str{idx_band}])

        subplot(2,nbands, idx_band + nbands)
        plot(squeeze(s_m3_logband(:,idx_band, t)))
        hold on;
        plot(squeeze(s_m3_hilbert(:,idx_band, t)))
        xline(minDurFix, '--r', 'Cue', 'LabelOrientation', 'horizontal');
        xline(minDurCue+minDurFix, '--r', 'Cf', 'LabelOrientation', 'horizontal');
        hold off;
        legend([{'band power'},{'hilbert'}])
        xticks(sampleRate:sampleRate:min_trial_data)
        xlim([0 min_trial_data])
        ylim([0, 1])
        xticklabels(string((sampleRate:sampleRate:min_trial_data) / sampleRate));
        title(['sparsity m3 (hoyer) | ' bands_str{idx_band}])
    end
end
set(handles, 'clim', [0, cl])



%%
[sparsityEntropy, spatialSpread, COM, s_m3_logband] = compute_sparsity(prova)

% sparsity and where
function [sparsityEntropy, m1, m2,m3] = compute_sparsity(signalVec)
    % Ensure column vector
    N = length(signalVec);
    % take the one which contribute at the 95% of the energy
    A2 = signalVec.^2;
    [A2_sorted, ~] = sort(A2, 'descend');
    cumulative = cumsum(A2_sorted); % Compute cumulative energy
    total = sum(A2_sorted);
    % Find how many components cover 95%
    threshold_index = find(cumulative >= 0.95 * total, 1, 'first');
    % Get the minimum value (in magnitude) that contributes to 95%
    value_cutoff = sqrt(A2_sorted(threshold_index));

    activeIdx = find(signalVec >= value_cutoff);

    % --- Entropy-Based Sparsity --- --> in general
    signalVec_entropy = signalVec(activeIdx).^2;
    total_entropy = sum(signalVec_entropy);
    if total_entropy == 0
        sparsityEntropy = NaN;
    else
        p = signalVec_entropy / total_entropy;
        H = -sum(p .* log2(p + eps));  % entropy
        Hmax = log2(N);
        sparsityEntropy = 1 - (H / Hmax);  % normalized: 0 to 1
    end
   
    % --- metodo 1 ---
    frontal = [1 2 3 4 5 20 21]; central = [6 7 8 9 10 11 12 22 23 24 25 26 27 28]; occipital = [13 14 15 16 17 18 29 30 31 32 33 34 35 36 37 38 39];
    c_l = [6 22 25 8 11 27]; c_r = [7 24 26 10 12 28]; o_l = [29 13 30 37 33 34 17]; o_r = [31 15 32 35 36 38 18];
    nFrontal = sum(ismember(activeIdx, frontal)); f_mean = mean(signalVec(ismember(activeIdx, frontal)));
    nCentral = sum(ismember(activeIdx, central)); c_mean = mean(signalVec(ismember(activeIdx, central)));
    nOccipital = sum(ismember(activeIdx, occipital)); o_mean = mean(signalVec(ismember(activeIdx, occipital)));
    total_mean = f_mean + c_mean + o_mean;
    n_total = nFrontal + nCentral + nOccipital;
    m1 = nFrontal/(nCentral + nOccipital) * f_mean/total_mean + nCentral/(nFrontal + nOccipital) * c_mean/total_mean + nOccipital/(nFrontal + nCentral) * o_mean/total_mean; 
    m1 = m1 * (nchoosek(length(frontal), nFrontal)*nchoosek(length(central), nCentral) * nchoosek(length(occipital), nOccipital)) / nchoosek(length(signalVec)-1, n_total);

    % --- metodo 2 ---
    mean_m2 = [mean(signalVec(activeIdx(ismember(activeIdx,frontal)))), mean(signalVec(activeIdx(ismember(activeIdx,central)))), mean(signalVec(activeIdx(ismember(activeIdx,occipital))))];
    mean_m2(isnan(mean_m2)) = 0;
    n = length(mean_m2);
    l1 = norm(mean_m2, 1);      % Norma L1
    l2 = norm(mean_m2, 2);      % Norma L2
    if l2 == 0
        m2 = 0;  % caso limite: vettore nullo
    else
        m2 = (sqrt(n) - (l1 / l2)) / (sqrt(n) - 1);
    end

%     signalVec_entropy = mean_m2.^2;
%     total_entropy = sum(signalVec_entropy);
%     if total_entropy == 0
%         m2 = NaN;
%     else
%         p = signalVec_entropy / total_entropy;
%         H = -sum(p .* log2(p + eps));  % entropy
%         Hmax = log2(n);
%         m2 = 1 - (H / Hmax);  % normalized: 0 to 1
%     end

    % --- metodo 3 ---
    m_f = max(signalVec(activeIdx(ismember(activeIdx,frontal))));
    m_c_l = max(signalVec(activeIdx(ismember(activeIdx,c_l))));
    m_c_r = max(signalVec(activeIdx(ismember(activeIdx,c_r))));
    m_o_l = max(signalVec(activeIdx(ismember(activeIdx,o_l))));
    m_o_r = max(signalVec(activeIdx(ismember(activeIdx,o_r))));
    s_m = sum([m_f, m_c_l, m_c_r, m_o_l, m_o_r]);
    mean_m3 = [mean(signalVec(activeIdx(ismember(activeIdx,frontal))))*m_f/s_m, mean(signalVec(activeIdx(ismember(activeIdx,c_l))))*m_c_l/s_m, ...
        mean(signalVec(activeIdx(ismember(activeIdx,c_r))))*m_c_r/s_m, mean(signalVec(activeIdx(ismember(activeIdx,o_l))))*m_o_l/s_m, ...
        mean(signalVec(activeIdx(ismember(activeIdx,o_r))))*m_o_r/s_m];
    mean_m3(isnan(mean_m3)) = 0;
    n = length(mean_m3);
    l1 = norm(mean_m3, 1);      % Norma L1
    l2 = norm(mean_m3, 2);      % Norma L2
    if l2 == 0
        m3 = 0;  % caso limite: vettore nullo
    else
        m3 = (sqrt(n) - (l1 / l2)) / (sqrt(n) - 1);
    end

%     signalVec_entropy = mean_m3.^2;
%     total_entropy = sum(signalVec_entropy);
%     if total_entropy == 0
%         m3 = NaN;
%     else
%         p = signalVec_entropy / total_entropy;
%         H = -sum(p .* log2(p + eps));  % entropy
%         Hmax = log2(n);
%         m3 = 1 - (H / Hmax);  % normalized: 0 to 1
%     end
end



