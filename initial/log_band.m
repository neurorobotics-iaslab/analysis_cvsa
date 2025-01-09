clc; clear all; close all force;
% compute the log band and the fisher score related for each
% sample computed on log band signals. In addition compute the PSD for all channels

%% files
subject = 'h8';
day = '/20241015';
path = ['/home/paolo/cvsa_ws/record/' subject day '/gdf/calibration'];
files = dir(fullfile(path, '*.gdf'));

channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
        '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

% channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

lap_path = '/home/paolo/laplacians/lap_39ch_CVSA.mat';
load(lap_path);

classes = [730 731];
sampleRate = 512;

normalization = false; % if false log band power otherwise ERD/ERS


%% log band
for idx_f = 1:length(files)
    disp(['[STARTING] file: ' files(idx_f).name])
    file = fullfile(path, files(idx_f).name);

    [signal, header] = sload(file);

    %% Extract infoo data
    nchannels = length(channels_label);
    signal = signal(:, 1:nchannels);
    [fixPOS, fixDUR, cuePOS, cueDUR, cueTYP, cfPOS, cfDUR, n_trial] = extract_info_data(header, 786, [730 731], 781);

    %% Extract trial infoo
    [trialStart, trialStop, fixStart, fixStop, ck, tk] = extract_trial_infoo(signal, fixPOS, ...
        fixDUR, cfPOS, cfDUR, cueTYP, n_trial);

    %% Processing
    band = [8 14];
    avg = 1;
    filtOrder = 4;
    signal_processed = processing_offline(signal, nchannels, sampleRate, band, filtOrder, avg);

    %% log band
    [ERD, minDur, minDurFix] = compute_ERDERS(signal_processed, trialStart, trialStop, fixStart, fixStop, nchannels, n_trial, ck, normalization);

    disp('[INFO] computing difference (bottom right - bottom left) in the ERD/ERS');
    img = (mean(ERD(:,:,cueTYP == classes(2)), 3) - mean(ERD(:,:,cueTYP == classes(1)), 3))';
    img = img(find(~strcmp(channels_label, '')),:);

    xcue = floor(mean(cuePOS - fixPOS));
    xcf  = floor(mean(cfPOS - fixPOS));

    figure();
    imagesc(img);
    colorbar;
    set(gca, 'CLim', [-0.6 0.6])
    line([xcue xcue], [0 40], 'Color', 'black', 'LineWidth', 2);
    line([xcf xcf], [0 40], 'Color', 'black', 'LineWidth', 2);
    xlabel('Time [s]')
    ylabel('Channels');
    xticks(0:sampleRate:minDur);
    xticklabels((0:sampleRate:minDur)/sampleRate);
    yticks(1:numel(channels_label(find(~strcmp(channels_label, '')))));
    yticklabels(channels_label(find(~strcmp(channels_label, ''))));
    
    if ~normalization
        title('Feature Map-- log band power -- bottom right - bottom left');
    else
        title('Feature Map-- ERD/ERS -- bottom right - bottom left');
    end
    

    %% fisher with log band
    cmu = zeros(nchannels, length(classes));
    csigma = zeros(nchannels, length(classes));
    fisher_ERD = zeros(nchannels, size(ERD, 1));

    for idx_t = 1:size(ERD, 1)

        for i= 1:length(classes)
            cl = classes(i);

            c_ERD = rearrange(ERD(idx_t,:,cueTYP==cl));
            cmu(:,i) = mean(c_ERD, 1);
            csigma(:,i) = std(c_ERD, [], 1);
        end


        fisher_ERD(:,idx_t) = abs(cmu(:, 2) - cmu(:, 1)) ./ sqrt(( csigma(:, 1).^2 + csigma(:, 2).^2 ));
    end

    fisher_ERD = fisher_ERD(find(~strcmp(channels_label, '')), :);

    xcue = floor(mean(cuePOS - fixPOS));
    xcf  = floor(mean(cfPOS - fixPOS));

    figure();
    imagesc(fisher_ERD);
    line([xcue xcue], [0 40], 'Color', 'black', 'LineWidth', 2);
    line([xcf xcf], [0 40], 'Color', 'black', 'LineWidth', 2);
    colorbar;
    %set(gca, 'CLim', [0 2])
    xlabel('Time [s]')
    ylabel('Channels');
    xticks(0:sampleRate:minDur);
    xticklabels((0:sampleRate:minDur)/sampleRate);
    yticks(1:numel(channels_label(find(~strcmp(channels_label, '')))));
    yticklabels(channels_label(find(~strcmp(channels_label, ''))));
    
    if ~normalization
        title('Fisher score -- log band power');
    else
        title('Fisher score -- ERD/ERS');
    end

    %% PSD
    %{
    psd_wlength = 0.5;
    psd_wshift = 0.0625;
    psd_pshift = 0.25;
    psd_mlength =  1;
    band = [8 14];

    signal = signal * lap;
    [features, f] = proc_spectrogram(signal, psd_wlength, psd_wshift, psd_pshift, sampleRate, psd_mlength);
    features = log(features);
    header.EVENT.POS = proc_pos2win(header.EVENT.POS, psd_wshift*sampleRate, 'backward', psd_mlength*sampleRate);
    header.EVENT.DUR = floor(header.EVENT.DUR/(psd_wshift*sampleRate)) + 1;

    features = squeeze(mean(features(:, find(f >= band(1) & f <= band(2)),:),2));

    fixPOS = header.EVENT.POS(header.EVENT.TYP == 786);
    cueTYP = header.EVENT.TYP(header.EVENT.TYP == 730 | header.EVENT.TYP == 731);
    cuePOS = header.EVENT.POS(header.EVENT.TYP == 730 | header.EVENT.TYP == 731);
    cueDUR = header.EVENT.DUR(header.EVENT.TYP == 730 | header.EVENT.TYP == 731);
    cfPOS  = header.EVENT.POS(header.EVENT.TYP == 781);
    cfDUR  = header.EVENT.DUR(header.EVENT.TYP == 781);
    startTrial = nan(n_trial, 1);
    stopTrial  = nan(n_trial, 1);
    tk = cueTYP;
    for idx_trial = 1:n_trial
        startTrial(idx_trial) = fixPOS(idx_trial);
        stopTrial(idx_trial) = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;
    end

    minDur = min(stopTrial - startTrial);
    dataTrial = nan(minDur, nchannels, n_trial);

    for idx_trial = 1:n_trial
        c_start = startTrial(idx_trial);
        dataTrial(:, :, idx_trial) = features(c_start:c_start+minDur-1,:);
    end

    psd_1 = mean(dataTrial(:,:,tk == classes(1)), 3);
    psd_2 = mean(dataTrial(:,:,tk == classes(2)), 3);

    psd = (psd_2 - psd_1);
    psd = psd';

    xcue = mean(cuePOS - fixPOS) +1;
    xcf = mean(cfPOS - fixPOS) +1;

    figure();
    imagesc(psd);
    line([xcue xcue], [1 39], 'Color', 'black', 'LineWidth', 2);
    line([xcf xcf], [1 39], 'Color', 'black', 'LineWidth', 2);
    colorbar;
    set(gca, 'CLim', [-0.6 0.6])
    xlabel('Time [s]')
    xticks(1:size(psd,2));
    xticklabels(0:1/16:size(psd,2)/16);
    ylabel('Channels');
    yticks(1:numel(channels_label));
    yticklabels(channels_label);
    title('Feature Map: psd');
    
    %}
end