%% check for the ERP in freq
clear all; % close all;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/EOG')

%% Initialization
signals = [];
headers.TYP = [];
    headers.DUR = [];
    headers.POS = [];
classes = [730 731];
cf_event = 781;
fix_event = 786;
nchannels = 39;
nclasses = length(classes);
filterOrder = 4;
avg = 1;
eog_threshold = 500;
roi = {{'F3', 'F1', 'Fz', 'F2', 'F4'}; {'FC3', 'FC1', 'FCz', 'FC2', 'FC4'}; {'C3', 'C1', 'Cz', 'C2', 'C4'}; {'CP3', 'CP1', 'CPz', 'CP2', 'CP4'};
    {'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6'}; {'PO7', 'PO5', 'PO3', 'POz', 'PO4', 'PO6', 'PO8'}; {'O1', 'Oz', 'O2'}};
roi_label = {{'F channels'}, {'FC channels'}, {'C channels'}, {'CP channels'}, {'P channels'}, {'PO channels'}, {'O channels'}};
roi_indices = cell(size(roi));
nroi = length(roi);

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);

%% concatenate the files
nFiles = length(filenames);
trial_with_eog = [];
for idx_file= 1: nFiles
    fullpath_file_shift = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_shift);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;
    for i = 1:nroi
        roi_indices{i} = find(ismember(channels_label, roi{i}));
    end

    c_trial_with_eog = eog_detection(c_signal, header, eog_threshold, {'FP1', 'FP2', 'EOG'});
    trial_with_eog = [trial_with_eog; c_trial_with_eog];

    signal_roi = nan(size(c_signal, 1), nroi); % for roi
    for i = 1:nroi
        signal_roi(:,i) = mean(c_signal(:, roi_indices{i}), 2);
    end

    disp('      [proc] applying filtering')
    [b, a] = butter(filterOrder, 40*(2/header.SampleRate),'low');
    s_low = filter(b,a,signal_roi);
    [b, a] = butter(filterOrder, 0.1*(2/header.SampleRate),'high');
    s_filt = filter(b,a,s_low);

    signal_roi_processed = s_filt;

    headers.sampleRate = header.SampleRate;
    headers.channels_labels = header.Label;
    if isempty(find(header.EVENT.TYP == 2, 1)) % no eye calibration
        headers.TYP = cat(1, headers.TYP, header.EVENT.TYP);
        headers.DUR = cat(1, headers.DUR, header.EVENT.DUR);
        headers.POS = cat(1, headers.POS, header.EVENT.POS + size(signals, 1));
    else
        k = find(header.EVENT.TYP == 1, 1);
        headers.TYP = cat(1, headers.TYP, header.EVENT.TYP(k:end));
        headers.DUR = cat(1, headers.DUR, header.EVENT.DUR(k:end));
        headers.POS = cat(1, headers.POS, header.EVENT.POS(k:end) + size(signals, 1));
    end

    signals = cat(1, signals, signal_roi_processed(:,:));
 
end

%% Labelling data 
sampleRate = headers.sampleRate;
cuePOS = headers.POS(ismember(headers.TYP, classes));
cueDUR = headers.DUR(ismember(headers.TYP, classes));
cueTYP = headers.TYP(ismember(headers.TYP, classes));

fixPOS = headers.POS(headers.TYP == 786);
fixDUR = headers.DUR(headers.TYP == 786);

cfPOS = headers.POS(headers.TYP == 781);
cfDUR = headers.DUR(headers.TYP == 781);

minDurCue = min(cueDUR);
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
trial_data = nan(min_trial_data, nroi, ntrial);
for trial = 1:ntrial
    c_start = trial_start(trial);
    c_end = trial_start(trial) + min_trial_data - 1;
    trial_data(:,:,trial) = signals(c_start:c_end,:);
end


%% refactoring the data
% now the data are placed alternatevely, so odd trial id is for 730, even for 731
balanced_trial_idx = balanced_data_afterEOG(trial_with_eog, trial_typ, classes);
balanced_trial_data = trial_data(:,:,logical(balanced_trial_idx));
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
        tmp(:,:,idx_trial_class + idx_class - 1) = balanced_trial_data(:,:,idx_classes_trial(i, idx_class));
        trial_typ(idx_trial_class + idx_class - 1) = classes(idx_class);
    end
    i = i + 1;
end
trial_data = tmp;

% trial_data(:,:,:,3) = [];

%% plot the frequency analysis for all the selected channels
for idx_roi = 1:nroi
    params.epochTime = (1:min_trial_data) / 512;
    params.fsamp = 512;
    eegEpochs_pre.data = squeeze(trial_data(:,:,:));
    eegEpochs_pre.label = ones(size(trial_data,3),1);
    compute_theta_peak_v3_1Dcursor(eegEpochs_pre, params, false, true, true, idx_roi, roi_label{idx_roi}{1});
end

function compute_theta_peak_v3_1Dcursor(eegEpochs_pre, params, baseline_subtract, baseline_normalise, error_correct_subtraction, chan, chan_name)

%close all

signal_range = dsearchn(params.epochTime', 0.0):dsearchn(params.epochTime', 1.0);

wavtime = -2:1/params.fsamp:2;
%cut the first half and the last half of the convolution (N + M - 1)
half_wave = (length(wavtime)-1)/2;
wavelet_num_cycles_list = [6];

plotting = true;

min_freq =  1;
max_freq = 30;
%increase frequency resolution
freq_res = 0.1;
frex = min_freq:freq_res:max_freq;
amplitude_threshold = 10000;
%Same baseline as Reinhart PNAS
baseline_range = [dsearchn(params.epochTime', -0.3):dsearchn(params.epochTime', -0.1)];

for wav_id = 1:length(wavelet_num_cycles_list)
    keep_index = not(squeeze(any(any(abs(eegEpochs_pre.data(signal_range, :, :)) > amplitude_threshold))));
    logical_index_negative = eegEpochs_pre.label == 1 & keep_index ;
    temp_negative = eegEpochs_pre.data(:, :, logical_index_negative);
    logical_index_positive = eegEpochs_pre.label == 2 & keep_index ;
    temp_positive = eegEpochs_pre.data(:, :, logical_index_positive);
    logical_index_neutral = eegEpochs_pre.label == 3 & keep_index ;
    temp_neutral = eegEpochs_pre.data(:, :, logical_index_neutral);

    %FFT parameters (block of trials -> computational efficiency)
    nWave = length(wavtime);
    nData = size(temp_negative, 1) * size(temp_negative, 3); %number of samples x number of trials
    nConv = nWave + nData - 1;
    
    % now compute the FFT of all trials concatenated
    %reshape to one vector
    alldata = reshape( temp_negative(:, chan, :) ,1,[]);
    dataX   = fft( alldata ,nConv );
    
    % initialize output time-frequency data
    tf_negative = zeros(length(frex),size(temp_negative,1));
    % loop over frequencies
    for fi=1:length(frex)
        wavelet_freq = frex(fi);
        s = wavelet_num_cycles_list(wav_id)/(2*pi*wavelet_freq);
        %morlet wavelet
        wavelet  = exp(2*1i*pi*wavelet_freq.*wavtime) .* exp(-wavtime.^2./(2*s^2));
        %mexican hat wavelet
        
        waveletX = fft(wavelet,nConv);
        %Normalise against amplitude
        waveletX = waveletX ./ max(waveletX);
        
        % now run convolution in one step
        as = ifft(waveletX .* dataX);
        % only select the part relevant to the signal
        as = as(half_wave+1:end-half_wave);
        
        % and reshape back to time X trials
        as = reshape( as, size(temp_negative, 1), size(temp_negative, 3) );
        
        % compute power and average over trials
        tf_negative(fi,:) = mean( abs(as).^2 ,2);
        
    end
    
    nWave = length(wavtime);
    nData = size(temp_positive, 1) * size(temp_positive, 3); %number of samples x number of trials
    nConv = nWave + nData - 1;
    
    % now compute the FFT of all trials concatenated
    %reshape to one vector
    alldata = reshape(temp_positive(:, chan, :) ,1,[]);
    dataX   = fft( alldata ,nConv );
    
    % initialize output time-frequency data
    tf_positive = zeros(length(frex),size(temp_positive,1));
    % loop over frequencies
    for fi=1:length(frex)
        wavelet_freq = frex(fi);
        s = wavelet_num_cycles_list(wav_id)/(2*pi*wavelet_freq);
        %morlet wavelet
        wavelet  = exp(2*1i*pi*wavelet_freq.*wavtime) .* exp(-wavtime.^2./(2*s^2));
        %mexican hat wavelet
        
        waveletX = fft(wavelet,nConv);
        %Normalise against amplitude
        waveletX = waveletX ./ max(waveletX);
        
        % now run convolution in one step
        as = ifft(waveletX .* dataX);
        % only select the part relevant to the signal
        as = as(half_wave+1:end-half_wave);
        
        % and reshape back to time X trials
        as = reshape( as, size(temp_positive, 1), size(temp_positive, 3) );
        
        % compute power and average over trials
        tf_positive(fi,:) = mean(abs(as).^2 ,2);
        
    end

    nWave = length(wavtime);
    nData = size(temp_neutral, 1) * size(temp_neutral, 3); %number of samples x number of trials
    nConv = nWave + nData - 1;
    
    % now compute the FFT of all trials concatenated
    %reshape to one vector
    alldata = reshape(temp_neutral(:, chan, :) ,1,[]);
    dataX   = fft( alldata ,nConv );
    
    % initialize output time-frequency data
    tf_neutral = zeros(length(frex),size(temp_neutral,1));
    % loop over frequencies
    for fi=1:length(frex)
        wavelet_freq = frex(fi);
        s = wavelet_num_cycles_list(wav_id)/(2*pi*wavelet_freq);
        %morlet wavelet
        wavelet  = exp(2*1i*pi*wavelet_freq.*wavtime) .* exp(-wavtime.^2./(2*s^2));
        %mexican hat wavelet
        
        waveletX = fft(wavelet,nConv);
        %Normalise against amplitude
        waveletX = waveletX ./ max(waveletX);
        
        % now run convolution in one step
        as = ifft(waveletX .* dataX);
        % only select the part relevant to the signal
        as = as(half_wave+1:end-half_wave);
        
        % and reshape back to time X trials
        as = reshape( as, size(temp_neutral, 1), size(temp_neutral, 3) );
        
        % compute power and average over trials
        tf_neutral(fi,:) = mean(abs(as).^2 ,2);
        
    end
    
    %Apply to trial average power
    % db conversion
    tf_negative_db_conversion = zeros(size(tf_negative,1), size(tf_negative,2));
    tf_negative_percentage_change = zeros(size(tf_negative,1), size(tf_negative,2));
    if baseline_normalise
        for freq_id = 1:size(tf_negative,1)
            tf_negative_db_conversion(freq_id,:) = 10*log10(tf_negative(freq_id,:) ./ mean(tf_negative(freq_id,baseline_range),2));
            tf_negative_percentage_change(freq_id,:) = (tf_negative(freq_id,:) - mean(tf_negative(freq_id,baseline_range),2)) ./ mean(tf_negative(freq_id,baseline_range),2);
        end
    end

    
    tf_positive_db_conversion = zeros(size(tf_positive,1), size(tf_positive,2));
    tf_positive_percentage_change = zeros(size(tf_positive,1), size(tf_positive,2));
    if baseline_normalise
        for freq_id = 1:size(tf_positive,1)
            tf_positive_db_conversion(freq_id,:) = 10*log10(tf_positive(freq_id,:) ./ mean(tf_positive(freq_id,baseline_range),2));
            tf_positive_percentage_change(freq_id,:) = (tf_positive(freq_id,:) - mean(tf_positive(freq_id,baseline_range),2)) ./ mean(tf_positive(freq_id,baseline_range),2);
        end
    end

    tf_neutral_db_conversion = zeros(size(tf_neutral,1), size(tf_neutral,2));
    tf_neutral_percentage_change = zeros(size(tf_neutral,1), size(tf_neutral,2));
    if baseline_normalise
        for freq_id = 1:size(tf_neutral,1)
            tf_neutral_db_conversion(freq_id,:) = 10*log10(tf_neutral(freq_id,:) ./ mean(tf_neutral(freq_id,baseline_range),2));
            tf_neutral_percentage_change(freq_id,:) = (tf_neutral(freq_id,:) - mean(tf_neutral(freq_id,baseline_range),2)) ./ mean(tf_neutral(freq_id,baseline_range),2);
        end
    end


    %Compute the difference
    
    if (error_correct_subtraction)
        tf_difference_neg_db_conversion = tf_negative_db_conversion - tf_neutral_db_conversion;
        tf_difference_pos_db_conversion = tf_positive_db_conversion - tf_neutral_db_conversion;
    end
  

    if plotting == true

        % Not plotting if NaN values
        isAllNaN = all(isnan(tf_negative_db_conversion), 'all');
        if ~isAllNaN
            figure; clf;
            contourf(params.epochTime(80:end-80),frex,tf_negative_db_conversion(:,(80:end-80)))
            title(['tf of negative valence trials (db conversion): ' num2str(wavelet_num_cycles_list(wav_id)) ' cycles | ' chan_name]);
            %contourf(params.epochTime,frex,tf .* 100,40,'linecolor','none')
            %set(gca,'clim',[0 400],'ydir','normal')
            xlabel('Time [s]', 'FontSize',10);
            ylabel('Frequency [Hz]', 'FontSize',10);
            colorbar
        end
        
        isAllNaN = all(isnan(tf_positive_db_conversion), 'all');
        if ~isAllNaN
            figure; clf;
            contourf(params.epochTime(80:end-80),frex,tf_positive_db_conversion(:,(80:end-80)))
            title(['tf of positive valence trials (db conversion): ' num2str(wavelet_num_cycles_list(wav_id)) ' cycles | ' chan_name]);
            %contourf(params.epochTime,frex,tf .* 100,40,'linecolor','none')
            %set(gca,'clim',[0 400],'ydir','normal')
            xlabel('Time [s]', 'FontSize',10);
            ylabel('Frequency [Hz]', 'FontSize',10);
            colorbar
        end

        isAllNaN = all(isnan(tf_neutral_db_conversion), 'all');
        if ~isAllNaN
            figure; clf;
            contourf(params.epochTime(80:end-80),frex,tf_neutral_db_conversion(:,(80:end-80)))
            title(['tf of neutral trials (db conversion): ' num2str(wavelet_num_cycles_list(wav_id)) ' cycles | ' chan_name]);
            %contourf(params.epochTime,frex,tf .* 100,40,'linecolor','none')
            %set(gca,'clim',[0 400],'ydir','normal')
            xlabel('Time [s]', 'FontSize',10);
            ylabel('Frequency [Hz]', 'FontSize',10);
            colorbar
        end

        isAllNaN = all(isnan(tf_difference_neg_db_conversion), 'all');
        if ~isAllNaN
            figure; clf;
            contourf(params.epochTime(80:end-80),frex,tf_difference_neg_db_conversion(:,(80:end-80)))
            title(['tf of negative - neutral trials (db conversion): ' num2str(wavelet_num_cycles_list(wav_id)) ' cycles | ' chan_name]);
            %contourf(params.epochTime,frex,tf .* 100,40,'linecolor','none')
            %set(gca,'clim',[0 400],'ydir','normal')
            xlabel('Time [s]', 'FontSize',10);
            ylabel('Frequency [Hz]', 'FontSize',10);
            colorbar
        end

        isAllNaN = all(isnan(tf_difference_pos_db_conversion), 'all');
        if ~isAllNaN
            figure; clf;
            contourf(params.epochTime(80:end-80),frex,tf_difference_pos_db_conversion(:,(80:end-80)))
            title(['tf of positive - neutral trials (db conversion): ' num2str(wavelet_num_cycles_list(wav_id)) ' cycles | ' chan_name]);
            %contourf(params.epochTime,frex,tf .* 100,40,'linecolor','none')
            %set(gca,'clim',[0 400],'ydir','normal')
            xlabel('Time [s]', 'FontSize',10);
            ylabel('Frequency [Hz]', 'FontSize',10);
            colorbar
        end
        
        % figure; clf;
        % contourf(params.epochTime(80:end-80),frex,tf_negative_percentage_change(:,(80:end-80)))
        % title(['tf of negative trials (percentage change): ' num2str(wavelet_num_cycles_list(wav_id)) 'cycles']);
        % %contourf(params.epochTime,frex,tf .* 100,40,'linecolor','none')
        % %set(gca,'clim',[0 400],'ydir','normal')
        % colorbar
        % 
        % figure; clf;
        % contourf(params.epochTime(80:end-80),frex,tf_positive_percentage_change(:,(80:end-80)))
        % title(['tf of positive trials (percentage change): ' num2str(wavelet_num_cycles_list(wav_id)) 'cycles']);
        % %contourf(params.epochTime,frex,tf .* 100,40,'linecolor','none')
        % %set(gca,'clim',[0 400],'ydir','normal')
        % colorbar
        % 
        % figure; clf;
        % contourf(params.epochTime(80:end-80),frex,tf_neutral_percentage_change(:,(80:end-80)))
        % title(['tf of neutral trials (percentage change): ' num2str(wavelet_num_cycles_list(wav_id)) 'cycles']);
        % %contourf(params.epochTime,frex,tf .* 100,40,'linecolor','none')
        % %set(gca,'clim',[0 400],'ydir','normal')
        % colorbar

    end
    
end
end