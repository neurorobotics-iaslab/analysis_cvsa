clc; clear all; close all;
% create the dataset log band power, selected channels and selected freq

%% informations
c_subject = 'c7';   %%%%%%%%%%%%% subject -------------
train_percentage = 0.75;
classes = [730 731];

filterOrder = 4;
% select the channels and bands according to the name of the subject
if c_subject == 'c7'
%     bands = {[15 18],  [18 21]};
%     selchs = {{'P5', 'PO5', 'PO7'}, {'P3', 'P5'}};
    bands = {[10 12],  [16 18]};
    selchs = {{'PO5', 'PO7'}, {'PO4', 'PO6'}};
elseif c_subject == 'g2'
    bands = {[6 9], [15 18], [18 21]};
    selchs = {{'PO8', 'P6'},  {'PO7'}, {'PO7', 'OZ'}};
elseif c_subject == 'h7'
    selchs = {{'P4', 'PO4', 'O2', 'PO6', 'PO8'}}; 
    bands = {[8 14]};
elseif c_subject == 'd61111'
    selchs = {{'P3', 'P1', 'P5'}};
    bands = {[8 14], [5, 12]};
elseif c_subject == 'c7_new'
    selchs = {{'P2', 'P4'},{'PO5', 'PO7'}};
    bands = {[12 14], [14 16]};
elseif c_subject == 'c7_vis'
    selchs = {{'P5', 'PO5', 'PO7'},{'P6', 'PO6', 'PO8'}};
    bands = {[10 12], [14 16]};
elseif c_subject == 'c7_aud'
    selchs = {{'P4', 'P5'},{'P5', 'P6'}};
    bands = {[8 10], [10 12]};
elseif c_subject == 'c72vis'
    selchs = {{'PO5', 'PO7'}, {'PO5'}};
    bands = {[10 12], [12 14]};
else
    disp('No such subject')
end

% not modification needed for these informations
sampleRate = 512;
sfile = ['/home/paolo/cvsa_ws/record/' c_subject '/dataset/logband_f_cf_selectedband1.mat'];

% path = ['/home/paolo/cvsa_ws/record/' c_subject '/mat_selectedTrials'];
% files = dir(fullfile(path, '*.mat'));

path = ['/home/paolo/cvsa_ws/record/' c_subject '/gdf'];
files = dir(fullfile(path, '*.gdf'));

channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
        'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};


%% initialization variable to save
X = [];
y = [];


% to have simple passage for qda
final_bands = [];
for i=1:size(selchs,2)
    c_selchs = selchs{i};
    c_band = bands{i};
    for j=1:size(c_selchs, 2)
        final_bands = cat(1, final_bands, c_band);
    end
end
info.classes = classes;
info.sampleRate = sampleRate;
info.selchs = selchs;
info.idx_selchs = [];
info.files = {};
info.band = final_bands;
info.trialStart = [];
info.trialDUR = [];
info.filterOrder = 4;
info.startNewFile = [];


%% take only interested data
for idx_f = 1:length(files)
    file = fullfile(path, files(idx_f).name);
    disp(['file (' num2str(idx_f) '/' num2str(length(files))  '): ', file])
    %file = '/home/paolo/prova32ch.gdf';
%     load(file);
    [signal,header] = sload(file);
    info.files = cat(1, info.files, files(idx_f).name);
    nchannels = length(channels_label);

    %% labeling
    disp('   Labelling')
    signal = signal(:,1:nchannels);
    events = header.EVENT;
    cuePOS = events.POS(events.TYP == 730 | events.TYP == 731);
    cueDUR = events.DUR(events.TYP == 730 | events.TYP == 731);
    cueTYP = events.TYP(events.TYP == 730 | events.TYP == 731);
    cfPOS  = events.POS(events.TYP == 781);
    cfDUR  = events.DUR(events.TYP == 781);
    nTrials = length(cueTYP);
    if(contains(file, 'calibration')) %% for gdf not mat
        cuePOS = cuePOS(3:end) - 1;
        cueTYP = cueTYP(3:end);
        cueDUR = cueDUR(3:end);
        cfPOS  = events.POS(events.TYP == 781)-1;
        nTrials = length(cueTYP);
    end
%     if(contains(file, 'evaluation')) %% for gdf not mat if in the evaluation the eye calibration is done
%         cuePOS = cuePOS(3:end) - 1;
%         cueTYP = cueTYP(3:end);
%         cueDUR = cueDUR(3:end);
%         cfPOS  = events.POS(events.TYP == 781)-1;
%         nTrials = length(cueTYP);
%     end

    %% Initialization variables
    disp('   Initialization variables')
    frameSize = 32;
    bufferSize = 512; 
    X_band = [];
%     s_band = [];
%     s_pow = [];
%     s_avg = [];
%     s_log = [];
    
    if idx_f==1
        prev_file = 0;
    end
    

    for idx_band = 1:length(bands)
        c_band = bands{idx_band};
        disp(['   band: ' num2str(c_band(1)) '-' num2str(c_band(2))]);
        [c_b_low, c_a_low] = butter(filterOrder, c_band(2)*(2/sampleRate),'low');
        [c_b_high, c_a_high] = butter(filterOrder, c_band(1)*(2/sampleRate),'high');
        zi_low = [];
        zi_high = [];
        X_temp = [];
        y_temp = [];

        %% Iterate over trials
        for i=1:nTrials
            disp(['      trial ' num2str(i) '/' num2str(nTrials)])
            % initialization variables
            buffer = nan(bufferSize, nchannels);
            start_trial = cuePOS(i);
            %end_trial = cfPOS(i) + cfDUR(i) - 1;
            end_trial = cfPOS(i) + 256 - 1; % take only first 0.5s
            % division for frameSize
            end_trial = int64(ceil(single(end_trial-start_trial)/32)*32) + start_trial;
            data = signal(start_trial:end_trial,:);

            % application of the buffer
            if idx_band == 1
                info.trialStart = cat(1, info.trialStart, size(X_temp,1)+prev_file);
            end
            nchunks = (end_trial-start_trial) / 32;
            for j = 1:nchunks
                frame = data((j-1)*frameSize+1:j*frameSize,:);
                buffer(1:end-frameSize,:) = buffer(frameSize+1:end,:);
                buffer(end-frameSize+1:end, :) = frame;

                % check
                if any(isnan(buffer))
                    continue;
                end

                % apply low and high pass filters
                [s_low, zi_low] = filter(c_b_low,c_a_low,buffer,zi_low);
                [tmp_data,zi_high] = filter(c_b_high,c_a_high,s_low,zi_high);
%                 s_band = cat(1, s_band, tmp_data);

                % apply pow
                tmp_data = power(tmp_data, 2);
%                 s_pow = cat(1, s_pow, tmp_data);

                % apply average
                tmp_data = mean(tmp_data, 1);
%                 s_avg = cat(1, s_avg, tmp_data);

                % apply log
                tmp_data = log(tmp_data);
%                 s_log = cat(1, s_log, tmp_data);

                % save in the dataset
                X_temp = cat(1, X_temp, tmp_data);
                y_temp = cat(1, y_temp, repmat(cueTYP(i), size(tmp_data,1), 1));

            end

            % save the dur of the trial only the first time a band is done
            if idx_band == 1
                info.trialDUR = cat(1, info.trialDUR, size(X_temp,1)+prev_file-info.trialStart(end));
            end
        end

        %% take only interested values
        disp('      Take only interested channels for that band')
        selch = selchs{idx_band};
        idx_interest_ch = zeros(1, numel(selch));
        for k=1:numel(selch)
            idx_interest_ch(k) = find(strcmp(channels_label, selch{k}));
        end

        % as before
        if idx_band == 1
            info.startNewFile = cat(1, info.startNewFile, size(X,1));
            y = cat(1, y, y_temp);
        end

        X_band = cat(2, X_band, X_temp(:,idx_interest_ch));

        if idx_f == 1
            info.idx_selchs = cat(2, info.idx_selchs, idx_interest_ch);
        end
    
    end
    X = cat(1, X, X_band);
    prev_file = size(X,1);
end
info.startTest = info.trialStart(floor(train_percentage * size(info.trialStart,1)));

%% save the values
%X = X(:, idx_interest_ch);
save(sfile, 'X', 'y', 'info');