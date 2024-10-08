%% Test the features extracted, the qda probs and the integrated probs only during CF
% works with only ne file
% clc; clear all; close all;

% loading data, gdf and classifier
subject = 'h8';
day = '/20240926'; % place / before the data
path_gdf = ['/home/paolo/cvsa_ws/record/' subject day '/gdf/calibration'];
gdf_files = dir(fullfile(path_gdf, '*.gdf'));

chanlog_path = '/home/paolo/new_chanlocs64.mat';

% variables for parameters given in the qda file + general parameters
channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
        '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
filterOrder = 4;
sampleRate = 512;
idx_selchs = find(~cellfun('isempty', channels_label));
% bands = [[6 8]; [8 10]; [10 12]; [12 14]; [14 16]; [8 14]];
bands = [8 14];
frameSize = 32;
nchannels = 39;
bufferSize_data = 512;
classes = [730 731];

% processing of gdf file
for idx_file=1:length(gdf_files)
    % load the gdf
    disp(['Loaded file: ' gdf_files(idx_file).name]);
    c_file = [path_gdf '/' gdf_files(idx_file).name];
    [s,header] = sload(c_file);
    s = s(:,1:nchannels);
    cueTYP_gdf = header.EVENT.TYP(ismember(header.EVENT.TYP, classes));
    fixPOS_gdf = header.EVENT.POS(header.EVENT.TYP == 786);
    cfPOS_gdf = header.EVENT.POS(ismember(header.EVENT.TYP, 781));
    cfDUR_gdf = header.EVENT.DUR(ismember(header.EVENT.TYP, 781));
    ntrials = size(fixPOS_gdf, 1);
    if size(cueTYP_gdf,1) ~= ntrials
        cueTYP_gdf = cueTYP_gdf(size(cueTYP_gdf, 1) - size(fixPOS_gdf,1) + 1:end);
    end

    % define for the filtering 
    zi_high = cell(1, size(bands, 1));
    zi_low  = cell(1, size(bands, 1));
    for i = 1:size(bands, 1)
        zi_high{i} = [];
        zi_low{i}  = [];
    end

    % vairables for keeping the features and the pointers
    all_features = cell(size(bands, 1), 1);
    trialStart_features = []; 
    trialEnd_features = [];
    in_trial = false;
    idx_trial = 1;

    % define variable for simulate the buffer used in rosneuro
    buffer = nan(bufferSize_data, nchannels);
    
    % apply processing on all the data
    disp('   [info] applying same processing used in ros')
    nchunk = floor(size(s,1)/frameSize);
    for idx_c=1:nchunk
        frame = s((idx_c-1)*frameSize+1:idx_c*frameSize,:);
        buffer(1:end-frameSize,:) = buffer(frameSize+1:end,:);
        buffer(end-frameSize+1:end, :) = frame;

        % check if the buffer is full
        if any(isnan(buffer))
            continue;
        end

        % iterate over bands
        for idx_band = 1:size(bands, 1)
            c_band = bands(idx_band,:);
            [c_b_low, c_a_low] = butter(filterOrder, c_band(2)*(2/sampleRate),'low');
            [c_b_high, c_a_high] = butter(filterOrder, c_band(1)*(2/sampleRate),'high');
            c_zi_low = zi_low{idx_band};
            c_zi_high = zi_high{idx_band};

            % apply low and high pass filters
            [s_low, c_zi_low] = filter(c_b_low,c_a_low,buffer,c_zi_low);
            [tmp_data,c_zi_high] = filter(c_b_high,c_a_high,s_low,c_zi_high);
            zi_low{idx_band} = c_zi_low;
            zi_high{idx_band} = c_zi_high;

            % apply pow
            tmp_data = power(tmp_data, 2);

            % apply average
            tmp_data = mean(tmp_data, 1);

            % apply log
            tmp_data = log(tmp_data);

            idx_feature = idx_band; % we repeat the bands to have same index for band and features
            features = tmp_data;
            all_features{idx_band} = cat(1, all_features{idx_band}, features);
        end

        if ~in_trial && idx_trial <= ntrials
            if (idx_c-1)*frameSize+1 >= fixPOS_gdf(idx_trial) - frameSize
                in_trial = true;
                trialStart_features = cat(1, trialStart_features, size(all_features{1}, 1)+1); % the -1 is to anlign with rosneuro
            end
        else
            if idx_trial <= ntrials
                if idx_c*frameSize > cfPOS_gdf(idx_trial) + cfDUR_gdf(idx_trial)
                    in_trial = false;
                    trialEnd_features = cat(1, trialEnd_features, size(all_features{1}, 1)-1);
                    idx_trial = idx_trial + 1;
                end
            end
        end
    end

    % create the ERD structure
    minDurTrial = min(trialEnd_features - trialStart_features);
    for idx_band = 1:size(bands, 1)
        data = all_features{idx_band};
        band = bands(idx_band,:);
        ERD = nan(minDurTrial, nchannels, ntrials);
        for idx_trial = 1:ntrials
            ERD(:,:,idx_trial) = data(trialStart_features(idx_trial):trialStart_features(idx_trial)+minDurTrial-1,:);
        end

        sampleRate_ros = 16;
        showTopoplot_ERDERS(chanlog_path, ERD, sampleRate_ros, 2, channels_label, band, minDurTrial, cueTYP_gdf, classes, nchannels, ['ros ' gdf_files(idx_file).name])
    end
end




