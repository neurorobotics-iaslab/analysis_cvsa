%% We compute the ERD for 8-14 and 1-30 to see the ERD in order to do the labelling of the data
%  we allow to take the calibration files but also the evaluation
%  we show each trial and ask to the user for the start and stop in time of the ERD

clc; clear all; close all;
%% Show the ERD/ERS for a file and for the channels selected
% compute and show the ERD/ERS as a single image for a selected channel

subject = 'c7';
day = '';

bands = [[8 14]; [1 30]]; 
plot_nrow = 1;
plot_ncol = size(bands, 1) / plot_nrow;

cal_eval = input('calibration (1) or evaluation (2): ');
if cal_eval == 1
    a = 'calibration';
elseif cal_eval == 2
    a = 'evaluation';
else
    disp('Error on the input, only 1 or 2 are allowd');
    return;
end
path = ['/home/paolo/cvsa_ws/record/' subject day '/gdf/' a];
files = dir(fullfile(path, '*.gdf'));

chanlog_path = '/home/paolo/new_chanlocs64.mat';

channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
        '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

% channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

eog_channels =  {'FP1', 'FP2', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'EOG', ...
        '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''};

classes = [730 731];
sampleRate = 512;
normalization = false; % if false compute the log power, if true the ERD/ERS
avg = 1;
filtOrder = 4;
th_eog = 3.0e4; %2.5e4;
channelSelected = find(~cellfun(@isempty, channels_label));

% check all files
for idx_f = 1:length(files)
    file = fullfile(path, files(idx_f).name);
    disp(['loading file: ' file])
    [signal, header] = sload(file);
    nchannels = length(channels_label);

    % Extract infoo data
    signal = signal(:, 1:nchannels);
    [fixPOS, fixDUR, cuePOS, cueDUR, cueTYP, cfPOS, cfDUR, n_trial] = extract_info_data(header, 786, [730 731], 781);

    % Extract trial infoo
    [trialStart, trialStop, fixStart, fixStop, ck, tk] = extract_trial_infoo(signal, fixPOS, fixDUR, cfPOS, cfDUR, cueTYP, n_trial);

    % EOG check
    prev_n_trial = n_trial;
    for i=1:n_trial
        c_start = cfPOS(i);
        c_end = cfPOS(i) + cfDUR(i) - 1;
        data_sum = signal(c_start:c_end,:);
        if eye_movement_check(data_sum, eog_channels, th_eog, sampleRate)
            trialStart(i) = -1;
            trialStop(i) = -1;
            fixStop(i) = -1;
            fixStart(i) = -1;
            cueTYP(i) = -1;
        end
    end
    trialStart(trialStart == - 1) = [];
    trialStop(trialStop == - 1) = [];
    fixStart(fixStart == - 1) = [];
    fixStop(fixStop == -1) = [];
    cueTYP(cueTYP == -1) = [];
    n_trial = size(trialStart, 1);
    disp(['   [INFO] from ' num2str(prev_n_trial) ' to ' num2str(n_trial)]);
    if n_trial == 0
        disp('ERROR no trials');
        return;
    end

    % variable for concatenation
    current_file.signal_processed = cell(1, size(bands, 1));
    current_file.ERD = cell(size(bands, 1), n_trial);
    current_file.minDur = Inf;

    % iterate over bands in order to apply the processing only once
    for idx_band=1:size(bands,1)
        band = bands(idx_band,:);
        disp(['   [INFO] band used: [' num2str(band(1)) ',' num2str(band(2)) ']']);
        % processing
        signal_processed = processing_offline(signal, nchannels, sampleRate, band, filtOrder, avg);
        
        current_file.signal_processed{idx_band} = signal_processed; % nsample x channels
    end

    % iterate over the trials --> to compute erd for single trial
    for idx_trial=1:n_trial

        figure;
        % iterate over the bands
        for idx_band = 1:size(bands, 1)
            band = bands(idx_band, :);

            c_signal_processed = current_file.signal_processed{idx_band};
            c_signal_trialSelected = c_signal_processed(trialStart(idx_trial):trialStop(idx_trial)-1,:);

            if normalization
                c_signal_fixation = c_signal_processed(fixStart(idx_trial):fixStop(idx_trial)-1,:);
                baseline = repmat(mean(c_signal_fixation), [size(c_signal_trialSelected, 1) 1 1]);
                ERD = log(c_signal_trialSelected ./ baseline);
                
            else
                ERD = log(c_signal_trialSelected);
            end
            ERD = ERD(:,channelSelected);

            % concatenation
            if current_file.minDur <= size(ERD, 1)
                current_file.ERD{idx_band, idx_trial} = ERD(1:current_file.minDur,:);
            else
                current_file.minDur = size(ERD, 1);
                current_file.ERD{idx_band, idx_trial} = ERD(1:current_file.minDur,:);
                % tutte le band prima
                for i=1:idx_band-1
                    a = current_file.ERD{i, idx_trial};
                    current_file.ERD{i, idx_trial} = a(1:current_file.minDur,:);
                end
                % tutti i trial prima
                for j = 1:idx_trial - 1
                    for i = 1: size(bands, 1)
                        a = current_file.ERD{i, j};
                        current_file.ERD{i, j} = a(1:current_file.minDur,:);
                    end
                end
            end

            % plot the erd with the name of the class asked
            band_strings = strings(size(bands, 1), 1);
            for i = 1:size(bands, 1)
                band_strings(i) = sprintf('%d-%d', bands(i, 1), bands(i, 2));
            end
            t = linspace(0, size(ERD, 1)/sampleRate, size(ERD, 1));
            subplot(plot_nrow, plot_ncol, idx_band);
            cdata = ERD;
            imagesc(t, 1:size(channelSelected, 1), cdata');
            set(gca,'YDir','normal');
            colormap(hot);
            colorbar;
            title(['Band ' band_strings(idx_band)]);
            xlabel('Time [s]');
            ylabel('Channels');
            yticks(1:size(channelSelected, 2));
            yticklabels(channels_label(channelSelected));
            ylim([0.5, size(channelSelected, 2) + 0.5])
            line([2 2],get(gca,'YLim'),'Color',[0 0 0])
            line([3 3],get(gca,'YLim'),'Color',[0 0 0])
            drawnow;
        end
        sgtitle(['Trial ' num2str(idx_trial) '/' num2str(n_trial) ' | cue: ' num2str(cueTYP(idx_trial))]);

        % ask to the user for this trial the start and the stop
        startTask = input('Input the start of ERD: ');
        endTask   = input('Input the end of ERD: ');
    end
end