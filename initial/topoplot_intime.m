% clc; clear all;
% compute and show the erd/ers of files in .mat show the chanlocs
% perform the laplacian over all 39 channels and in the chanlocs shows only
% the interest channels

subject = 'g2';

path = ['/home/paolo/cvsa_ws/record/' subject '/mat_selectedTrials'];
% path = ['/home/paolo/cvsa_ws/record/' subject '/gdf'];

lap_path39 = '/home/paolo/laplacians/lap_39ch_CVSA.mat';
chanlog_path = '/home/paolo/new_chanlocs64.mat';

channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
        '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

% channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

classes = [730 731];
sampleRate = 512;
normalization = false; % if false compute the log power, if true the ERD/ERS


%% with all channel for the laplacian
files = dir(fullfile(path, '*.mat'));
total_signal = [];
total_header.EVENT.TYP = [];
total_header.EVENT.DUR = [];
total_header.EVENT.POS = [];
for idx_f = 1:length(files)
    file = fullfile(path, files(idx_f).name);

    disp(['loading file: ' file])
    load(file);
%     [signal, header] = sload(file);
    
    %% ERD_ERS
    [ERD, ~, ~] = ERDERS_CVSA(signal, header, channels_label, lap_path39, sampleRate, chanlog_path, classes, normalization, files(idx_f).name);

    header.EVENT.POS = header.EVENT.POS + length(total_signal);
    total_signal = cat(1, total_signal, signal);
    total_header.EVENT.POS = cat(1, total_header.EVENT.POS, header.EVENT.POS);
    total_header.EVENT.DUR = cat(1, total_header.EVENT.DUR, header.EVENT.DUR);
    total_header.EVENT.TYP = cat(1, total_header.EVENT.TYP, header.EVENT.TYP);
end

%%
disp('All files' )
ERDERS_CVSA(total_signal, total_header, channels_label, lap_path39, sampleRate, chanlog_path, classes, normalization, 'all files merged');
