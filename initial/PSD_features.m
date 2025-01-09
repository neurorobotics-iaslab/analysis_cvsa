% clc; clear all; close all;
% show fisher and cva of the PSD signal of the channels selected

%% initial informations
subject = 'h8';
day = '/20241015';
% path = ['/home/paolo/cvsa_ws/record/' subject '/mat_selectedTrials'];
path = ['/home/paolo/cvsa_ws/recordings/' subject day '/gdf/calibration'];
files = dir(fullfile(path, '*.gdf'));

lap_path39 = '/home/paolo/laplacians/lap_39ch_CVSA.mat';

channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
        '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

% channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

classes = [730 731];
sampleRate = 512;
% selFreqs = 4:2:48;
selFreqs = 8:2:14;

%% for total file merged
total_features = [];
total_header.EVENT.POS = [];
total_header.EVENT.DUR = [];
total_header.EVENT.TYP = [];

for idx_file=1:length(files)
    disp(['file ' num2str(idx_file) '/' num2str(length(files))])
    
    file = fullfile(path, files(idx_file).name);
    disp(['   [info] file: ' file])
%     load(file);
    [signal, header] = sload(file);

    %% processing online
    disp('      [PROC] Apply lap filter');
    load(lap_path39)
    signal = signal(:,1:length(channels_label));
    signal = signal - mean(signal, 2);
    signal = signal * lap;
    disp('      [PROC] Apply PSD');
    psd_wlength = 0.5;
    psd_wshift = 0.0625;
    psd_pshift = 0.25;
    psd_mlength =  1;
    [features, f] = proc_spectrogram(signal, psd_wlength, psd_wshift, psd_pshift, sampleRate, psd_mlength);
    disp('      [PROC] Apply log');
    features = log(features);
    header.EVENT.POS = proc_pos2win(header.EVENT.POS, psd_wshift*sampleRate, 'backward', psd_mlength*sampleRate);
    header.EVENT.DUR = floor(header.EVENT.DUR/(psd_wshift*sampleRate)) + 1;

    total_header.EVENT.POS = cat(1, total_header.EVENT.POS, header.EVENT.POS + size(total_features, 1));
    total_features = cat(1, total_features, features);
    total_header.EVENT.DUR = cat(1, total_header.EVENT.DUR, header.EVENT.DUR);
    total_header.EVENT.TYP = cat(1, total_header.EVENT.TYP, header.EVENT.TYP);

    %% fisherscore
    idx_freqs = find(ismember(f,selFreqs));
    idx_interest_ch = find(~strcmp(channels_label, ''));
    interval_step = 1; %in sec
    frameRate = 1/psd_wshift;
    [fisher, cva] = fisherAndCVA_CVSA(features, header, [730 731],  idx_freqs, idx_interest_ch, interval_step, frameRate);

    %% show fisher
    %{
    disp('[proc] |- Visualizing fisher score for offline runs');
    fig1 = figure;
    n = ceil(size(fisher,1)/2);
    climits = [];
    handles = nan(size(fisher,1), 1);
    
    for idx_inter=1:size(fisher,1)
        subplot(2, n, idx_inter)
        colormap('jet');
        imagesc(squeeze(fisher(idx_inter,:,:))');
        axis square;
        colorbar;
        set(gca, 'XTick', 1:length(selFreqs));
        set(gca, 'XTickLabel', selFreqs);
        set(gca, 'YTick', 1:length(idx_interest_ch));
        set(gca, 'YTickLabel', channels_label(idx_interest_ch));
        %xtickangle(90);
        xlabel('Hz');
        ylabel('channel');
        if idx_inter==size(fisher,1)
            c_title = 'Fisher score: cue + cf';
        elseif idx_inter==size(fisher,1)-1
            c_title = 'Fisher score: cf';
        elseif idx_inter==size(cva,1)-2
            c_title = ['CVA: from ' num2str(interval_step*(idx_inter-1)) 's to the finish of cf'];
        else
            c_title = ['Fisher score: from ' num2str(idx_inter-1) 's to ' num2str(idx_inter) 's of cf'];
        end
        title(c_title); 

        climits = cat(2, climits, get(gca, 'CLim'));
        handles(idx_inter) = gca;
    end
    set(handles, 'clim', [0 max(max(climits))]);
    sgtitle(['file: ' files(idx_file).name]);   
    
    %}

    %% show cva
    disp('[proc] |- Visualizing cva for offline runs');
    fig1 = figure;
    n = ceil(size(cva,1)/2);
    climits = [];
    handles = nan(size(cva,1), 1);
    
    for idx_inter=1:size(cva,1)
        subplot(2, n, idx_inter)
        colormap('jet');
        imagesc(squeeze(cva(idx_inter,:,:))');
        axis square;
        colorbar;
        set(gca, 'XTick', 1:length(selFreqs));
        set(gca, 'XTickLabel', selFreqs);
        set(gca, 'YTick', 1:length(idx_interest_ch));
        set(gca, 'YTickLabel', channels_label(idx_interest_ch));
        %xtickangle(90);
        xlabel('Hz');
        ylabel('channel');
        if idx_inter==size(cva,1)
            c_title = 'CVA: cue + cf';
        elseif idx_inter==size(cva,1)-1
            c_title = 'CVA: cf';
        elseif idx_inter==size(cva,1)-2
            c_title = ['CVA: from ' num2str(interval_step*(idx_inter-1)) 's to the finish of cf'];
        else
            c_title = ['CVA: from ' num2str(interval_step*(idx_inter-1)) 's to ' num2str(interval_step*idx_inter) 's of cf'];
        end
        title(c_title); 

        climits = cat(2, climits, get(gca, 'CLim'));
        handles(idx_inter) = gca;
    end
    set(handles, 'clim', [0 max(max(climits))]);
    sgtitle(['file: ' files(idx_file).name]);

end

%% show the merging of the two files for cva and fisher
[fisher, cva] = fisherAndCVA_CVSA(total_features, total_header, [730 731],  idx_freqs, idx_interest_ch, interval_step, frameRate);
disp('[proc] |- Visualizing cva for offline runs');
fig1 = figure;
n = ceil(size(cva,1)/2);
climits = [];
handles = nan(size(cva,1), 1);

for idx_inter=1:size(cva,1)
    subplot(2, n, idx_inter)
    colormap('jet');
    imagesc(squeeze(cva(idx_inter,:,:))');
    axis square;
    colorbar;
    set(gca, 'XTick', 1:length(selFreqs));
    set(gca, 'XTickLabel', selFreqs);
    set(gca, 'YTick', 1:length(idx_interest_ch));
    set(gca, 'YTickLabel', channels_label(idx_interest_ch));
    %xtickangle(90);
    xlabel('Hz');
    ylabel('channel');
    if idx_inter==size(cva,1)
        c_title = 'CVA: cue + cf';
    elseif idx_inter==size(cva,1)-1
        c_title = 'CVA: cf';
    elseif idx_inter==size(cva,1)-2
        c_title = ['CVA: from ' num2str(interval_step*(idx_inter-1)) 's to the finish of cf'];
    else
        c_title = ['CVA: from ' num2str(interval_step*(idx_inter-1)) 's to ' num2str(interval_step*idx_inter) 's of cf'];
    end
    title(c_title);

    climits = cat(2, climits, get(gca, 'CLim'));
    handles(idx_inter) = gca;
end
set(handles, 'clim', [0 max(max(climits))]);
sgtitle('merge of all files');

% fisher
%{
disp('[proc] |- Visualizing fisher for offline runs');
fig1 = figure;
n = ceil(size(fisher,1)/2);
climits = [];
handles = nan(size(fisher,1), 1);

for idx_inter=1:size(fisher,1)
    subplot(2, n, idx_inter)
    colormap('jet');
    imagesc(squeeze(fisher(idx_inter,:,:))');
    axis square;
    colorbar;
    set(gca, 'XTick', 1:length(selFreqs));
    set(gca, 'XTickLabel', selFreqs);
    set(gca, 'YTick', 1:length(idx_interest_ch));
    set(gca, 'YTickLabel', channels_label(idx_interest_ch));
    %xtickangle(90);
    xlabel('Hz');
    ylabel('channel');
    if idx_inter==size(fisher,1)
        c_title = 'fisher: cue + cf';
    elseif idx_inter==size(fisher,1)-1
        c_title = 'fisher: cf';
    elseif idx_inter==size(cva,1)-2
            c_title = ['CVA: from ' num2str(interval_step*(idx_inter-1)) 's to the finish of cf'];
    else
        c_title = ['fisher: from ' num2str(interval_step*(idx_inter-1)) 's to ' num2str(interval_step*idx_inter) 's of cf'];
    end
    title(c_title);

    climits = cat(2, climits, get(gca, 'CLim'));
    handles(idx_inter) = gca;
end
set(handles, 'clim', [0 max(max(climits))]);
sgtitle('merge of all files');
%}