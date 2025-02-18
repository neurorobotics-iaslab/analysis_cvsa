%% for each file compute the fisher for a band and the topoplot. Obv we use the min dur of the cf/cue/fix
clc; clear all; %close all;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/utils');

%% files variables
subject  = 'c7';
day      = '20250217';
path_gdf = ['/home/paolo/cvsa_ws/record/' subject '/' day '/calibration/gdf'];
files = dir(fullfile(path_gdf, '*.gdf')); 

%% processing variables
avg = 1;
filterOrder = 4;
nchannels = 39;
windows_size = 512;
classes = [730 731];
fix_event = 786;
cf_event = 781;
startTrial_event = 786; %% in future check for fake_rest
nclasses = length(classes);
bands = {[8 10] [10 12] [12 14] [14 16] [8 14]};
% bands = {[ 8 14]};
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
% channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
channels_select = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
        'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
nselected_channels = size(channels_select, 2);
normalize_std = false;
normalization_baseline = false;

%% variables chancslog
chanlocs_path = '/home/paolo/chanlocs39.mat';
load(chanlocs_path);

%% Concatenate all + processing + extract selected channels
for idx_file= 1: length(files)
    file = fullfile(path_gdf, files(idx_file).name);
    disp(['file (' num2str(idx_file) '/' num2str(length(files))  '): ', file]);
    [signal,header] = sload(file);
    signal = signal(:,1:nchannels);
    channels_label = header.Label;
    [~, idx_channels_select] = ismember(channels_select, channels_label);

    for idx_band = 1:nbands
        band = bands{idx_band};
        [signal_processed, header_ros] = processing_offline(signal, header, nchannels, header.SampleRate,filterOrder, band, 32);
        c_header = headers{idx_band};
        c_header.sampleRate = header_ros.SampleRate;
        c_header.channels_labels = header_ros.Label(idx_channels_select);
        c_header.TYP = cat(1, c_header.TYP, header_ros.EVENT.TYP);
        c_header.DUR = cat(1, c_header.DUR, header_ros.EVENT.DUR);
        c_header.POS = cat(1, c_header.POS, header_ros.EVENT.POS + size(signals{idx_band}, 1));

        signals{idx_band} = cat(1, signals{idx_band}, signal_processed(:,idx_channels_select));
        headers{idx_band} = c_header;
    end
end

%% extract all data
% for each band the signal are: samples x channels x trial
%   in data the signal keep this division for data.trial, data.cue data.cf data.fix
data = cell(1, nbands); 
for idx_band = 1:nbands
    data{idx_band} = extract_all(signals{idx_band}, headers{idx_band}, classes, fix_event, cf_event, startTrial_event);
end

%% compute fisher score for all the files
fisher_band = nan(nbands, nselected_channels);
for idx_band = 1:nbands
    fisher_band(idx_band, :) = compute_fiser(classes, data{idx_band}.cf, data{idx_band}.typ, normalize_std);
end

% show fisher score
figure;
colormap('jet');
imagesc(fisher_band');
axis square;
colorbar;
set(gca, 'XTick', 1:nbands);
set(gca, 'XTickLabel', bands_str);
set(gca, 'YTick', 1:nselected_channels);
set(gca, 'YTickLabel', headers{1}.channels_labels);
xtickangle(90);
xlabel('Hz');
ylabel('channel');
title(['Total FS Subj: ' subject]);


%% Compute cva for all the files
cva = nan(nbands, nselected_channels);
for idx_band=1:nbands
    cva(idx_band, :) = compute_cva(data{idx_band}.cf, data{idx_band}.typ);
end

% show fisher score
figure;
colormap('jet');
imagesc(cva');
axis square;
colorbar;
set(gca, 'XTick', 1:nbands);
set(gca, 'XTickLabel', bands_str);
set(gca, 'YTick', 1:nselected_channels);
set(gca, 'YTickLabel', headers{1}.channels_labels);
xtickangle(90);
xlabel('Hz');
ylabel('channel');
title(['Total CVA Subj: ' subject]);

%% compute the topoplot
topo = cell(1, nbands);
for idx_band = 1:nbands
    if normalization_baseline
        baseline = mean(data{1}.fix, 1)
        topo{idx_band}.cf = data{idx_band}.cf ./ repmat(baseline, [size(data{idx_band}.cf, 1) 1 1]);
    else
        topo{idx_band}.cf = data{idx_band}.cf;
    end
end

% show the plotting
figure();
handles = [];
cl = -inf;
for idx_band = 1:nbands
    data_1 = mean(mean(topo{idx_band}.cf(:,:,data{idx_band}.typ == classes(1)), 3), 1);
    data_2 = mean(mean(topo{idx_band}.cf(:,:,data{idx_band}.typ == classes(2)), 3), 1);
    diff = data_2 - data_1;
    chanlocs_data = zeros(nchannels, 1);
    chanlocs_data(idx_channels_select) = diff(idx_channels_select);
    subplot(2, ceil(nbands/2), idx_band)
    topoplot(squeeze(chanlocs_data), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(chanlocs_data)) max(abs(chanlocs_data))], 'electrodes', 'labelpoint');
    cl = max(cl, max(abs(chanlocs_data)));
    handles = [handles gca];
    axis image;
    colorbar;
    title(['band: ' bands_str{idx_band}])
end


set(handles, 'clim', [-cl cl])
all_title = 'br-bl | only cf | all files concatenated';
sgtitle(all_title)

%% show topoplot in time
step_time_toplot = 0.25; %sec
rows_plot = 3;
handles = [];
cl = -inf;
for idx_band=1:nbands
    sampleRate = headers{idx_band}.sampleRate;
    nplotFIX = floor(size(data{idx_band}.fix, 1) / (sampleRate * step_time_toplot));
    nplotCUE = floor(size(data{idx_band}.cue, 1) / (sampleRate * step_time_toplot));
    nplotCF  = floor(size(data{idx_band}.cf,  1) / (sampleRate * step_time_toplot));
    nplot = nplotFIX + nplotCUE + nplotCF;

    figure();
    
    % plot the fixation
    for idx_nplot = 1:nplotFIX
        start_topo = (idx_nplot - 1) * (sampleRate * step_time_toplot) + 1;
        end_topo   = start_topo + (sampleRate * step_time_toplot) - 1;
        data_1 = mean(mean(data{idx_band}.fix(start_topo:end_topo,:, data{idx_band}.typ == classes(1)), 3), 1);
        data_2 = mean(mean(data{idx_band}.fix(start_topo:end_topo,:, data{idx_band}.typ == classes(2)), 3), 1);
        diff = data_2 - data_1;
        chanlocs_data = zeros(nchannels, 1);
        chanlocs_data(idx_channels_select) = diff(idx_channels_select);
        subplot(rows_plot, ceil(nplot/rows_plot), idx_nplot)
        topoplot(squeeze(chanlocs_data), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(chanlocs_data)) max(abs(chanlocs_data))]);
        cl = max(cl, max(abs(chanlocs_data)));
        handles = [handles gca];
        axis image;
        colorbar;
        title(['FIX: ' num2str(((start_topo -1 ) / sampleRate)) 's - ' num2str((end_topo/sampleRate)) 's'])
    end

    % plot the cue
    for idx_nplot = 1:nplotCUE
        start_topo = (idx_nplot - 1) * (sampleRate * step_time_toplot) + 1;
        end_topo   = start_topo + (sampleRate * step_time_toplot) - 1;
        data_1 = mean(mean(data{idx_band}.cue(start_topo:end_topo,:, data{idx_band}.typ == classes(1)), 3), 1);
        data_2 = mean(mean(data{idx_band}.cue(start_topo:end_topo,:, data{idx_band}.typ == classes(2)), 3), 1);
        diff = data_2 - data_1;
        chanlocs_data = zeros(nchannels, 1);
        chanlocs_data(idx_channels_select) = diff(idx_channels_select);
        subplot(rows_plot, ceil(nplot/rows_plot), idx_nplot + nplotFIX)
        topoplot(squeeze(chanlocs_data), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(chanlocs_data)) max(abs(chanlocs_data))]);
        cl = max(cl, max(abs(chanlocs_data)));
        handles = [handles gca];
        axis image;
        colorbar;
        title(['CUE: ' num2str(((start_topo -1)/ sampleRate)) 's - ' num2str((end_topo/sampleRate)) 's'])
    end

    % plot the cf
    for idx_nplot = 1:nplotCF
        start_topo = (idx_nplot - 1) * (sampleRate * step_time_toplot) + 1;
        end_topo   = start_topo + (sampleRate * step_time_toplot) - 1;
        data_1 = mean(mean(data{idx_band}.cf(start_topo:end_topo,:, data{idx_band}.typ == classes(1)), 3), 1);
        data_2 = mean(mean(data{idx_band}.cf(start_topo:end_topo,:, data{idx_band}.typ == classes(2)), 3), 1);
        diff = data_2 - data_1;
        chanlocs_data = zeros(nchannels, 1);
        chanlocs_data(idx_channels_select) = diff(idx_channels_select);
        subplot(rows_plot, ceil(nplot/rows_plot), idx_nplot + nplotFIX + nplotCUE)
        topoplot(squeeze(chanlocs_data), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(chanlocs_data)) max(abs(chanlocs_data))]);
        cl = max(cl, max(abs(chanlocs_data)));
        handles = [handles gca];
        axis image;
        colorbar;
        title(['CF: ' num2str(((start_topo -1) / sampleRate)) 's - ' num2str((end_topo/sampleRate)) 's'])
    end

    set(handles, 'clim', [-cl cl])
    all_title = ['br-bl | band: ' bands_str{idx_band} ' | all files concatenated'];
    sgtitle(all_title)
end

