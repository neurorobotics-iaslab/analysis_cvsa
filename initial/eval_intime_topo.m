%% We use only evaluation file and show the topoplot in time for each cf
clc; clearvars; close all;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/utils');

%% variables
bands = {[8 10], [10 12], [12 14], [14 16], [8 14]};
avg = 1;
bufferSize_integrator = 48;
chunkSize = 32;
nchannels = 39;
bufferSize_processing = 512;
filterOrder = 4;
classes = [730 731];
fix_event = 786;
cf_event = 781;
hit_event = 897;
miss_event = 898;
timeout_event = 899;
startTrial_event = 786; %% we use only data inside fix to boom
nclasses = length(classes);
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
nbands = length(bands);
signals = cell(1, nbands);
headers = cell(1, nbands);
for idx_band = 1:nbands
    headers{idx_band}.TYP = [];
    headers{idx_band}.DUR = [];
    headers{idx_band}.POS = [];
    headers{idx_band}.startNewFile = [];
    headers{idx_band}.ths_rejection = [];
    headers{idx_band}.ths = [];
    headers{idx_band}.bufferSize_integrator = [];
    headers{idx_band}.alpha = [];
    signals{idx_band} = [];
end

%% load the gdfs
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF evaluation Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end

% Concatenate all + processing + extract selected channels
nFiles = length(filenames);
hit = zeros(nFiles, 1); miss = zeros(nFiles, 1); timeout = zeros(nFiles, 1);
for idx_file= 1: nFiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [signal,header] = sload(fullpath_file);
    signal = signal(:,1:nchannels);
    channels_label = header.Label;
    [~, idx_channels_select] = ismember(channels_select, channels_label);
    % load the parameters used in ros
    parameters_file = [pathname(1:end-4) 'parameters/' filenames{idx_file}(1:end-3) 'yaml'];
    a = ReadYaml(parameters_file);
    ths_rejection = a.integrator.thresholds_rejection;
    ths = a.trainingCVSA_node.thresholds;
    if contains(filenames{idx_file}, 'expo')
        alpha = a.integrator.alpha;
    else
        bufferSize = a.integrator.buffer_size;
    end

    for idx_band = 1:nbands
        band = bands{idx_band};
        [signal_processed, header_processed] = processing_offline(signal, header, nchannels, bufferSize_processing, filterOrder, band, chunkSize);
        c_header = headers{idx_band};
        c_header.sampleRate = header_processed.SampleRate;
        c_header.channels_labels = header.Label;
        c_header.TYP = cat(1, c_header.TYP, header_processed.EVENT.TYP);
        c_header.DUR = cat(1, c_header.DUR, header_processed.EVENT.DUR);
        c_header.POS = cat(1, c_header.POS, header_processed.EVENT.POS + size(signals{idx_band}, 1));
        c_header.startNewFile = cat(1, c_header.startNewFile, size(signals{idx_band}, 1) + 1);

        signals{idx_band} = cat(1, signals{idx_band}, signal_processed(:,:));
        c_header.ths = cat(1, c_header.ths, repmat(ths,size(signal_processed, 1), 1));
        c_header.ths_rejection = cat(1, c_header.ths_rejection, repmat(ths_rejection, size(signal_processed, 1), 1));
        if contains(filenames{idx_file}, 'expo')
            c_header.alpha = cat(1, c_header.alpha, repmat(alpha, size(signal_processed, 1), 1));
        else
            c_header.bufferSize_integrator = cat(1, c_header.bufferSize_integrator, repmat(bufferSize, size(signal_processed, 1), 1));
        end
        headers{idx_band} = c_header;
    end

    % for the accuracy
    hit(idx_file) = sum(header.EVENT.TYP == hit_event);
    miss(idx_file) = sum(header.EVENT.TYP == miss_event);
    timeout(idx_file) = sum(header.EVENT.TYP == timeout_event);
end

%% extract all data
data = cell(1, nbands); 
for idx_band = 1:nbands
    [data{idx_band}, info] = extract_cf(signals{idx_band}, headers{idx_band}, classes, cf_event);
end

% reshape the data
maxDUR_cf = max(info.endCf - info.startCf);
ntrials = size(info.startCf, 1);
data_reshape = cell(1, nbands);
for idx_band=1:nbands
    data_reshape{idx_band} = nan(maxDUR_cf, nchannels, ntrials);

    for idx_trial=1:ntrials
        c_start = info.startCf(idx_trial);
        c_end = info.endCf(idx_trial);
        c_dur = c_end - c_start;
        tmp = data{idx_band}.data(c_start:c_end,:);
        data_reshape{idx_band}(1:c_dur+1,:,idx_trial) = tmp;
    end
end

%% variables chancslog
chanlocs_path = '/home/paolo/chanlocs39.mat';
load(chanlocs_path);

%% show topoplot in time
step_time_toplot = 0.25; %sec
rows_plot = 3;
handles = [];
cl = -inf;
for idx_band=1:nbands
    sampleRate = headers{idx_band}.sampleRate;
    nplot  = floor(size(data_reshape{idx_band},  3) / (sampleRate * step_time_toplot));

    figure();
    
    % plot the data
    for idx_nplot = 1:nplot
        start_topo = (idx_nplot - 1) * (sampleRate * step_time_toplot) + 1;
        end_topo   = start_topo + (sampleRate * step_time_toplot) - 1;
        data_1 = nanmean(nanmean(data_reshape{idx_band}(start_topo:end_topo,:, data{idx_band}.typ == classes(1)), 3), 1);
        data_1(isnan(data_1)) = 0;
        nan_trials_1 = any(any(isnan(data_reshape{idx_band}(start_topo:end_topo,:, data{idx_band}.typ == classes(1))), 1), 2);
        data_2 = nanmean(nanmean(data_reshape{idx_band}(start_topo:end_topo,:, data{idx_band}.typ == classes(2)), 3), 1);
        data_2(isnan(data_2)) = 0;
        nan_trials_2 = any(any(isnan(data_reshape{idx_band}(start_topo:end_topo,:, data{idx_band}.typ == classes(2))), 1), 2);
        diff = data_2 - data_1;
        chanlocs_data = zeros(nchannels, 1);
        chanlocs_data(idx_channels_select) = diff(idx_channels_select);
        subplot(rows_plot, ceil(nplot/rows_plot), idx_nplot)
        topoplot(squeeze(chanlocs_data), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(chanlocs_data)) max(abs(chanlocs_data))]);
        cl = max(cl, max(abs(chanlocs_data)));
        handles = [handles gca];
        axis image;
        colorbar;
        title({['CF: ' num2str(((start_topo -1) / sampleRate)) 's - ' num2str((end_topo/sampleRate)) 's'] ...
            ['Trials: ' num2str(sum(data{idx_band}.typ == classes(1)) - sum(nan_trials_1)) '/' num2str(sum(data{idx_band}.typ == classes(1))) ...
            ' ' num2str(sum(data{idx_band}.typ == classes(2)) - sum(nan_trials_2)) '/' num2str(sum(data{idx_band}.typ == classes(2)))]})
    end

    set(handles, 'clim', [-cl cl])
    all_title = ['br-bl | band: ' bands_str{idx_band} ' | all files concatenated'];
    sgtitle(all_title)
end



