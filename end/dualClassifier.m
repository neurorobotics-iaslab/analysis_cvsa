%% file which shows the features map for each subject the fischer for each 100ms and the topoplot
clc; clear all; close all;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/utils')

%% Variables
bands = {[8 10], [10 12], [12 14], [14 16], [8 14]};
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
classes = [730 731];
cf_event = 781;
fix_event = 786;
nchannels = 39;
nclasses = length(classes);
% channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
channels_select = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
        'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
filterOrder = 4;
avg = 1;
ths_rejection = [0.5 0.5];
ths = [1.0 1.0];
bufferSize_integrator = 64;
alpha = 0.97;

%% load files
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end

%% concatenate the files
nFiles = length(filenames);
for idx_file= 1: nFiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [signal,header] = sload(fullpath_file);
    signal = signal(:,1:nchannels);
    channels_label = header.Label;
    [~, idx_channels_select] = ismember(channels_select, channels_label);

    for idx_band = 1:nbands
        band = bands{idx_band};

        signal_processed = proc_512hz(signal, header.SampleRate, band, filterOrder, avg);
        
        c_header = headers{idx_band};
        c_header.sampleRate = header.SampleRate;
        c_header.channels_labels = header.Label;
        c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP);
        c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR);
        c_header.POS = cat(1, c_header.POS, header.EVENT.POS + size(signals{idx_band}, 1));
        c_header.startNewFile = cat(1, c_header.startNewFile, size(signals{idx_band}, 1) + 1);
        c_header.ths = cat(1, c_header.ths, repmat({ths},size(signal_processed, 1), 1));
        c_header.ths_rejection = cat(1, c_header.ths_rejection, repmat({ths_rejection}, size(signal_processed, 1), 1));

        signals{idx_band} = cat(1, signals{idx_band}, signal_processed(:,:));
        c_header.ths = cat(1, c_header.ths, repmat({ths},size(signal_processed, 1), 1));
        c_header.ths_rejection = cat(1, c_header.ths_rejection, repmat({ths_rejection}, size(signal_processed, 1), 1));
        c_header.bufferSize_integrator = cat(1, c_header.bufferSize_integrator, repmat(bufferSize_integrator, size(signal_processed, 1), 1));
        c_header.alpha = cat(1, c_header.alpha, repmat(alpha, size(signal_processed, 1), 1));
        headers{idx_band} = c_header;
    end
end
subject = filenames{1}(1:2);

%% extract 
data = cell(1, nbands); 
for idx_band = 1:nbands
    data{idx_band} = extract_all(signals{idx_band}, headers{idx_band}, classes, fix_event, cf_event, fix_event);
end

%% Topoplot for 100ms for each runs
chanlocs_path = '/home/paolo/chanlocs39.mat';
load(chanlocs_path);
ntrial_file = size(data{1}.cue, 3) / nFiles;
sampleRate = headers{1}.sampleRate;
time_step = 0.2; %s
period_analysis = 1.0; %s
ncols = period_analysis / time_step;
nsamples_4topoplot = round(sampleRate*time_step);
handles = [];
cl = -inf;
for idx_band=nbands:nbands
    figure();
    for idx_file=1:nFiles
        tmp_data = data{idx_band}.cf(:,:, (idx_file-1)*ntrial_file +1:idx_file*ntrial_file);
        tmp_typ = data{idx_band}.typ((idx_file-1)*ntrial_file +1:idx_file*ntrial_file);
        for idx_col = 1:ncols
            start_period = (idx_col - 1)*nsamples_4topoplot + 1;
            end_period = idx_col * nsamples_4topoplot;
            data_1 = mean(mean(tmp_data(start_period:end_period,:,tmp_typ == classes(1)), 3), 1);
            data_2 = mean(mean(tmp_data(start_period:end_period,:,tmp_typ == classes(2)), 3), 1);

            diff = data_2 - data_1;
            chanlocs_data = zeros(nchannels, 1);
            chanlocs_data(idx_channels_select) = diff(idx_channels_select);
            subplot(nFiles, ncols, (idx_file-1)*ncols + idx_col)
            topoplot(squeeze(chanlocs_data), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(chanlocs_data)) max(abs(chanlocs_data))]);
            cl = max(cl, max(abs(chanlocs_data)));
            handles = [handles gca];
            axis image;
            colorbar;
            title([num2str((idx_col-1)*time_step) '-' num2str(idx_col*time_step)])
            drawnow;
        end
    end
    set(handles, 'clim', [-cl cl])
    sgtitle(['subject: ' subject ' | band: ' bands_str{idx_band}])
end

%% show fischer score
normalize_std = true;
figure();
nfeatures = size(data{idx_band}.cf, 2);
ntrial_file = size(data{1}.cue, 3) / nFiles;
sampleRate = headers{1}.sampleRate;
time_step = 0.5; %s
period_analysis = 1.0; %s
ncols = period_analysis / time_step;
nsamples_4topoplot = round(sampleRate*time_step);
handles = [];
cl = -inf;
for idx_file=1:nFiles
    for idx_col = 1:ncols
        start_period = (idx_col - 1)*nsamples_4topoplot + 1;
        end_period = idx_col * nsamples_4topoplot;
        fischer = nan(nfeatures, nbands);
        for idx_band = 1:nbands
            tmp_data = data{idx_band}.cf(:,:, (idx_file-1)*ntrial_file +1:idx_file*ntrial_file);
            tmp_typ = data{idx_band}.typ((idx_file-1)*ntrial_file +1:idx_file*ntrial_file);

            c_data = tmp_data(start_period:end_period,:,:);
            fischer(:,idx_band) = compute_fiser(classes, c_data, tmp_typ, normalize_std);
        end
        subplot(nFiles, ncols, (idx_file-1)*ncols + idx_col)
        imagesc(fischer')
        cl = max(cl, max(abs(chanlocs_data)));
        handles = [handles gca];
        xticks(1:nchannels);
        xticklabels(channels_label);
        yticks(1:nbands)
        yticklabels(bands_str);
        axis image;
        title([num2str((idx_col-1)*time_step) '-' num2str(idx_col*time_step)])
        drawnow;
    end

    set(handles, 'clim', [-cl cl])
    sgtitle(['subject: ' subject])
end

%% Topoplot and fischer first x ms and from x to end trial
chanlocs_path = '/home/paolo/chanlocs39.mat';
load(chanlocs_path);
ntrial_file = size(data{1}.cue, 3) / nFiles;
sampleRate = headers{1}.sampleRate;
time_division = [0 0.5]*sampleRate;
time_division = cat(2, time_division, size(data{1}.cf,1));
ncols = 2;
handles = [];
cl = -inf;
for idx_band=nbands:nbands
    figure();
    for idx_file=1:nFiles
        tmp_data = data{idx_band}.cf(:,:, (idx_file-1)*ntrial_file +1:idx_file*ntrial_file);
        tmp_typ = data{idx_band}.typ((idx_file-1)*ntrial_file +1:idx_file*ntrial_file);
        for idx_col = 1:ncols
            start_period = time_division(idx_col) + 1;
            end_period = time_division(idx_col + 1);
            data_1 = mean(mean(tmp_data(start_period:end_period,:,tmp_typ == classes(1)), 3), 1);
            data_2 = mean(mean(tmp_data(start_period:end_period,:,tmp_typ == classes(2)), 3), 1);

            diff = data_2 - data_1;
            chanlocs_data = zeros(nchannels, 1);
            chanlocs_data(idx_channels_select) = diff(idx_channels_select);
            subplot(nFiles, ncols, (idx_file-1)*ncols + idx_col)
            topoplot(squeeze(chanlocs_data), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(chanlocs_data)) max(abs(chanlocs_data))]);
            cl = max(cl, max(abs(chanlocs_data)));
            handles = [handles gca];
            axis image;
            colorbar;
            title([num2str(time_division(idx_col)/sampleRate) '-' num2str(time_division(idx_col+1)/sampleRate)])
            drawnow;
        end
    end
    set(handles, 'clim', [-cl cl])
    sgtitle(['subject: ' subject ' | band: ' bands_str{idx_band}])
end

% show fischer score
figure();
nfeatures = size(data{idx_band}.cf, 2);
ncols = 2;
handles = [];
cl = -inf;
for idx_file=1:nFiles
    for idx_col = 1:ncols
        start_period = time_division(idx_col) + 1;
        end_period = time_division(idx_col + 1);
        fischer = nan(nfeatures, nbands);
        for idx_band = 1:nbands
            tmp_data = data{idx_band}.cf(:,:, (idx_file-1)*ntrial_file +1:idx_file*ntrial_file);
            tmp_typ = data{idx_band}.typ((idx_file-1)*ntrial_file +1:idx_file*ntrial_file);

            c_data = tmp_data(start_period:end_period,:,:);
            fischer(:,idx_band) = compute_fiser(classes, c_data, tmp_typ, normalize_std);
        end
        subplot(nFiles, ncols, (idx_file-1)*ncols + idx_col)
        imagesc(fischer')
        cl = max(cl, max(abs(chanlocs_data)));
        handles = [handles gca];
        xticks(1:nchannels);
        xticklabels(channels_label);
        yticks(1:nbands)
        yticklabels(bands_str);
        axis image;
        title([num2str(time_division(idx_col)/sampleRate) '-' num2str(time_division(idx_col+1)/sampleRate)])
        drawnow;
    end

    set(handles, 'clim', [-cl cl])
    sgtitle(['subject: ' subject])
end

