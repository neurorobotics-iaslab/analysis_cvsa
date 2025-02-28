%% file which shows the features map for each subject the fischer for each 100ms and the topoplot
clc; clear all; close all;

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
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
filterOrder = 4;
avg = 1;
ths_rejection = [0.5 0.5];
ths = [1.0 1.0];
bufferSize_integrator = 64;
alpha = 0.97;

%% load files
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF evaluation Files', 'MultiSelect', 'on');
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

%% extract 
data = cell(1, nbands); 
for idx_band = 1:nbands
    data{idx_band} = extract_all(signals{idx_band}, headers{idx_band}, classes, fix_event, cf_event, fix_event);
end

%% Topoplot for 100ms for each runs
ntrial_file = size(data{1}.cue, 3) / nFiles;
sampleRate = headers{1}.sampleRate;
time_step = 0.1; %ms
period_analysis = 0.5; %ms
ncols = period_analysis / time_step;
figure();
for idx_file=1:nFiles
    for idx_col = 1:ncols
        
    end
end

