addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils');
%% from gdf to mat
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
%% concatenate the files
nFiles = length(filenames);
nchannels = 39;
filterOrder = 4;
avg = 1;
bands = {[8 10], [10 12], [12 14], [14 16], [8 14]};
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
info.bands = bands;
nbands = size(bands, 2);
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
info.channel_labels = channels_select;
signals = cell(nbands, 1);

for idx_band = 1:nbands
    signals{idx_band} = [];
end
events.POS = [];
events.DUR = [];
events.TYP = [];
for idx_band = 1:nbands
    for idx_file= 1: nFiles
        fullpath_file = fullfile(pathname, filenames{idx_file});
        disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
        [c_signal,header] = sload(fullpath_file);
        c_signal = c_signal(:,1:nchannels);
        [~, idx_channels_select] = ismember(channels_select, info.channel_labels);
        band = bands{idx_band};

        info.sampleRate = header.SampleRate;

        if idx_file == 1
            POS = header.EVENT.POS + size(signals{idx_band}, 1);
            events.POS = [events.POS; POS];
            events.DUR = [events.DUR; header.EVENT.DUR];
            events.TYP = [events.TYP; header.EVENT.TYP];
        end
        signal_processed = proc_512hz(c_signal(:,idx_channels_select), header.SampleRate, band, filterOrder, avg);
        signals{idx_band} = [signals{idx_band}; signal_processed];
    end
end
%% save
path_to_save = [pathname(1:end-4) 'mat/' pathname(end-15:end-12) '_proc512.mat'];
save(path_to_save, 'events', 'signals', "filenames", 'info');