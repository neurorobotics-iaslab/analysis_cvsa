%% processing without using the function equal to ros
clc; clearvars;% close all;
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
ths_rejection = [0.5 0.5];
bufferSize_integrator = 64;
alpha = 0.97;

%% load the gdfs
[filenames, pathname] = uigetfile('*.gdf', 'Select gdf CALIBRATION files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
nfiles = size(filenames, 2);

%% iterate over files
for idx_file=1:nfiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nfiles)  '): ', filenames{idx_file}]);
    [signal, header] = sload(fullpath_file);
    signal = signal(:, 1:nchannels);
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

        signals{idx_band} = cat(1, signals{idx_band}, signal_processed(:,:));
        c_header.ths = cat(1, c_header.ths, repmat({[1.0 1.0]},size(signal_processed, 1), 1));
        c_header.ths_rejection = cat(1, c_header.ths_rejection, repmat({ths_rejection}, size(signal_processed, 1), 1));
        c_header.bufferSize_integrator = cat(1, c_header.bufferSize_integrator, repmat(bufferSize_integrator, size(signal_processed, 1), 1));
        c_header.alpha = cat(1, c_header.alpha, repmat(alpha, size(signal_processed, 1), 1));
        headers{idx_band} = c_header;
    end
end
subject = filenames{1}(1:2);

%% topoplot for each band in time
chanlocs_path = '/home/paolo/chanlocs39.mat';
load(chanlocs_path);
data = cell(1, nbands); 
for idx_band = 1:nbands
    data{idx_band} = extract_all(signals{idx_band}, headers{idx_band}, classes, fix_event, cf_event, startTrial_event);
end
step_time_toplot = 0.5; %sec
rows_plot = 3;
handles = [];
cl = -inf;
for idx_band=nbands:nbands
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

    drawnow;
    set(handles, 'clim', [-cl cl])
    all_title = ['subject: ' subject ' | br-bl | band: ' bands_str{idx_band} ' | all files concatenated'];
    sgtitle(all_title)
end


% %% diff features plot in mean (br - bl)
% time_plot = 0.5;
% cl = - inf;
% for idx_band=nbands:nbands
%     handles = [];
%     figure();
%     c_fix = nan(size(data{idx_band}.fix, 1), size(data{idx_band}.fix, 2), nclasses);
%     c_cue = nan(size(data{idx_band}.cue, 1), size(data{idx_band}.cue, 2), nclasses);
%     c_cf = nan(size(data{idx_band}.cf, 1), size(data{idx_band}.cf, 2), nclasses);
%     label_ch = [headers{1}.channels_labels(idx_channels_select)];
%     for idx_class=1:nclasses
%         c_fix(:,:,idx_class) = mean(data{idx_band}.fix(:,:, data{idx_band}.typ == classes(idx_class)), 3);
%         c_cue(:,:,idx_class) = mean(data{idx_band}.cue(:,:, data{idx_band}.typ == classes(idx_class)), 3);
%         c_cf(:,:,idx_class)  = mean(data{idx_band}.cf(:,:, data{idx_band}.typ == classes(idx_class)), 3);
%     end
% 
%     subplot(1,3,1)
%     diff_fix = squeeze(c_fix(:,:,2) - c_fix(:,:,1));
%     diff_fix = diff_fix(:, idx_channels_select);
%     imagesc(diff_fix')
%     xlim([0 size(diff_fix, 1)])
%     xticks_ = [(1:sampleRate*time_plot:size(diff_fix, 1))-1 size(diff_fix, 1)];
%     xticks(xticks_)
%     xticklabels(string((xticks_)/sampleRate))
%     ylim([1, size(diff_fix,2)])
%     yticks(1:size(diff_fix,2))
%     yticklabels(label_ch)
%     colorbar;
%     title('fixation')
%     handles = [handles gca];
%     cl = max(cl, max(abs(diff_fix), [], 'all'));
% 
%     subplot(1,3,2)
%     diff_cue = squeeze(c_cue(:,:,2) - c_cue(:,:,1));
%     diff_cue = diff_cue(:, idx_channels_select);
%     imagesc(diff_cue')
%     xlim([0 size(diff_cue, 1)])
%     xticks_ = [(1:sampleRate*time_plot:size(diff_cue, 1))-1 size(diff_cue, 1)];
%     xticks(xticks_)
%     xticklabels(string((xticks_)/sampleRate))
%     ylim([1, size(diff_cue,2)])
%     yticks(1:size(diff_cue,2))
%     yticklabels(label_ch)
%     title('cue')
%     colorbar;
%     handles = [handles gca];
%     cl = max(cl, max(abs(diff_cue), [], 'all'));
% 
%     subplot(1,3,3)
%     diff_cf = squeeze(c_cf(:,:,2) - c_cf(:,:,1));
%     diff_cf = diff_cf(:, idx_channels_select);
%     imagesc(diff_cf')
%     xlim([0 size(diff_cf, 1)])
%     xticks_ = [(1:sampleRate*time_plot:size(diff_cf, 1))-1 size(diff_cf, 1)];
%     xticks(xticks_)
%     xticklabels(string((xticks_)/sampleRate))
%     ylim([1, size(diff_cf,2)])
%     yticks(1:size(diff_cf,2))
%     yticklabels(label_ch)
%     title('cf')
%     handles = [handles gca];
%     cl = max(cl, max(abs(diff_cue), [], 'all'));
%     sgtitle(['subject: ' subject ' | band: ' bands_str{idx_band} ' | diff mean features (br - bl)'])
%     set(handles, 'clim', [-cl cl])
%     colorbar;
%     drawnow;
% end

%% the log band of the trials for each file
fix_data = cell(1, nbands);
cue_data = cell(1, nbands);
cf_data = cell(1, nbands); 
time_plot = 0.5;

for idx_band = 1:nbands
    [cf_data{idx_band}, cf_info] = extract_event(signals{idx_band}, headers{idx_band}, classes, [cf_event]);
    [cue_data{idx_band}, cue_info] = extract_event(signals{idx_band}, headers{idx_band}, classes, classes);
    [fix_data{idx_band}, fix_info] = extract_event(signals{idx_band}, headers{idx_band}, classes, [fix_event]);
end


trials4file = 20;
ntrial = size(cf_data{1}.typ,1);
nfeatures = size(idx_channels_select, 2);
nrows = 2;
for idx_file = 1:nfiles
    for idx_band =nbands:nbands
        figure();
        handles = [];
        cl = -inf;
        for idx_trial=1:trials4file
            c_idx_trial = (idx_file - 1)*trials4file + idx_trial;
            start_cf = cf_info.startEvent(c_idx_trial);
            end_cf = cf_info.endEvent(c_idx_trial);
            start_cue = cue_info.startEvent(c_idx_trial);
            end_cue = cue_info.endEvent(c_idx_trial);
            start_fix = fix_info.startEvent(c_idx_trial);
            end_fix = fix_info.endEvent(c_idx_trial);

            fix_features = fix_data{idx_band}.data(start_fix:end_fix, idx_channels_select);
            cue_features = cue_data{idx_band}.data(start_cue:end_cue, idx_channels_select);
            cf_features = cf_data{idx_band}.data(start_cf:end_cf, idx_channels_select);
            c_features = [fix_features; cue_features; cf_features];

            x = 1:size(c_features, 1) + 1;
            subplot(nrows, ceil(trials4file/nrows), idx_trial)
            hold on
            imagesc(c_features')
            plot([size(fix_features,1) size(fix_features, 1)], ylim, 'Color', 'k');
            plot([size(fix_features,1) + size(cue_features, 1) size(fix_features, 1)+ size(cue_features, 1)], ylim, 'Color', 'k');
            hold off
            xlim([0 size(c_features, 1)])
            xticks_ = [(1:cf_info.sampleRate*time_plot:max(x))-1 max(x)];
            xticks(xticks_)
            xticklabels(string((xticks_)/cf_info.sampleRate))
            ylim([0.5, nfeatures+0.5])
            yticks(1:nfeatures)
            yticklabels([headers{1}.channels_labels(idx_channels_select)])
            title({['cue: ' num2str(cf_data{1}.typ(c_idx_trial)) ' | trial ' num2str(idx_trial)]});
            drawnow;
            handles = [handles gca];
            cl = max(cl, max(abs(c_features), [], 'all'));
        end
        set(handles, 'clim', [0 cl])
        sgtitle(['subject: ' subject ' file: ' filenames{idx_file} ' | band: ' bands_str{idx_band}]);
        colorbar;
    end
end