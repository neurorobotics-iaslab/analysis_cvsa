clear all; % close all;

addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/equal_ros')
addpath('/home/paolo/cvsa/ic_cvsa_ws/src/analysis_bci/utils')

%% Initialization
DATAPAH = '/home/paolo/cvsa/ic_cvsa_ws/src/';
classes = [769 770];
nchannels = 16;
nclasses = length(classes);
filterOrder = 4;
avg = 1;
threshold_gmm_ic = 0.7;
channels_label = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'Fp1', 'CP1', 'Pz', 'CP2', 'Fp2'};


%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);

%% understand the band
nFiles = length(filenames);
peaks = zeros(1, nFiles);
for idx_file = 1:nFiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    peaks(idx_file) = analyze_alpha_peak(fullpath_file, 'RestTrigger', 786, 'band', [8 14], ...
        'target_regions', {'C1', 'C3', 'C2', 'C4'});
end

%% start processing data
bands = [{[8 13]} {[18 24]}];
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
nbands = length(bands);
signals = cell(1, nbands);
artifacts = [];
headers = cell(1, nbands);
for idx_band = 1:nbands
    headers{idx_band}.TYP = [];
    headers{idx_band}.DUR = [];
    headers{idx_band}.POS = [];
    signals{idx_band} = [];
end

for idx_file= 1: nFiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file);
    c_signal = c_signal(:,1:nchannels);
    sampleRate = header.SampleRate;

    excl_ch = {'FP1', 'FP2', 'EOG'};
    [found, indices] = ismember(excl_ch, channels_label);
    excl_chs = indices(found);

    % for power band using hilbert transformation and artefact remotion -----------------------------------------------
    bufferSize = floor(avg*sampleRate);
    chunkSize = 32;
    eog.filterOrder = 4;
    eog.band = [1 10];
    eog.label = {'FP1', 'FP2'};
    eog.h_threshold = 75; %60;
    eog.v_threshold = 75; %60;
    picks.filterOrder = 4;
    picks.freq = 1; % remove antneuro problems
    picks.threshold = 120; %100;
    artifact = artifact_rejection(c_signal, header, nchannels, bufferSize, chunkSize, eog, picks);
    artifacts = cat(1, artifacts, artifact(:,:));

    disp('   [proc] power band');
    for idx_band = 1:nbands
        band = bands{idx_band};

        [signal_processed, header_processed] = processing_onlineROS_CAR_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, band, chunkSize, excl_chs);
        
        c_header = headers{1, idx_band};
        c_header.sampleRate = header_processed.SampleRate/chunkSize;
        c_header.channels_labels = header_processed.Label;
        if isempty(find(header_processed.EVENT.TYP == 2, 1)) % no eye calibration
            c_header.TYP = cat(1, c_header.TYP, header_processed.EVENT.TYP);
            c_header.DUR = cat(1, c_header.DUR, header_processed.EVENT.DUR);
            c_header.POS = cat(1, c_header.POS, header_processed.EVENT.POS + size(signals{1, idx_band}, 1));
        else
            k = find(header_processed.EVENT.TYP == 1, 1);
            c_header.TYP = cat(1, c_header.TYP, header_processed.EVENT.TYP(k:end));
            c_header.DUR = cat(1, c_header.DUR, header_processed.EVENT.DUR(k:end));
            c_header.POS = cat(1, c_header.POS, header_processed.EVENT.POS(k:end) + size(signals{1, idx_band}, 1));
        end
        signals{1, idx_band} = cat(1, signals{1, idx_band}, signal_processed(:,:));
        headers{1, idx_band} = c_header;
    end
end

events = headers{1,1};

%%
nsparsity = 2;
% o_l_ch = {'P3', 'O1', 'P5', 'P1', 'PO5', 'PO3', 'PO7'};
% o_r_ch = {'P4', 'O2', 'P2', 'P6', 'PO4', 'PO6', 'PO8'};
% c_l_ch = {'FC1', 'C3', 'CP1', 'FC3', 'C1', 'CP3'};
% c_r_ch = {'FC2', 'C4', 'CP2', 'FC4', 'C2', 'CP4'};

o_l_ch = {};
o_r_ch = {};
c_l_ch = {'C3', 'CP1', 'C1'};
c_r_ch = {'C4', 'CP2', 'C2'};

[~, o_l] = ismember(o_l_ch, channels_label);
[~, o_r] = ismember(o_r_ch, channels_label);
[~, c_l] = ismember(c_l_ch, channels_label);
[~, c_r] = ismember(c_r_ch, channels_label);

type = 'mi';

sparsity = nan(size(signals{1},1), nsparsity*nbands); % samples x nsparsity*nbands
label_plot = [];

for idx_band = 1:nbands
    c_signal = signals{idx_band};
    for idx_sample = 1:size(c_signal, 1)
        [tmp, label_plot_tmp] =  compute_features_icnic(c_signal(idx_sample,:), type, o_l, o_r, c_l, c_r, nsparsity);
        sparsity(idx_sample, (idx_band-1)*nsparsity+1:(idx_band-1)*nsparsity+nsparsity) = tmp;
    end
    label_plot = [label_plot, label_plot_tmp];
end

for idx_band = 1:nbands
    for idx_s = 1:nsparsity
        idx = (idx_band-1)*nbands+idx_s;
        label_plot{idx} = [label_plot{idx}, ' ', bands_str{idx_band}];
    end
end

%% PLOT INTERATTIVO DEFINITIVO
close all;

% --- CONFIGURAZIONE ---
channels_label = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'Fp1', 'CP1', 'Pz', 'CP2', 'Fp2'};
n_ch_labels = length(channels_label);

samples_window = 1000; % Campioni visibili nella finestra
tot_samples = size(signals{1}, 1);

% Codici ed Eventi
code_fix = 786;         
code_cue = [769 770];   
code_cf  = 781;         

idx_fix = find(events.TYP == code_fix);
idx_cue = find(ismember(events.TYP, code_cue));
idx_cf  = find(events.TYP == code_cf);

% Preparazione Nomi Legenda
legend_str = label_plot;

% Creazione Figura
f = figure('Name', 'Analisi BCI: Log Power & Sparsity', 'Color', 'w');

% =========================================================================
% SUBPLOT 1: Banda 1
% =========================================================================
ax1 = subplot(3,1,1);
imagesc(log(signals{1})'); 
title(['Log Power: ' bands_str{1} ' Hz']);
ylabel('Canali');
set(gca, 'YTick', 1:n_ch_labels, 'YTickLabel', channels_label); 
axis tight; 
% colorbar;
hold on;

% Eventi (Con HandleVisibility = off)
yl = ylim;
for i = 1:length(idx_fix), line([events.POS(idx_fix(i)) events.POS(idx_fix(i))], yl, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); end
for i = 1:length(idx_cue), line([events.POS(idx_cue(i)) events.POS(idx_cue(i))], yl, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off'); end
for i = 1:length(idx_cf),  line([events.POS(idx_cf(i))  events.POS(idx_cf(i))],  yl, 'Color', 'b', 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off'); end
hold off;


% =========================================================================
% SUBPLOT 2: Banda 2
% =========================================================================
ax2 = subplot(3,1,2);
imagesc(log(signals{2})');
title(['Log Power: ' bands_str{2} ' Hz']);
ylabel('Canali');
set(gca, 'YTick', 1:n_ch_labels, 'YTickLabel', channels_label);
axis tight;
% colorbar;
hold on;

% Eventi (Con HandleVisibility = off)
yl = ylim;
for i = 1:length(idx_fix), line([events.POS(idx_fix(i)) events.POS(idx_fix(i))], yl, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); end
for i = 1:length(idx_cue), line([events.POS(idx_cue(i)) events.POS(idx_cue(i))], yl, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off'); end
for i = 1:length(idx_cf),  line([events.POS(idx_cf(i))  events.POS(idx_cf(i))],  yl, 'Color', 'b', 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off'); end
hold off;


% =========================================================================
% SUBPLOT 3: Sparsity Indices
% =========================================================================
ax3 = subplot(3,1,3);

% Plot Sparsity
h_sparsity = plot(sparsity, 'LineWidth', 1.5);

ylim([0 4]); 
title('Indici di Sparsity');
ylabel('Ampiezza');
xlabel('Campioni');
grid on;

hold on;
% Eventi (Con HandleVisibility = off -> IMPOSSIBILE che finiscano in legenda)
yl = ylim; 
for i = 1:length(idx_fix), line([events.POS(idx_fix(i)) events.POS(idx_fix(i))], yl, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); end
for i = 1:length(idx_cue), line([events.POS(idx_cue(i)) events.POS(idx_cue(i))], yl, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off'); end
for i = 1:length(idx_cf),  line([events.POS(idx_cf(i))  events.POS(idx_cf(i))],  yl, 'Color', 'b', 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off'); end
hold off;

% Legenda applicata alla fine (troverà solo h_sparsity perché il resto è invisibile)
legend(h_sparsity, legend_str, 'Location', 'BestOutside', 'Orientation', 'horizontal');


% =========================================================================
% SLIDER
% =========================================================================
linkaxes([ax1, ax2, ax3], 'x');
xlim([1, min(samples_window, tot_samples)]);

sliderPos = [0.15 0.01 0.7 0.04]; 
max_val = max(1, tot_samples - samples_window);
if max_val <= 1, max_val = 2; end 

hSlider = uicontrol('Style', 'slider', ...
    'Parent', f, ...
    'Units', 'normalized', ...
    'Position', sliderPos, ...
    'Min', 1, ...
    'Max', max_val, ...
    'Value', 1, ...
    'Callback', @(src, event) updateXLim(src, [ax1, ax2, ax3], samples_window));

function updateXLim(sliderHandle, axesHandles, windowSize)
    startVal = round(get(sliderHandle, 'Value'));
    endVal = startVal + windowSize;
    set(axesHandles, 'XLim', [startVal, endVal]);
end
