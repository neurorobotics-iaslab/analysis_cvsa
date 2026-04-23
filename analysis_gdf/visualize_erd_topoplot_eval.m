% script that emulates the offline processing to compute ERD/ERS and R^2
% Topoplots for EVALUATION FILES
clear all; close all;

addpath(genpath('/home/paolo/bci_vr_ws/src/analysis_bci/equal_ros'));
addpath(genpath('/home/paolo/bci_vr_ws/src/analysis_bci/utils'));

%% Initialization
nchannels = 32;

[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if isequal(filenames, 0)
    disp('No files selected.');
    return;
end
if ischar(filenames)
    filenames = {filenames};
end
nFiles = length(filenames);

for idx_file = 1:nFiles
    fullpath_file_gdf = fullfile(pathname, filenames{idx_file});
    disp(['File ' num2str(idx_file) '/' num2str(nFiles)]);
    disp(['   Loading gdf file : ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file_gdf);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label(1:nchannels);

    %% ----------------- load the file parameters -----------------
    disp(['   Loading parameters file: ', filenames{idx_file}(1:end-3) 'yaml'])
    fullpath_file_parameters = [pathname(1:end-4) '/parameters/' filenames{idx_file}(1:end-3) 'yaml'];
    [ringBufferCfg, artifactCfg, processingCfg, qda_mi, qda_cvsa, integratorCfg, paradigm] = loadParameters(fullpath_file_parameters);

    if strcmp(paradigm, 'hybrid')
        classes = [750 751]; 
    elseif strcmp(paradigm, 'mi_lhrh')
        qda_path = qda_mi.path_to_model;
        qda_mi.model = loadQDA(qda_path);
        classes = qda_mi.model.classes;
    elseif strcmp(paradigm, 'cvsa_blbr')
        qda_path = qda_cvsa.path_to_model;
        qda_cvsa.model = loadQDA(qda_path);
        classes = qda_cvsa.model.classes;
    end

    %% ----------------- Artifact (just for eog_channels) -----------------
    bufferSize = ringBufferCfg.size;
    chunkSize = processingCfg.chunkSize;
    eog.label = channels_label(cell2mat(artifactCfg.EOG_ch));

    %% ----------------- data processing -----------------
    disp('   processing EEG data');
    [~, indici] = ismember(eog.label, channels_label);
    eog_channels = indici(indici > 0);
    filterOrder = processingCfg.filterOrder;
    
    bands = processingCfg.bands; 
    nbands = size(bands, 1);
    do_hann = processingCfg.do_hann;
    
    % Use the first band for visualization
    idx_bands = 1;
    c_band = bands(idx_bands,:);
    [signal_processed, header_processed] = processing_onlineROS_CAR_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, c_band, chunkSize, eog_channels, do_hann);
    
    %% ----------------- labels for the data -----------------
    disp('   extracting epochs for ERD/ERS and Topoplot') 
    events = header_processed.EVENT;
    sampleRate = header.SampleRate; 
    fs_processed = sampleRate / chunkSize; 
    
    cueDUR = events.DUR(ismember(events.TYP, classes));
    cueTYP = events.TYP(ismember(events.TYP, classes));
    cuePOS = events.POS(ismember(events.TYP, classes));
    
    if isempty(cuePOS)
        warning('Nessun evento della classe trovata. Skipping file.');
        continue;
    end
    
    baseline_sec = 2; % 2 secondi pre-stimolo
    trial_sec = 4;    % 4 secondi post-stimolo
    
    baseline_samples = floor(baseline_sec * fs_processed);
    trial_samples = floor(trial_sec * fs_processed);
    t_axis = linspace(-baseline_sec, trial_sec, baseline_samples + trial_samples + 1);
    
    unique_classes = unique(cueTYP);
    epochs_erd = cell(length(unique_classes), 1);
    r2_data = [];
    r2_labels = [];
    
    for c_idx = 1:length(unique_classes)
        c_class = unique_classes(c_idx);
        c_cues = cuePOS(cueTYP == c_class);
        
        c_epochs = [];
        for i = 1:length(c_cues)
            start_samp = c_cues(i) - baseline_samples;
            end_samp = c_cues(i) + trial_samples;
            
            if start_samp > 0 && end_samp <= size(signal_processed, 1)
                c_pow = signal_processed(start_samp:end_samp, :);
                
                % Baseline Correction (Relative %)
                base_pow = mean(c_pow(1:baseline_samples, :), 1);
                c_erd = (c_pow - repmat(base_pow, size(c_pow, 1), 1)) ./ repmat(base_pow, size(c_pow, 1), 1) * 100;
                
                c_epochs(end+1, :, :) = c_erd;
                
                % Dati per R2 Topoplot: periodo attivo [+0.5 a +3.0]s
                r2_start = c_cues(i) + floor(0.5 * fs_processed);
                r2_end   = c_cues(i) + floor(3.0 * fs_processed);
                if r2_start > 0 && r2_end <= size(signal_processed, 1)
                    r2_data = [r2_data; log(signal_processed(r2_start:r2_end, :))];
                    r2_labels = [r2_labels; repmat(c_class, r2_end - r2_start + 1, 1)];
                end
            end
        end
        epochs_erd{c_idx} = c_epochs;
    end
    
    %% ----------------- Plot ERD/ERS Timecourse -----------------
    chans_to_plot = {'C3', 'C4', 'Cz'}; % standard channels
    figure('Color','w', 'Name', ['ERD Timecourse | Band: ' num2str(c_band(1)) '-' num2str(c_band(2)) ' Hz']);
    colors = {'b', 'r', 'g', 'm'};
    
    for p = 1:length(chans_to_plot)
        ch_name = chans_to_plot{p};
        ch_idx = find(strcmpi(channels_label, ch_name));
        
        if ~isempty(ch_idx)
            subplot(1, length(chans_to_plot), p);
            hold on;
            for c_idx = 1:length(unique_classes)
                if ~isempty(epochs_erd{c_idx})
                    mean_erd = squeeze(mean(epochs_erd{c_idx}(:, :, ch_idx), 1));
                    plot(t_axis, mean_erd, 'Color', colors{c_idx}, 'LineWidth', 2, 'DisplayName', ['Class ' num2str(unique_classes(c_idx))]);
                end
            end
            xline(0, 'k--', 'Cue');
            yline(0, 'k-');
            title(['ERD/ERS - ' ch_name]);
            xlabel('Time [s]');
            ylabel('% ERD/ERS');
            legend('Location', 'best');
            grid on;
        end
    end
    
    %% ----------------- Calculate R2 & Plot Topoplot -----------------
    if length(unique_classes) == 2
        disp('   calculating R2 for topoplot...');
        [r2_values] = calc_r2_from_data(r2_data, r2_labels, 'Plot', false);
        
        figure('Color', 'w', 'Name', ['Signed R2 Topoplot | ' filenames{idx_file}]);
        try
            % Se topoplot e' compatibile, prova a disegnarlo
            disp('   Attempting to plot topoplot...');
            if exist('topoplot', 'file')
                topoplot(r2_values, channels_label);
                title(['Signed R^2: Class ' num2str(unique_classes(1)) ' vs ' num2str(unique_classes(2))]);
                colorbar;
            else
                error('topoplot non presente nel workspace MATLAB.');
            end
        catch ME
            disp('   [Warning] Error plotting topoplot or not configured properly. Using fallback bar chart.');
            disp(['   Info: ' ME.message]);
            
            b = bar(r2_values);
            b.FaceColor = 'flat';
            b.CData(r2_values > 0, :) = repmat([0.8 0.2 0.2], sum(r2_values > 0), 1);
            b.CData(r2_values < 0, :) = repmat([0.2 0.4 0.8], sum(r2_values < 0), 1);
            grid on;
            ylabel('Signed r^2');
            xticks(1:nchannels);
            xticklabels(strrep(channels_label, 'EEG ', ''));
            xtickangle(90);
            title(['Fallback R^2 | Band: ' num2str(c_band(1)) '-' num2str(c_band(2)) ' Hz']);
        end
    else
        disp('   R2 calcolato solo se ci sono esattamente 2 classi!');
    end
    
    disp('   Finished visualization for current file.');
end
