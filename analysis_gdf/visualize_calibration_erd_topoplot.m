% OFFLINE SIMULATION: Topoplot and ERD/ERS visualization for CALIBRATION
% script that emulates the offline processing for Calibration GDF files.
clear all; close all;

addpath(genpath('/home/paolo/bci_vr_ws/src/analysis_bci/equal_ros'));
addpath(genpath('/home/paolo/bci_vr_ws/src/analysis_bci/utils'));

%% Initialization
nchannels = 32; % Numero di canali standard

[filenames, pathname] = uigetfile('*.gdf', 'Select Calibration GDF Files', 'MultiSelect', 'on');
if isequal(filenames, 0)
    disp('No files selected.');
    return;
end
if ischar(filenames)
    filenames = {filenames};
end
nFiles = length(filenames);

freq_band = [8, 14];

for idx_file = 1:nFiles
    fullpath_file_gdf = fullfile(pathname, filenames{idx_file});
    disp(['File ' num2str(idx_file) '/' num2str(nFiles)]);
    disp(['   Loading gdf file : ', filenames{idx_file}]);
    [c_signal, header] = sload(fullpath_file_gdf);
    c_signal = c_signal(:, 1:nchannels);
    channels_label = header.Label(1:nchannels);

    %% ----------------- load the file parameters -----------------
    disp(['   Loading parameters file: ', filenames{idx_file}(1:end-3) 'yaml'])
    fullpath_file_parameters = [pathname(1:end-4) 'parameters/' filenames{idx_file}(1:end-3) 'yaml'];
    
    % La funzione loadParameters è stata modificata per non crashare con i file di calibrazione
    [ringBufferCfg, artifactCfg, processingCfg, qda_mi, qda_cvsa, integratorCfg, paradigm] = loadParameters(fullpath_file_parameters);

    % Determinare le classi in base al paradigm trovato nel yaml (protocol.task)
    if contains(paradigm, 'mi')
        classes = [769 770];
        disp('   Detected MI Paradigm -> Classes: [769, 770]');
    elseif contains(paradigm, 'cvsa')
        classes = [730 731];
        disp('   Detected CVSA Paradigm -> Classes: [730, 731]');
    elseif contains(paradigm, 'hybrid')
        classes = [750 751];
        disp('   Detected HYBRID Paradigm -> Classes: [750, 751]');
    else
        % Fallback se il paradigm non è molto chiaro
        classes = [730 731];
        warning('Paradigm not strictly identified. Using default fallback classes [730, 731]');
    end

    %% ----------------- data processing -----------------
    disp(['   processing EEG data on band ' num2str(freq_band(1)) '-' num2str(freq_band(2)) ' Hz']);
    
    bufferSize = ringBufferCfg.size;
    chunkSize = processingCfg.chunkSize;
    filterOrder = processingCfg.filterOrder;
    do_hann = processingCfg.do_hann;
    
    % Nessun canale EOG da rigettare per default durante l'analisi della calibrazione
    eog_channels = []; 
    
    [signal_processed, header_processed] = processing_onlineROS_CAR_hilbert(c_signal, header, nchannels, bufferSize, filterOrder, freq_band, chunkSize, eog_channels, do_hann);
    
    %% ----------------- epochs extraction -----------------
    disp('   extracting epochs for ERD/ERS and Topoplot') 
    events = header_processed.EVENT;
    sampleRate = header.SampleRate; 
    fs_processed = sampleRate / chunkSize; 
    
    cuePOS = events.POS(ismember(events.TYP, classes));
    cueTYP = events.TYP(ismember(events.TYP, classes));
    
    if isempty(cuePOS)
        warning('Nessun evento della classe trovata nel file di calibrazione. Skipping file.');
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
    figure('Color','w', 'Name', ['Calibration ERD | Band: ' num2str(freq_band(1)) '-' num2str(freq_band(2)) ' Hz']);
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
    
    %% ----------------- Calculate R2, Mean Diff & Topoplots -----------------
    if length(unique_classes) == 2
        disp('   calculating metrics for topoplots...');
        
        % 1. Pulizia dei nomi dei canali (FONDAMENTALE)
        clean_labels = strrep(channels_label, 'EEG ', '');
        clean_labels = strrep(clean_labels, '-Ref', '');
        
        % 2. Generazione sicura delle coordinate (Allineamento BrainProducts -> EEGLAB)
        clear chanlocs;
        for ch = 1:nchannels
            chanlocs(ch).labels = clean_labels{ch};
            
            % Se abbiamo le coordinate, le ruotiamo di 90 gradi per far felice EEGLAB
            if isfield(header, 'ELEC') && ~isempty(header.ELEC) && size(header.ELEC.XYZ, 1) >= ch
                % L'asse Y (Naso in BP) diventa l'asse X di EEGLAB
                chanlocs(ch).X =  header.ELEC.XYZ(ch, 2); 
                
                % L'asse X (Orecchio Dx in BP) diventa l'asse Y invertito (Orecchio Sx in EEGLAB)
                chanlocs(ch).Y = -header.ELEC.XYZ(ch, 1); 
                
                % L'asse Z (Alto) rimane invariato
                chanlocs(ch).Z =  header.ELEC.XYZ(ch, 3);
            end
        end
        
        % Diciamo a EEGLAB di calcolare gli angoli polari basandosi SULLE TUE coordinate appena ruotate
        try
            chanlocs = pop_chanedit(chanlocs, 'convert', {'cart2all'});
        catch
            warning('Conversione coordinate EEG fallita. Il topoplot userà una disposizione approssimata.');
        end

        % 3. Controllo dati estratti
        if isempty(r2_data)
            warning('Nessun dato valido estratto. Impossibile fare plot.');
        end

        % 4. Calcolo dei dati per le due classi
        data_class1 = r2_data(r2_labels == unique_classes(1), :);
        data_class2 = r2_data(r2_labels == unique_classes(2), :);
        
        N1 = size(data_class1, 1);
        N2 = size(data_class2, 1);
        mean1 = mean(data_class1, 1);
        mean2 = mean(data_class2, 1);
        
        % 5. Calcolo Signed R-square
        var1 = var(data_class1, 1);
        var2 = var(data_class2, 1);
        var_tot = (var1 * (N1 - 1) + var2 * (N2 - 1)) / (N1 + N2 - 2); 
        var_tot(var_tot == 0) = eps; % Previene divisioni per zero (crash)
        
        t_values = (mean1 - mean2) ./ sqrt(var_tot .* (1/N1 + 1/N2));
        r2_values = (t_values.^2) ./ (t_values.^2 + N1 + N2 - 2);
        signed_r2_values = r2_values .* sign(mean1 - mean2);
        
        % 6. Calcolo Mean Difference
        mean_diff = mean1 - mean2;

        % --- MODIFICA CRITICA ---
        % Trasformiamo i risultati in VETTORI COLONNA per evitare controlli buggati in topoplot
        signed_r2_values = signed_r2_values(:);
        mean_diff = mean_diff(:);

        % --- PLOT 1: SIGNED R^2 ---
        figure('Color', 'w', 'Name', ['Signed R2 | ' filenames{idx_file}]);
        try
            topoplot(signed_r2_values, chanlocs, 'style', 'both', 'electrodes', 'ptslabels');
            title(['Signed R^2: Class ' num2str(unique_classes(1)) ' vs ' num2str(unique_classes(2))]);
            colorbar;
            limit_r2 = max(abs(signed_r2_values));
            if limit_r2 > 0, caxis([-limit_r2, limit_r2]); end
            colormap jet;
        catch ME
            disp('   [Warning] Topoplot fallito per R2. Uso barre.');
            disp(['   Errore: ' ME.message]);
            b = bar(signed_r2_values); b.FaceColor = 'flat';
            b.CData(signed_r2_values > 0, :) = repmat([0.8 0.2 0.2], sum(signed_r2_values > 0), 1);
            b.CData(signed_r2_values < 0, :) = repmat([0.2 0.4 0.8], sum(signed_r2_values < 0), 1);
            grid on; ylabel('Signed R^2'); xticks(1:nchannels); xticklabels(clean_labels); xtickangle(90);
        end

        % --- PLOT 2: MEAN DIFFERENCE ---
        figure('Color', 'w', 'Name', ['Mean Diff | ' filenames{idx_file}]);
        try
            topoplot(mean_diff, chanlocs, 'style', 'both', 'electrodes', 'ptslabels');
            title(['Mean Diff (Log Power): ' num2str(unique_classes(1)) ' - ' num2str(unique_classes(2))]);
            colorbar;
            limit_diff = max(abs(mean_diff));
            if limit_diff > 0, caxis([-limit_diff, limit_diff]); end
            colormap jet;
        catch ME
            disp('   [Warning] Topoplot fallito per Mean Diff. Uso barre.');
            b = bar(mean_diff); b.FaceColor = 'flat';
            b.CData(mean_diff > 0, :) = repmat([0.8 0.2 0.2], sum(mean_diff > 0), 1);
            b.CData(mean_diff < 0, :) = repmat([0.2 0.4 0.8], sum(mean_diff < 0), 1);
            grid on; ylabel('Mean Diff (Log Power)'); xticks(1:nchannels); xticklabels(clean_labels); xtickangle(90);
        end

    else
        disp('   R2 e Differenza calcolati solo se ci sono esattamente 2 classi!');
    end
    
    disp('   Finished calibration visualization for current file.');
end
