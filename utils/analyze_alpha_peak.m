function [peak_freq] = analyze_alpha_peak(filename, varargin)
% ANALYZE_ALPHA_PEAK Calcola PSD via Spettrogramma per evitare salti temporali.
%
% LOGICA AGGIORNATA:
% 1. Filtro 2Hz su tutto il file.
% 2. Calcolo dello SPETTROGRAMMA su tutto il file (slide continuo).
% 3. Selezione delle sole colonne temporali che cadono nei trigger di riposo.
% 4. Media di queste colonne.
%
% Questo evita gli artefatti da concatenazione temporale.

    %% 1. Parsing Input
    p = inputParser;
    addRequired(p, 'filename', @ischar);
    addParameter(p, 'RestTrigger', [], @isnumeric);
    parse(p, filename, varargin{:});
    
    rest_typ = p.Results.RestTrigger;
    target_regions = {'O1', 'O2', 'OZ', 'PO7', 'PO8', 'PO3', 'PO4', 'PZ', 'POZ'};
    
    %% 2. Caricamento e Pre-processing
    disp('--------------------------------------------------');
    disp(['[INFO] Analisi IAF (No-Concatenation) su: ' filename]);
    
    try
        [s, h] = sload(filename);
    catch
        error('Funzione sload non trovata (BIOSIG).');
    end
    
    sampleRate = h.SampleRate;
    signal_raw = s(:, 1:end-1); % Via canale eventi
    
    % --- Filtro Passa-Alto 2Hz (ANT Fix) ---
    disp('[INFO] Applicazione Filtro 2Hz su tutto il flusso...');
    signal_raw = signal_raw - mean(signal_raw, 1); 
    cutoff_hz = 2;
    [b_high, a_high] = butter(4, cutoff_hz / (sampleRate/2), 'high');
    signal_clean = filter(b_high, a_high, signal_raw);
    
    % --- Selezione Canali ---
    ch_indices = [];
    if isfield(h, 'Label')
        all_labels = upper(h.Label);
        for i = 1:length(target_regions)
            idx = find(contains(all_labels, target_regions{i}));
            ch_indices = [ch_indices; idx];
        end
        if isempty(ch_indices)
            warning('Canali Occipitali non trovati. Uso tutti.');
            ch_indices = 1:size(signal_clean,2);
        else
            ch_indices = unique(ch_indices);
            disp(['[INFO] Canali: ' strjoin(h.Label(ch_indices)', ', ')]);
        end
    else
        ch_indices = 1:size(signal_clean,2);
    end

    %% 3. Setup Parametri Spettrali (Simili a Pwelch)
    % Usiamo finestre di 2 secondi per avere risoluzione 0.5 Hz
    window_sec = 2; 
    window_samples = round(window_sec * sampleRate);
    noverlap = round(window_samples / 2); % 50% overlap
    nfft = 4 * sampleRate;
    
    %% 4. Mappatura Temporale (Dove sono i riposi?)
    disp('[INFO] Mappatura finestre temporali...');
    
    % Eseguiamo una chiamata dummy a spectrogram su un canale per ottenere il vettore Tempo T
    % Questo ci dice esattamente a che secondo corrisponde ogni colonna dello spettrogramma
    [~, ~, T_vector] = spectrogram(signal_clean(:,1), window_samples, noverlap, nfft, sampleRate);
    
    % Creiamo una maschera logica per le colonne da tenere
    % T_vector contiene il centro (in secondi) di ogni finestra calcolata
    mask_keep = false(size(T_vector));
    
    if ~isempty(rest_typ) && isfield(h, 'EVENT')
        idx_ev = find(h.EVENT.TYP == rest_typ);
        if isempty(idx_ev), error('Trigger riposo non trovato.'); end
        
        fprintf('[INFO] Individuazione segmenti su %d eventi...\n', length(idx_ev));
        
        for k = 1:length(idx_ev)
            % Tempi dell'evento (in secondi)
            t_start = h.EVENT.POS(idx_ev(k)) / sampleRate;
            t_dur   = h.EVENT.DUR(idx_ev(k)) / sampleRate;
            t_end   = t_start + t_dur;
            
            % Scarta il primo secondo (transiente)
            t_start_valid = t_start + 1.0; 
            
            if t_start_valid < t_end
                % Trova quali finestre dello spettrogramma cadono in questo intervallo
                % Accettiamo finestre il cui centro è dentro l'intervallo valido
                idx_in_window = (T_vector >= t_start_valid) & (T_vector <= t_end);
                mask_keep = mask_keep | idx_in_window;
            end
        end
    else
        warning('Nessun trigger. Uso tutto il tempo.');
        mask_keep = true(size(T_vector));
    end
    
    if sum(mask_keep) == 0
        error('Nessuna finestra temporale valida trovata (durata riposo troppo breve?).');
    end
    
    fprintf('[INFO] Selezionate %d finestre temporali valide su %d totali.\n', sum(mask_keep), length(T_vector));

    %% 5. Calcolo PSD ed Estrazione (Loop sui Canali)
    disp('[INFO] Calcolo Spettrogramma ed Estrazione colonne...');
    
    % Pre-allocazione per la media dei canali
    % spectrogram restituisce P di dimensione [Freq x Time]
    % Noi vogliamo accumulare la media PSD di ogni canale
    
    % Calcoliamo le frequenze una volta
    [~, F_vector] = spectrogram(signal_clean(:,1), window_samples, noverlap, nfft, sampleRate);
    
    all_ch_psds = zeros(length(F_vector), length(ch_indices));
    
    for i = 1:length(ch_indices)
        ch_idx = ch_indices(i);
        
        % A. Calcola Spettrogramma su TUTTO il segnale (no salti)
        % [~, ~, ~, P] restituisce la Power Spectral Density direttamente
        [~, ~, ~, P_full] = spectrogram(signal_clean(:, ch_idx), window_samples, noverlap, nfft, sampleRate);
        
        % B. Estrai SOLO le colonne (tempo) che corrispondono al riposo
        P_rest = P_full(:, mask_keep);
        
        % C. Fai la media temporale subito (otteniamo la PSD media di questo canale)
        all_ch_psds(:, i) = mean(P_rest, 2);
    end
    
    %% 6. Statistiche Spaziali
    % Ora abbiamo all_ch_psds [Freq x Canali]
    
    mean_pxx = mean(all_ch_psds, 2); % Media tra i canali
    std_pxx  = std(all_ch_psds, 0, 2); % Varianza tra i canali
    
    % Conversioni dB e Plotting bounds
    mean_db = 10*log10(mean_pxx);
    upper_db = 10*log10(mean_pxx + std_pxx);
    lower_val = mean_pxx - std_pxx; 
    lower_val(lower_val <= 0) = 1e-10;
    lower_db = 10*log10(lower_val);
    
    %% 7. Ricerca Picco (8-14 Hz)
    search_mask = (F_vector >= 8 & F_vector <= 14);
    f_band = F_vector(search_mask);
    p_band = mean_pxx(search_mask);
    
    if isempty(f_band)
        peak_freq = NaN; max_p_db = 0;
    else
        [max_p, idx_max] = max(p_band);
        peak_freq = f_band(idx_max);
        max_p_db = 10*log10(max_p);
    end
    
    %% 8. Visualizzazione
    figure('Color', 'w', 'Name', 'IAF No-Jumps Analysis'); hold on;
    
    % Plot range 4-30 Hz
    plot_mask = (F_vector >= 4 & F_vector <= 30);
    f_plot = F_vector(plot_mask);
    
    % Area Varianza
    x_poly = [f_plot; flipud(f_plot)];
    y_poly = [upper_db(plot_mask); flipud(lower_db(plot_mask))];
    fill(x_poly, y_poly, [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', 'Std Dev (Canali)');
    
    % Banda Alpha
    y_area = mean_db(search_mask);
    if ~isempty(f_band)
        area(f_band, y_area, 'FaceColor', [0.2 0.8 1], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'DisplayName', 'Alpha (8-14 Hz)');
    end
    
    % Media
    plot(f_plot, mean_db(plot_mask), 'k', 'LineWidth', 2, 'DisplayName', 'PSD Media');
    
    % Picco
    if ~isnan(peak_freq)
        plot(peak_freq, max_p_db, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'IAF Peak');
        xline(peak_freq, 'r--');
        text(peak_freq + 0.5, max_p_db, sprintf('IAF: %.1f Hz', peak_freq), 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    title(['IAF (Analisi Continua) - Peak: ' num2str(peak_freq) ' Hz']);
    xlabel('Frequenza (Hz)'); ylabel('Potenza (dB/Hz)');
    legend('Location', 'northeast'); grid on; xlim([4 30]);
    
    fprintf('\n   >>> PICCO ALPHA (8-14 Hz): %.2f Hz <<<\n', peak_freq);
end