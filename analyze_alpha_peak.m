function [peak_freq] = analyze_alpha_peak(filename, varargin)
% ANALYZE_ALPHA_PEAK_ANT Calcola PSD, IAF e Varianza (Specifico per ANT Neuro).
%
% Novità:
% - Mostra la variabilità (Deviazione Standard) tra i canali come area grigia.
% - Evidenzia la banda 8-14 Hz.
% - Filtro Passa-Alto 2Hz incluso.
%
% Usage:
%   peak = analyze_alpha_peak_ant('soggetto01.gdf'); 
%   peak = analyze_alpha_peak_ant('soggetto01.gdf', 'RestTrigger', 786);

    %% 1. Parsing Input
    p = inputParser;
    addRequired(p, 'filename', @ischar);
    addParameter(p, 'RestTrigger', [], @isnumeric);
    parse(p, filename, varargin{:});
    
    rest_typ = p.Results.RestTrigger;
    
    % Canali target (Occipitali/Parietali)
    target_regions = {'O1', 'O2', 'Oz', 'PO7', 'PO8', 'PO3', 'PO4', 'Pz', 'POz'};
    
    %% 2. Caricamento Dati
    disp('--------------------------------------------------');
    disp(['[INFO] Analisi IAF + Varianza su: ' filename]);
    
    try
        [s, h] = sload(filename);
    catch
        error('Funzione sload non trovata. Installa BIOSIG.');
    end
    
    fs = h.SampleRate;
    
    % Rimuovi canale eventi
    n_ch_total = size(s, 2);
    if n_ch_total > 1 && max(s(:,end)) > 100 
        signal_raw = s(:, 1:end-1);
    else
        signal_raw = s;
    end

    %% 3. Filtro Passa-Alto 2Hz (ANT Fix)
    disp('[INFO] Applicazione Filtro Passa-Alto 2Hz...');
    signal_raw = signal_raw - mean(signal_raw, 1); % Rimuovi DC
    
    cutoff_hz = 2;
    [b_high, a_high] = butter(4, cutoff_hz / (fs/2), 'high');
    signal_clean = filtfilt(b_high, a_high, signal_raw);
    
    %% 4. Selezione Canali
    ch_indices = [];
    if isfield(h, 'Label')
        all_labels = upper(h.Label);
        for i = 1:length(target_regions)
            idx = find(contains(all_labels, target_regions{i}));
            ch_indices = [ch_indices; idx];
        end
        if isempty(ch_indices)
            warning('Canali Occipitali non trovati. Uso tutti i canali.');
            ch_indices = 1:size(signal_clean,2);
        else
            ch_indices = unique(ch_indices);
            ch_indices = ch_indices(ch_indices <= size(signal_clean,2));
            disp(['[INFO] Canali analizzati: ' strjoin(h.Label(ch_indices)', ', ')]);
        end
    else
        warning('Nessuna Label. Uso tutti i canali.');
        ch_indices = 1:size(signal_clean,2);
    end
    
    %% 5. Estrazione Segmenti (Riposo/Fixation)
    data_segments = [];
    
    if ~isempty(rest_typ) && isfield(h, 'EVENT')
        idx_ev = find(h.EVENT.TYP == rest_typ);
        if isempty(idx_ev)
            warning(['Trigger ' num2str(rest_typ) ' non trovato. Uso tutto il file.']);
            data_segments = signal_clean(:, ch_indices);
        else
            fprintf('[INFO] Trovati %d segmenti di riposo.\n', length(idx_ev));
            for k = 1:length(idx_ev)
                pos = h.EVENT.POS(idx_ev(k));
                dur = h.EVENT.DUR(idx_ev(k));
                
                % Scarta il primo secondo (transiente)
                start_s = pos + fs; 
                stop_s = pos + dur - 1;
                
                if stop_s <= size(signal_clean,1) && start_s < stop_s
                    seg = signal_clean(start_s:stop_s, ch_indices);
                    data_segments = [data_segments; seg];
                end
            end
        end
    else
        disp('[INFO] Analisi sull''intero file.');
        data_segments = signal_clean(:, ch_indices);
    end
    
    if isempty(data_segments), error('Nessun dato estratto.'); end

    %% 6. Calcolo PSD e Statistiche
    disp('[INFO] Calcolo PSD e Varianza spaziale...');
    
    % Parametri Welch
    window = 2 * fs;
    noverlap = window / 2; 
    nfft = 4 * fs;          
    
    % pxx sarà [Freq x Canali]
    [pxx, f] = pwelch(data_segments, window, noverlap, nfft, fs);
    
    % Statistiche TRA i canali (Spatial Variability)
    mean_pxx = mean(pxx, 2);       % Media
    std_pxx  = std(pxx, 0, 2);     % Deviazione Standard
    
    % Calcoliamo i bordi per l'area grigia (Varianza) in dB
    % Usiamo la formula dei dB: 10*log10(valore)
    mean_db = 10*log10(mean_pxx);
    upper_db = 10*log10(mean_pxx + std_pxx);
    % Protezione per log(numeri negativi/zero) nel lower bound
    lower_val = mean_pxx - std_pxx;
    lower_val(lower_val <= 0) = 1e-10; % piccolo epsilon
    lower_db = 10*log10(lower_val);
    
    %% 7. Ricerca Picco (Banda 8-14 Hz)
    % Aggiornato a 8-14 come richiesto
    search_idx = (f >= 8 & f <= 14);
    
    f_band = f(search_idx);
    p_band = mean_pxx(search_idx);
    
    [max_p, idx_max] = max(p_band);
    peak_freq = f_band(idx_max);
    max_p_db = 10*log10(max_p);
    
    %% 8. Visualizzazione
    figure('Color', 'w', 'Name', 'IAF con Varianza');
    hold on;
    
    % Plot range
    plot_idx = (f >= 2 & f <= 30);
    f_plot = f(plot_idx);
    
    % A. Disegna l'area della Varianza (Grigio)
    % Creiamo il poligono per 'fill'
    x_poly = [f_plot; flipud(f_plot)];
    y_poly = [upper_db(plot_idx); flipud(lower_db(plot_idx))];
    
    fill(x_poly, y_poly, [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
        'DisplayName', 'Deviazione Standard (Canali)');
    
    % B. Evidenzia la banda Alpha 8-14 Hz (Azzurrino)
    % Per colorare solo sotto la curva media
    y_area = mean_db(search_idx);
    area(f_band, y_area, 'FaceColor', [0.2 0.8 1], 'FaceAlpha', 0.4, ...
        'EdgeColor', 'none', 'DisplayName', 'Banda Alpha (8-14 Hz)');
    
    % C. Plot della Media (Linea Nera)
    plot(f_plot, mean_db(plot_idx), 'k', 'LineWidth', 2, 'DisplayName', 'PSD Media');
    
    % D. Segna il Picco
    plot(peak_freq, max_p_db, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'IAF Peak');
    
    % E. Cosmetica
    xline(peak_freq, 'r--');
    text(peak_freq + 0.5, max_p_db, sprintf('IAF: %.1f Hz', peak_freq), ...
        'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
    
    title(['Spettro di Potenza (HP 2Hz) - IAF: ' num2str(peak_freq) ' Hz']);
    xlabel('Frequenza (Hz)');
    ylabel('Potenza (dB/Hz)');
    legend('Location', 'northeast');
    grid on;
    xlim([2 30]);
    
    fprintf('\n   >>> PICCO ALPHA (8-14 Hz): %.2f Hz <<<\n', peak_freq);
end