function plot_erd_timecourse(filename)
    % Carica file
    [s, h] = sload(filename);
    fs = h.SampleRate;
    
    % Filtra nella banda del soggetto (es. 8-12 Hz)
    % IMPORTANTE: Usa la banda trovata con la PSD!
    band = [8 12]; 
    [b,a] = butter(4, band/(fs/2), 'bandpass');
    s_filt = filtfilt(b,a, s(:, 1:end-1));
    
    % Inviluppo
    pow = abs(hilbert(s_filt)).^2;
    
    % Canali da analizzare (es. PO8 - Emisfero Destro)
    % Ci aspettiamo che PO8 salga per "Left" (Ipsilaterale) e scenda per "Right" (Contralaterale)
    ch_name = 'PO8'; 
    ch_idx = find(contains(h.Label, ch_name, 'IgnoreCase',true), 1);
    
    % Estrai epoche (-1s prima del cue, +4s dopo)
    events = h.EVENT;
    cues = find(events.TYP == 730 | events.TYP == 731);
    
    epoch_len = 4 * fs; 
    baseline_len = 1 * fs;
    
    trials_730 = [];
    trials_731 = [];
    t_axis = linspace(-1, 4, epoch_len + baseline_len);
    
    for i = 1:length(cues)
        pos = events.POS(cues(i));
        typ = events.TYP(cues(i));
        
        start = pos - baseline_len;
        stop = pos + epoch_len - 1;
        
        if start < 1 || stop > size(pow,1), continue; end
        
        % Estrai segmento
        segment = pow(start:stop, ch_idx);
        
        % Baseline Correction (Diviso per la media pre-stimolo)
        base_val = mean(segment(1:baseline_len));
        segment_rel = (segment - base_val) / base_val * 100; % % Change
        
        if typ == 730
            trials_730 = [trials_730, segment_rel];
        else
            trials_731 = [trials_731, segment_rel];
        end
    end
    
    % Media
    mean_730 = mean(trials_730, 2);
    mean_731 = mean(trials_731, 2);
    
    % Plot
    figure('Color','w'); hold on;
    plot(t_axis, mean_730, 'b', 'LineWidth', 2, 'DisplayName', 'Class 730 (Left)');
    plot(t_axis, mean_731, 'r', 'LineWidth', 2, 'DisplayName', 'Class 731 (Right)');
    xline(0, 'k--');
    
    title(['ERD Time Course - Canale ' ch_name]);
    xlabel('Tempo dal Cue (s)');
    ylabel('% Variazione Potenza (ERD/ERS)');
    legend; grid on;
    
    disp('Interpretazione:');
    disp('Se le linee si separano (una su, una giù), il soggetto controlla.');
    disp('Se vanno insieme, non lateralizza.');
end