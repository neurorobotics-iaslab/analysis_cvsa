%% --- 1. CARICAMENTO DATI ---
% Utilizzo il percorso fornito nel tuo file check_psd.m
full_path = '/home/paolo/bci_vr_ws/recordings/c7/20260422/calibration/gdf/c7.20260422.131344.calibration.mi_lhrh.gdf';
[s, h] = sload(full_path);

setting.samplerate = h.SampleRate;
setting.labels = h.Label;
setting.elec = h.ELEC; 

errp_comput_f(s, setting)

function [] = errp_comput_f(T, setting)
    % --- Parametri ---
    fs = setting.samplerate;
    f_max = 40;
    labels = setting.labels;
    num_channels = size(T, 2);
    T = double(T); % Conversione necessaria per i calcoli

    % --- PRE-PROCESSING (Butterworth 2-40 Hz) ---
    [b, a] = butter(4, [2 40]/(fs/2), 'bandpass');
    T_filt = filtfilt(b, a, T); 

    % --- CALCOLO PSD (Welch) ---
    win_len = 2 * fs; 
    noverlap = win_len * 0.5;
    [pxx, f] = pwelch(T_filt, hanning(win_len), noverlap, win_len, fs);
    
    idx = (f >= 1) & (f <= f_max);
    f_p = f(idx);
    psd_db = 10 * log10(pxx(idx, :)); % [Frequenze x Canali]

    % --- PLOT PSD OVERLAY ---
    figure('Color', 'w', 'Name', 'PSD Overlay (1-40 Hz)');
    hold on;
    colors = lines(num_channels); 
    for ch = 1:num_channels
        plot(f_p, psd_db(:, ch), 'Color', [colors(ch,:), 0.5]);
    end
    plot(f_p, mean(psd_db, 2), 'k', 'LineWidth', 2.5, 'DisplayName', 'Media');
    grid on; xlabel('Frequenza (Hz)'); ylabel('Sottopotenza (dB/Hz)');
    title('PSD - Matching EEGLAB Style');

    % --- SEZIONE TOPOPLOT ---
    % Creazione struttura chanlocs
    chanlocs = struct('labels', labels);
    if isfield(setting.elec, 'Pos')
        for c = 1:num_channels
            chanlocs(c).X = setting.elec.Pos(c,1);
            chanlocs(c).Y = setting.elec.Pos(c,2);
            chanlocs(c).Z = setting.elec.Pos(c,3);
        end
        % Conversione coordinate (risolve l'errore plotrad)
        chanlocs = pop_chanedit(chanlocs, 'convert', 'chancenter', [], 'optimize', 'on');
    else
        % Fallback su locazioni standard se Pos è vuoto
        chanlocs = pop_chanedit(chanlocs, 'lookup', 'standard_1005.elc');
    end

    % Frequenze target per i topoplot
    target_freqs = [6 10 22];
    figure('Color', 'w', 'Name', 'Topographic Maps (Power Distribution)');
    
    for i = 1:length(target_freqs)
        subplot(1, 3, i);
        
        % Trova l'indice della frequenza target
        [~, f_idx] = min(abs(f_p - target_freqs(i)));
        topodata = psd_db(f_idx, :); 
        
        % Topoplot senza parametri sconosciuti
        % 'maplimits', 'absmax' aiuta a rendere i colori comparabili
        topoplot(topodata, chanlocs, 'style', 'both', 'electrodes', 'on');
        title([num2str(target_freqs(i)) ' Hz']);
        colorbar;
    end
end