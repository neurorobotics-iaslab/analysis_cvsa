%% --- 1. SETUP AMBIENTE ---
clear; clc;
eeglab_repo = '/home/paolo/Local/Matlab/eeglab'; 

if ~exist('eeglab', 'file'), addpath(eeglab_repo); end
[ALLEEG, EEG, CURRENTSET] = eeglab;

%% --- 2. SELEZIONE E CARICAMENTO FILE ---
[files, folder] = uigetfile('*.gdf', 'Seleziona uno o più file GDF', 'MultiSelect', 'on');
if isequal(files,0), disp('Operazione annullata'); return; end
if ischar(files), files = {files}; end 

ALLEEG = []; 
all_indices = [];

for f = 1:length(files)
    fprintf('Caricamento file %d/%d: %s\n', f, length(files), files{f});
    full_path = fullfile(folder, files{f});
    
    % Caricamento dati
    [s, h] = sload(full_path);
    
    % --- RIMOZIONE CANALI TESTA ---
    % Supponiamo che s abbia dimensioni [campioni x canali]
    % Escludiamo gli ultimi 3 canali
    total_chans = size(s, 2);
    eeg_chans_idx = 1:(total_chans - 3);
    s_eeg = s(:, eeg_chans_idx);
    labels_eeg = h.Label(eeg_chans_idx);
    
    % Creazione struttura EEG
    EEG = pop_importdata('dataformat', 'array', 'nbchan', length(eeg_chans_idx), ...
        'data', s_eeg', 'srate', h.SampleRate, 'setname', files{f});
    
    % Coordinate elettrodi
    EEG.chanlocs = struct('labels', labels_eeg);
    EEG = pop_chanedit(EEG, 'lookup', 'standard_1005.elc');
    
    % Importazione Eventi
    if isfield(h, 'EVENT')
        for i = 1:length(h.EVENT.POS)
            EEG.event(i).type = num2str(h.EVENT.TYP(i));
            EEG.event(i).latency = h.EVENT.POS(i);
            EEG.event(i).duration = h.EVENT.DUR(i);
        end
    end
    
    EEG = eeg_checkset(EEG);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
    all_indices = [all_indices, CURRENTSET];
end

%% --- 3. MERGE (Dataset: "merge") ---
if length(all_indices) > 1
    EEG = pop_mergeset(ALLEEG, all_indices);
else
    EEG = ALLEEG(1);
end

EEG.setname = 'merge';
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
fprintf('Dataset "merge" creato con %d canali (head position rimossa).\n', EEG.nbchan);

%% --- 4. FILTRAGGIO (Dataset: "merge_1_40") ---
% Creiamo una copia per il filtraggio partendo dal dataset 'merge'
EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'hicutoff', 40);
EEG.setname = 'merge_1_40';

% Memorizziamo come nuovo dataset
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
eeglab redraw;

%% --- 5. VISUALIZZAZIONE ---
fprintf('Visualizzazione del dataset: %s\n', EEG.setname);
pop_eegplot(EEG, 1, 1, 1);

%% --- 6. ANALISI PSD (Power Spectral Density) ---
% Selezioniamo il dataset filtrato (l'ultimo creato)
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

% Definiamo i parametri della PSD
freq_range = [2 30]; % Range di frequenze da visualizzare

figure('Name', 'Analisi PSD - Ricerca Ritmo Mu (10Hz)', 'Color', 'w');

% pop_spectopo calcola la PSD per tutti i canali
% 'percent', 15: usa il 15% dei dati per velocizzare (puoi mettere 100 per precisione massima)
% 'freqrange': limita l'asse X del grafico
% 'electrodes', 'on': mostra i nomi dei canali nel grafico
pop_spectopo(EEG, 1, [0  EEG.pnts], 'EEG' , ...
    'percent', 20, ...
    'freqrange', freq_range, ...
    'electrodes', 'on');

% Istruzioni per l'interpretazione
fprintf('\n--- ANALISI PSD ---\n');
fprintf('Cerca la "collinetta" (peak) tra 8-13 Hz.\n');
fprintf('Controlla specialmente i canali C3, Cz, C4 (motori).\n');

%% --- 7. EPOCHING ---
% Logica di selezione dei trigger in base al nome del file
if any(contains(files{1}, 'cvsa_blbr'))
    event_types = {'730', '731'};
    fprintf('Pattern "cvsa_blbr" rilevato. Uso trigger: 730, 731\n');
elseif any(contains(files{1}, 'mi_lhrh'))
    event_types = {'769', '770'};
    fprintf('Pattern "mi_lhrh" rilevato. Uso trigger: 769, 770\n');
elseif any(contains(files{1}, 'hybrid'))
    event_types = {'750', '751'};
    fprintf('Pattern "hybrid" rilevato. Uso trigger: 750, 751\n');
else
    event_types = {'769', '770'}; % Default
    warning('Nessun pattern riconosciuto. Uso i trigger di default (769, 770).');
end

% Definiamo i limiti temporali dell'epoca (in secondi)
epoch_limits = [-2  5]; 

% Estraiamo le epoche basandoci sui trigger selezionati
EEG_ep = pop_epoch(EEG, event_types, epoch_limits, 'newname', [EEG.setname '_epoched'], 'epochinfo', 'yes');

% Rimozione della Baseline
% Sottraiamo la media del periodo di fixation (-2000ms a 0ms)
EEG_ep = pop_rmbase(EEG_ep, [-2000 0]);

% Salviamo il dataset in EEGLAB
[ALLEEG, EEG_ep, CURRENTSET] = eeg_store(ALLEEG, EEG_ep, 0);
eeglab redraw;

%% --- 8. ANALISI IN FREQUENZA ERD/ERS singoli canali ---
target_labels = {'C3', 'C4', 'Cz', 'FC1', 'FC5', 'FC2', 'FC6', 'CP1', 'CP5', 'CP2', 'CP6', ...
    'O1', 'O2', 'Oz', 'P8', 'P7'};
num_chans = length(target_labels);

% Identificazione indici epoche (Robust version)
idx_left = []; idx_right = [];
for i = 1:length(EEG_ep.epoch)
    ev_latencies = [EEG_ep.epoch(i).eventlatency{:}];
    zero_idx = find(ev_latencies == 0, 1);
    if ~isempty(zero_idx)
        t = EEG_ep.epoch(i).eventtype;
        if iscell(t), t = t{zero_idx}; end
        if isequal(t, event_types{1}) || isequal(t, str2num(event_types{1})), idx_left = [idx_left i];
        elseif isequal(t, event_types{2}) || isequal(t, str2num(event_types{2})), idx_right = [idx_right i]; end
    end
end

class_idx = {idx_left, idx_right};
class_lbl = {event_types{1}, event_types{2}};

for c = 1:2
    h_fig = figure('Name', ['ERSP Grid: ' class_lbl{c}], 'Color', 'w');
    set(h_fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    
    % Definiamo la griglia del subplot
    rows = 4; cols = 4; 
    
    for ch_i = 1:num_chans
        label = target_labels{ch_i};
        idx_in_eeg = find(strcmp({EEG_ep.chanlocs.labels}, label));
        
        if isempty(idx_in_eeg), continue; end
        
        subplot(rows, cols, ch_i);
        
        % Chiamata a newtimef: 'plotitc','off' e 'plotersp','off' per gestire noi il plot
        % Usiamo 'nwin', 50 per forzare una dimensione della finestra compatibile
        [ersp, itc, powbase, times, freqs] = newtimef(EEG_ep.data(idx_in_eeg, :, class_idx{c}), ...
            EEG_ep.pnts, [EEG_ep.times(1) EEG_ep.times(end)], EEG_ep.srate, 0, ...
            'baseline', [-1500 -500], 'freqs', [6 35], 'winsize', 256, ...
            'plotersp', 'off', 'plotitc', 'off', 'verbose', 'off');
        
        % Plot manuale della heatmap per avere pieno controllo
        imagesc(times, freqs, ersp);
        set(gca, 'YDir', 'normal');
        colormap(jet);
        clim([-4 4]); % Scala dB fissa per confronto
        
        title(label, 'FontSize', 10);
        if mod(ch_i-1, cols) == 0, ylabel('Hz'); end
        if ch_i > (num_chans - cols), xlabel('ms'); end
    end
    
    % Aggiungiamo una colorbar comune a destra
    cb = colorbar('Position', [0.93 0.11 0.02 0.81]);
    ylabel(cb, 'ERSP (dB)');
end


%% --- 9. DIFFERENZA DI POTENZA NEL TEMPO (HEATMAP L-R) ---
% 1. Filtriamo nella banda Alpha/Mu (8-14 Hz)
low_freq = 8;
high_freq = 14;
% low_freq = input('low freq: ');
% high_freq = input('high freq: ');
EEG_mu = pop_eegfiltnew(EEG_ep, low_freq, high_freq);

% 2. Identificazione epoche
idx_left = []; idx_right = [];
for i = 1:length(EEG_mu.epoch)
    if any(strcmp(EEG_mu.epoch(i).eventtype, event_types{1})), idx_left = [idx_left i];
    elseif any(strcmp(EEG_mu.epoch(i).eventtype, event_types{2})), idx_right = [idx_right i]; end
end

% 3. Calcolo Potenza (Hilbert)
hilb_data = abs(hilbert(reshape(EEG_mu.data, EEG_mu.nbchan, [])'))';
pow_data_abs = reshape(hilb_data.^2, EEG_mu.nbchan, EEG_mu.pnts, []);

% Trasformata Logaritmica
pow_data_log = 10 * log10(pow_data_abs + eps); 

% --- NOVITÀ: NORMALIZZAZIONE BASELINE TRIAL-BY-TRIAL (Come in ERSP) ---
% Troviamo gli indici del periodo di fixation (es. da -1500 a -500 ms)
base_idx = find(EEG_mu.times >= -1500 & EEG_mu.times <= -500);

% Calcoliamo la media della baseline per ogni canale e per ogni trial
base_pow = mean(pow_data_log(:, base_idx, :), 2);

% Sottraiamo la baseline (in logaritmo equivale al calcolo dell'ERD%)
pow_data = pow_data_log - base_pow; 

% Calcolo medie per il plot
mean_all_left  = mean(pow_data(:, :, idx_left), 3);
mean_all_right = mean(pow_data(:, :, idx_right), 3);

% 4. Calcolo Differenza (LEFT - RIGHT)
mean_diff = mean_all_left - mean_all_right;

% 5. PLOT UNICO (Heatmap)
figure('Name', sprintf('Differenza ERD (dB) [%d-%d Hz]: Left - Right', low_freq, high_freq), 'Color', 'w');

lim = max(abs(mean_diff(:))); 
imagesc(EEG_mu.times, 1:EEG_mu.nbchan, mean_diff, [-lim lim]);
colormap('jet'); 
h = colorbar;
ylabel(h, '\Delta ERD (dB)'); 

set(gca, 'YTick', 1:EEG_mu.nbchan, 'YTickLabel', {EEG_mu.chanlocs.labels}, 'FontSize', 9);
xlabel('Tempo (ms)');
ylabel('Elettrodi');
title(sprintf('Evoluzione Differenza (L-R) normalizzata | Banda: %d-%d Hz', low_freq, high_freq));

line([0 0], [0.5 EEG_mu.nbchan+0.5], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2); 
line([1000 1000], [0.5 EEG_mu.nbchan+0.5], 'Color', [0 1 0], 'LineStyle', ':', 'LineWidth', 2); 
text(0, -0.5, 'CUE', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(1000, -0.5, 'CF START', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', [0 0.5 0]);
grid on;


%% --- 10. MATRICE 3x7: DIFF, R^2 E FISHER SCORE ---
intervals = [0, 5; 1, 5; 0, 1; 1, 2; 2, 3; 3, 4; 4, 5];
titles_cols = {'Cue + CF (0-5s)', 'Solo CF (1-5s)', '0-1s', '1-2s', '2-3s', '3-4s', '4-5s'};
num_intervals = size(intervals, 1);

all_diff_z = cell(1, num_intervals);
all_r2 = cell(1, num_intervals);
all_fisher = cell(1, num_intervals);
max_r2_val = 0; max_f_val = 0;

labels_vec = [ones(1, length(idx_left)), zeros(1, length(idx_right))];

for i = 1:num_intervals
    t_idx = find(EEG_mu.times >= intervals(i,1)*1000 & EEG_mu.times < intervals(i,2)*1000);
    pow_trial_avg = squeeze(mean(pow_data(:, t_idx, :), 2));
    
    % 1. Diff Z-Score
    m1 = mean(pow_trial_avg(:, idx_left), 2);   v1 = var(pow_trial_avg(:, idx_left), 0, 2);
    m2 = mean(pow_trial_avg(:, idx_right), 2);  v2 = var(pow_trial_avg(:, idx_right), 0, 2);
    std_tot = std(pow_trial_avg, 0, 2);
    all_diff_z{i} = (m1 - m2) ./ std_tot;
    
    % 2. R^2
    r2_v = zeros(EEG_mu.nbchan, 1);
    for ch = 1:EEG_mu.nbchan
        r = corrcoef(pow_trial_avg(ch, [idx_left, idx_right]), labels_vec);
        if size(r,1) > 1, r2_v(ch) = r(1,2)^2; end
    end
    all_r2{i} = r2_v;
    max_r2_val = max(max_r2_val, max(r2_v));
    
    % 3. Fisher Score: (m1-m2)^2 / (v1+v2)
    f_score = ((m1 - m2).^2) ./ (v1 + v2);
    all_fisher{i} = f_score;
    max_f_val = max(max_f_val, max(f_score));
end

% PLOT
h_fig = figure('Name', 'Analisi BCI 3x7', 'Color', 'w');
set(h_fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
tlo = tiledlayout(3, num_intervals, 'TileSpacing', 'compact', 'Padding', 'tight');

% Disattiva warning per evitare spam in console durante topoplot
warning('off', 'all');

for r = 1:3
    for c = 1:num_intervals
        nexttile;
        if r == 1 % Riga Diff
            data = all_diff_z{c}; lim = [-max(abs(data)) max(abs(data))]; label = 'DIFF (Z)';
        elseif r == 2 % Riga R2
            data = all_r2{c}; lim = [0 max_r2_val]; label = 'R^2';
        else % Riga Fisher
            data = all_fisher{c}; lim = [0 max_f_val]; label = 'FISHER';
        end
        
        topoplot(data, EEG_mu.chanlocs, 'style', 'both', 'maplimits', lim, 'electrodes', 'labels', 'whitebk', 'on');
        
        if r == 1, title(titles_cols{c}); end
        if c == 1, ylabel(label, 'FontSize', 12, 'FontWeight', 'bold', 'Visible', 'on'); end
        colorbar;
    end
end

%warning('on', 'all'); % Riattiva i warning

%% --- 11. COMMON SPATIAL PATTERN (CSP) ---
fprintf('Calcolo Common Spatial Pattern (CSP)...\n');

% Usiamo i dati EEG filtrati nella banda Mu (8-14 Hz) nell'intervallo di task
t_fb = find(EEG_mu.times >= 1000 & EEG_mu.times <= 4000);

% Inizializziamo le matrici di covarianza
cov_left = zeros(EEG_mu.nbchan, EEG_mu.nbchan);
cov_right = zeros(EEG_mu.nbchan, EEG_mu.nbchan);

% 1. Calcolo covarianza media per la classe Sinistra
for i = 1:length(idx_left)
    tr_data = EEG_mu.data(:, t_fb, idx_left(i));
    tr_data = tr_data - mean(tr_data, 2); % Centratura dati
    c = tr_data * tr_data';
    cov_left = cov_left + (c / trace(c)); % Normalizzazione traccia
end
cov_left = cov_left / length(idx_left);

% 2. Calcolo covarianza media per la classe Destra
for i = 1:length(idx_right)
    tr_data = EEG_mu.data(:, t_fb, idx_right(i));
    tr_data = tr_data - mean(tr_data, 2);
    c = tr_data * tr_data';
    cov_right = cov_right + (c / trace(c));
end
cov_right = cov_right / length(idx_right);

% 3. Problema generalizzato degli autovalori (Equazione centrale del CSP)
[V, D] = eig(cov_left, cov_left + cov_right);

% Ordiniamo gli autovalori in modo decrescente
[~, sort_idx] = sort(diag(D), 'descend');
W = V(:, sort_idx); % Matrice dei Filtri CSP

% 4. Calcolo Spatial Patterns (come il cervello genera i segnali)
% I pattern spaziali si ottengono invertendo la matrice dei filtri
CSP_patterns = inv(W); 

% 5. PLOT Topografico dei Pattern CSP
figure('Name', 'Common Spatial Patterns (CSP)', 'Color', 'w');
set(gcf, 'Units', 'normalized', 'OuterPosition', [0.1 0.2 0.8 0.4]);
sgtitle('Filtri Spaziali Ottimi (CSP) - 8-14 Hz', 'FontSize', 14, 'FontWeight', 'bold');

% Le prime due componenti massimizzano la Sinistra, le ultime due la Destra
comps_to_plot = [1, 2, EEG_mu.nbchan-1, EEG_mu.nbchan];
titles_csp = {'Comp 1 (Max Left)', 'Comp 2 (Max Left)', 'Comp N-1 (Max Right)', 'Comp N (Max Right)'};

for i = 1:4
    subplot(1, 4, i);
    % Plottiamo la riga corrispondente della matrice dei pattern
    topoplot(CSP_patterns(comps_to_plot(i), :), EEG_mu.chanlocs, 'style', 'both');
    title(titles_csp{i});
end
colormap(jet);
