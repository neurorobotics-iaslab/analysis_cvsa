%% --- 1. SETUP AMBIENTE ---
clear; clc; close all;
eeglab_repo = '/home/paolo/Local/Matlab/eeglab'; 

if ~exist('eeglab', 'file'), addpath(eeglab_repo); end
[ALLEEG, EEG, CURRENTSET] = eeglab('nogui');

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
    total_chans = size(s, 2);
    eeg_chans_idx = 1:(total_chans - 3);
%     eeg_chans_idx = 1:total_chans-1;
    s_eeg = s(:, eeg_chans_idx);
    labels_eeg = h.Label(eeg_chans_idx);
%     labels_eeg = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C3', 'C1', 'Cz', 'C2', 'C4', 'Fp1', 'CP1', 'CPz', 'CP2', 'Fp2'}';
%     labels_eeg = {'P5', 'P3', 'P1', 'P2', 'P4', 'P6', 'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'FP1', 'O1', 'Oz', 'O2', 'FP2'};
    
    % Creazione struttura EEG
    EEG = pop_importdata('dataformat', 'array', 'nbchan', length(eeg_chans_idx), ...
        'data', s_eeg', 'srate', h.SampleRate, 'setname', files{f});
    
    EEG.chanlocs = struct('labels', labels_eeg);
    EEG = pop_chanedit(EEG, 'lookup', 'standard_1005.elc');
    
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

%% --- 4. downsampling and filtering (Dataset: "merge_1_40") ---
fprintf('\n--- DOWNSAMPLING A 128 HZ ---\n');
EEG = pop_resample(EEG, 128);
EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'hicutoff', 40);
EEG.setname = 'merge_1_40';

[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
eeglab redraw;

%% --- 5. VISUALIZZAZIONE ---
fprintf('Visualization dataset: %s\n', EEG.setname);
pop_eegplot(EEG, 1, 1, 1);

%% --- 6. ANALISI PSD PRE-PROCESAMENTO (Dati Continui) ---
fprintf('\n--- STARTING PSD ---\n');
freq_range = [2 30]; 
figure('Name', 'PSD 1: Data Pre-ICA and CAR', 'Color', 'w');
pop_spectopo(EEG, 1, [0  EEG.pnts], 'EEG' ,  'freq', [8 12 14], 'freqrange', freq_range, 'electrodes', 'on');

%% --- 6.5 ICA ---
fprintf('\n--- COMPUTING ICA ON ALL CONTINUOUS DATA ---\n');
EEG = pop_runica(EEG, 'icatype', 'runica');
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% --- 7. EPOCHING ---
if any([contains(files{1}, 'cvsa_blbr'), contains(files{1}, 'cvsa_lbrb')])
    event_types = {'730', '731'};
    fprintf('Pattern "cvsa_blbr" rilevato. Trigger: 730, 731\n');
elseif any(contains(files{1}, 'mi_lhrh'))
    event_types = {'769', '770'};
    fprintf('Pattern "mi_lhrh" rilevato. Trigger: 769, 770\n');
elseif any(contains(files{1}, 'hybrid'))
    event_types = {'750', '751'};
    fprintf('Pattern "hybrid" rilevato. Trigger: 750, 751\n');
else
    event_types = {'771', '773'}; 
    warning('Pattern non riconosciuto. Uso default: 771, 773.');
end

epoch_limits = [-2  5]; 
EEG_ep = pop_epoch(EEG, event_types, epoch_limits, 'newname', [EEG.setname '_epoched'], 'epochinfo', 'yes');
EEG_ep = pop_rmbase(EEG_ep, [-2000 0]);

%% --- 8. ISPEZIONE E RIMOZIONE ICA INTERATTIVA (Sulle Epoche) ---
fprintf('\n--- APPLY ICA ---\n');
EEG_ep = iclabel(EEG_ep);

thresh = [NaN NaN; 0.75 1; 0.75 1; 0.75 1; NaN NaN; 0.75 1; 0.75 1];
EEG_ep = pop_icflag(EEG_ep, thresh);

bad_comps = find(EEG_ep.reject.gcompreject);
if ~isempty(bad_comps)
    fprintf('ICLabel SUGGEST remotion of %d components: [%s]\n', length(bad_comps), num2str(bad_comps));
end

assignin('base', 'EEG', EEG_ep);
[ALLEEG, EEG_ep, CURRENTSET] = eeg_store(ALLEEG, EEG_ep, 0);

% Apre l'interfaccia
pop_selectcomps(EEG_ep, 1:size(EEG_ep.icaweights,1));

h_gui = gcf; 

if isgraphics(h_gui, 'figure')
    uiwait(h_gui);
end

EEG_ep = evalin('base', 'EEG');
final_bad_comps = find(EEG_ep.reject.gcompreject);

if ~isempty(final_bad_comps)
    fprintf('\nRemotion components: [%s]...\n', num2str(final_bad_comps));
    EEG_ep = pop_subcomp(EEG_ep, final_bad_comps, 0);
    
    assignin('base', 'EEG', EEG_ep);
    [ALLEEG, EEG_ep, CURRENTSET] = eeg_store(ALLEEG, EEG_ep, CURRENTSET);
    
    fprintf('Remotion ended.\n');
else
    fprintf('\nNothing removed.\n');
end

%% --- 8.5 ANALISI PSD POST-ICA ---
figure('Name', 'PSD 2: After ICA', 'Color', 'w');
pop_spectopo(EEG_ep, 1, [EEG_ep.times(1) EEG_ep.times(end)], 'EEG' ,  'freq', [8 12 14], 'freqrange', freq_range, 'electrodes', 'on');

%% --- 9. APPLICAZIONE CAR (Common Average Reference) ---
fprintf('\n--- APPLY CAR ---\n');
EEG_ep = pop_reref(EEG_ep, []);
[ALLEEG, EEG_ep, CURRENTSET] = eeg_store(ALLEEG, EEG_ep, CURRENTSET);

%% --- 9.5 ANALISI PSD POST-CAR ---
figure('Name', 'PSD 3: After ICA and CAR', 'Color', 'w');
pop_spectopo(EEG_ep, 1, [EEG_ep.times(1) EEG_ep.times(end)], 'EEG' , 'freq', [8 12 14], 'freqrange', freq_range, 'electrodes', 'on');

%% --- 10. SALVATAGGIO MAPPE ERD/ERS SINGOLI CANALI ---
fprintf('\n--- ERD/ERS MAPS ---\n');
output_dir = fullfile(folder, 'results_eeglab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

target_labels = {EEG_ep.chanlocs.labels}; 
num_chans = length(target_labels);

idx_left = []; idx_right = [];
for i = 1:length(EEG_ep.epoch)
    ev_latencies = [EEG_ep.epoch(i).eventlatency{:}];
    zero_idx = find(ev_latencies == 0, 1);
    if ~isempty(zero_idx)
        t = EEG_ep.epoch(i).eventtype;
        if iscell(t), t = t{zero_idx}; end
        if any([isequal(t, event_types{1}), isequal(t, str2num(event_types{1}))]), idx_left = [idx_left i];
        elseif any([isequal(t, event_types{2}), isequal(t, str2num(event_types{2}))]), idx_right = [idx_right i]; end
    else
        disp('empty epoch');
    end
end

class_idx = {idx_left, idx_right};
class_lbl = {event_types{1}, event_types{2}};

fprintf('Starting save maps ERD/ERS per channel...\n');

for ch_i = 1:num_chans
    label = target_labels{ch_i};
    h_fig = figure('Name', ['ERSP: ' label], 'Color', 'w', 'Visible', 'off');

    for c = 1:length(class_lbl)
        
        idx_in_eeg = find(strcmp({EEG_ep.chanlocs.labels}, label));
        if isempty(idx_in_eeg)
            warning('Channel %s not founded, skip saving.', label);
            continue; 
        end
        
        % Chiamata a newtimef
        [ersp, itc, powbase, times, freqs]  = newtimef(EEG_ep.data(idx_in_eeg, :, class_idx{c}), ...
            EEG_ep.pnts, [EEG_ep.times(1) EEG_ep.times(end)], EEG_ep.srate, 0, ...
            'baseline', 0, 'freqs', [4 30], 'winsize', 128, ...
            'plotersp', 'off', 'plotitc', 'off', 'verbose', 'off');
        
        % Plot della heatmap
        subplot(1,2,c)
        imagesc(times, freqs, ersp);
        set(gca, 'YDir', 'normal');
        colormap(jet);
        clim([-4 4]);
        title(sprintf(' Classe %s', class_lbl{c}), 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Frequenza (Hz)');
        xlabel('Tempo (ms)');
        cb = colorbar;
        ylabel(cb, 'ERSP (dB)');
        hold on;
        line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
        hold off;
    end

    sgtitle(sprintf('ERSP Canale %s', label), 'FontSize', 12, 'FontWeight', 'bold')
    
    % Save
    save_name = sprintf('erd_%s_%s_%s.png', label, class_lbl{1}, class_lbl{2});
    save_path = fullfile(output_dir, save_name);
    saveas(h_fig, save_path);
    close(h_fig);

    fprintf('Saved: %s\n', save_path);
end

fprintf('Saved everything on: %s\n', output_dir);


%% --- 10. COMPUTE ERSP TOPOPLOT ---
fprintf('\n--- COMPUTE ERD/ERS TOPOPLOT ---\n');

bands = [4 8; 8 14; 14 24; 14 30];
band_names = arrayfun(@(i) sprintf('%d_%d', bands(i,1), bands(i,2)), 1:size(bands,1), 'UniformOutput', false);
intervals = [0, 5; 1, 5; 0, 1; 1, 2; 2, 3; 3, 4; 4, 5];
titles_cols = {'Cue + CF (0-5s)', 'Solo CF (1-5s)', '0-1s', '1-2s', '2-3s', '3-4s', '4-5s'};
num_intervals = size(intervals, 1);

all_p_db = nan(EEG.nbchan, length(class_lbl), length(freqs), length(times));

for ch = 1:EEG_ep.nbchan
    for idx_cl = 1:length(class_lbl)
        [ersp, itc, powbase, times, freqs] = newtimef(EEG_ep.data(ch, :, class_idx{idx_cl}), ...
            EEG_ep.pnts, [EEG_ep.times(1) EEG_ep.times(end)], EEG_ep.srate, 0, ...
            'baseline', 0, 'freqs', [4 30], 'winsize', 128, ...
            'plotersp', 'off', 'plotitc', 'off', 'verbose', 'off');

        P_percentual = (10.^(ersp/10)-1)*100;% db: 10 * log10(abs(ersp).^2 + eps);

        all_p_db(ch, idx_cl,:,:) = P_percentual;
    end
end

for idx_b = 1:size(bands, 1)
    low_f = bands(idx_b, 1); high_f = bands(idx_b, 2);
    freq_idx = find(freqs >= low_f & freqs <= high_f);
    
    pow_band_1 = squeeze(mean(all_p_db(:,1,freq_idx,:), 3));
    pow_band_2 = squeeze(mean(all_p_db(:,2,freq_idx,:), 3));

    h_topos = figure('Name', sprintf('Topoplots %d-%d Hz', low_f, high_f), 'Color', 'w', 'Visible','off');
    set(h_topos, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    tlo = tiledlayout(2, num_intervals, 'TileSpacing', 'compact', 'Padding', 'tight');

    all_m1 = cell(1, num_intervals);
    all_m2 = cell(1, num_intervals);
    max_p = 0;

    for i = 1:num_intervals
        t_idx = find(times >= intervals(i,1)*1000 & times < intervals(i,2)*1000);
        
        all_m1{i} = squeeze(mean(pow_band_1(:, t_idx), 2));
        all_m2{i} = squeeze(mean(pow_band_2(:, t_idx), 2));
        
        max_p = max([max_p; abs(all_m1{i}); abs(all_m2{i})]);
    end

    % Plotting
    for r = 1:2
        for c = 1:num_intervals
            nexttile;
            if r == 1, data = all_m1{c}; lim = [-20 0]; lbl = ['Cl ' class_lbl{1}];
            elseif r == 2, data = all_m2{c}; lim = [-20 0]; lbl = ['Cl ' class_lbl{2}]; 
            end
            
            topoplot(data, EEG_ep.chanlocs, 'style', 'both', 'maplimits', lim, 'electrodes', 'labels', 'whitebk', 'on');
            if r == 1, title(titles_cols{c}); end
            if c == 1, ylabel(lbl, 'Visible', 'on', 'FontWeight', 'bold'); end
            if c == num_intervals, colorbar; end
        end
    end
    
    % Salvataggio
    saveas(h_topos, fullfile(output_dir, sprintf('topo_%s__%s_%s.png', band_names{idx_b}, class_lbl{1}, class_lbl{2})));
    close(h_topos);
end


%% --- 11. TIME-COURSE ERD/ERS ALL CHANNELS ---
fprintf('\n--- ERD/ERS IN TIME ---\n');

mu_band = [8 14]; 
freq_idx = find(freqs >= mu_band(1) & freqs <= mu_band(2));

pow_tc_cls1 = squeeze(mean(all_p_db(:, 1, freq_idx, :), 3));
pow_tc_cls2 = squeeze(mean(all_p_db(:, 2, freq_idx, :), 3));

% SMOOTHING
win_size = round((1000 / (times(2) - times(1))) * 0.5); 
pow_tc_cls1_sm = movmean(pow_tc_cls1, win_size, 2);
pow_tc_cls2_sm = movmean(pow_tc_cls2, win_size, 2);

h_tc_all = figure('Name', sprintf('Time-Course ERD/ERS (%d-%d Hz)', mu_band(1), mu_band(2)), 'Color', 'w', 'Visible', 'off');
set(h_tc_all, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

num_cols = ceil(sqrt(EEG_ep.nbchan));
num_rows = ceil(EEG_ep.nbchan / num_cols);
tlo = tiledlayout(num_rows, num_cols, 'TileSpacing', 'compact', 'Padding', 'tight');

y_lims = [min(min([pow_tc_cls1_sm, pow_tc_cls2_sm])) - 0.5, max(max([pow_tc_cls1_sm, pow_tc_cls2_sm])) + 0.5];

for ch = 1:EEG_ep.nbchan
    nexttile;
    plot(times, pow_tc_cls1_sm(ch, :), 'b', 'LineWidth', 2); hold on;
    plot(times, pow_tc_cls2_sm(ch, :), 'r', 'LineWidth', 2);
    
    line([0 0], y_lims, 'Color', 'k', 'LineStyle', '--');
    title(EEG_ep.chanlocs(ch).labels, 'FontSize', 10, 'FontWeight', 'bold');
    ylim(y_lims); xlim([times(1) times(end)]); grid on;
    
    if mod(ch-1, num_cols) == 0, ylabel('dB'); else, yticklabels(''); end
    if ch > EEG_ep.nbchan - num_cols, xlabel('ms'); else, xticklabels(''); end
    if ch == 1, legend(['Class ' class_lbl{1}], ['Class ' class_lbl{2}], 'Location', 'best'); end
end

save_path_tc = fullfile(output_dir, sprintf('TimeCourse_AllChans_%d_%dHz.png', mu_band(1), mu_band(2)));
saveas(h_tc_all, save_path_tc);
fprintf('Graph saved in: %s\n', save_path_tc);


%% --- 12. ERSP DIFF ---
fprintf('\n--- ERSP DIFF ---\n');

num_bands = size(bands, 1);
num_interv = size(intervals, 1);
num_chans = EEG_ep.nbchan;

ERSP_mean_Cls = zeros(length(class_lbl), num_chans, num_bands, num_interv);

for ch = 1:num_chans
    for c = 1:2

        [ersp, ~, ~, times_tf, freqs_tf, ~, ~, ~] = newtimef(EEG_ep.data(ch, :, class_idx{c}), ...
            EEG_ep.pnts, [EEG_ep.times(1) EEG_ep.times(end)], EEG_ep.srate, 0, ...
            'baseline', 0, 'freqs', [4 30], 'winsize', 128, 'plotersp', 'off', 'plotitc', 'off', 'verbose', 'off');
        
        %ersp = ersp <= 0;
        P_percent = (10.^(ersp/10)-1)*100; 
        
        for idx_b = 1:num_bands
            f_idx = find(freqs_tf >= bands(idx_b,1) & freqs_tf <= bands(idx_b,2));
            for idx_t = 1:num_interv
                t_idx = find(times_tf >= intervals(idx_t,1)*1000 & times_tf < intervals(idx_t,2)*1000);
                
                val_mean = mean(mean(P_percent(f_idx, t_idx), 1), 2);
                
                ERSP_mean_Cls(c, ch, idx_b, idx_t) = val_mean;
            end
        end
    end
    if mod(ch,10)==0, fprintf('Processing ch %d...\n', ch); end
end

h_ersp_heat = figure('Name', 'Analisi Reattività ERSP %', 'Color', 'w', 'Visible','off');
set(h_ersp_heat, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
tlo = tiledlayout(2,ceil(num_interv/2),  'TileSpacing', 'compact');

handles = []; max_val = 0;
for i = 1:num_interv
    nexttile;
    imagesc(squeeze(abs(ERSP_mean_Cls(1,:,:,i) - ERSP_mean_Cls(2,:,:,i)))); 

    max_val = max([max_val, max(squeeze(abs(ERSP_mean_Cls(1,:,:,i) - ERSP_mean_Cls(2,:,:,i))))]);
    
    colormap(jet); colorbar;
    title(titles_cols{i});
    set(gca, 'XTick', 1:num_bands, 'XTickLabel', band_names);
    set(gca, 'YTick', 1:num_chans, 'YTickLabel', {EEG_ep.chanlocs.labels});
    handles = [handles, gca];
end
set(handles, 'clim', [0 max_val]);
sgtitle('ERSP ABS DIFF')

save_path_tc = fullfile(output_dir, sprintf('diff_ERSP_%s_%s.png', class_lbl{1}, class_lbl{2}));
saveas(h_ersp_heat, save_path_tc);
fprintf('Grafico salvato in: %s\n', save_path_tc);

%% --- 13. FISHER SCORE HEATMAP ---
fprintf('\n--- COMPUTING FISHER SCORE ---\n');

Fisher_Matrix = zeros(num_chans, num_bands, num_interv);

for ch = 1:num_chans
    [~, ~, ~, times_tf, freqs_tf, ~, ~, tfdata] = newtimef(EEG_ep.data(ch, :, [idx_left, idx_right]), ...
        EEG_ep.pnts, [EEG_ep.times(1) EEG_ep.times(end)], EEG_ep.srate, 0, ...
        'baseline', 0, 'freqs', [4 30], 'winsize', 128, 'plotersp', 'off', 'plotitc', 'off', 'verbose', 'off');
    
    P_trial = abs(tfdata).^2; 
    
    base_idx = find(times_tf >= -2000 & times_tf <= 0);
    for tr = 1:size(P_trial, 3)
        mean_base = mean(P_trial(:, base_idx, tr), 2);
        P_trial(:, :, tr) = ((P_trial(:, :, tr) ./ mean_base) - 1) * 100;
    end
    
    trials_cls1 = P_trial(:, :, 1:length(idx_left));
    trials_cls2 = P_trial(:, :, length(idx_left)+1:end);

    for idx_b = 1:num_bands
        f_idx = find(freqs_tf >= bands(idx_b,1) & freqs_tf <= bands(idx_b,2));
        for idx_t = 1:num_interv
            t_idx = find(times_tf >= intervals(idx_t,1)*1000 & times_tf < intervals(idx_t,2)*1000);
            
            m_trials_1 = squeeze(mean(mean(trials_cls1(f_idx, t_idx, :), 1), 2));
            m_trials_2 = squeeze(mean(mean(trials_cls2(f_idx, t_idx, :), 1), 2));
            
            mean1 = mean(m_trials_1); mean2 = mean(m_trials_2);
            var1 = var(m_trials_1);   var2 = var(m_trials_2);
            
            Fisher_Matrix(ch, idx_b, idx_t) = (mean1 - mean2)^2 / (var1 + var2 + eps);
        end
    end
    if mod(ch,10)==0, fprintf('Fisher Score: ch %d processed\n', ch); end
end

h_fisher = figure('Name', 'Fisher Score Analysis', 'Color', 'w', 'Visible','off');
set(h_fisher, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
tlo = tiledlayout(2, ceil(num_interv/2), 'TileSpacing', 'compact');

for i = 1:num_interv
    nexttile;
    imagesc(Fisher_Matrix(:, :, i)); 
    colormap(jet); colorbar; 
    title(titles_cols{i});
    set(gca, 'XTick', 1:num_bands, 'XTickLabel', band_names);
    set(gca, 'YTick', 1:num_chans, 'YTickLabel', {EEG_ep.chanlocs.labels});
end
sgtitle('Fisher Score Heatmap (Feature Importance)');
saveas(h_fisher, fullfile(output_dir, 'Fisher_Score_Heatmap.png'));

% --- TOPOPLOT FISHER SCORE ---
fprintf('\n--- GENERATING FISHER SCORE TOPOPLOTS ---\n');

for idx_b = 1:num_bands
    b_name = band_names{idx_b};
    
    h_topo_fish = figure('Name', sprintf('Fisher Topo: %s', b_name), 'Color', 'w', 'Visible', 'off');
    set(h_topo_fish, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    tlo_f = tiledlayout(1, num_interv, 'TileSpacing', 'compact', 'Padding', 'tight');
    
    max_fish_band = max(max(Fisher_Matrix(:, idx_b, :)));
    
    for idx_t = 1:num_interv
        nexttile;
        data_topo = Fisher_Matrix(:, idx_b, idx_t);
        
        topoplot(data_topo, EEG_ep.chanlocs, 'style', 'both', ...
            'maplimits', [0 max_fish_band + eps], 'electrodes', 'labels', 'whitebk', 'on');
        
        title(titles_cols{idx_t}, 'FontSize', 10);
        
        if idx_t == num_interv
            cb = colorbar;
            ylabel(cb, 'Fisher Score');
        end
    end
    
    sgtitle(sprintf('Fisher Score Topographic Map - Band: %s Hz', strrep(b_name, '_', '-')), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    save_topo_name = sprintf('Fisher_Score_topo_%s.png', b_name);
    saveas(h_topo_fish, fullfile(output_dir, save_topo_name));
    close(h_topo_fish);
    
    fprintf('Saved Fisher Topoplot for band %s: %s\n', b_name, save_topo_name);
end