clear all; close all; clc;

% --- Loading data ---
file_pattern = '/home/paolo/cvsa/ic_cvsa_ws/record_mi/results_SMC/results_*.mat';
results_files = dir(file_pattern);
if isempty(results_files), error('No files found.'); end

FullDatabase = [];
for i = 1:length(results_files)
    data = load([results_files(1).folder '/' results_files(i).name]);
    for j = 1:length(data.Database)
        data.Database(j).SubjectID = results_files(i).name(9:10);
    end
    FullDatabase = [FullDatabase; data.Database];
end

subjects = unique({FullDatabase.SubjectID});
n_sub = length(subjects);
metrics = {'Act_Acc_NoTo', 'Rest_Acc', 'Act_Num_Timeout'};
m_labels = {'Trial Accuracy (no timeout)', 'Rest Accuracy', 'Timeout Rate'};
metric_colors = [0.2, 0.4, 0.6; 0.4, 0.6, 0.2; 0.8, 0.6, 0.2];

figure('Color', 'w', 'Position', [50 50 1500 700]); hold on;

% Parametri spaziali
bw = 0.10; 
gap_methods = 0.04; 
gap_subs = 0.6;

for s = 1:n_sub
    sub_id = subjects{s};
    center_sub = (s-1) * (6*bw + gap_methods + gap_subs);
    
    data_gmm = FullDatabase(strcmpi({FullDatabase.SubjectID}, sub_id) & strcmpi({FullDatabase.Method}, 'gmm'));
    data_trad = FullDatabase(strcmpi({FullDatabase.SubjectID}, sub_id) & strcmpi({FullDatabase.Method}, 'traditional'));
    
    % Blocchi individuali (Our e Trad)
    renderBlock(center_sub, data_gmm, data_trad, bw, gap_methods, metrics, metric_colors, sub_id);
end

%% BLOCCO GRAND AVERAGE (AVG)
center_avg = n_sub * (6*bw + gap_methods + gap_subs);

all_gmm_data = FullDatabase(strcmpi({FullDatabase.Method}, 'gmm'));
all_trad_data = FullDatabase(strcmpi({FullDatabase.Method}, 'traditional'));

% Renderizza il blocco media con asterisco per Rest Accuracy (p=0.0020)
renderBlock(center_avg, all_gmm_data, all_trad_data, bw, gap_methods, metrics, metric_colors, 'AVG');

h = []; for m=1:3, h(m) = bar(nan, nan, 'FaceColor', metric_colors(m,:)); end
legend(h, m_labels, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 12);
ylabel('Percentage (%)', 'FontSize', 14, 'FontWeight', 'bold');
title('\bf Online MI-BCI Performance: Individual Subjects and Grand Average', 'FontSize', 16);
grid on; ylim([0 125]); set(gca, 'XTick', [], 'XColor', 'none'); box off;

%
fprintf('\n--- STATISTICAL ANALYSIS (Grand Average) ---\n');
for m = 1:3
    field = metrics{m};
    % Aggrega tutti i valori di tutti i soggetti per il metodo GMM e TRAD
    all_gmm = [FullDatabase(strcmpi({FullDatabase.Method}, 'gmm')).(field)];
    all_trad = [FullDatabase(strcmpi({FullDatabase.Method}, 'traditional')).(field)];
    
    % Se è la metrica dei timeout, converti in %
    if m == 3, all_gmm = all_gmm/20*100; all_trad = all_trad/20*100; end
    
    [p, h_stat] = signrank(all_gmm, all_trad);
    fprintf('%s: p-value = %.4f | Significant: %s\n', m_labels{m}, p, mat2str(h_stat));
end


%% TIMING ANALYSIS
% Metriche: Hit (Comando Corretto), Err (Falso Positivo in Rest), Miss (Comando Errato)
metrics_t = {'Act_Time_Hit', 'Rest_Time_Err', 'Act_Time_Miss'};
m_labels_t = {'Active: Time to Hit', 'Rest: Time to Miss', 'Active: Time to Miss'};

colors_t = [0.2, 0.4, 0.6; 0.4, 0.6, 0.2; 0.8, 0.6, 0.2];

figure('Color', 'w', 'Position', [50 50 1500 700]); hold on;
bw = 0.08; gap_m = 0.04; gap_s = 0.7;

for s = 1:n_sub
    sub_id = subjects{s};
    center_sub = (s-1) * (6*bw + gap_m + gap_s);
    
    d_gmm = FullDatabase(strcmpi({FullDatabase.SubjectID}, sub_id) & strcmpi({FullDatabase.Method}, 'gmm'));
    d_trad = FullDatabase(strcmpi({FullDatabase.SubjectID}, sub_id) & strcmpi({FullDatabase.Method}, 'traditional'));
    
    % Renderizziamo i blocchi Our e Trad
    renderTimingBlock(center_sub, d_gmm, d_trad, bw, gap_m, metrics_t, colors_t, sub_id);
end

center_avg_t = n_sub * (6*bw + gap_m + gap_s);
renderTimingBlock(center_avg_t, all_gmm_data, all_trad_data, bw, gap_m, metrics_t, colors_t, 'AVG');

h_t = []; for m=1:3, h_t(m) = bar(nan, nan, 'FaceColor', colors_t(m,:)); end
legend(h_t, m_labels_t, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 12);
ylabel('Duration (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
title('\bf Control Dynamics: Timing Comparison across Subjects', 'FontSize', 16);
grid on; set(gca, 'XTick', [], 'XColor', 'none'); box off;

fprintf('\n--- STATISTICAL ANALYSIS: TIMING (Grand Average) ---\n');
time_metrics = {'Act_Time_Hit', 'Rest_Time_Err', 'Act_Time_Miss'};
for m = 1:3
    all_o = [FullDatabase(strcmpi({FullDatabase.Method}, 'gmm')).(time_metrics{m})];
    all_t = [FullDatabase(strcmpi({FullDatabase.Method}, 'traditional')).(time_metrics{m})];
    p = signrank(all_o, all_t);
    fprintf('%s: p-value = %.4f\n', time_metrics{m}, p);
end

%% CROSS-VALIDATION ANALYSIS
subjects = unique({FullDatabase.SubjectID});
n_sub = length(subjects);
metrics_sim = {'Act_Acc_NoTo', 'Rest_Acc', 'Act_Num_Timeout'};
m_labels = {'Trial Accuracy (no timeout)', 'Rest Accuracy','Timeout Rate'};

figure('Color', 'w', 'Position', [50 50 1400 700]); hold on;
bw = 0.12; gap_m = 0.05; gap_s = 0.7;

for s = 1:n_sub
    sub_id = subjects{s};
    center_sub = (s-1) * (6*bw + gap_m + gap_s);
    
    % Filtriamo le run registrate con il metodo TRADITIONAL
    db_trad = FullDatabase(strcmpi({FullDatabase.SubjectID}, sub_id) & strcmpi({FullDatabase.Method}, 'gmm'));
    
    x_real_start = center_sub - (3*bw + gap_m/2);
    x_sim_start  = center_sub + gap_m/2;
    
    for m = 1:3
        % --- REAL ---
        v_real = [db_trad.(metrics_sim{m})];
        if m == 3, v_real = (v_real / 20) * 100; end
        m_r = mean(v_real, 'omitnan');
        sem_r = std(v_real, 'omitnan') / sqrt(sum(~isnan(v_real)));
        
        pos_r = x_real_start + (m-1)*bw + bw/2;
        bar(pos_r, m_r, bw, 'FaceColor', metric_colors(m,:), 'EdgeColor', 'k', 'FaceAlpha', 0.3);
        errorbar(pos_r, m_r, sem_r, 'k.', 'LineWidth', 1);

        % --- SIMULATED ---
        v_sim = arrayfun(@(x) x.simulation.(metrics_sim{m}), db_trad);
        if m == 3, v_sim = (v_sim / 20) * 100; end
        m_s = mean(v_sim, 'omitnan');
        sem_s = std(v_sim, 'omitnan') / sqrt(sum(~isnan(v_sim)));
        
        pos_s = x_sim_start + (m-1)*bw + bw/2;
        bar(pos_s, m_s, bw, 'FaceColor', metric_colors(m,:), 'EdgeColor', 'k');
        errorbar(pos_s, m_s, sem_s, 'k.', 'LineWidth', 1);
    end
    
    text(x_real_start + 1.5*bw, -2, 'real', 'HorizontalAlignment', 'center', 'FontSize', 9);
    text(x_sim_start + 1.5*bw, -2, 'sim', 'HorizontalAlignment', 'center', 'FontSize', 9);
    plot([x_real_start, x_sim_start + 3*bw], [-5 -5], 'k-', 'LineWidth', 1.2);
    text(center_sub, -8, sub_id, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

center_avg = n_sub * (6*bw + gap_m + gap_s);
db_all_trad = FullDatabase(strcmpi({FullDatabase.Method}, 'gmm'));

x_avg_r = center_avg - (3*bw + gap_m/2);
x_avg_s = center_avg + gap_m/2;

for m = 1:3
    % Real AVG
    v_r = [db_all_trad.(metrics_sim{m})];
    if m == 3, v_r = (v_r / 20) * 100; end
    bar(x_avg_r + (m-1)*bw + bw/2, mean(v_r, 'omitnan'), bw, 'FaceColor', metric_colors(m,:), 'EdgeColor', 'k', 'FaceAlpha', 0.3);
    errorbar(x_avg_r + (m-1)*bw + bw/2, mean(v_r, 'omitnan'), std(v_r, 'omitnan')/sqrt(sum(~isnan(v_r))), 'k.', 'LineWidth', 1);
    
    % Sim AVG
    v_s = arrayfun(@(x) x.simulation.(metrics_sim{m}), db_all_trad);
    if m == 3, v_s = (v_s / 20) * 100; end
    bar(x_avg_s + (m-1)*bw + bw/2, mean(v_s, 'omitnan'), bw, 'FaceColor', metric_colors(m,:), 'EdgeColor', 'k');
    errorbar(x_avg_s + (m-1)*bw + bw/2, mean(v_s, 'omitnan'), std(v_s, 'omitnan')/sqrt(sum(~isnan(v_s))), 'k.', 'LineWidth', 1);
end
text(x_avg_r + 1.5*bw, -2, 'real', 'HorizontalAlignment', 'center', 'FontSize', 9);
text(x_avg_s + 1.5*bw, -2, 'sim', 'HorizontalAlignment', 'center', 'FontSize', 9);
text(center_avg, -8, 'AVG', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
plot([x_avg_r, x_avg_s + 3*bw], [-5 -5], 'k-', 'LineWidth', 1.5);

% Legenda e Labeling
h_sim = []; for m=1:3, h_sim(m) = bar(nan, nan, 'FaceColor', metric_colors(m,:)); end
legend(h_sim, m_labels, 'Location', 'southoutside', 'Orientation', 'horizontal');
ylabel('Percentage (%)', 'FontSize', 14, 'FontWeight', 'bold');
title('\bf Cross-Validation: GMM Gating Recovery on Traditional Baseline Sessions', 'FontSize', 16);
grid on; ylim([0 120]); set(gca, 'XTick', [], 'XColor', 'none'); box off;

fprintf('\n--- STATISTICAL ANALYSIS: CROSS-VALIDATION (Real Trad vs Sim GMM) ---\n');

% Metriche da testare
metrics_sim = {'Rest_Acc', 'Act_Acc_NoTo', 'Act_Num_Timeout'};
m_labels = {'Rest Accuracy', 'Effective Acc', 'Timeout Rate'};

db_only = FullDatabase(strcmpi({FullDatabase.Method}, 'gmm'));

for m = 1:3
    field = metrics_sim{m};
    v_real = [db_only.(field)];
    v_sim = arrayfun(@(x) x.simulation.(field), db_only);
    if m == 3
        v_real = (v_real / 20) * 100;
        v_sim = (v_sim / 20) * 100;
    end
    
    % Test di Wilcoxon per campioni appaiati
    [p, h] = signrank(v_real, v_sim);
    
    fprintf('%s:\n', m_labels{m});
    fprintf('  - Mean Real GMM: %.2f%%\n', mean(v_real, 'omitnan'));
    fprintf('  - Mean Sim Trad:  %.2f%%\n', mean(v_sim, 'omitnan'));
    fprintf('  - p-value: %.4f %s\n', p, char(60*(p<0.05) + 32));
end

%% --- FUNCTIONS ---
function renderTimingBlock(center, d_gmm, d_trad, bw, gap, metrics, colors, label)
    x_our = center - (3*bw + gap/2);
    x_tr  = center + gap/2;
    for m = 1:3
        % Our
        val_o = mean([d_gmm.(metrics{m})], 'omitnan');
        sem_o = std([d_gmm.(metrics{m})], 'omitnan') / sqrt(sum(~isnan([d_gmm.(metrics{m})])));
        pos_o = x_our + (m-1)*bw + bw/2;
        bar(pos_o, val_o, bw, 'FaceColor', colors(m,:), 'EdgeColor', 'k');
        errorbar(pos_o, val_o, sem_o, 'k.', 'LineWidth', 1.1);
        
        % Trad
        val_t = mean([d_trad.(metrics{m})], 'omitnan');
        sem_t = std([d_trad.(metrics{m})], 'omitnan') / sqrt(sum(~isnan([d_trad.(metrics{m})])));
        pos_t = x_tr + (m-1)*bw + bw/2;
        bar(pos_t, val_t, bw, 'FaceColor', colors(m,:), 'EdgeColor', 'k', 'FaceAlpha', 0.4);
        errorbar(pos_t, val_t, sem_t, 'k.', 'LineWidth', 1.1);
    end
    text(x_our + 1.5*bw, -0.2, 'our', 'HorizontalAlignment', 'center', 'FontSize', 10);
    text(x_tr + 1.5*bw, -0.2, 'trad', 'HorizontalAlignment', 'center', 'FontSize', 10);
    plot([x_our, x_tr + 3*bw], [-0.5 -0.5], 'k-', 'LineWidth', 1.5);
    text(center, -0.8, label, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12);
end

function renderBlock(center, d_gmm, d_trad, bw, gap, metrics, colors, label)
    x_our = center - (3*bw + gap/2);
    x_tr  = center + gap/2;
    
    for m = 1:3
        % --- GMM (OUR) ---
        vals_gmm = [d_gmm.(metrics{m})];
        if m == 3, vals_gmm = vals_gmm / 20 * 100; end
        m_gmm = mean(vals_gmm, 'omitnan');
        
        % CALCOLO SEM: SD diviso radice di N
        sem_gmm = std(vals_gmm, 'omitnan') / sqrt(sum(~isnan(vals_gmm)));
        
        p_our = x_our + (m-1)*bw + bw/2;
        bar(p_our, m_gmm, bw, 'FaceColor', colors(m,:), 'EdgeColor', 'k');
        errorbar(p_our, m_gmm, sem_gmm, 'k.', 'LineWidth', 1.1);
        
        % --- TRADITIONAL ---
        vals_tr = [d_trad.(metrics{m})];
        if m == 3, vals_tr = vals_tr / 20 * 100; end
        m_tr = mean(vals_tr, 'omitnan');
        
        % CALCOLO SEM: SD diviso radice di N
        sem_tr = std(vals_tr, 'omitnan') / sqrt(sum(~isnan(vals_tr)));
        
        p_tr = x_tr + (m-1)*bw + bw/2;
        bar(p_tr, m_tr, bw, 'FaceColor', colors(m,:), 'EdgeColor', 'k', 'FaceAlpha', 0.4);
        errorbar(p_tr, m_tr, sem_tr, 'k.', 'LineWidth', 1.1);
        
    end
    
    % Etichette asse X (stessa logica del codice precedente)
    text(x_our + 1.5*bw, -2, 'our', 'HorizontalAlignment', 'center', 'FontSize', 10);
    text(x_tr + 1.5*bw, -2, 'trad', 'HorizontalAlignment', 'center', 'FontSize', 10);
    plot([x_our, x_tr + 3*bw], [-5 -5], 'k-', 'LineWidth', 1.5);
    text(center, -8, label, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12);
end

