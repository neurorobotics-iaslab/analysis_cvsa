clear; clc; close all;

%% CARICAMENTO
data_path = '/home/paolo/cvsa/ic_cvsa_ws/record_mi/results_graz/fisher';
files = dir(fullfile(data_path, 'Metrics_*.mat'));
n_subjects = length(files);
if n_subjects == 0, error('Nessun file trovato.'); end


%% INIZIALIZZAZIONE STRUTTURE
% Indici: 1 = Mu (8-13Hz), 2 = Beta (18-24Hz)
fisher_IC = cell(1,2); fisher_All = cell(1,2); fisher_NIC = cell(1,2);
r2_IC = cell(1,2);     r2_All = cell(1,2);     r2_NIC = cell(1,2);
topo_IC = cell(1,2);   topo_All = cell(1,2);   topo_NIC = cell(1,2);

band_names = {'\mu-band (8-13 Hz)', '\beta-band (18-24 Hz)'};

% Carica tutti i soggetti per entrambe le bande
for i = 1:n_subjects
    load(fullfile(files(i).folder, files(i).name));
    
    for b = 1:2
        offset_fisher = (b-1)*3; 
        
        % Fisher
        f_ic  = fisher_score_subj(1 + offset_fisher, :);
        f_all = fisher_score_subj(2 + offset_fisher, :);
        f_nic = fisher_score_subj(3 + offset_fisher, :);
        max_f_subj = max([f_ic, f_all, f_nic]); 
        if max_f_subj == 0, max_f_subj = 1; end 
        fisher_IC{b}(i, :)  = f_ic / max_f_subj;
        fisher_All{b}(i, :) = f_all / max_f_subj;
        fisher_NIC{b}(i, :) = f_nic / max_f_subj;
        
        % R2
        r_ic  = r2_struct.ic{b};
        r_all = r2_struct.all{b};
        r_nic = r2_struct.nic{b};
        max_r_subj = max(abs([r_ic, r_all, r_nic]));
        if max_r_subj == 0, max_r_subj = 1; end
        r2_IC{b}(i, :)  = r_ic / max_r_subj;
        r2_All{b}(i, :) = r_all / max_r_subj;
        r2_NIC{b}(i, :) = r_nic / max_r_subj;
        
        % Topoplot
        topo_IC{b}(i, :) = topo_diff_struct.ic{b};
        topo_All{b}(i, :) = topo_diff_struct.all{b};
        topo_NIC{b}(i, :) = topo_diff_struct.nic{b};
    end
end

% Parametri comuni
n_ch = size(fisher_IC{1}, 2);
x = 1:n_ch;
width = 0.25;
colors = {[0.8 0.2 0.2], [0.6 0.6 0.6], [0.2 0.2 0.8]};


%% --- plot ---
figure('Color', 'w', 'Name', 'Fig 2: Class Separability Metrics', 'Position', [50, 50, 1400, 800]);

ax_fisher = zeros(1,2);
ax_r2 = zeros(1,2);

for b = 1:2
    % --- FISHER SCORE ---
    ax_fisher(b) = subplot(2, 2, b); hold on;
    
    m_F_IC = mean(fisher_IC{b}, 1);  s_F_IC = std(fisher_IC{b}, 0, 1)/sqrt(n_subjects);
    m_F_All = mean(fisher_All{b}, 1); s_F_All = std(fisher_All{b}, 0, 1)/sqrt(n_subjects);
    m_F_NIC = mean(fisher_NIC{b}, 1); s_F_NIC = std(fisher_NIC{b}, 0, 1)/sqrt(n_subjects);
    
    b1 = bar(x - width, m_F_IC, width, 'FaceColor', colors{1}, 'EdgeColor', 'none');
    b2 = bar(x, m_F_All, width, 'FaceColor', colors{2}, 'EdgeColor', 'none');
    b3 = bar(x + width, m_F_NIC, width, 'FaceColor', colors{3}, 'EdgeColor', 'none');
    errorbar(x - width, m_F_IC, s_F_IC, 'k.', 'LineWidth', 1);
    errorbar(x, m_F_All, s_F_All, 'k.', 'LineWidth', 1);
    errorbar(x + width, m_F_NIC, s_F_NIC, 'k.', 'LineWidth', 1);
    ylabel('Fisher Score');
    xticks(x); xticklabels(channels_label); xtickangle(45);
    title(['Fisher Score: ' band_names{b}], 'FontSize', 12, 'FontWeight', 'bold');
    grid on; xlim([0.5 n_ch+0.5]); 
    if b == 1
        legend([b1 b2 b3], {'GMM Selected (IC)', 'Standard (All Data)', 'Rejected (INC)'}, 'Location', 'northeast');
    end

    % --- SIGNED R2 ---
    ax_r2(b) = subplot(2, 2, b+2); hold on;
    m_R_IC = mean(r2_IC{b}, 1);  s_R_IC = std(r2_IC{b}, 0, 1)/sqrt(n_subjects);
    m_R_All = mean(r2_All{b}, 1); s_R_All = std(r2_All{b}, 0, 1)/sqrt(n_subjects);
    m_R_NIC = mean(r2_NIC{b}, 1); s_R_NIC = std(r2_NIC{b}, 0, 1)/sqrt(n_subjects);

    bar(x - width, m_R_IC, width, 'FaceColor', colors{1}, 'EdgeColor', 'none');
    bar(x, m_R_All, width, 'FaceColor', colors{2}, 'EdgeColor', 'none');
    bar(x + width, m_R_NIC, width, 'FaceColor', colors{3}, 'EdgeColor', 'none');
    errorbar(x - width, m_R_IC, s_R_IC, 'k.', 'LineWidth', 1);
    errorbar(x, m_R_All, s_R_All, 'k.', 'LineWidth', 1);
    errorbar(x + width, m_R_NIC, s_R_NIC, 'k.', 'LineWidth', 1);
    yline(0, 'k-', 'LineWidth', 1); 
    ylabel('Signed R^2');
    xticks(x); xticklabels(channels_label); xtickangle(45);
    title(['Signed R^2: ' band_names{b}], 'FontSize', 12, 'FontWeight', 'bold');
    grid on; xlim([0.5 n_ch+0.5]);
end

linkaxes(ax_fisher, 'y'); 
linkaxes(ax_r2, 'y');     

sgtitle(['Calibration Quality Assessment (Grand Average, N=' num2str(n_subjects) ')'], 'FontSize', 14, 'FontWeight', 'bold');

%% Friedman + Wilcoxon
disp('===================================================================');
disp('          CHANNEL-BY-CHANNEL STATISTICAL ANALYSIS                  ');
disp('===================================================================');

p_friedman_fisher = zeros(2, n_ch);
p_friedman_r2 = zeros(2, n_ch);

for b = 1:2
    fprintf('\n---> %s <---\n', band_names{b});
    
    for ch = 1:n_ch
        % --- FISHER SCORE STATS ---
        data_F_ch = [fisher_IC{b}(:, ch), fisher_All{b}(:, ch), fisher_NIC{b}(:, ch)]; % [subj x 3]
        p_friedman_fisher(b, ch) = friedman(data_F_ch, 1, 'off');
        
        if p_friedman_fisher(b, ch) < 0.05
            p_IC_vs_All_F = signrank(data_F_ch(:,1), data_F_ch(:,2), 'tail', 'right', 'method', 'exact');
            p_IC_vs_NIC_F = signrank(data_F_ch(:,1), data_F_ch(:,3), 'tail', 'right', 'method', 'exact');
            
            fprintf('[Fisher] %s: P-Friedman = %.4f | Wilcoxon (IC vs All) = %.4f | Wilcoxon (IC vs INC) = %.4f\n', ...
                channels_label{ch}, p_friedman_fisher(b, ch), p_IC_vs_All_F, p_IC_vs_NIC_F);
        end
        
        % --- SIGNED R2 STATS ---
        data_R_ch = abs([r2_IC{b}(:, ch), r2_All{b}(:, ch), r2_NIC{b}(:, ch)] );
        p_friedman_r2(b, ch) = friedman(data_R_ch, 1, 'off');
        
        if p_friedman_r2(b, ch) < 0.05
            p_IC_vs_All_R = signrank(data_R_ch(:,1), data_R_ch(:,2), 'tail', 'right', 'method', 'exact');
            p_IC_vs_NIC_R = signrank(data_R_ch(:,1), data_R_ch(:,3), 'tail', 'right', 'method', 'exact');
            
            fprintf('[Abs R2] %s: P-Friedman = %.4f | Wilcoxon (IC vs All) = %.4f | Wilcoxon (IC vs INC) = %.4f\n', ...
                channels_label{ch}, p_friedman_r2(b, ch), p_IC_vs_All_R, p_IC_vs_NIC_R);
        end
    end
end
disp('===================================================================');

%% --- TOPOPLOTS ---
if exist('chanlocs_subset', 'var')
    figure('Color', 'w', 'Name', 'Fig 3: Topoplots Comparison', 'Position', [100, 50, 800, 1000]);
    
    for b = 1:2
        grand_topo_IC = mean(topo_IC{b}, 1);
        grand_topo_All = mean(topo_All{b}, 1);
        grand_topo_NIC = mean(topo_NIC{b}, 1);
        
        max_val = max(abs([grand_topo_IC, grand_topo_All, grand_topo_NIC]));
        clim_range = [-max_val, max_val];
        if max_val == 0, clim_range = [-1 1]; end 
        
        % --- GMM Selected ---
        subplot(3, 2, b);
        topoplot(grand_topo_IC, chanlocs_subset, 'electrodes', 'on', 'style', 'map', 'shading', 'interp', 'conv', 'on');
        colorbar; caxis(clim_range);
        title(['GMM Selected (' band_names{b} ')'], 'FontSize', 12, 'FontWeight', 'bold');
        
        % --- Standard ---
        subplot(3, 2, b + 2);
        topoplot(grand_topo_All, chanlocs_subset, 'electrodes', 'on', 'style', 'map', 'shading', 'interp', 'conv', 'on');
        colorbar; caxis(clim_range);
        title(['Standard (' band_names{b} ')'], 'FontSize', 12, 'FontWeight', 'bold');
        
        % --- Rejected ---
        subplot(3, 2, b + 4);
        topoplot(grand_topo_NIC, chanlocs_subset, 'electrodes', 'on', 'style', 'map', 'shading', 'interp', 'conv', 'on');
        colorbar; caxis(clim_range);
        title(['Rejected Data (' band_names{b} ')'], 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    colormap(jet);
    % sgtitle('Neurophysiological Validation (Mean Difference Class 1 - Class 2)', 'FontSize', 14, 'FontWeight', 'bold');
else
    warning('Variabile "chanlocs_subset" non trovata. Impossibile fare i topoplot.');
end