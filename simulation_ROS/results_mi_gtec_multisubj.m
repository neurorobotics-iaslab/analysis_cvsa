%% PLOT GLOBAL ANALYSIS (MULTI-SUBJECT)
clear; clc; close all;

% 1. Caricamento Dati
files = dir('*_results_database.mat'); % Trova tutti i file salvati prima
if isempty(files)
    error('Nessun file _results_database.mat trovato!');
end

GlobalDatabase = [];
subjects_list = {};

for i = 1:length(files)
    fprintf('Caricamento %s...\n', files(i).name);
    loaded = load(files(i).name);
    
    % Aggiungi una colonna "Subject" al database per distinguere
    temp_db = loaded.Database;
    sub_name = extractBefore(files(i).name, '_results');
    [temp_db(:).Subject] = deal(sub_name);
    subjects_list{end+1} = sub_name;
    
    % Unisci al database globale
    GlobalDatabase = [GlobalDatabase; temp_db];
end

% 2. Setup Plotting
methods = {GlobalDatabase.Method};
is_gmm = strcmpi(methods, 'gmm');
is_trad = strcmpi(methods, 'traditional') | strcmpi(methods, 'trad');
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980]; % Blu/Arancio

% Funzioni helper (copiate dal tuo script precedente)
get_d = @(field) extract_data_nan_zero(GlobalDatabase, is_gmm, is_trad, field);

%% FIGURA 1: GLOBAL SAFETY & PERFORMANCE
figure('Name', 'Global Analysis: Safety vs Perf', 'Color', 'w', 'Position', [100 100 1200 500]);

% Rest Accuracy (Safety)
subplot(1, 3, 1);
[d, g] = get_d('Rest_Acc');
custom_boxplot(d, g, colors, 'Rest Accuracy (%)', 'Global Safety', [-5 105]);
subtitle(['Mean GMM: ' num2str(mean(d(strcmp(g,'GMM'))), '%.1f') '%']);

% Active Accuracy (No Timeout)
subplot(1, 3, 2);
[d, g] = get_d('Act_Acc_NoTo');
custom_boxplot(d, g, colors, 'Accuracy (%)', 'Active Precision (Decision)', [40 105]);

% Control SNR (Quality)
subplot(1, 3, 3);
[d, g] = get_d('SNR_Ratio');
custom_boxplot(d, g, colors, 'SNR Ratio', 'Control Quality', []);
yline(1, 'k--');
subtitle('Values > 1 mean Signal > Noise');

sgtitle(['Analysis on ' num2str(length(files)) ' Subjects: ' strjoin(unique(subjects_list), ', ')]);

%% FIGURA 2: PARETO PROOF (Se hai implementato la metrica Pareto nel DB)
% Nota: Questo plot ha senso farlo "per soggetto" o "medio". 
% Qui facciamo un esempio di scatter plot: Safety vs Performance per ogni file

figure('Name', 'Efficiency Frontier Scatter', 'Color', 'w');
hold on; grid on;

% Dati GMM
gmm_rest = [GlobalDatabase(is_gmm).Rest_Acc];
gmm_act  = [GlobalDatabase(is_gmm).Act_Trial_Acc];
scatter(gmm_rest, gmm_act, 80, colors(1,:), 'filled', 'o', 'DisplayName', 'GMM Runs');

% Dati Traditional
trad_rest = [GlobalDatabase(is_trad).Rest_Acc];
trad_act  = [GlobalDatabase(is_trad).Act_Trial_Acc];
scatter(trad_rest, trad_act, 80, colors(2,:), 'filled', 's', 'DisplayName', 'Traditional Runs');

% Estetica
xlabel('Safety (Rest Accuracy) [%]');
ylabel('Performance (Raw Trial Accuracy) [%]');
title('Safety vs Performance Trade-off');
legend('Location', 'SouthEast');
xlim([0 105]); ylim([0 105]);

% Zona Ottimale (Alto Destra)
annotation('arrow',[0.2 0.8],[0.2 0.8], 'Color', 'k', 'LineStyle', '--');
text(50, 50, 'Better System Direction', 'Rotation', 45);

%% FUNZIONI DI SUPPORTO (Da includere in fondo allo script)
function [data, groups] = extract_data_nan_zero(db, idx_gmm, idx_trad, field)
    if ~isfield(db, field), error(['Campo mancante: ' field]); end
    d_gmm = [db(idx_gmm).(field)]';
    d_trad = [db(idx_trad).(field)]';
    d_gmm(isnan(d_gmm)) = 0; d_trad(isnan(d_trad)) = 0;
    if isempty(d_gmm), d_gmm=[]; end; if isempty(d_trad), d_trad=[]; end
    data = [d_gmm; d_trad];
    groups = [repmat({'GMM'}, length(d_gmm), 1); repmat({'Traditional'}, length(d_trad), 1)];
end

function custom_boxplot(data, groups, colors, y_label, t_title, y_lims)
    if isempty(data), text(0.5,0.5,'No Data'); title(t_title); axis off; return; end
    h = boxplot(data, groups, 'Colors', 'k', 'Symbol', 'o', 'Widths', 0.5);
    set(h, 'LineWidth', 1.2); ylabel(y_label, 'FontWeight', 'bold'); title(t_title); grid on;
    if ~isempty(y_lims), ylim(y_lims); end
    h_box = findobj(gca, 'Tag', 'Box'); idx = length(h_box);
    % Colorazione (Trick per colorare i box inversi)
    if idx >= 2
        patch(get(h_box(2),'XData'), get(h_box(2),'YData'), colors(1,:), 'FaceAlpha', 0.5); % GMM
        patch(get(h_box(1),'XData'), get(h_box(1),'YData'), colors(2,:), 'FaceAlpha', 0.5); % Trad
    end
end