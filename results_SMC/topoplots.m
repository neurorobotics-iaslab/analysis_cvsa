clear; clc; close all;

%% CARICAMENTO
data_path = '/home/paolo/cvsa/ic_cvsa_ws/recordings/record_mi/results_graz/fisher';
files = dir(fullfile(data_path, 'Metrics_*.mat'));
n_subjects = length(files);
if n_subjects == 0, error('Nessun file trovato.'); end


%% INIZIALIZZAZIONE STRUTTURE
topo_IC = cell(1,2);   topo_All = cell(1,2);   topo_NIC = cell(1,2);

band_names = {'\mu-band (8-13 Hz)', '\beta-band (18-24 Hz)'};

% Carica tutti i soggetti per entrambe le bande
for i = 1:n_subjects
    load(fullfile(files(i).folder, files(i).name));
    
    for b = 1:2
        % Topoplot
        topo_IC{b}(i, :) = topo_diff_struct.ic{b};
        topo_All{b}(i, :) = topo_diff_struct.all{b};
        topo_NIC{b}(i, :) = topo_diff_struct.nic{b};
    end
end

%% --- TOPOPLOTS SUBJ ---
if exist('chanlocs_subset', 'var')

    for c_subject = 1:n_subjects
        figure('Color', 'w', 'Name', 'Fig 3: Topoplots SUBJ', 'Position', [100, 50, 800, 1000]);
        c_subj_label = files(c_subject).name(9:10);
        for b = 1:2
            grand_topo_IC = topo_IC{b}(c_subject,:);
            grand_topo_All = topo_All{b}(c_subject,:);
            grand_topo_NIC = topo_NIC{b}(c_subject,:);

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

%         sgtitle(['subject: ' c_subj_label]);
    end

    colormap(jet);
    % sgtitle('Neurophysiological Validation (Mean Difference Class 1 - Class 2)', 'FontSize', 14, 'FontWeight', 'bold');
else
    warning('Variabile "chanlocs_subset" non trovata. Impossibile fare i topoplot.');
end

%% --- TOPOPLOTS AVG ---
if exist('chanlocs_subset', 'var')
    figure('Color', 'w', 'Name', 'Fig 3: Topoplots AVG', 'Position', [100, 50, 800, 1000]);
    
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