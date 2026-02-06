function [m] = computeMetrics_signal(gmm_prob, qda_prob, artifact, signals, event, event_start, task_classes, gmm_classes, integratorCfg)
% COMPUTEMETRICS_SIGNAL - Focus su Topomappe per bin di confidenza e Accuratezza QDA.

% --- Configurazione ---
CODE_SX = task_classes(1); CODE_DX = task_classes(2);
ic_index = find(gmm_classes == integratorCfg.ic_class_label);

% Definiamo i bin di "messa a fuoco" (GMM Probability)
m.bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0];
nbins = length(m.bins) - 1;
nbands = length(signals);
nchans = size(signals{1}, 2);

% Accumulatori
m.topo_power = zeros(nbins, nbands, nchans);
m.topo_count = zeros(nbins, 1);
m.error_conf = [];    % Lista confidenze GMM
m.error_correct = []; % Lista 1/0 (QDA indovina?)

% Estrazione Trial
cfPOS = event.POS(event.TYP == event_start);
cfDUR = event.DUR(event.TYP == event_start);
all_cues = event.TYP(ismember(event.TYP, task_classes));

for i = 1:length(cfPOS)
    idx_range = cfPOS(i) : (cfPOS(i) + cfDUR(i) - 1);
    cue = all_cues(i);
    
    % Dati trial correnti
    s_gmm = gmm_prob(idx_range, ic_index);
    s_qda = qda_prob(idx_range, :);
    s_art = artifact(idx_range);
    
    for t = 1:length(s_gmm)
        if s_art(t) == 1, continue; end % Salta artefatti
        
        conf = s_gmm(t);
        
        % 1. METRICA: TOPOMAPPE (Accumulo potenza per bin)
        for b = 1:nbins
            if conf >= m.bins(b) && conf <= m.bins(b+1)
                m.topo_count(b) = m.topo_count(b) + 1;
                for band_idx = 1:nbands
                    m.topo_power(b, band_idx, :) = squeeze(m.topo_power(b, band_idx, :))' + signals{band_idx}(idx_range(t), :);
                end
            end
        end
        
        % 2. METRICA: ERRORE vs CONFIDENZA (Solo Trial Active)
        if cue == CODE_SX || cue == CODE_DX
            is_sx = (cue == CODE_SX);
            qda_pred_sx = (s_qda(t, 1) >= 0.5);
            
            m.error_conf = [m.error_conf; conf];
            m.error_correct = [m.error_correct; (qda_pred_sx == is_sx)];
        end
    end
end

% Media finale per le mappe
for b = 1:nbins
    if m.topo_count(b) > 0
        m.topo_power(b, :, :) = m.topo_power(b, :, :) / m.topo_count(b);
    end
end

% Calcolo Accuratezza Media per Bin (per il plot finale)
m.acc_per_bin = zeros(nbins, 1);
for b = 1:nbins
    idx = m.error_conf >= m.bins(b) & m.error_conf < m.bins(b+1);
    if any(idx)
        m.acc_per_bin(b) = mean(m.error_correct(idx)) * 100;
    else
        m.acc_per_bin(b) = NaN;
    end
end

end