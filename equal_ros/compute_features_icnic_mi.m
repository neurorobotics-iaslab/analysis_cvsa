function [sparsity, label_sparsity] = compute_features_icnic_mi(c_signal, c_l, c_r, c_c)

    % 1. Calcola potenze medie per le 3 ROI
    P_left  = mean(c_signal(c_l));
    P_right = mean(c_signal(c_r));
    P_foot  = mean(c_signal(c_c));
    
    % Vettore delle 3 macro-aree
    motors_roi = [P_left, P_right, P_foot];
    
    % --- FEATURE: FOCALITY (Focus Ratio) ---
    % Logica: Contrast between Active Area (Min) and Idle Background (Median)
    % Aggiungiamo eps per evitare divisioni per zero
    baseline_locale = median(motors_roi); 
    active_area = min(motors_roi);
    
    val_focus = log((baseline_locale + eps) / (active_area + eps));
    
    % --- OUTPUT ---
    % Restituisce SOLO il Focus Ratio.
    sparsity = val_focus;
    
    % Etichetta corretta (Importante per i plot!)
    label_sparsity = {'LogMedianRatio (Focus)'};
end