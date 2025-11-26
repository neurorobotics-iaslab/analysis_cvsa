function [sparsity, label_sparsity] = compute_features_icnic(c_signal, type)
% signal: signal 1 x channels, according to the notion of my 39 channels
    % take the one which contribute at the 95% of the energy
    sparsity = nan(3,1);
    label_sparsity = [{'LI'},{'GI'},{'GB'}];
    o_l = sort([29 13 30 37 33 34 17]); o_r = sort([31 15 32 35 36 38 18]);
    frontal = sort([3 4 5 20 21]); c_l = sort([6 22 25 8 11 27]); c_r = sort([7 24 26 10 12 28]);
    excluded_chs = [1,2,19];

    % --- LAP --- Calcola il LAP per tutti i punti nella finestra passata (window_signal)
    % show the occipital lateralization that is strongand present during the CVSA
    P_left_window  = mean(c_signal(o_l));
    P_right_window = mean(c_signal(o_r));
    LAP_history = (P_right_window - P_left_window) ./ (P_right_window + P_left_window + eps);
    sparsity(1) = abs(LAP_history); % LAP_Mean

    % --- Gini Index + Occipital Power ---  -> when high there is a zone stronger, so IC
    % show the focusing is weighted in with the power in the occipital part, in this way ig CVSA then strong value 
    non_zeros_chs = setdiff(1:size(c_signal,2), excluded_chs);
    global_mean = mean(c_signal(non_zeros_chs)); % car filter
    current_signal_normalized = c_signal - global_mean; % remove the global energy
    mean_roi_raw = [mean(current_signal_normalized(c_l)), ...
        mean(current_signal_normalized(c_r)), mean(current_signal_normalized(o_l)), ...
        mean(current_signal_normalized(o_r))];
    mean_roi = abs(mean_roi_raw); % make sure the energy is positive--> we are using peak and valli with same significance
    mean_roi_ordered = sort(mean_roi);
    n = length(mean_roi_ordered);
    sum_roi_p = 0;
    for i = 1:n
        sum_roi_p = sum_roi_p + (n+1-i) * mean_roi_ordered(i);
    end
    total_sum = sum(mean_roi_ordered);
    if total_sum > 0
        gi = (1/n) * (n+1-2*sum_roi_p/total_sum);
    else
        gi = 0;
    end
    sparsity(2) = gi;

    % --- GB ---
    % return the global power mean, 
    sparsity(3) = global_mean;
end