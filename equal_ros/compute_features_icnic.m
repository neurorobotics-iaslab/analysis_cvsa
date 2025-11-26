function [sparsity, label_sparsity] = compute_features_icnic(c_signal, type, o_l, o_r, c_l, c_r, excl_chs)
% signal: signal 1 x channels, according to the notion of my 39 channels
    
    sparsity = nan(3,1);
    label_sparsity = [{'LI'},{'GI'},{'GB'}];

    % --- LI --- Lateralization index
    if all(type == 'cvsa')
        % show the occipital lateralization --> CVSA
        P_left_window  = mean(c_signal(o_l));
        P_right_window = mean(c_signal(o_r));
    elseif all(type == 'mi')
        % show the central lateralization --> MI
        P_left_window  = mean(c_signal(c_l));
        P_right_window = mean(c_signal(c_r));
    else
        disp('ERROR')
    end
    LI = (P_right_window - P_left_window) ./ (P_right_window + P_left_window + eps);
    sparsity(1) = abs(LI);

    % --- Gini Index  ---  -> IC means focus on a specific zone
    non_zeros_chs = setdiff(1:size(c_signal,2), excl_chs);
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