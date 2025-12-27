function [sparsity, label_sparsity] = compute_features_icnic(c_signal, type, o_l, o_r, c_l, c_r, nsparsity)
% signal: signal 1 x channels, according to the notion of my 39 channels
    
    sparsity = nan(nsparsity,1);
    label_sparsity = [{'LI'}, {'Gini'}];

    % --- LI --- Lateralization index
    if strcmp(type, 'cvsa')
        % show the occipital lateralization --> CVSA
        P_left_window  = mean(c_signal(o_l));
        P_right_window = mean(c_signal(o_r));
    elseif strcmp(type, 'mi')
        % show the central lateralization --> MI
        P_left_window  = mean(c_signal(c_l));
        P_right_window = mean(c_signal(c_r));
    else
        disp('ERROR')
    end
    LI = (P_right_window - P_left_window) ./ (P_right_window + P_left_window + eps);
    sparsity(1) = abs(LI);

    % --- max desyncronization ---
    Act_L = -log(P_left_window + eps);
    Act_R = -log(P_right_window + eps);
    Max_Activation = max(Act_L, Act_R);
    sparsity(1) = Max_Activation;

    % --- Gini Index  ---  -> IC means focus on a specific zone
    if isempty(o_l) | isempty(o_r) | isempty(c_l) | isempty(c_r)
        mean_roi_raw = c_signal;
    else
        mean_roi_raw = [mean(c_signal(c_l)), ...
                        mean(c_signal(c_r)), ...
                        mean(c_signal(o_l)), ...
                        mean(c_signal(o_r))];
    end
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

end