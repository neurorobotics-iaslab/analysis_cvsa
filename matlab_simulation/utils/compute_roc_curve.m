function [fpr, tpr, thr, auc_val] = compute_roc_curve(scores, labels)
% COMPUTE_ROC_CURVE  Non-parametric ROC (no Statistics Toolbox required):
%   sweep every distinct score as a threshold, count cumulative TP/FP.
%   labels: logical, true = positive (class 1).
%
%   Shared by main_roc_analysis.m (per-file/per-paradigm classifier ROC) and
%   main_group_analysis.m (per-subject re-pooled ROC, from raw scores/labels
%   pooled across that subject's sessions -- valid because the same
%   deployed classifier/calibration produced all of a subject's frames,
%   unlike pooling raw scores ACROSS subjects/classifiers).
    scores = scores(:);
    labels = logical(labels(:));
    P = sum(labels);
    N = sum(~labels);
    if P == 0 || N == 0
        fpr = [0; 1]; tpr = [0; 1]; thr = [Inf; -Inf]; auc_val = NaN;
        return;
    end
    [sorted_scores, order] = sort(scores, 'descend');
    sorted_labels = labels(order);
    tp_cum = cumsum(sorted_labels);
    fp_cum = cumsum(~sorted_labels);
    is_last_of_tie = [diff(sorted_scores) ~= 0; true];
    tp_cum = tp_cum(is_last_of_tie);
    fp_cum = fp_cum(is_last_of_tie);
    thr_pts = sorted_scores(is_last_of_tie);
    tpr = [0; tp_cum(:) / P];
    fpr = [0; fp_cum(:) / N];
    thr = [Inf; thr_pts(:)];
    if fpr(end) < 1 || tpr(end) < 1
        fpr(end+1) = 1; tpr(end+1) = 1; thr(end+1) = -Inf; %#ok<AGROW>
    end
    auc_val = trapz(fpr, tpr);
end
