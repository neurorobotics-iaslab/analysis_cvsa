%% MAIN_COMPARE_PARADIGMS  Side-by-side comparison of MI, CVSA, and Hybrid.
%
%   Loads the three .mat files produced by main_batch_evaluate (one per
%   paradigm), prints a comparison table, and plots grouped bar charts.
%
%   USAGE
%     Run after main_batch_evaluate for each paradigm (mi, cvsa, hybrid).
%     Three successive file dialogs will ask for each .mat.

clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
addpath(fullfile(this_dir, 'utils'));

% ── Select one .mat per paradigm ─────────────────────────────────────────────
PARADIGM_ORDER = {'mi', 'cvsa', 'hybrid'};
N_PAR          = numel(PARADIGM_ORDER);

default_dir = '/home/paolo/bci_vr_ws/recordings';
if ~isfolder(default_dir), default_dir = pwd; end

mat_paths = cell(N_PAR, 1);
for p = 1:N_PAR
    [fname, fdir] = uigetfile('*.mat', ...
        sprintf('Select eval .mat for paradigm: %s  (%d/%d)', ...
                upper(PARADIGM_ORDER{p}), p, N_PAR), default_dir);
    if isequal(fname, 0)
        error('main_compare_paradigms:cancel', ...
              'Cancelled at paradigm %s.', PARADIGM_ORDER{p});
    end
    mat_paths{p} = fullfile(fdir, fname);
    default_dir  = fdir;
end

% ── Load ──────────────────────────────────────────────────────────────────────
DATA   = cell(N_PAR, 1);
LABELS = cell(N_PAR, 1);
for p = 1:N_PAR
    d        = load(mat_paths{p});
    DATA{p}  = d;
    LABELS{p} = upper(d.paradigm_label);
end

% ── Print comparison table ────────────────────────────────────────────────────
METRICS = { ...
    'trial_acc',           'Trial accuracy (%)';         ...
    'trial_acc_no_reject', 'Trial acc no-reject (%)';    ...
    'mean_tth',            'Mean time-to-hit (s)';       ...
    'mean_sample_acc',     'Sample accuracy (%)';        ...
    'mean_confidence',     'Mean confidence (%)';        ...
    'mean_art_rate',       'Artifact rate (%)';          ...
    'mean_peak_norm',      'Mean peak norm';             ...
};
n_metrics = size(METRICS, 1);

fprintf('\n');
hdr = sprintf('%-28s', 'Metric');
for p = 1:N_PAR
    hdr = [hdr, sprintf('  %-18s', LABELS{p})]; %#ok<AGROW>
end
sep = repmat('-', 1, numel(hdr));
fprintf('%s\n%s\n', hdr, sep);

for m = 1:n_metrics
    fld   = METRICS{m, 1};
    label = METRICS{m, 2};
    row   = sprintf('%-28s', label);
    for p = 1:N_PAR
        v   = DATA{p}.aggregate.(fld);
        row = [row, sprintf('  %6.2f +/- %-7.2f', v.mean, v.std)]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end

% ── Per-class breakdown ───────────────────────────────────────────────────────
fprintf('\n%s\nPer-class metrics:\n%s\n', sep, sep);
for p = 1:N_PAR
    agg = DATA{p}.aggregate;
    fprintf('\n  %s  (%d files):\n', LABELS{p}, agg.n_files);
    for c = 1:numel(agg.per_class)
        pc = agg.per_class(c);
        fprintf('    Class %d (code %d)  hit=%.1f+/-%.1f%%  prec=%.1f+/-%.1f%%  ' ...
                'recall=%.1f+/-%.1f%%  tth=%.2f+/-%.2fs  peak=%.3f+/-%.3f\n', ...
                c, pc.class_code, ...
                pc.hit_rate.mean,   pc.hit_rate.std, ...
                pc.precision.mean,  pc.precision.std, ...
                pc.recall.mean,     pc.recall.std, ...
                pc.mean_tth.mean,   pc.mean_tth.std, ...
                pc.peak_norm.mean,  pc.peak_norm.std);
    end
end
fprintf('\n%s\n', sep);

% ── Bar chart: scalar metrics ─────────────────────────────────────────────────
PLOT_METRICS = { ...
    'trial_acc',           'Trial accuracy (%)';       ...
    'trial_acc_no_reject', 'Trial acc no-reject (%)';  ...
    'mean_tth',            'Mean TTH (s)';             ...
    'mean_sample_acc',     'Sample accuracy (%)';      ...
    'mean_confidence',     'Mean confidence (%)';      ...
    'mean_art_rate',       'Artifact rate (%)';        ...
    'mean_peak_norm',      'Mean peak norm';           ...
};
n_plot = size(PLOT_METRICS, 1);
n_col  = 4;
n_row  = ceil(n_plot / n_col);

COLORS = [0.22 0.48 0.76;   % MI    (blue)
           0.87 0.47 0.12;   % CVSA  (orange)
           0.17 0.63 0.30];  % Hybrid (green)

figure('Name','Paradigm comparison — scalar metrics', 'Color','w', ...
       'NumberTitle','off', 'Position',[50 50 1400 680]);

for m = 1:n_plot
    fld   = PLOT_METRICS{m, 1};
    label = PLOT_METRICS{m, 2};
    ax    = subplot(n_row, n_col, m);
    hold(ax, 'on');

    means = zeros(N_PAR, 1);
    stds  = zeros(N_PAR, 1);
    for p = 1:N_PAR
        v       = DATA{p}.aggregate.(fld);
        means(p) = v.mean;
        stds(p)  = v.std;
    end

    b = bar(ax, 1:N_PAR, means, 'FaceColor','flat');
    for p = 1:N_PAR, b.CData(p,:) = COLORS(p,:); end
    errorbar(ax, 1:N_PAR, means, stds, 'k.', 'LineWidth', 1.2, 'CapSize', 8);

    % Annotate bar values
    for p = 1:N_PAR
        text(ax, p, means(p) + stds(p) * 0.15 + (max(means)*0.03), ...
             sprintf('%.1f', means(p)), 'HorizontalAlignment','center', ...
             'FontSize', 7.5, 'Color','k');
    end

    ax.XTick      = 1:N_PAR;
    ax.XTickLabel = LABELS;
    ylabel(ax, label, 'FontSize', 9);
    title(ax, label, 'FontSize', 9);
    grid(ax, 'on');
end

sgtitle('Paradigm comparison: MI vs CVSA vs Hybrid', 'FontSize', 13, 'FontWeight','bold');

% ── Grouped bar: per-class hit rate ───────────────────────────────────────────
n_cls = numel(DATA{1}.aggregate.per_class);
figure('Name','Per-class hit rate by paradigm', 'Color','w', 'NumberTitle','off', ...
       'Position',[50 50 700 420]);
ax2 = axes();
hold(ax2, 'on');

% Matrix: rows = paradigms, cols = classes (bar groups by paradigm)
hit_mat  = zeros(N_PAR, n_cls);
std_mat  = zeros(N_PAR, n_cls);
for p = 1:N_PAR
    for c = 1:n_cls
        hit_mat(p,c) = DATA{p}.aggregate.per_class(c).hit_rate.mean;
        std_mat(p,c) = DATA{p}.aggregate.per_class(c).hit_rate.std;
    end
end
b2 = bar(ax2, 1:N_PAR, hit_mat);
for c = 1:n_cls
    % Lighten the colour slightly for class 2
    b2(c).FaceColor    = COLORS(min(c,size(COLORS,1)),:) .* (0.6 + 0.4*c/n_cls);
    b2(c).DisplayName  = sprintf('Class %d (code %d)', c, ...
                                  DATA{1}.aggregate.per_class(c).class_code);
end
ax2.XTick      = 1:N_PAR;
ax2.XTickLabel = LABELS;
ylabel(ax2, 'Hit rate (%)');
title(ax2, 'Per-class hit rate by paradigm', 'FontSize', 11);
legend(ax2, 'Location','best');
grid(ax2, 'on');
ylim(ax2, [0, 110]);

% ── Per-class precision vs recall scatter ─────────────────────────────────────
figure('Name','Precision vs Recall by paradigm', 'Color','w', 'NumberTitle','off', ...
       'Position',[50 50 600 500]);
ax3 = axes();
hold(ax3, 'on');
markers = {'o','s','^','d'};
for p = 1:N_PAR
    for c = 1:n_cls
        pc = DATA{p}.aggregate.per_class(c);
        errorbar(ax3, pc.precision.mean, pc.recall.mean, ...
                 pc.recall.std, pc.recall.std, ...
                 pc.precision.std, pc.precision.std, ...
                 markers{mod(c-1, numel(markers))+1}, ...
                 'Color', COLORS(p,:), 'MarkerFaceColor', COLORS(p,:), ...
                 'MarkerSize', 10, 'LineWidth', 1.2, 'CapSize', 6, ...
                 'DisplayName', sprintf('%s cls%d', LABELS{p}, c));
    end
end
plot(ax3, [0 100], [0 100], 'k:', 'HandleVisibility','off');
xlabel(ax3, 'Precision (%)');
ylabel(ax3, 'Recall / Hit rate (%)');
title(ax3, 'Precision vs Recall by paradigm and class', 'FontSize', 11);
xlim(ax3, [0 105]); ylim(ax3, [0 105]);
legend(ax3, 'Location','best', 'FontSize', 8);
grid(ax3, 'on');

fprintf('Done.\n');
