%% MAIN_COMPARE_SESSIONS  Cross-session metric comparison grouped by paradigm.
%
%   Selects a root folder via GUI. Recursively finds all
%   eval_single_*.mat files produced by main_evaluate_single, loads them,
%   and for each paradigm found (mi / cvsa / hybrid):
%     1. Prints a per-file metric table.
%     2. Prints mean ± std for every scalar metric.
%     3. Produces a summary figure with bar charts per session.
%
%   To skip the GUI, set  root_dir = '/your/path'  in the workspace
%   before running this script.

clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir,'utils'));

% ── Folder picker ─────────────────────────────────────────────────────────────
if ~exist('root_dir','var') || isempty(root_dir)
    root_dir = uigetdir('/home/paolo/bci_vr_ws/recordings', ...
                        'Select root folder containing eval_single_*.mat files');
    if isequal(root_dir, 0)
        error('main_compare_sessions:cancel', 'No folder selected.');
    end
end

% ── Recursive mat search ──────────────────────────────────────────────────────
hits = dir(fullfile(root_dir, '**', 'eval_single_*.mat'));
if isempty(hits)
    fprintf('No eval_single_*.mat files found under:\n  %s\n', root_dir);
    return
end
fprintf('Found %d file(s) under  %s\n\n', numel(hits), root_dir);

% ── Load all ──────────────────────────────────────────────────────────────────
all_res = [];
for i = 1:numel(hits)
    d = load(fullfile(hits(i).folder, hits(i).name), 'results');
    if ~isfield(d, 'results'), continue; end
    all_res = [all_res; d.results]; %#ok<AGROW>
end
if isempty(all_res)
    fprintf('No valid results structs found.\n'); return
end

paradigms = unique({all_res.paradigm});

% ── Per-file table ────────────────────────────────────────────────────────────
fprintf('%-40s  %-7s  %5s  %5s  %5s  %7s  %7s\n', ...
        'File', 'Paradigm', 'Acc%', 'Ev%', 'SA%', 'TTH_C1', 'TTH_C2');
fprintf('%s\n', repmat('-', 1, 85));
for i = 1:numel(all_res)
    r = all_res(i);
    fprintf('%-40s  %-7s  %5.1f  %5.1f  %5.1f  %7.2f  %7.2f\n', ...
            r.basename, r.paradigm, ...
            100*r.trial_acc, 100*r.trial_acc_events, 100*r.sample_acc_total, ...
            safe_val(r.tth_per_class, 1), safe_val(r.tth_per_class, 2));
end
fprintf('\n');

% ── Per-paradigm aggregate ────────────────────────────────────────────────────
for pi = 1:numel(paradigms)
    par  = paradigms{pi};
    mask = strcmp({all_res.paradigm}, par);
    sub  = all_res(mask);
    n    = numel(sub);

    fprintf('══  %s  (%d session(s))  ══════════════════════════════════\n', upper(par), n);

    scalar_fields  = {'trial_acc','trial_acc_events','trial_acc_no_rej', ...
                      'sample_acc_total','sample_acc_mean_trial','art_rate'};
    scalar_labels  = {'trial acc','acc events (897/898)','acc no-art-rej', ...
                      'sample acc total','sample acc mean/trial','artifact rate'};
    for f = 1:numel(scalar_fields)
        if ~isfield(sub(1), scalar_fields{f}), continue; end
        vals = 100 * [sub.(scalar_fields{f})];
        fprintf('  %-26s  %.1f ± %.1f %%\n', scalar_labels{f}, ...
                mean(vals,'omitnan'), std(vals,'omitnan'));
    end

    n_cls = numel(sub(1).tth_per_class);
    for c = 1:n_cls
        tth_v  = arrayfun(@(r) safe_val(r.tth_per_class,   c), sub);
        hr_v   = arrayfun(@(r) safe_val(r.hit_rate_cls,    c), sub) * 100;
        miss_v = arrayfun(@(r) safe_val(r.t_miss_peak_cls, c), sub);
        fprintf('  class %d  hit=%.1f±%.1f%%  tth=%.2f±%.2fs  t_miss_peak=%.2f±%.2fs\n', ...
                c, mean(hr_v,'omitnan'),   std(hr_v,'omitnan'), ...
                   mean(tth_v,'omitnan'),  std(tth_v,'omitnan'), ...
                   mean(miss_v,'omitnan'), std(miss_v,'omitnan'));
    end
    fprintf('\n');

    % ── Summary figure ────────────────────────────────────────────────────────
    fig = figure('Name', sprintf('Sessions — %s', upper(par)), ...
                 'Color','w','NumberTitle','off','Position',[60 60 1000 450]);
    subj_labels = {sub.basename};

    subplot(1, 3, 1);
    acc_mat = 100 * [[sub.trial_acc]; [sub.trial_acc_events]; [sub.trial_acc_no_rej]];
    bar(acc_mat.');
    legend({'Simulated','Events (897/898)','No-art-rej'}, 'FontSize',7,'Location','south');
    set(gca,'XTick',1:n,'XTickLabel',subj_labels,'XTickLabelRotation',30,'YLim',[0,105]);
    ylabel('%'); grid on; title(sprintf('%s — Trial accuracy', upper(par)), 'FontSize',9);
    yline(50,'--k');

    subplot(1, 3, 2);
    sa_mat = 100 * [[sub.sample_acc_total]; [sub.sample_acc_mean_trial]];
    bar(sa_mat.');
    legend({'Total','Mean/trial'}, 'FontSize',7,'Location','south');
    set(gca,'XTick',1:n,'XTickLabel',subj_labels,'XTickLabelRotation',30,'YLim',[0,105]);
    ylabel('%'); grid on; title('Frame-level accuracy', 'FontSize',9);
    yline(50,'--k');

    subplot(1, 3, 3);
    tth_mat = nan(n_cls, n);
    for c = 1:n_cls
        tth_mat(c,:) = arrayfun(@(r) safe_val(r.tth_per_class, c), sub);
    end
    bar(tth_mat.');
    legend(arrayfun(@(c) sprintf('C%d',c), 1:n_cls,'un',0), 'FontSize',7,'Location','north');
    set(gca,'XTick',1:n,'XTickLabel',subj_labels,'XTickLabelRotation',30);
    ylabel('s'); grid on; title('Mean time to hit', 'FontSize',9);

    sgtitle(fig, sprintf('%s — %d session(s)', upper(par), n), 'FontSize',11);
end

% ── Local helper ──────────────────────────────────────────────────────────────
function v = safe_val(arr, idx)
    if idx <= numel(arr), v = arr(idx); else, v = NaN; end
end
