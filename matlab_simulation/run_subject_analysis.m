%% RUN_SUBJECT_ANALYSIS  Batch launcher: analyses for one subject's recording data.
%
%   Expected folder layout:
%     recordings/
%       <subject>/
%         <day>/
%           calibration/   ← topo_erders only (ERD/ERS topoplots)
%           evaluation/    ← all 7 analysis steps
%         <day>/
%           evaluation/
%
%   Usage (GUI — user picks a folder and the script figures out the rest):
%       run_subject_analysis
%
%   The GUI prompt accepts any of:
%     (a) an evaluation/ folder directly  → analyses that single session
%     (b) a <day>/ folder with evaluation/ inside → analyses that session
%     (c) a <subject>/ folder (name must start with 'a') → lists available days,
%         prompts user to pick which day(s) to process
%     (d) the recordings/ root → finds all subject folders starting with 'a',
%         processes ALL their days automatically (no day-selection prompt)
%
%   Programmatic equivalents:
%       run_subject_analysis('/path/recordings/S01/20250101/evaluation')
%       run_subject_analysis('/path/recordings/S01/20250101')
%       run_subject_analysis('/path/recordings/a001')
%       run_subject_analysis('/path/recordings')   % all 'a*' subjects
%
%   For each evaluation/ folder found, the following analyses run in order:
%     1. main_simulate          — all GDFs: per-trial pipeline replay + figure
%     2. main_session_overview  — all GDFs: 10-figure session summary (incl.
%                                 per-class confusion matrix, console ITR +
%                                 chance-level test)
%     3. main_trial_dynamics    — all GDFs: learning/fatigue trend figures
%     4. main_threshold_sweep   — all GDFs: (th1, th2) hit-rate heatmaps
%     5. main_roc_analysis      — all GDFs: classifier-level ROC/AUC
%                                 (pre-integrator); saves raw pooled
%                                 scores/labels per paradigm so
%                                 main_group_analysis can build a per-subject
%                                 and cross-subject averaged ROC
%     6. main_hybrid_advantage_probs   — hybrid GDFs only (skipped if none)
%     7. main_hybrid_advantage_integ   — hybrid GDFs only (skipped if none):
%                                 6-figure counterfactual analysis (incl.
%                                 per-stream confusion matrix, console ITR +
%                                 chance-level test)
%     8. main_validate_counterfactual  — needs outputs of steps 2+7 (auto-skipped)
%
%   For the sibling calibration/ folder (when present):
%     C. topo_erders            — all GDFs: ERD/ERS topoplots  (EEGLAB required)
%
%   All figures saved hidden (SHOW_FIGURES = false).
%   main_group_analysis.m is a multi-subject aggregator — run it separately.

function run_subject_analysis(root_dir)

this_dir   = fileparts(mfilename('fullpath'));
topo_dir   = fullfile(this_dir, '..', 'analysis_gdf');
addpath(this_dir, fullfile(this_dir,'io'), fullfile(this_dir,'utils'), topo_dir);

% ── Pick root folder (GUI or argument) ───────────────────────────────────
if nargin < 1 || isempty(root_dir)
    root_dir = uigetdir('/home/paolo/bci_vr_ws/recordings', ...
        'Select subject, day, or evaluation folder');
    if isequal(root_dir, 0)
        error('run_subject_analysis:cancel', 'No folder selected.');
    end
end
if ~isfolder(root_dir)
    error('run_subject_analysis:notfound', 'Folder not found: %s', root_dir);
end

% ── Resolve to a list of evaluation/ folders ─────────────────────────────
eval_dirs = resolve_eval_dirs(root_dir);
if isempty(eval_dirs)
    error('run_subject_analysis:noeval', ...
          'No evaluation/ folder found under %s', root_dir);
end

fprintf('\n══════════════════════════════════════════════════════\n');
fprintf('  run_subject_analysis\n');
fprintf('  Root   : %s\n', root_dir);
fprintf('  Sessions found: %d\n', numel(eval_dirs));
for k = 1:numel(eval_dirs)
    fprintf('    [%d] %s\n', k, eval_dirs{k});
end
fprintf('══════════════════════════════════════════════════════\n');

% ── Process each evaluation folder (+ sibling calibration/) ──────────────
for si = 1:numel(eval_dirs)
    eval_dir = eval_dirs{si};
    fprintf('\n\n══════ Session %d/%d : %s ══════\n', si, numel(eval_dirs), eval_dir);
    process_eval_folder(eval_dir);

    % calibration/ sits next to evaluation/ (same day folder)
    day_dir  = fileparts(eval_dir);
    calib_dir = fullfile(day_dir, 'calibration');
    if isfolder(calib_dir) && ~isempty(dir(fullfile(calib_dir, '*.gdf')))
        fprintf('\n── Calibration : %s ──\n', calib_dir);
        process_calibration_folder(calib_dir);
    end
end

fprintf('\n══════ run_subject_analysis complete ══════\n\n');

end % run_subject_analysis


% ── Local helpers ─────────────────────────────────────────────────────────────

function process_eval_folder(eval_dir)
% Run the full analysis pipeline on a single evaluation/ folder.

    f = dir(fullfile(eval_dir, '*.gdf'));
    if isempty(f)
        fprintf('  [WARN] No GDF files found — skipping %s\n', eval_dir);
        return;
    end
    all_names = {f.name};

    mi_names = {}; cvsa_names = {}; hybrid_names = {};
    for k = 1:numel(all_names)
        par = detect_paradigm_local(all_names{k});
        switch par
            case 'hybrid', hybrid_names{end+1} = all_names{k}; %#ok<AGROW>
            case 'cvsa',   cvsa_names{end+1}   = all_names{k}; %#ok<AGROW>
            case 'mi',     mi_names{end+1}      = all_names{k}; %#ok<AGROW>
            otherwise, fprintf('  [WARN] paradigm unknown for %s — skipped\n', all_names{k});
        end
    end
    fprintf('  GDFs: %d total  (MI:%d  CVSA:%d  Hybrid:%d)\n', ...
            numel(all_names), numel(mi_names), numel(cvsa_names), numel(hybrid_names));

    run_step(1, 'main_simulate',         @() main_simulate(eval_dir, all_names));
    run_step(2, 'main_session_overview', @() main_session_overview(eval_dir, all_names));
    run_step(3, 'main_trial_dynamics',   @() main_trial_dynamics(eval_dir, all_names));
    run_step(4, 'main_threshold_sweep',  @() main_threshold_sweep(eval_dir, all_names));
    run_step(5, 'main_roc_analysis',     @() main_roc_analysis(eval_dir, all_names));

    if isempty(hybrid_names)
        fprintf('[Step 6/9] main_hybrid_advantage_probs — SKIPPED (no hybrid GDFs)\n');
        fprintf('[Step 7/9] main_hybrid_advantage_integ — SKIPPED (no hybrid GDFs)\n');
    else
        run_step(6, 'main_hybrid_advantage_probs', ...
                 @() main_hybrid_advantage_probs(eval_dir, hybrid_names));
        run_step(7, 'main_hybrid_advantage_integ', ...
                 @() main_hybrid_advantage_integ(eval_dir, hybrid_names));
    end

    session_mat = fullfile(eval_dir, 'analysis_results', 'session_overview',       'session_summary.mat');
    cf_mat      = fullfile(eval_dir, 'analysis_results', 'hybrid_advantage_integ', 'counterfactual_summary.mat');
    if ~isfile(session_mat)
        fprintf('[Step 8/9] main_validate_counterfactual — SKIPPED (session_summary.mat missing)\n');
    elseif ~isfile(cf_mat)
        fprintf('[Step 8/9] main_validate_counterfactual — SKIPPED (counterfactual_summary.mat missing)\n');
    else
        run_step(8, 'main_validate_counterfactual', ...
                 @() main_validate_counterfactual(session_mat, cf_mat));
    end

    run_step(9, 'topo_erders',            @() topo_erders(eval_dir, all_names));

    fprintf('  Results: %s\n', fullfile(eval_dir, 'analysis_results'));
end

function process_calibration_folder(calib_dir)
% Run topo_erders on all GDFs in a calibration/ folder.
    f = dir(fullfile(calib_dir, '*.gdf'));
    if isempty(f)
        fprintf('  [WARN] No GDF files found — skipping %s\n', calib_dir);
        return;
    end
    all_names = {f.name};
    fprintf('  GDFs: %d (calibration)\n', numel(all_names));
    run_step_calib('topo_erders', @() topo_erders(calib_dir, all_names));
    fprintf('  Results: %s\n', fullfile(calib_dir, 'analysis_results'));
end

function eval_dirs = resolve_eval_dirs(root_dir)
% Return list of evaluation/ folders to process.
%
%   Case 1: root_dir is an evaluation/ folder itself → use directly.
%   Case 2: root_dir is a <day>/ folder            → use its evaluation/ child.
%   Case 3: root_dir is a <subject>/ folder         → list days, ask which to process.
%   Case 4: root_dir is the recordings root          → scan 'a*' subject folders,
%                                                      process ALL their days (batch).

    SUBJECT_PREFIX = 'a';   % only subject folders starting with this letter

    % Case 1: root_dir IS already an evaluation folder
    if is_eval_folder(root_dir)
        eval_dirs = {root_dir};
        return;
    end

    % Case 2: root_dir is a <day>/ folder (has evaluation/ directly inside)
    candidate = fullfile(root_dir, 'evaluation');
    if isfolder(candidate) && is_eval_folder(candidate)
        eval_dirs = {candidate};
        return;
    end

    % Case 4: root_dir is the recordings root — look for 'a*' subject subfolders
    %   that themselves contain <day>/evaluation/ structure.
    subj_entries = dir(root_dir);
    subj_entries = subj_entries([subj_entries.isdir]);
    subj_names   = {subj_entries.name};
    a_mask       = startsWith(lower(subj_names), SUBJECT_PREFIX) & ...
                   ~startsWith(subj_names, '.');
    a_subjects   = sort(subj_names(a_mask));

    recordings_root_eval = {};
    for s = 1:numel(a_subjects)
        subj_path = fullfile(root_dir, a_subjects{s});
        recordings_root_eval = [recordings_root_eval, ...
                                collect_eval_dirs_for_subject(subj_path)]; %#ok<AGROW>
    end

    if ~isempty(recordings_root_eval)
        % This is the recordings root — batch mode, process everything found
        eval_dirs = sort(recordings_root_eval);
        fprintf('  Found %d subject(s) starting with "%s": %s\n', ...
                numel(a_subjects), SUBJECT_PREFIX, strjoin(a_subjects, ', '));
        return;
    end

    % Case 3: root_dir is a <subject>/ folder — find days, ask which to process
    eval_dirs = collect_eval_dirs_for_subject(root_dir);
    if isempty(eval_dirs)
        return;
    end
    if numel(eval_dirs) > 1
        eval_dirs = ask_day_selection(eval_dirs, root_dir);
    end
end

function evals = collect_eval_dirs_for_subject(subj_path)
% Return all <day>/evaluation/ paths under subj_path, sorted chronologically.
    evals = {};
    day_entries = dir(subj_path);
    day_entries = day_entries([day_entries.isdir] & ~startsWith({day_entries.name}, '.'));
    for k = 1:numel(day_entries)
        candidate = fullfile(subj_path, day_entries(k).name, 'evaluation');
        if isfolder(candidate) && is_eval_folder(candidate)
            evals{end+1} = candidate; %#ok<AGROW>
        end
    end
    evals = sort(evals);
end

function eval_dirs = ask_day_selection(all_eval, subj_path)
% Show a listdlg so the user can pick which days to process.
    day_labels = cellfun(@(p) strsplit(p, filesep), all_eval, 'UniformOutput', false);
    day_labels = cellfun(@(parts) parts{end-1}, day_labels, 'UniformOutput', false);

    [sel, ok] = listdlg( ...
        'ListString',    day_labels, ...
        'SelectionMode', 'multiple', ...
        'Name',          'Select recording day(s)', ...
        'PromptString',  {sprintf('Subject: %s', subj_path), ...
                          'Select day(s) to analyse:'}, ...
        'ListSize',      [340, 200]);
    if ~ok || isempty(sel)
        error('run_subject_analysis:cancel', 'No day selected.');
    end
    eval_dirs = all_eval(sel);
end

function tf = is_eval_folder(d)
% True if d contains at least one GDF directly (no subdirectory scan).
    tf = ~isempty(dir(fullfile(d, '*.gdf')));
end

function run_step(n, name, fn)
    fprintf('[Step %d/9] %-35s ...  ', n, name);
    t0 = tic;
    try
        fn();
        fprintf('done (%.1f s)\n', toc(t0));
    catch ME
        fprintf('ERROR\n');
        fprintf('           %s: %s\n', ME.identifier, ME.message);
    end
end

function run_step_calib(name, fn)
    fprintf('[Calib]    %-35s ...  ', name);
    t0 = tic;
    try
        fn();
        fprintf('done (%.1f s)\n', toc(t0));
    catch ME
        fprintf('ERROR\n');
        fprintf('           %s: %s\n', ME.identifier, ME.message);
    end
end

function par = detect_paradigm_local(fname)
    s = lower(fname);
    if contains(s, 'hybrid'),   par = 'hybrid';
    elseif contains(s, 'cvsa'), par = 'cvsa';
    elseif contains(s, 'mi'),   par = 'mi';
    else,                       par = 'unknown';
    end
end
