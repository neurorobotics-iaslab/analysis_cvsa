%% BATCH_RUN_SUBJECTS  Per-subject analysis pipeline for a hardcoded list of
%   subjects, no GUI prompts.
%
%   Edit RECORDINGS_ROOT and SUBJECTS below, then run this script directly.
%
%   Calls run_subject_analysis(<day>/evaluation) once per (subject, day)
%   found under each listed subject's folder. Passing a direct evaluation/
%   folder always resolves as run_subject_analysis's Case 1 (single
%   session) -- this is what avoids its own multi-day listdlg GUI picker,
%   which only fires when a whole <subject>/ folder with more than one day
%   is passed directly. For each session this also runs topo_erders on the
%   sibling calibration/ folder, same as run_subject_analysis does normally.
%
%   SAVE_PER_FILE_TOPOPLOTS controls topo_erders' per-file topoplot/heatmap/
%   ROI figures (grand-average figures + both .mat summaries are always
%   produced regardless). Each per-file figure is an EEGLAB topoplot() call
%   plus an SVG export, which profiling shows dominates a batch run's total
%   time (not the FBCSP/CSP signal pipeline) — set to false for a much
%   faster run when the per-file figures aren't needed.

RECORDINGS_ROOT = '/home/paolo/bci_vr_ws/recordings';
SUBJECTS = {'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', ...
            'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', ...
            'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c8'};   % <-- edit this list
SUBJECTS = {'a5', 'c6', 'c8'};
SAVE_PER_FILE_TOPOPLOTS = false;   % <-- set true to also generate per-file topoplots

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir,'io'), fullfile(this_dir,'utils'));

fprintf('\n══════════════════════════════════════════════════════\n');
fprintf('  batch_run_subjects\n');
fprintf('  Root     : %s\n', RECORDINGS_ROOT);
fprintf('  Subjects : %s\n', strjoin(SUBJECTS, ', '));
fprintf('══════════════════════════════════════════════════════\n');

for si = 1:numel(SUBJECTS)
    subj_dir = fullfile(RECORDINGS_ROOT, SUBJECTS{si});
    if ~isfolder(subj_dir)
        fprintf('\n[batch_run_subjects] subject folder not found, skipping: %s\n', subj_dir);
        continue;
    end

    day_entries = dir(subj_dir);
    day_entries = day_entries([day_entries.isdir] & ~startsWith({day_entries.name}, '.'));

    n_done = 0;
    for di = 1:numel(day_entries)
        eval_dir = fullfile(subj_dir, day_entries(di).name, 'evaluation');
        if ~isfolder(eval_dir), continue; end
        fprintf('\n\n########## %s / %s ##########\n', SUBJECTS{si}, day_entries(di).name);
        run_subject_analysis(eval_dir, SAVE_PER_FILE_TOPOPLOTS);
        n_done = n_done + 1;
    end

    if n_done == 0
        fprintf('\n[batch_run_subjects] no evaluation/ day found for subject %s\n', SUBJECTS{si});
    end
end

fprintf('\n\n══════ batch_run_subjects complete (%d subjects) ══════\n\n', numel(SUBJECTS));
