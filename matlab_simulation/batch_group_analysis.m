%% BATCH_GROUP_ANALYSIS  Group-level (cross-subject) aggregation for a
%   hardcoded list of subjects, no GUI prompts.
%
%   Edit RECORDINGS_ROOT and SUBJECTS below, then run this script directly.
%   Keep SUBJECTS in sync with batch_run_subjects.m if you want the same
%   cohort in both the per-subject and group-level analyses.
%
%   Calls, both restricted to SUBJECTS via their subject-filter argument:
%     main_group_analysis  — accuracy/counterfactual/fusion/ROC/cluster-test
%                             group figures + group_summary.mat
%     group_topo_erders     — cross-subject averaged ERD/ERS topoplots

RECORDINGS_ROOT = '/home/paolo/bci_vr_ws/recordings';
SUBJECTS = {'a1', 'a2', 'a3'};   % <-- edit this list
SHOW_FIGURES = false;

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir, 'utils'));
addpath(fullfile(this_dir, '..', 'analysis_gdf'));

fprintf('\n══════════════════════════════════════════════════════\n');
fprintf('  batch_group_analysis\n');
fprintf('  Root     : %s\n', RECORDINGS_ROOT);
fprintf('  Subjects : %s\n', strjoin(SUBJECTS, ', '));
fprintf('══════════════════════════════════════════════════════\n');

main_group_analysis(RECORDINGS_ROOT, SUBJECTS, SHOW_FIGURES);
group_topo_erders(RECORDINGS_ROOT, SUBJECTS, SHOW_FIGURES);

fprintf('\n\n══════ batch_group_analysis complete (%d subjects) ══════\n\n', numel(SUBJECTS));
