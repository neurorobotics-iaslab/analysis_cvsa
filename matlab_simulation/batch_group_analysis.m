%% BATCH_GROUP_ANALYSIS  Group-level (cross-subject) aggregation for a
%   hardcoded list of subjects, no GUI prompts.
%
%   Edit RECORDINGS_ROOT and SUBJECTS below, then run this script directly.
%   Keep SUBJECTS in sync with batch_run_subjects.m if you want the same
%   cohort in both the per-subject and group-level analyses.
%
%   Also edit WELL_SUBJECTS/BAD_SUBJECTS (cellstr, e.g. {'a1','a3'}) to run
%   two EXTRA group-level reports restricted to those manually-assigned
%   "strong"/"weak" performer subsets, on top of the ALL-subjects report --
%   a subject left out of both is still included in the ALL report. Leave
%   either list empty ({}) to skip that report. main_group_analysis tags
%   every output filename with the group ('01_all_...', '01_well_...',
%   '01_bad_...', 'group_summary_all.mat', etc.) so all three coexist under
%   the same <root>/group_analysis/ folder. If both WELL_SUBJECTS and
%   BAD_SUBJECTS are non-empty, well_vs_bad_comparison.m additionally runs
%   an unpaired comparison (16_well_vs_bad_comparison.svg) testing whether
%   real accuracy / counterfactual advantage / fusion advantage / ROC AUC
%   actually differ between the two groups.
%
%   Calls:
%     main_group_analysis     — accuracy/counterfactual/fusion/ROC/cluster-
%                                test group figures + group_summary_<tag>.mat,
%                                once per non-empty group (ALL/well/bad)
%     well_vs_bad_comparison   — unpaired well-vs-bad comparison figure,
%                                only if both well and bad groups ran
%     group_topo_erders        — cross-subject averaged ERD/ERS topoplots
%                                (ALL subjects only)
%     group_csp_slda_importance — cross-subject averaged CSP/sLDA channel
%                                importance topoplots + feature-selection
%                                consistency, per band (ALL subjects only)

RECORDINGS_ROOT = '/home/paolo/bci_vr_ws/recordings';
SUBJECTS = {'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', ...
            'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', ...
            'c1', 'c2', 'c3', 'c4', 'c5', 'c6'};   % <-- edit this list

WELL_SUBJECTS = {'a1', 'a2', 'a3', 'a5', 'a6', 'a7', 'a8', 'a9', ...
                 'b3', 'b4', 'b6', 'b7', 'b8', ...
                 'c2', 'c4 ', 'c5', 'c6'};   % e.g. {'a1', 'a3', 'b2'} -- "strong" performers
BAD_SUBJECTS  = {'a4',  ...
                 'b1', 'b2', 'b5', 'b9', ...
                 'c1', 'c3'};   % e.g. {'a2', 'a4', 'c1'} -- "weak" performers

SHOW_FIGURES = false;

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir, fullfile(this_dir, 'utils'));
addpath(fullfile(this_dir, '..', 'analysis_gdf'));

fprintf('\n══════════════════════════════════════════════════════\n');
fprintf('  batch_group_analysis\n');
fprintf('  Root     : %s\n', RECORDINGS_ROOT);
fprintf('  Subjects : %s\n', strjoin(SUBJECTS, ', '));
fprintf('══════════════════════════════════════════════════════\n');

main_group_analysis(RECORDINGS_ROOT, SUBJECTS, SHOW_FIGURES, 'all');

gs_well = [];
if ~isempty(WELL_SUBJECTS)
    try
        gs_well = main_group_analysis(RECORDINGS_ROOT, WELL_SUBJECTS, SHOW_FIGURES, 'well');
    catch ME
        if strcmp(ME.identifier, 'main_group_analysis:nodata')
            fprintf('[batch_group_analysis] Skipping "well" group: %s\n', ME.message);
        else
            rethrow(ME);
        end
    end
end

gs_bad = [];
if ~isempty(BAD_SUBJECTS)
    try
        gs_bad = main_group_analysis(RECORDINGS_ROOT, BAD_SUBJECTS, SHOW_FIGURES, 'bad');
    catch ME
        if strcmp(ME.identifier, 'main_group_analysis:nodata')
            fprintf('[batch_group_analysis] Skipping "bad" group: %s\n', ME.message);
        else
            rethrow(ME);
        end
    end
end

if ~isempty(gs_well) && ~isempty(gs_bad)
    well_vs_bad_comparison(gs_well, gs_bad, fullfile(RECORDINGS_ROOT, 'group_analysis'), SHOW_FIGURES);
end

group_topo_erders(RECORDINGS_ROOT, SUBJECTS, SHOW_FIGURES);
group_csp_slda_importance(RECORDINGS_ROOT, SUBJECTS, SHOW_FIGURES);

fprintf('\n\n══════ batch_group_analysis complete (%d subjects) ══════\n\n', numel(SUBJECTS));
