clear; clc; close all;

% Data
data_path = '/home/paolo/cvsa/ic_cvsa_ws/record_mi/results_graz/qda';
files = dir(fullfile(data_path, '*.mat'));
n_subj = length(files);

if n_subj == 0
    error('Nessun file .mat trovato nella cartella selezionata.');
end
% --- initialization --
sub_s_trad = zeros(n_subj, 1);
sub_s_my   = zeros(n_subj, 1);
sub_t_act  = zeros(n_subj, 1);
sub_t_noto = zeros(n_subj, 1);
sub_t_rest = zeros(n_subj, 1);
sub_t_timeout = zeros(n_subj, 1);
std_s_trad = zeros(n_subj, 1);
std_s_my   = zeros(n_subj, 1);
std_t_act  = zeros(n_subj, 1);
std_t_noto = zeros(n_subj, 1);
std_t_rest = zeros(n_subj, 1);
std_t_timeout = zeros(n_subj, 1);

for i = 1:n_subj
    load(fullfile(files(i).folder, files(i).name), 'Database');
    
    run_s_trad = [Database.ActAccQDATrad];
    run_s_my   = [Database.ActAccQDAMy];
    run_t_act  = [Database.ActAccTrial];
    run_t_noto = [Database.ActAccNo_timeout];
    run_t_rest = [Database.RestAccTrial];
    run_t_timouet = [Database.PercTimeout];
    
    % mean
    sub_s_trad(i) = mean(run_s_trad, 'omitnan');
    sub_s_my(i)   = mean(run_s_my, 'omitnan');
    sub_t_act(i)  = mean(run_t_act, 'omitnan');
    sub_t_noto(i) = mean(run_t_noto, 'omitnan');
    sub_t_rest(i) = mean(run_t_rest, 'omitnan');
    sub_t_timeout(i) = mean(run_t_timouet, 'omitnan');
    
    % std
    std_s_trad(i) = std(run_s_trad, 0, 'omitnan');
    std_s_my(i)   = std(run_s_my, 0, 'omitnan');
    std_t_act(i)  = std(run_t_act, 0, 'omitnan');
    std_t_noto(i) = std(run_t_noto, 0, 'omitnan');
    std_t_rest(i) = std(run_t_rest, 0, 'omitnan');
    std_t_timeout(i) = std(run_t_timouet, 0, 'omitnan');
end

disp('===================================================================');
disp('                  ONLINE PERFORMANCE METRICS                       ');
disp('===================================================================');

for i = 1:n_subj
    fprintf('\n--- SUBJECT %d ---\n', i);
    fprintf('  Sample Acc (Standard) : %5.1f +/- %4.1f %%\n', sub_s_trad(i), std_s_trad(i));
    fprintf('  Sample Acc (Gated)    : %5.1f +/- %4.1f %%\n', sub_s_my(i), std_s_my(i));
    fprintf('  Trial Acc (Active)    : %5.1f +/- %4.1f %%\n', sub_t_act(i), std_t_act(i));
    fprintf('  Trial Acc (No TimeOut): %5.1f +/- %4.1f %%\n', sub_t_noto(i), std_t_noto(i));
    fprintf('  Trial Acc (Rest)      : %5.1f +/- %4.1f %%\n', sub_t_rest(i), std_t_rest(i));
    fprintf('  Timeout (percentual)  : %5.1f +/- %4.1f %%\n', sub_t_timeout(i), std_t_timeout(i));
end

disp(' ');
disp('===================================================================');
fprintf('                GRAND AVERAGE (N = %d SUBJECTS)                    \n', n_subj);
disp('===================================================================');
fprintf('  Sample Acc (Standard) : %5.1f +/- %4.1f %%\n', mean(sub_s_trad), std(sub_s_trad));
fprintf('  Sample Acc (Gated)    : %5.1f +/- %4.1f %%\n', mean(sub_s_my), std(sub_s_my));
fprintf('  Trial Acc (Active)    : %5.1f +/- %4.1f %%\n', mean(sub_t_act), std(sub_t_act));
fprintf('  Trial Acc (No TimeOut): %5.1f +/- %4.1f %%\n', mean(sub_t_noto), std(sub_t_noto));
fprintf('  Trial Acc (Rest)      : %5.1f +/- %4.1f %%\n', mean(sub_t_rest), std(sub_t_rest));
fprintf('  Timeout (percentual)  : %5.1f +/- %4.1f %%\n', mean(sub_t_timeout), std(std_t_timeout));
disp('===================================================================');
