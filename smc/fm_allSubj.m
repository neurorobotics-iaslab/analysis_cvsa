%% results over all subjects
clc; clearvars;
% load the files
classes = [730 731];
threshold_pseudo_online = 0.7;
[filenames, pathname] = uigetfile('*.mat', 'Load RAW probabilities', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subjects = cellfun(@(x) erase(x, '.mat'), filenames, 'UniformOutput', false);
nfiles = length(filenames);

raw = cell(1, nfiles);
parameters = cell(1, nfiles);
min_trial_dur = inf;
for idx_file = 1:nfiles
    c_val = load([pathname filenames{idx_file}]);
    raw{idx_file} = c_val.raw;
    parameters{idx_file} = c_val.parameters;

    if size(c_val.raw.test_cf, 2) < min_trial_dur
        min_trial_dur = size(c_val.raw.test_cf, 2);
    end
end
mincueDur = 512;

%% for each subject perform: the integration and the AUC
integrated = cell(1, nfiles);
start_integrator = 0.5;
prob_integ = nan(nfiles +1, min_trial_dur + 1, 2); % contains all the data integrated mean across trials
auc = nan(nfiles+1, 2); % contains auc for merge - classical
acc = nan(nfiles + 1, 2);
nfeatures = nan(nfiles +1, 4);
intersect_point = nan(nfiles+1, 1);
pseudo_online_acc = nan(nfiles +1, 2);
pseudo_online_time = nan(nfiles+1, 2);
for idx_file = 1:nfiles
    c_raw = raw{idx_file};
    c_parameters = parameters{idx_file};
    ntrial = size(c_raw.test_cf, 1);
    intersect_point(idx_file) = c_parameters.intersect;

    nfeatures(idx_file, 1) = length(c_parameters.features.shift.channels);
    nfeatures(idx_file, 2) = length(c_parameters.features.sus.channels);
    nfeatures(idx_file, 3) = length(c_parameters.features.cf.channels);
    count_commonFeatures = 0;
    for idx_feature = 1:length(c_parameters.features.cf.channels)
        c_feature_chs = c_parameters.features.cf.channels(idx_feature);
        c_feature_band = c_parameters.features.cf.bands(idx_feature,:);
        if ismember(c_feature_chs, c_parameters.features.shift.channels)
            idx = find(strcmp(c_parameters.features.shift.channels, c_feature_chs{1}));
            bands_subset = c_parameters.features.shift.bands(idx, :);
            if ismember(c_feature_band, bands_subset, 'rows')
                count_commonFeatures = count_commonFeatures + 1;
            end
        elseif ismember(c_feature_chs, c_parameters.features.sus.channels)
            idx = find(strcmp(c_parameters.features.sus.channels, c_feature_chs{1}));
            bands_subset = c_parameters.features.sus.bands(idx, :);
            if ismember(c_feature_band, bands_subset, 'rows')
                count_commonFeatures = count_commonFeatures + 1;
            end
        end
    end
    nfeatures(idx_file, 4) = count_commonFeatures;
    
    % compute the raw prob for each trial -> merge case
    % save values for AUC
    y_pred_merge = nan(ntrial * min_trial_dur, 1);
    y_pred_classic = nan(ntrial * min_trial_dur, 1);
    y_true = nan(ntrial * min_trial_dur, 1);
    raw_merge = nan(size(c_raw.test_cf));
    integ_merge = ones(ntrial, min_trial_dur) * start_integrator;
    integ_cf = ones(ntrial, min_trial_dur) * start_integrator;
    integ_merge_soft = ones(ntrial, min_trial_dur) * start_integrator;
    integ_cf_soft = ones(ntrial, min_trial_dur - mincueDur) * start_integrator;
%     alpha = c_parameters.alpha;
    alpha = 0.98;
    for idx_trial = 1:ntrial
        % compute the merge
        c_shift = squeeze(c_raw.test_shift(idx_trial,:,:));
        c_sus = squeeze(c_raw.test_sus(idx_trial,:,:));
        c_g = c_parameters.g; c_f = c_parameters.f;
        c_raw_merge = c_shift .* c_f + c_sus .* c_g;
        raw_merge(idx_trial,:,:) = c_raw_merge;

        % for AUC
        c_true = repmat(c_parameters.typ(idx_trial), min_trial_dur, 1);
        c_true(c_true == classes(1)) = 0;
        c_true(c_true == classes(2)) = 1;
        y_true((idx_trial-1)*min_trial_dur+1:idx_trial*min_trial_dur) = c_true;
        y_pred_merge((idx_trial-1)*min_trial_dur+1:idx_trial*min_trial_dur) = raw_merge(idx_trial, 1:min_trial_dur, 2);
        y_pred_classic((idx_trial-1)*min_trial_dur+1:idx_trial*min_trial_dur) = c_raw.test_cf(idx_trial,1:min_trial_dur,2);

        % for integration
        idx_class = find(c_parameters.typ(idx_trial) == classes);
        hard_merge = ones(min_trial_dur, 1);
        hard_merge(c_raw_merge(:,idx_class) < 0.5) = 0;
        hard_merge(c_raw_merge(:,idx_class) == 0.5) = 0.5;
        hard_cf = ones(min_trial_dur, 1);
        hard_cf(c_raw.test_cf(idx_trial,:,idx_class) < 0.5) = 0;
        hard_cf(c_raw.test_cf(idx_trial,:,idx_class) == 0.5) = 0.5;
        for idx_sample = 1:min_trial_dur
            integ_merge(idx_trial, idx_sample + 1) = integ_merge(idx_trial, idx_sample) * alpha + ...
                (1-alpha) * hard_merge(idx_sample);
            integ_cf(idx_trial, idx_sample + 1) = integ_cf(idx_trial, idx_sample) * alpha + ...
                (1-alpha) * hard_cf(idx_sample);
            integ_merge_soft(idx_trial, idx_sample + 1) = integ_merge_soft(idx_trial, idx_sample) * alpha + ...
                (1-alpha) * c_raw_merge(idx_sample,idx_class);
            if idx_sample > mincueDur
                c_idx_sample = idx_sample - mincueDur;
                integ_cf_soft(idx_trial, c_idx_sample + 1) = integ_cf_soft(idx_trial, c_idx_sample) * alpha + ...
                    (1-alpha) * c_raw.test_cf(idx_trial,idx_sample,idx_class);
            end
%             integ_cf_soft(idx_trial, idx_sample + 1) = integ_cf_soft(idx_trial, idx_sample) * alpha + ...
%                     (1-alpha) * c_raw.test_cf(idx_trial,idx_sample,idx_class);
        end
    end
    c_raw.test_merge = raw_merge;
    raw{idx_file} = c_raw;

    prob_integ(idx_file, :, 1) = mean(integ_merge, 1);
    prob_integ(idx_file,:,2) = mean(integ_cf, 1);

    time_method = nan(ntrial, 2);
    acc_method = zeros(ntrial, 2);
    for i = 1:ntrial
        t_merge_correct = find(integ_merge_soft(i,:) >= threshold_pseudo_online, 1, 'first');
        t_merge_wrong = find(integ_merge_soft(i,:) <= 1 -threshold_pseudo_online, 1, 'first');
        t_cf_correct = find(integ_cf_soft(i,:) >= threshold_pseudo_online, 1, 'first');
        t_cf_wrong = find(integ_cf_soft(i,:) <= 1- threshold_pseudo_online, 1, 'first');

        if ~isempty(t_merge_correct)
            if ~isempty(t_merge_wrong)
                if t_merge_wrong > t_merge_correct
                    time_method(i, 1) = t_merge_correct / 512;
                    acc_method(i, 1) = 1;
                end
            else
                time_method(i, 1) = t_merge_correct / 512;
                acc_method(i, 1) = 1;
            end
        end

        if ~isempty(t_cf_correct)
            if ~isempty(t_cf_wrong)
                if t_cf_wrong > t_cf_correct
                    time_method(i, 2) = t_cf_correct / 512;
                    acc_method(i, 2) = 1;
                end
            else
                time_method(i, 2) = t_cf_correct / 512;
                acc_method(i, 2) = 1;
            end
        end

    end
    % pseudo online
    pseudo_online_acc(idx_file, 1) = mean(acc_method(:,1));
    pseudo_online_time(idx_file, 1) = nanmean(time_method(:,1));
    pseudo_online_acc(idx_file, 2) = mean(acc_method(:,2));
    pseudo_online_time(idx_file, 2) = nanmean(time_method(:,2));

    % compute the AUC
    [~, ~, ~, auc(idx_file, 1)] = perfcurve(y_true, y_pred_merge, 1);
    [~, ~, ~, auc(idx_file, 2)] = perfcurve(y_true, y_pred_classic, 1);

    tmp = y_pred_merge;
    tmp(tmp < 0.5) = 0; tmp(tmp> 0.5) = 1;
    acc(idx_file, 1) = sum(tmp == y_true) / length(y_true);
    tmp = y_pred_classic;
    tmp(tmp < 0.5) = 0; tmp(tmp> 0.5) = 1;
    acc(idx_file, 2) = sum(tmp == y_true) / length(y_true);

end
pseudo_online_acc(end, :) = mean(pseudo_online_acc(1:end-1,:), 1);
pseudo_online_time(end, :) = nanmean(pseudo_online_time(1:end-1,:), 1);
nfeatures(end, :) = mean(nfeatures(1:end-1,:), 1);
auc(end,:) = mean(auc(1:end-1,:), 1);
prob_integ(end,:,1) = mean(prob_integ(1:end-1,:,1), 1);
prob_integ(end,:,2) = mean(prob_integ(1:end-1,:,2), 1);
intersect_point(end) = mean(intersect_point(1:end-1));
acc(end,:) = mean(acc(1:end-1,:), 1);

% acc pseudo online
figure();
bar(pseudo_online_acc);
x_label = subjects; x_label{end+1} = 'average';
xlabel('Subjects');
ylabel('Accuracy');
legend({'Shift-sustained', 'Traditional'}, 'Location', 'northeast');
xticklabels(x_label);  % Optional: custom group labels
grid on;

% acc pseudo online
figure();
bar(pseudo_online_time);
x_label = subjects; x_label{end+1} = 'average';
xlabel('Subjects');
ylabel('Time');
legend({'Shift-sustained', 'Traditional'}, 'Location', 'northeast');
xticklabels(x_label);  % Optional: custom group labels
grid on;

% show the auc as bar plot (each subject 2 auc) + mean over subjects
figure();
bar(auc);
x_label = subjects; x_label{end+1} = 'average';
xlabel('Subjects');
ylabel('AUC');
legend({'Shift-sustained', 'Traditional'}, 'Location', 'northeast');
xticklabels(x_label);  % Optional: custom group labels
grid on;

%% plot as difference
diff_auc = nan(nfiles+1,1);
diff_acc = nan(nfiles+1,1);
for i=1:nfiles+1
    diff_acc(i) = pseudo_online_acc(i, 1) - pseudo_online_acc(i,2);
    diff_auc(i) = auc(i, 1) - auc(i,2);
end
auc1 = auc(1:9,1); auc2 = auc(1:9,2);
[~, p_auc, ~, ~] = ttest(auc2, auc1);

pseudo_online_acc1 = pseudo_online_acc(1:9,1); pseudo_online_acc2 = pseudo_online_acc(1:9,2);
[~, p_pseudo_online_acc, ci, ~] = ttest(pseudo_online_acc1, pseudo_online_acc2);


figure();
subplot(1,2,1)
bar(diff_auc)
x_label = subjects; x_label{end+1} = 'average';
xlabel('Subjects');
ylabel('AUC');
xticklabels(x_label);  % Optional: custom group labels

subplot(1,2,2)
bar(diff_acc)
x_label = subjects; x_label{end+1} = 'average';
xlabel('Subjects');
ylabel('AUC');
xticklabels(x_label);  % Optional: custom group labels
