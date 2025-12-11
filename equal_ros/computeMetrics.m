function [accuracy, number, time] = computeMetrics(integratorCfg, artifact, gmm_prob, qda_prob, event, event_start, gmm_classes, task_classes)
cfPOS = event.POS(event.TYP == event_start);
cfDUR = event.DUR(event.TYP == event_start);
cueTYP = event.TYP(ismember(event.TYP, task_classes));
hitmissTYP = event.TYP(ismember(event.TYP, [897, 898, 899]));
ntrial = length(cfPOS);
ic_index = find(gmm_classes == integratorCfg.ic_class_label);

% variables for counting
n_correct_gmm = 0;         n_total_gmm = 0;
n_correct_traditional = 0; n_total_traditional = 0;
n_discarded_artifact = 0;  n_total_all = 0;

% variables for time
time_hit = 0;
time_miss = 0;
time_timeout = 0;

for idx_trial = 1:ntrial
    start_trial = cfPOS(idx_trial);
    end_trial = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;

    nsamples_trial = end_trial - start_trial + 1;
    c_gmm_prob = gmm_prob(start_trial:end_trial,:);
    c_qda_prob = qda_prob(start_trial:end_trial,:);
    c_artifact = artifact(start_trial:end_trial,:);

    prob_index = find(task_classes == cueTYP(idx_trial));

    % counting for accuracy and number variable
    for idx_sample = 1:nsamples_trial
        if c_artifact(idx_sample) == 0
            if c_qda_prob(idx_sample, prob_index) >= 0.5
                if c_gmm_prob(idx_sample,ic_index) >= integratorCfg.ic_threshold
                    n_correct_gmm = n_correct_gmm + 1;
                    n_total_gmm = n_total_gmm + 1;
                end
                n_correct_traditional = n_correct_traditional + 1;
                n_total_traditional = n_total_traditional + 1;
            else
                if c_gmm_prob(idx_sample,ic_index) >= integratorCfg.ic_threshold
                    n_total_gmm = n_total_gmm + 1;
                end
                n_total_traditional = n_total_traditional + 1;
            end
        else
            n_discarded_artifact = n_discarded_artifact + 1;
        end
    end
    n_total_all = n_total_all + nsamples_trial;

    % time variables
    if hitmissTYP(idx_trial) == 897
        time_hit = time_hit + nsamples_trial;
    elseif hitmissTYP(idx_trial) == 898
        time_miss = time_miss + nsamples_trial;
    elseif hitmissTYP(idx_trial) == 899
        time_timeout = time_timeout + nsamples_trial;
    else
        disp('ERROR in type for time hit-miss-timeout')
    end
end

accuracy.sample.icnic = n_correct_gmm/n_total_gmm;
accuracy.sample.traditional   = n_correct_traditional/n_total_traditional;
accuracy.trial.hit = sum(event.TYP == 897)/ntrial;
accuracy.trial.miss = sum(event.TYP == 898)/ntrial;
accuracy.trial.timeout = sum(event.TYP == 899)/ntrial;

number.sample.icnic.correct = n_correct_gmm;
number.sample.icnic.total = n_total_gmm;
number.sample.traditional.correct = n_correct_traditional;
number.sample.traditional.total = n_total_traditional;
number.sample.discrded = n_discarded_artifact;
number.sample.all = n_total_all;
number.trial.hit = sum(event.TYP == 897);
number.trial.miss = sum(event.TYP == 898);
number.trial.timeout = sum(event.TYP == 899);
number.trial.ntrial = ntrial;

time.hit = time_hit;
time.miss = time_miss;
time.timeout = time_timeout;

end