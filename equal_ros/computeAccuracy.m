function [accuracy, number] = computeAccuracy(integratorCfg, artifact, gmm_prob, qda_prob, event, event_start, gmm_classes, task_classes)
cfPOS = event.POS(event.TYP == event_start);
cfDUR = event.DUR(event.TYP == event_start);
cueTYP = event.TYP(ismember(event.TYP, task_classes));
ntrial = length(cfPOS);
ic_index = find(gmm_classes == integratorCfg.ic_class_label);

for idx_trial =1:ntrial
    start_trial = cfPOS(idx_trial);
    end_trial = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;

    nsamples_trial = end_trial - start_trial + 1;
    c_gmm_prob = gmm_prob(start_trial:end_trial,:);
    c_qda_prob = qda_prob(start_trial:end_trial,:);
    c_artifact = artifact(start_trial:end_trial,:);

    prob_index = find(task_classes == cueTYP(idx_trial));

    correct_gmm = 0; total_gmm = 0;
    correct_all = 0; total_all = 0;

    for idx_sample = 1:nsamples_trial
        if c_qda_prob(idx_sample, prob_index) >= 0.5 
            if c_gmm_prob(idx_sample,ic_index) >= integratorCfg.ic_threshold && c_artifact(idx_sample) == 0
                correct_gmm = correct_gmm + 1;
            end
            correct_all = correct_all + 1;
        else
            total_all = total_all + 1;
            total_gmm = total_gmm + 1;
        end
    end
end

accuracy.sample.icnic = correct_gmm/total_gmm;
accuracy.sample.all   = correct_all/total_all;
accuracy.trial.hit = length(event.TYP == 897)/ntrial;
accuracy.trial.miss = length(event.TYP == 898)/ntrial;
accuracy.trial.timeout = length(event.TYP == 899)/ntrial;

number.sample.icnic.correct = correct_gmm;
number.sample.icnic.total = total_gmm;
number.sample.all.correct = correct_all;
number.sample.all.total = total_all;
number.trial.hit = length(event.TYP == 897);
number.trial.miss = length(event.TYP == 898);
number.trial.timeout = length(event.TYP == 899);
number.trial.ntrial = ntrial;

end