function [integrated_prob, mask] = applyIntegration(integratorCfg, artifact, qda_prob_cvsa, qda_prob_mi, event, event_start, rejection, paradigm, fs)
cfPOS = event.POS(event.TYP == event_start);
cfDUR = event.DUR(event.TYP == event_start);
ntrial = length(cfPOS);

integrated_prob = ones(size(qda_prob_cvsa, 1), 2) .* cell2mat(integratorCfg.init_val);
mask = nan(size(qda_prob_cvsa, 1), 1);

if all(integratorCfg.type == 'Buffer')
    bufferSize = integratorCfg.bufferSize;
elseif all(integratorCfg.type == 'Exponential')
    alpha = integratorCfg.alpha;
end

for idx_trial =1:ntrial
    start_trial = cfPOS(idx_trial);
    end_trial = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;

    nsamples_trial = end_trial - start_trial + 1;
    if strcmp(paradigm, 'mi')
        merged_prob = qda_prob_mi(start_trial:end_trial,:);
    elseif strcmp(paradigm, 'cvsa')
        merged_prob = qda_prob_cvsa(start_trial:end_trial,:);
    elseif strcmp(paradigm, 'hybrid')
        c_cvsa = qda_prob_cvsa(start_trial:end_trial,:);
        c_mi = qda_prob_mi(start_trial:end_trial,:);

        t = (0 : nsamples_trial-1)' / fs;
        t_alpha = zeros(nsamples_trial, 1);
        idx_decay = t <= 2.5; 
        t_alpha(idx_decay) = 0.5 * (1.0 + cos(pi * t(idx_decay) / 2.5));
        tempered_priors = c_cvsa .^ t_alpha;
        tempered_priors = tempered_priors ./ sum(tempered_priors, 2);

        merged_prob = c_mi .* tempered_priors;
        merged_prob = merged_prob ./ sum(merged_prob, 2);
    end
    c_artifact = artifact(start_trial:end_trial,:);
    c_mask = [];
    c_integrated = ones(nsamples_trial + 1, 2) .* cell2mat(integratorCfg.init_val);

    for idx_sample = 1:nsamples_trial
        if c_artifact(idx_sample) == 0 % no artifact
            % integration
            if all(integratorCfg.type == 'Buffer')
                if merged_prob(idx_sample, 1) >= rejection 
                    if integratorCfg.increment_type == 0 % HARDINCREMENT
                        inc = 1/bufferSize;
                    elseif integratorCfg.increment_type == 1 % SOFTINCREMENT
                        vel = (merged_prob(idx_sample, 1) - 0.5) * 2.0 * integratorCfg.k_gain;
                        vel = min(vel, 1);
                        inc = 1/bufferSize * vel;
                    end
                    c_integrated(idx_sample+1, 1)=min(c_integrated(idx_sample, 1) + inc, 1);
                    c_integrated(idx_sample+1, 2)=max(c_integrated(idx_sample, 2) - inc, 0);
                    
                else
                    if integratorCfg.increment_type == 0
                        inc = -1/bufferSize;
                    elseif integratorCfg.increment_type == 1
                        vel = (merged_prob(idx_sample, 2) - 0.5) * 2.0 * integratorCfg.k_gain;
                        vel = min(vel, 1);
                        inc = -1/bufferSize * vel;
                    end
                    c_integrated(idx_sample+1, 1)=max(c_integrated(idx_sample, 1) + inc,0);
                    c_integrated(idx_sample+1, 2)=max(c_integrated(idx_sample, 2) - inc,0);
                end
            elseif all(integratorCfg.type == 'Exponential')
                if c_qda_prob(idx_sample, 1) >= rejection
                    tmp = [1 0];
                else
                    tmp = [0 1];
                end
                c_integrated(idx_sample+1,:) = c_integrated(idx_sample,:) * alpha + (1-alpha) * tmp;
            else
                disp('   ERROR')
            end
            c_mask = cat(1, c_mask, 1);
        else
            % no integration
            c_mask = cat(1, c_mask, 0);
            c_integrated(idx_sample+1,:) = c_integrated(idx_sample,:);

        end
    end

    integrated_prob(start_trial:end_trial,:) = c_integrated(2:end,:);
    mask(start_trial:end_trial,:) = c_mask;
end
end