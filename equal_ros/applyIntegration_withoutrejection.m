function [integrated_prob, mask] = applyIntegration_withoutrejection(integratorCfg, artifact, gmm_prob, qda_prob, event, event_start, gmm_classes)
cfPOS = event.POS(event.TYP == event_start);
cfDUR = event.DUR(event.TYP == event_start);
ntrial = length(cfPOS);
ic_index = find(gmm_classes == integratorCfg.ic_class_label);

integrated_prob = ones(size(gmm_prob, 1), 2) .* cell2mat(integratorCfg.init_val);
mask = zeros(size(gmm_prob, 1), 1);

if all(integratorCfg.type == 'Buffer')
    bufferSize = integratorCfg.bufferSize;
elseif all(integratorCfg.type == 'Exponential')
    alpha = integratorCfg.alpha;
end

for idx_trial =1:ntrial
    start_trial = cfPOS(idx_trial);
    end_trial = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;

    nsamples_trial = end_trial - start_trial + 1;
    c_gmm_prob = gmm_prob(start_trial:end_trial,:);
    c_qda_prob = qda_prob(start_trial:end_trial,:);
    c_artifact = artifact(start_trial:end_trial,:);
    c_mask = [];

    c_integrated = ones(nsamples_trial + 1, 2) .* cell2mat(integratorCfg.init_val);
    for idx_sample = 1:nsamples_trial
        if c_artifact(idx_sample) == 0
            % integration
            if all(integratorCfg.type == 'Buffer')
                merged_prob = (1-c_gmm_prob(idx_sample,ic_index)) * 0.5 + c_gmm_prob(idx_sample,ic_index) * c_qda_prob(idx_sample, :);
                if merged_prob(1) >= 0.5
                    vel = (merged_prob(1) - 0.5) * 2.0 * integratorCfg.k_gain;
                    vel = min(vel, 1);
                    inc = 1/bufferSize * vel;
                    c_integrated(idx_sample+1, 1)=min(c_integrated(idx_sample, 1) + inc, 1);
                    c_integrated(idx_sample+1, 2)=max(c_integrated(idx_sample, 2) - inc, 0);
                else
                    vel = (merged_prob(2) - 0.5) * 2.0 * integratorCfg.k_gain;
                    vel = min(vel, 1);
                    inc = -1/bufferSize * vel;
                    c_integrated(idx_sample+1, 1)=max(c_integrated(idx_sample, 1) + inc,0);
                    c_integrated(idx_sample+1, 2)=max(c_integrated(idx_sample, 2) - inc,0);
                end
            elseif all(integratorCfg.type == 'Exponential')
                if c_qda_prob(idx_sample, 1) >= 0.5
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