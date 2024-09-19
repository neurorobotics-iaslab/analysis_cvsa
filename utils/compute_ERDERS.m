function [ERD, minDur, minDurFix] = compute_ERDERS(signal, trialStart, trialStop, fixStart, fixStop, nchannels, n_trial, ck, normalization)
disp('      [INFO] computing ERD/ERS')
% extract useful data
minDur = min(trialStop-trialStart);
trial_data = nan(minDur, nchannels, n_trial);
tCk = zeros(n_trial, 1);
for idx_tr = 1:n_trial
    cstart = trialStart(idx_tr);
    cstop = trialStart(idx_tr) + minDur - 1;
    trial_data(:,:,idx_tr) = signal(cstart:cstop,:);
    tCk(idx_tr) = unique(ck(cstart:cstop));
end

% extract baseline
minDurFix = min(fixStop - fixStart);
fix_data = nan(minDurFix, nchannels, n_trial);
for idx_tr = 1:n_trial
    cstart = fixStart(idx_tr);
    cstop = fixStart(idx_tr) + minDurFix - 1;
    fix_data(:,:,idx_tr) = signal(cstart:cstop,:);
end


%% ERD/ERS
if normalization
    disp('         [proc] normalize it with the baseline (fixation)')
    baseline = repmat(mean(fix_data), [size(trial_data, 1) 1 1]);
    ERD = log(trial_data ./ baseline);
else
    ERD = log(trial_data);
end
end