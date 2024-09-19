function [trialStart, trialStop, fixStart, fixStop, ck, tk] = extract_trial_infoo(signal, fixPOS, fixDUR, cfPOS, cfDUR, cueTYP, n_trial)
disp('   [INFO] extracting trials')
ck = zeros(size(signal,1), 1);
tk = zeros(size(signal,1), 1);
trialStart = nan(n_trial, 1);
trialStop  = nan(n_trial, 1);

fixStart = nan(n_trial, 1);
fixStop  = nan(n_trial, 1);

for idx_tr = 1:n_trial
    cstart = fixPOS(idx_tr);
    cstop = cfPOS(idx_tr) + cfDUR(idx_tr) - 1;

    ck(cstart:cstop) = cueTYP(idx_tr);
    tk(cstart:cstop) = idx_tr;

    trialStart(idx_tr) = cstart;
    trialStop(idx_tr) = cstop;

    fixStart(idx_tr) = fixPOS(idx_tr);
    fixStop(idx_tr) = fixPOS(idx_tr) + fixDUR(idx_tr) - 1;
end
end