%% extract cf, fix, cue and all the trial
% OUTPUT: 
%        - signal struct contains all the data in samples x channels x trial:
%               signal.trial: signal of each trial
%               signal.fix: data during only fization
%               signal.cue: data during only cue
%               signal.cf:  data during only cf
%        - typ is a vector 1 x trial contain the value of the cue
function signal = extract_all(s, header, classes, fix_event, cf_event, startTrial_event)
cueDUR = header.DUR(ismember(header.TYP, classes));
cueTYP = header.TYP(ismember(header.TYP, classes));
cuePOS = header.POS(ismember(header.TYP, classes));

fixPOS = header.POS(header.TYP == fix_event);
fixDUR = header.DUR(header.TYP == fix_event);

cfPOS = header.POS(header.TYP == cf_event);
cfDUR = header.DUR(header.TYP == cf_event);

trialPOS = header.POS(header.TYP == startTrial_event);
trialDUR = header.DUR(header.TYP == startTrial_event);

minDURfix = min(fixDUR);
minDURcue = min(cueDUR);
minDURcf = min(cfDUR);
minDURtrial = min(trialDUR);
ntrial = length(cueTYP);
nchannels = size(s, 2);

signal.trial = nan(minDURtrial, nchannels, ntrial);
signal.cue = nan(minDURcue, nchannels, ntrial);
signal.fix = nan(minDURfix, nchannels, ntrial);
signal.cf  = nan(minDURcf, nchannels, ntrial);
signal.typ = cueTYP;

for idx_trial = 1:ntrial
    fixStart = fixPOS(idx_trial);
    fixEnd = fixStart + minDURfix - 1;
    signal.fix(:,:,idx_trial) = s(fixStart:fixEnd,:);

    cueStart = cuePOS(idx_trial);
    cueEnd = cueStart + minDURcue - 1;
    signal.cue(:,:,idx_trial) = s(cueStart:cueEnd,:);

    cfStart = cfPOS(idx_trial);
    cfEnd = cfStart + minDURcf - 1;
    signal.cf(:,:,idx_trial) = s(cfStart:cfEnd,:);

    trialStart = trialPOS(idx_trial);
    trialEnd = trialStart + minDURtrial - 1;
    signal.trial(:,:,idx_trial) = s(trialStart:trialEnd,:);
end


end