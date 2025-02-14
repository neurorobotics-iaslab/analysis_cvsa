%% extract only the cf
% INPUT:
%       - s: signal (samples x channels)
%       - header: header with TYP, POS and DUR
%       - classes: classes values to extract
%       - cf_event: event of the CF
% OUTPUT: 
%        - signal struct contains all the data in samples x channels x trial:
%               signal: signal of each trial
%        - typ is a vector 1 x trial contain the value of the cue
function [signal, info] = extract_cf(s, header, classes, cf_event)
hitMissTimeout = [897, 898, 899];
nclasses = size(classes, 2);

cueTYP = header.TYP(ismember(header.TYP, classes));
info.hit = header.TYP(ismember(header.TYP, hitMissTimeout));

cfPOS = header.POS(header.TYP == cf_event);
cfDUR = header.DUR(header.TYP == cf_event);

ntrial = length(cfPOS);

signal.data = [];
signal.typ = cueTYP;

info.ths = nan(ntrial, nclasses);
info.ths_rejection = nan(ntrial, nclasses);
info.bufferSize_integrator = nan(ntrial, 1);
info.startCf = [];
info.endCf = [];
info.startNewFile =  [];
idx_file = 0;
nfiles = size(header.startNewFile, 1);

info.sampleRate = header.sampleRate;

for idx_trial = 1:ntrial

    cfStart = cfPOS(idx_trial);
    info.startCf = cat(1, info.startCf, size(signal.data, 1)+1);
    cfEnd = cfStart + cfDUR(idx_trial) - 1;
    if ~(idx_file == nfiles) 
        if cfStart >= header.startNewFile(idx_file+1)
            info.startNewFile = cat(1, info.startNewFile, size(signal.data, 1)+1);
            idx_file = idx_file + 1;
        end
    end

    signal.data = cat(1, signal.data, s(cfStart:cfEnd,:));
    
    info.ths(idx_trial,:) = cell2mat(header.ths(cfStart,:));
    info.ths_rejection(idx_trial,:) = cell2mat(header.ths_rejection(cfStart,:));
    info.bufferSize_integrator(idx_trial) = header.bufferSize_integrator(cfStart);
    info.endCf = cat(1, info.endCf, size(signal.data, 1));
end


end