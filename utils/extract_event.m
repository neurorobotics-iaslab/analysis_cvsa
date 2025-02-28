%% extract only the cf
% INPUT:
%       - s: signal (samples x channels)
%       - header: header with TYP, POS and DUR
%       - classes: classes values to extract
%       - event: events to extract (example 781, 786, ..) as a list
% OUTPUT: 
%        - signal struct contains all the data in:
%               data: signal of each trial (samples x channels) x trial
%               typ: is a vector 1 x trial contain the value of the cue
function [signal, info] = extract_event(s, header, classes, event)
hitMissTimeout = [897, 898, 899];
nclasses = size(classes, 2);

cueTYP = header.TYP(ismember(header.TYP, classes));
info.hit = header.TYP(ismember(header.TYP, hitMissTimeout));
info.TYP = header.TYP(ismember(header.TYP, event));

POS = header.POS(ismember(header.TYP, event));
DUR = header.DUR(ismember(header.TYP, event));

ntrial = length(POS);

signal.data = [];
signal.typ = cueTYP;

info.ths = nan(ntrial, nclasses);
info.ths_rejection = nan(ntrial, nclasses);
info.bufferSize_integrator = nan(ntrial, 1);
info.alpha = nan(ntrial, 1);
info.startEvent = [];
info.endEvent = [];
info.startNewFile =  [];
idx_file = 0;
nfiles = size(header.startNewFile, 1);

info.sampleRate = header.sampleRate;

for idx_trial = 1:ntrial

    cfStart = POS(idx_trial);
    info.startEvent = cat(1, info.startEvent, size(signal.data, 1)+1);
    cfEnd = cfStart + DUR(idx_trial) - 1;
    if ~(idx_file == nfiles) 
        if cfStart >= header.startNewFile(idx_file+1)
            info.startNewFile = cat(1, info.startNewFile, size(signal.data, 1)+1);
            idx_file = idx_file + 1;
        end
    end

    signal.data = cat(1, signal.data, s(cfStart:cfEnd,:));
    
    info.ths(idx_trial,:) = cell2mat(header.ths(cfStart,:));
    info.ths_rejection(idx_trial,:) = cell2mat(header.ths_rejection(cfStart,:));
    info.alpha(idx_trial) = header.alpha(cfStart);
    info.bufferSize_integrator(idx_trial) = header.bufferSize_integrator(cfStart);
    
    info.endEvent = cat(1, info.endEvent, size(signal.data, 1));
end


end