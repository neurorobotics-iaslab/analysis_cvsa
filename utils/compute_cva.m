%% function which compute the cva using cva_opt by luke
%       INPUT:
%           - classes: vector with the classes values
%           - signal: signal to compute fisher score (samples x channels x trial)
%           - typ: vector with the cue type (1 x trial)
%       OUTPUT:
%           - c: computed cva (1 x channels)
function c = compute_cva(signal, typ)
addpath('/home/paolo/Local/cnbi-smrtrain/toolboxes/cva')
nchannels = size(signal, 2);
ntrial = length(typ);
nsamples = size(signal, 1);

% reshape the data to have samples*trial x channels
reshaped_signal = nan(nsamples * ntrial, nchannels);
ck = [];
for idx_trial = 1:ntrial
    reshaped_signal((idx_trial - 1)*nsamples + 1 : idx_trial*nsamples, :) = signal(:,:,idx_trial);
    ck = cat(1, ck, repmat(typ(idx_trial), [nsamples,1]));
end

c = cva_tun_opt(reshaped_signal, ck);

end