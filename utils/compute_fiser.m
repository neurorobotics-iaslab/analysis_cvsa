%% function that compute fisher score:
%       INPUT:
%               - classes: list of classes
%               - signal: signal to compute fisher score (samples x channels x trial)
%               - typ: vector with the cue type (1 x trial)
%       OUTPUT:
%               - fisher: fisher score computed (1 x channels)

function fisher = compute_fiser(classes, signal, typ, normalize_std)
nclasses = length(classes);
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

% comput std and mean
mu = nan(nclasses, nchannels);
sigma = nan(nclasses, nchannels);
for idx_class = 1:nclasses
    c_class = classes(idx_class);
    mu(idx_class,:) = squeeze(mean(reshaped_signal(ck == c_class,:), 1)); % mean first on samples then on trials
    sigma(idx_class,:) = std(reshaped_signal(ck == c_class, :));
end

% compute fisher
if normalize_std
    fisher = abs(mu(1,:) - mu(2,:)) ./ sqrt(sigma(1,:).^2 + sigma(2,:).^2);
else
    fisher = abs(mu(1,:) - mu(2,:));
end

end