%% this function extract the features giving the signal and the features selected
%   INPUT:
%       - signal: it is a struct of nband element
%           signal{band}: contains only cf values (samples*trial x channels)
%       - bands: vector of bands (1 x nbands)
%       - bands_features: bands asked for the features (nfeatures x 2)
%       - idxchans_features: channels for the featus (1 x nfeature)
%   OUTPUT:
%       - features: matrix with all the features (samples*trial x nfeatures)
function features = features_extraction(signal, bands, bands_features, idxchans_features)
    nfeatures = length(idxchans_features);
    nbands = size(bands, 2);

    % extract only the interested features
    features = nan(size(signal{1}.data, 1), nfeatures);
    for idx_feature=1:nfeatures
        for idx_band=1:nbands
            band = bands{idx_band};
            if isequal(band, bands_features(idx_feature,:))
                features(:,idx_feature) = signal{idx_band}.data(:,idxchans_features(idx_feature));
            end
        end
    end
end