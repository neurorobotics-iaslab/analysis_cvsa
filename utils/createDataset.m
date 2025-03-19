%% create the dataset giving the data and the features selected.
%       INPUT:
%           - files: name of the files uesed for the creation of the dataset
%           - data: struct composed by band x (samples x channels x trial)
%               - data{}.cf: contains only cf values (samples x channels x trial)
%               - data{}.fix:
%               - data{}.cue:
%               - data{}.typ: contains the label of the trial (1 x trial)
%           - selectedFeatures:  (features x 2) first element channel, second the band
%           - bands: list of all the bands (1 x bands)
%       OUTPUT:
%           - dataset:
%               - X: contain samples
%               - y: contains labels

% the dataset is created starting form the cue to the end of cf
function [X, y, info] = createDataset(files, data, selectedFeatures, bands, channels_label)

info.files = files;
% count occurencies to have a balanced dataset
[classes, ~, idx_classes] = unique(data{1}.typ);
count_trials4class = accumarray(idx_classes, 1);
trial4class = min(count_trials4class);

% reshape the data
X = [];
y = [];
info.startTrial    = [];
info.bandSelected  = [];
info.chSelected    = [];
info.channelsLabel = channels_label;
nbands = size(bands, 2);
ntrial = size(data{1}.typ);
features = cell(1, nbands);
for idx_band=1:nbands
    count_trials4class = [0, 0];
    for idx_trial=1:ntrial
        idx_class = find(classes == data{idx_band}.typ(idx_trial));
        if ~(count_trials4class(idx_class) >= trial4class)
            if idx_band == 1 % we need it only once
                info.startTrial = cat(1, info.startTrial, size(features{idx_band}, 1));
                y = cat(1, y, repmat(classes(idx_class), size(data{idx_band}.cf, 1), 1));
            end
            features{idx_band} = cat(1, features{idx_band}, data{idx_band}.cf(:,:,idx_trial));
            count_trials4class(idx_class) = count_trials4class(idx_class) + 1;
        end
    end
end

% extract only useful data
nfeatures = size(selectedFeatures, 1);

for idx_feature = 1:nfeatures
    tmp_features = features{selectedFeatures(idx_feature, 2)}; % take the correct band
    X = cat(2, X, tmp_features(:,selectedFeatures(idx_feature, 1))); % take the chns
    info.chSelected = cat(1, info.chSelected, selectedFeatures(idx_feature, 1));
    info.bandSelected = cat(1, info.bandSelected, bands{selectedFeatures(idx_feature, 2)});
end
end




