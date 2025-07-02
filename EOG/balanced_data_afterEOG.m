%% balance the number of typ in the data where the eog is not present
% INPUT:
%       - trial_with_eog: vector with 1 and 0 if present or not the eog
%       - trial_typ: trial type
% OUTPUT:
%       - balanced_trial_idx: index of the trial in which eog is not
%       present and it is balanced in number according to the class asked
function balanced_trial_idx = balanced_data_afterEOG(trial_with_eog, trial_typ, classes)
trial_start = trial_typ(~trial_with_eog);
nclasses = length(classes);
count_classes = nan(nclasses, 1);
for idx_class = 1:nclasses
    count_classes(idx_class) = sum(trial_start == classes(idx_class));
end

if any(count_classes ~= count_classes(1))
    min_ntrial = min(count_classes);
    balanced_trial_idx = zeros(size(trial_with_eog));
    count = zeros(nclasses, 1);
    for idx_trial = 1:length(trial_with_eog)
        for idx_class = 1:nclasses
            if trial_typ(idx_trial) == classes(idx_class) && count(idx_class) < min_ntrial && ~trial_with_eog(idx_trial)
                count(idx_class) = count(idx_class) + 1;
                balanced_trial_idx(idx_trial) = 1;
            end
        end
    end
else
    balanced_trial_idx = ~trial_with_eog;
end

end