function idx = resolve_channels(names_wanted, all_labels)
% RESOLVE_CHANNELS  Case-insensitive map from channel names to indices in
%   eeg.info.labels. Errors if any name is not found.
    all_labels = to_strcell(all_labels);
    names_wanted = to_strcell(names_wanted);
    idx = zeros(numel(names_wanted), 1);
    lower_labels = lower(all_labels);
    for k = 1:numel(names_wanted)
        hit = find(strcmp(lower_labels, lower(names_wanted{k})), 1);
        if isempty(hit)
            error('resolve_channels:missing', ...
                  'Channel "%s" not found in GDF labels.', names_wanted{k});
        end
        idx(k) = hit;
    end
end
