function dfet = compute_vectorFeature(channelsSelected, freqsSelected, all_freqs, psd_matrix)
    % inputs:
    %   - channelsSelected      all channel needed for the decoder
    %   - freqsSelected         all freqs needed for the decoder
    %   - all_freqs             all the freq used into the calculatio of psd
    %   - psd_matrix            psd values
    % output:
    %   - df                    vector with all the features

    dfet = [];
    for ch=channelsSelected
        freqs_idxs = [];
        freqs = cell2mat(freqsSelected(ch));
        for c_freq= freqs
            temp_idx = [find(all_freqs==c_freq)];
            freqs_idxs = cat(2, freqs_idxs, temp_idx);
        end
        c_dfet = psd_matrix(1,freqs_idxs,ch);
        dfet = cat(2, dfet, c_dfet);
     end