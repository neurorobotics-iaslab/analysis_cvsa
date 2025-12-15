%% function that processed the signal data like what is performed in ros
%   INPUT:
%       - signal: matrix of signals (samples x channels)
%       - header: header of the gdf
%       - nchannels: numer of channels
%       - bufferSize: size of the buffer (rosneuro uses 512)
%       - band: band to perform low and high filter
%       - filterOrder: int with th efilter order
%       - chunkSize: size of chunk (rosneuro uses 32)
%   OUTPUT:
%       - signal_processed: signal processed
%       - header: modification in the POS and DUR of the gdf header
function [signal_processed, header] = processing_onlineROS_CSD_hilbert(signal, header, nchannels, bufferSize, filterOrder, band, chunkSize)
disp(['   [proc] start processing like ros for band ' num2str(band(1)) '-' num2str(band(2))]);

persistent M_CSD valid_idx_eeg
if isempty(M_CSD)
    try
        load('csd_transform_matrix.mat', 'M_CSD', 'valid_idx'); 
        % Rinominato per chiarezza
        valid_idx_eeg = valid_idx; 
        disp('      Matrice CSD caricata correttamente.');
    catch
        error('      File csd_transform_matrix.mat non trovato! Esegui prima setup_csd.m');
    end
end

nchunks = floor(size(signal, 1)/chunkSize);
buffer = nan(bufferSize, nchannels);

[b_low, a_low] = butter(filterOrder, band(2)*(2/header.SampleRate),'low');
[b_high, a_high] = butter(filterOrder, band(1)*(2/header.SampleRate),'high');
zi_low = [];
zi_high = [];

signal_processed = nan(nchunks, nchannels);

for i=1:nchunks

    frame = signal((i-1)*chunkSize+1:i*chunkSize,:);
    frame_eeg_raw = frame(:, valid_idx_eeg);
    
    % Application CSD
    frame_eeg_csd = frame_eeg_raw * M_CSD;
    frame(:, valid_idx_eeg) = frame_eeg_csd;
    
    % apply low and high pass filters
    [tmp_data, zi_low] = filter(b_low,a_low,frame,zi_low);
    [tmp_data,zi_high] = filter(b_high,a_high,tmp_data,zi_high);

    buffer(1:end-chunkSize,:) = buffer(chunkSize+1:end,:);
    buffer(end-chunkSize+1:end, :) = tmp_data;

    % check
    if any(isnan(buffer))
        continue;
    end

    % apply power with hilbert
    analytic = hilbert(buffer);
    tmp_data = abs(analytic).^2;

    % apply average
    tmp_data = mean(tmp_data, 1);

    signal_processed(i,:) = tmp_data;
end

header.EVENT.DUR = round(header.EVENT.DUR/chunkSize);
header.EVENT.POS = round(header.EVENT.POS/chunkSize); 
end