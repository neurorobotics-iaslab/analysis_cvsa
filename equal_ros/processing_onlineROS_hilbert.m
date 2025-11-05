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
function [signal_processed, header] = processing_onlineROS_hilbert(signal, header, nchannels, bufferSize, filterOrder, band, chunkSize)
disp(['      [INFO] start processing like ros for band ' num2str(band(1)) '-' num2str(band(2))]);

nchunks = floor(size(signal, 1)/chunkSize);
buffer = nan(bufferSize, nchannels);

[b_low, a_low] = butter(filterOrder, band(2)*(2/header.SampleRate),'low');
[b_high, a_high] = butter(filterOrder, band(1)*(2/header.SampleRate),'high');
zi_low = [];
zi_high = [];

signal_processed = nan(nchunks, nchannels);

for i=1:nchunks
    % add
    frame = signal((i-1)*chunkSize+1:i*chunkSize,:);
    buffer(1:end-chunkSize,:) = buffer(chunkSize+1:end,:);
    buffer(end-chunkSize+1:end, :) = frame;

    % check
    if any(isnan(buffer))
        continue;
    end

    % apply low and high pass filters
    [tmp_data, zi_low] = filter(b_low,a_low,buffer,zi_low);
    [tmp_data,zi_high] = filter(b_high,a_high,tmp_data,zi_high);

    % apply power with hilbert
    analytic = hilbert(tmp_data);
    tmp_data = abs(analytic).^2;

    % apply average
    tmp_data = mean(tmp_data, 1);

    signal_processed(i,:) = tmp_data;
end

header.EVENT.DUR = round(header.EVENT.DUR/chunkSize);
header.EVENT.POS = round(header.EVENT.POS/chunkSize); % not this part since signal_processed has nan: - (bufferSize/chunkSize) + 1;  
% header.SampleRate = header.SampleRate/chunkSize;