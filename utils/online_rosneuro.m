function [s_buffer, s_band, s_pow, s_avg, s_log] = online_rosneuro(signal, bufferSize, frameSize, band, filterOrder, sampleRate)
    nchannels = size(signal, 2);
    nsamples = size(signal,1);
    nchunk = floor(nsamples/frameSize);
    buffer = nan(bufferSize, nchannels);
    s_buffer = [];
    s_band = [];
    s_pow = [];
    s_avg = [];
    s_log = [];

    [b_low, a_low] = butter(filterOrder, band(2)*(2/sampleRate),'low');
    [b_high, a_high] = butter(filterOrder, band(1)*(2/sampleRate),'high');

    zi_low = [];
    zi_high = [];
    for i=1:nchunk
        % add
        frame = signal((i-1)*frameSize+1:i*frameSize,:);
        buffer(1:end-frameSize,:) = buffer(frameSize+1:end,:);
        buffer(end-frameSize+1:end, :) = frame;

        % check
        if any(isnan(buffer))
            continue;
        end

        s_buffer = cat(1, s_buffer, buffer);

        % apply low and high pass filters
        [s_low, zi_low] = filter(b_low,a_low,buffer,zi_low);
        [tmp_data,zi_high] = filter(b_high,a_high,s_low,zi_high);
        s_band = cat(1, s_band, tmp_data);

        % apply pow
        tmp_data = power(tmp_data, 2);
        s_pow = cat(1, s_pow, tmp_data);

        % apply average
        tmp_data = mean(tmp_data, 1);
        s_avg = cat(1, s_avg, tmp_data);

        % apply log
        tmp_data = log(tmp_data);
        s_log = cat(1, s_log, tmp_data);
    end
end