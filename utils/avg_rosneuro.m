function avg = avg_rosneuro(data, bufferSize, frameSize)
    nchannels = size(data, 2);
    nsamples = size(data,1);
    nchunk = floor(nsamples/frameSize);
    buffer = nan(bufferSize, nchannels);
    idx_circle = 1;
    max_idx_circle = bufferSize/frameSize;
    avg = [];

    for i=1:nchunk
        % add
        buffer((idx_circle-1)*frameSize+1:idx_circle*frameSize,:) = data((i-1)*frameSize+1:i*frameSize,:);

        % check
        if isnan(any(buffer))
            continue;
        end

        % apply average
        avg = cat(1, avg, mean(buffer, 1));

        % update idx_circle
        idx_circle = idx_circle + 1;
        if idx_circle == max_idx_circle
            idx_circle = 1;
        end
    end
end