function signal = processing_offline(signal, nchannels, sampleRate, band, filtOrder, avg)
disp('      [INFO] start processing');

% filter alpha band
disp('         [proc] applying filtering')
[b, a] = butter(filtOrder, band(2)*(2/sampleRate),'low');
s_low = filter(b,a,signal);
[b, a] = butter(filtOrder, band(1)*(2/sampleRate),'high');
s_filt = filter(b,a,s_low);

% squaring
disp('         [proc] applying power')
s_rect = power(s_filt, 2);

% average windows 1 sec
disp('         [proc] applying average window')
s_moveavg = zeros(size(signal));
for idx_ch=1:nchannels
    s_moveavg(:, idx_ch) = (filter(ones(1,avg*sampleRate)/avg/sampleRate, 1, s_rect(:, idx_ch)));
end

% log
disp('         [proc] applying log')
signal = log(s_moveavg);

end