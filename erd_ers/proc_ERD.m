%% function to compute the ERD/ERS
function s_out = proc_ERD(signal, sampleRate, band, filtOrder)
disp(['   [INFO] band: ' num2str(band(1)) '-' num2str(band(2))])

% filter alpha band
disp('      [proc] applying filtering')
[b, a] = butter(filtOrder, band(2)*(2/sampleRate),'low');
s_low = filter(b,a,signal);
[b, a] = butter(filtOrder, band(1)*(2/sampleRate),'high');
s_filt = filter(b,a,s_low);

% squaring
disp('      [proc] applying power')
s_rect = power(s_filt, 2);

% average windows 1 sec
disp('      [proc] applying average window')
s_out = zeros(size(signal));
nchannels = size(signal, 2);
avg = 1;
for idx_ch=1:nchannels
    s_out(:, idx_ch) = (filter(ones(1,avg*sampleRate)/avg/sampleRate, 1, s_rect(:, idx_ch)));
end


end