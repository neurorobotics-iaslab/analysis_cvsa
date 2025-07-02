%% This function compute the EOG to check if in hte cf there is a eye artifact
% INPUT:
%       - signal: matrix with the eeg signal (samples x channels)
%       - header: struct whihc contains the position of the cf and the
%       channels labels
%       - threshold: value used as threshold for eog detection
%       - eog_channel: channels of the eog --> in this way 'FP1', 'FP2', 'EOG'
% OUTPUT:
%       - trial_with_eog: vector of 1 and 0. 1=trial has eog, 0 otherwise
function trial_with_eog = eog_detection(signal, header, threshold, eog_channel)

if size(eog_channel, 2) < 2 || size(eog_channel, 2) > 3
    disp('ERROR: this code works only having two/three electrodes (example FP1 and FP2 or FP1, FP2 and EOG)');
    return;
end

% filtering the signal
band = [1 7];
filtOrder = 3;
sampleRate = header.SampleRate;
[b, a] = butter(filtOrder, band(2)*(2/sampleRate),'low');
s_low = filter(b,a,signal);
[b, a] = butter(filtOrder, band(1)*(2/sampleRate),'high');
s_filt = filter(b,a,s_low);

% take the signal of the electrodes and compute the horizontal and vertical motion
channels_label = header.Label;
idx_eog = find(ismember(channels_label, eog_channel));
heog = s_filt(:, idx_eog(1)) - s_filt(:, idx_eog(2));
if size(eog_channel, 2) == 2
    % we have only FP1 and FP2
    veog = (s_filt(:, idx_eog(1)) + s_filt(:, idx_eog(2))) / 2;
else
    veog = ((s_filt(:, idx_eog(1)) + s_filt(:, idx_eog(2))) / 2) - s_filt(:, idx_eog(3));
end


% check for each trial during the cf
ntrial = sum(header.EVENT.TYP == 1);
trial_with_eog = zeros(ntrial, 1);
pos = header.EVENT.POS(header.EVENT.TYP == 781); 
dur = header.EVENT.DUR(header.EVENT.TYP == 781); 
for idx_trial = 1:ntrial
    c_start = pos(idx_trial);
    c_end = c_start + dur(idx_trial) - 1;
    c_veog = abs(veog(c_start:c_end,:));
    c_heog = abs(heog(c_start:c_end,:));

    if any(c_heog > threshold)
        trial_with_eog(idx_trial) = 1;
    elseif any(c_veog > threshold)
        trial_with_eog(idx_trial) = 1;
    end
end

end