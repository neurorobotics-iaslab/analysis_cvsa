function artifact = artifact_rejection(signal, header, nchannels, bufferSize, chunkSize, eog, picks)
% ARTIFACT_REJECTION Detects EOG and Muscle (EMG) artifacts in an EEG signal.
%
%   This function simulates an online "chunk-wise" processing using a
%   sliding buffer. For each data window (buffer), it performs two 
%   parallel checks:
%
%   1.  Picks (picks) Check: Applies a high-pass filter to the non-EOG 
%       channels and checks if the amplitude exceeds a threshold 
%       ('muscle.threshold').
%
%   2.  Ocular (EOG) Check: Applies a band-pass filter to the EOG channels,
%       calculates the HEOG (horizontal) and VEOG (vertical) signals, and
%       checks if they exceed their respective thresholds ('eog.h_threshold' 
%       and 'eog.v_threshold').
%
%   The function returns a binary vector indicating, for each processed 
%   chunk, whether an artifact was detected ('1') or not ('0').
%
%   Inputs:
%     signal      - Matrix (samples x channels) of the raw EEG signal.
%     header      - Header structure (from GDF/file loading) containing 
%                   at least .SampleRate and .Label (channel labels).
%     nchannels   - Number of EEG channels to use (excludes extra channels).
%     bufferSize  - Size of the buffer (in samples) on which to calculate 
%                   artifacts (e.g., 512 for 1 sec at 512Hz).
%     chunkSize   - Size of the chunk (in samples) simulating 
%                   data arrival (e.g., 32).
%     eog         - Configuration structure for EOG artifacts.
%       .label      - Cell array with EOG channel labels.
%                     IMPORTANT: Order is critical. 
%                     For HEOG: {'LeftChannel', 'RightChannel'} (e.g., {'Fp1', 'Fp2'}).
%                     For VEOG: {'LeftChannel', 'RightChannel', 'BottomChannel'} 
%                               (e.g., {'Fp1', 'Fp2', 'EOGv'}).
%       .band       - Vector (1x2) with the EOG frequency band [low high] 
%                     (e.g., [1 10]).
%       .filterOrder- Order of the Butterworth filter for EOG (e.g., 3).
%       .h_threshold- Threshold (in µV) for the HEOG artifact.
%       .v_threshold- Threshold (in µV) for the VEOG artifact.
%     muscle      - Configuration structure for muscle artifacts.
%       .freq       - Cutoff frequency (in Hz) of the high-pass filter 
%                     (e.g., 30). Because ANTNEURO CAP
%       .filterOrder- Order of the Butterworth filter for EMG (e.g., 5).
%       .threshold  - Threshold (in nV or µV) for the EMG artifact on non-EOG channels.
%
%   Outputs:
%     artifact    - Binary column vector (nchunks x 1). '1' if an artifact 
%                   was detected in that chunk, '0' otherwise.
disp('   [proc] start artifact rejection: EOG and picks' );

% initial variables
nchunks = floor(size(signal, 1)/chunkSize);
buffer_peak = nan(bufferSize, nchannels);
buffer_eog = nan(bufferSize, nchannels);
sampleRate = header.SampleRate;
if ~isempty(eog.label)
    [found, indices] = ismember(eog.label, header.Label);
    eog_idx = indices(found);
    non_eog = setdiff(1:nchannels, eog_idx);
else
    non_eog = 1:nchannels;
end

% correct the signal according to the nchannesl
signal = signal(:,1:nchannels);

% prepare the filters for both artefact removal
[b_high, a_high] = butter(picks.filterOrder, picks.freq*(2/sampleRate),'high');
zi_high = [];

if ~isempty(eog.label)
    eog_band = eog.band;
    [b_low_eog, a_low_eog] = butter(eog.filterOrder, eog_band(2)*(2/sampleRate),'low');
    z_low_eog = [];
    [b_high_eog, a_high_eog] = butter(eog.filterOrder, eog_band(1)*(2/sampleRate),'high');
    z_high_eog = [];
end

artifact = zeros(nchunks, 1);

for i=1:nchunks
    % add
    frame = signal((i-1)*chunkSize+1:i*chunkSize,:);
    frame_no_eog = frame(:,non_eog); 
    frame = frame - mean(frame_no_eog, 2);

    % --- muscle artefact part buffer ---
    [frame_peak,zi_high] = filter(b_high,a_high,frame,zi_high);

    buffer_peak(1:end-chunkSize,:) = buffer_peak(chunkSize+1:end,:);
    buffer_peak(end-chunkSize+1:end, :) = frame_peak;

    % --- eog artefact part buffer ---
    if ~isempty(eog.label)
        % compte eog with horizontal and vertical movement
        [data_eog,z_low_eog] = filter(b_low_eog,a_low_eog,frame,z_low_eog);
        [frame_eog,z_high_eog] = filter(b_high_eog,a_high_eog,data_eog,z_high_eog);

        buffer_eog(1:end-chunkSize,:) = buffer_eog(chunkSize+1:end,:);
        buffer_eog(end-chunkSize+1:end, :) = frame_eog;
    end
    
    % buffer pick
    data_non_eog = buffer_peak(:,non_eog);
    if ~any(isnan(buffer_peak))
        if any(abs(data_non_eog(:)) > picks.threshold)
            artifact(i) = 1;
        end
    end

    % vertical and horizontal movements
    if ~isempty(eog.label) && size(eog_idx, 2) > 0
        heog = buffer_eog(:, eog_idx(1)) - buffer_eog(:, eog_idx(2));
        if size(eog_idx, 2) == 2
            % we have only FP1 and FP2
            veog = (buffer_eog(:, eog_idx(1)) + buffer_eog(:, eog_idx(2))) / 2;
        else
            veog = ((buffer_eog(:, eog_idx(1)) + buffer_eog(:, eog_idx(2))) / 2) - buffer_eog(:, eog_idx(3));
        end
        if any(abs(heog(:)) > eog.h_threshold)
            artifact(i) = 1;
        end
        if any(abs(veog(:)) > eog.v_threshold)
            artifact(i) = 1;
        end
    end
end

end