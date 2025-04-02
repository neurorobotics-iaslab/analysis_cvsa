
ntrial = 75;
Fs = 512;
num_channels = size(channelsSelected, 2);

c_signal = squeeze(data{1}.cue(:,channelsSelected,1));
c_signal = [c_signal; squeeze(data{1}.cf(:,channelsSelected,1))];

window_size = 100;  % 1 sec window (adjust based on sampling rate)
step_size = 50; 

num_windows = floor((size(c_signal, 1) - window_size) / step_size) + 1;
entropy_matrix = zeros(num_windows, ntrial); % Matrix to store average entropy per trial

for t = 1:ntrial
    c_signal = squeeze(data{1}.cue(:,channelsSelected,t));
    c_signal = [c_signal; squeeze(data{1}.cf(:,channelsSelected,t))];

    trial_entropy = eeg_entropy_windowed(c_signal, 'shannon', window_size, step_size);
    entropy_matrix(:, t) = mean(trial_entropy, 2); % Average entropy over channels
end

entropy_diff = diff(entropy_matrix, 1, 1); % Compute entropy change per trial (time-wise diff)
shift_threshold = 0.5; %mean(entropy_diff, 2) + 2 * std(entropy_diff, 0, 2); % Compute threshold per time window
shift_indices = abs(entropy_diff) > shift_threshold; % Binary matrix of detected shifts

mean_shift = mean(shift_indices, 2); % Average over trials

time_vector = (1:size(mean_shift,1)) * step_size / Fs; % Convert to seconds

figure;
plot(time_vector, mean_shift, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Shift Detection Probability');
title('Attention Shift Detection Across Trials');




function entropy_matrix = eeg_entropy_windowed(eeg_signal, method, window_size, step_size)
    % EEG_ENTROPY_WINDOWED Computes entropy in sliding windows
    %
    % entropy_matrix = eeg_entropy_windowed(eeg_signal, method, window_size, step_size)
    %
    % eeg_signal: matrix (time x channels) of EEG data
    % method: 'shannon' or 'approx'
    % window_size: number of samples per window
    % step_size: step between consecutive windows
    %
    % Returns: entropy_matrix (windows x channels)

    [num_samples, num_channels] = size(eeg_signal);
    num_windows = floor((num_samples - window_size) / step_size) + 1;
    entropy_matrix = zeros(num_windows, num_channels);

    for win = 1:num_windows
        start_idx = (win - 1) * step_size + 1;
        end_idx = start_idx + window_size - 1;
        window_data = eeg_signal(start_idx:end_idx, :);

        entropy_matrix(win, :) = eeg_entropy(window_data, method);
    end
end
function entropy_values = eeg_entropy(eeg_signal, method)
    % EEG_ENTROPY Computes entropy for each channel of an EEG signal
    %
    % entropy_values = eeg_entropy(eeg_signal, method)
    %
    % eeg_signal: matrix (time x channels) of EEG data
    % method: 'shannon' or 'approx'
    %
    % Returns a vector of entropy values for each channel
    
    if nargin < 2
        method = 'shannon'; % Default method
    end
    
    [num_samples, num_channels] = size(eeg_signal);
    entropy_values = zeros(1, num_channels);
    
    for ch = 1:num_channels
        signal = eeg_signal(:, ch); % Extract individual channel
        
        switch lower(method)
            case 'shannon'
                % Normalize the signal between 0 and 1
                signal = signal - min(signal);
                signal = signal / max(signal);
                
                % Create probability histograms
                nbins = 50;
                p = histcounts(signal, nbins, 'Normalization', 'probability');
                
                % Remove zeros to avoid log(0)
                p(p == 0) = [];
                
                % Compute Shannon entropy
                entropy_values(ch) = -sum(p .* log2(p));
            
            case 'approx'
                % Typical parameters for approximate entropy
                m = 2;
                r = 0.2 * std(signal);
                
                entropy_values(ch) = approx_entropy(signal, m, r);
            
            otherwise
                error('Unrecognized method. Use "shannon" or "approx".');
        end
    end
end

function ApEn = approx_entropy(signal, m, r)
    % APPROX_ENTROPY Computes approximate entropy of a signal
    N = length(signal);
    r = r * std(signal); % Tolerance factor
    
    C = zeros(1, N - m + 1);
    for i = 1:N - m + 1
        count = 0;
        for j = 1:N - m + 1
            if max(abs(signal(i:i+m-1) - signal(j:j+m-1))) <= r
                count = count + 1;
            end
        end
        C(i) = count / (N - m + 1);
    end
    
    phi_m = sum(log(C)) / (N - m + 1);
    
    % Repeat for m+1
    C = zeros(1, N - m);
    for i = 1:N - m
        count = 0;
        for j = 1:N - m
            if max(abs(signal(i:i+m) - signal(j:j+m))) <= r
                count = count + 1;
            end
        end
        C(i) = count / (N - m);
    end
    
    phi_m1 = sum(log(C)) / (N - m);
    
    % Compute approximate entropy
    ApEn = phi_m - phi_m1;
end