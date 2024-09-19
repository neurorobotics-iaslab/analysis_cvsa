function result = eye_movement_check(data,channels_label,threshold,sampleRate)
    %Fai il check di Fp1,Fp2,EOG per valori > 25mv
    % Se sono maggiori il trial è da scartare, altrimenti tieni il trial e vai
    % avanti. I dati sono nell'ordine dei microVolt: 2.02*10^4*10^-6

    % INPUT:
    %   data: eeg signal per trial
    %   channels_label: cell vector with all the channels used for
    %   recordings
    %   threshold: value to define if the trial should be discarded

    % OUTPUT:
    %   result: returns true if any value of the specific channels is > threshold, otherwise return false
    result = false;
    target_electrodes = {'FP1','FP2'};
    eog_ch= find(ismember(channels_label, target_electrodes));

    %% horizontal and vertical movement
    h_mov = data(:,eog_ch(1)) - data(:,eog_ch(2));
    v_mov = (data(:,eog_ch(1)) + data(:,eog_ch(2)))/2;
    eye_mov = [h_mov v_mov];

    %% filtraggio
    % Filter Parameters
    filtOrder = 2;
    band = [1 10];
    [b_low, a_low] = butter(filtOrder, band(2)/(sampleRate/2), 'low');
    [b_high, a_high] = butter(filtOrder, band(1)/(sampleRate/2), 'high');
    % Apply filters
    sfilt = zeros(size(eye_mov));
    for i=1:length(target_electrodes)
    sfilt(:,i) = filter(b_low,a_low,eye_mov(:,i));
    sfilt(:,i) = filter(b_high,a_high,sfilt(:,i));
    end

    % Aggiungi filtraggio su movimento orizzontale (differenza FP1-FP2) e
    % verticale (media) di ordine 2 e freq. [1-10]). Di questi prendere gli
    % abs e fare confronto con soglia
        for j=1:length(eog_ch)
            if any(abs(sfilt(:,j))>threshold) %threshold a 30mV
                result = true;
                return;
            end
        end       
end
