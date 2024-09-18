% works only for psd
function [fisher, cva] = fisherAndCVA_CVSA(signal_psd, header, cueTYPs, idx_selFreqs, idx_interest_ch, interval_step, frameRate)

    %% extract information
    nchannels = length(idx_interest_ch);
    [fixPOS, fixDUR, cuePOS, cueDUR, cueTYP, cfPOS, cfDUR, n_trial] = extract_info_data(header, 786, cueTYPs, 781);

    %% Extract trial infoo -> a noi interessa solo il continuous feedback con la divisione in secondi
    disp('[INFO] extracting trials and fixation')

    minDur = min(cfDUR);
    intervals = double(1):double(interval_step):double(minDur/frameRate); % in sec
    ck = nan(size(signal_psd,1),length(intervals)+2);

    for idx_inter = 1:length(intervals)+2
        for idx_tr=1:n_trial
            if idx_inter == length(intervals) + 2
                % cue + cf
                c_start = cuePOS(idx_tr);
                c_stop = cfPOS(idx_tr) + cfDUR(idx_tr) - 1;

                ck(c_start:c_stop,idx_inter) = cueTYP(idx_tr);
            elseif idx_inter == length(intervals) + 1
                % cf
                c_start = cfPOS(idx_tr);
                c_stop = cfPOS(idx_tr) + cfDUR(idx_tr) - 1;

                ck(c_start:c_stop,idx_inter) = cueTYP(idx_tr);
            elseif idx_inter == length(intervals)
                % from 3 to end
                c_start = cfPOS(idx_tr) + (intervals(idx_inter)-1)*frameRate;
%                 c_start = cfPOS(idx_tr);
                c_stop = cfPOS(idx_tr) + cfDUR(idx_tr) - 1;

                ck(c_start:c_stop,idx_inter) = cueTYP(idx_tr);
            else
                % intervals
                c_start = cfPOS(idx_tr) + (intervals(idx_inter)-1)*frameRate;
%                 c_start = cfPOS(idx_tr);
                c_stop = cfPOS(idx_tr) + intervals(idx_inter)*frameRate;
                ck(c_start:c_stop,idx_inter) = cueTYP(idx_tr);
            end

        end
    end

    %% keep only required signal
    signal_fisher = signal_psd(:,idx_selFreqs,idx_interest_ch);
    signal_cva = signal_fisher(:,:);

    fisher = nan(length(intervals)+2, length(idx_selFreqs), nchannels);
    cva = nan(length(intervals)+2, length(idx_selFreqs), nchannels);

    %% compute fisher score
    for idx_inter=1:length(intervals) + 2
        cmu = zeros(length(idx_selFreqs), nchannels, length(cueTYPs));
        csigma = zeros(length(idx_selFreqs), nchannels, length(cueTYPs));

        for idx_class=1:length(cueTYPs)
            s = signal_fisher(ck(:,idx_inter)==cueTYPs(idx_class),:,:);
            cmu(:, :, idx_class) = squeeze(mean(s, 1));
            csigma(:,:, idx_class) = squeeze(std(s,1));
        end

        fisher(idx_inter,:,:) = abs(cmu(:,:,2) - cmu(:,:,1)) ./ sqrt((csigma(:,:,1).^2 + csigma(:,:,2).^2));
    end

    %% compute cva
    for idx_inter=1:length(intervals) + 2
        c_ck = ck(:,idx_inter);
        idx = ~isnan(c_ck);
        c = cva_tun_opt(signal_cva(idx,:), c_ck(idx));
        cva(idx_inter,:,:) = reshape(c, length(idx_selFreqs), nchannels);
    end
end