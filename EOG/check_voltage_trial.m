function trial_with_high_voltage = check_voltage_trial(signal, header, threshold)
    ntrial = sum(header.EVENT.TYP == 1);
    nchannels = size(signal, 2);
    trial_with_high_voltage = zeros(ntrial, 1);
    pos = header.EVENT.POS(header.EVENT.TYP == 781);
    dur = header.EVENT.DUR(header.EVENT.TYP == 781);
    for idx_trial = 1:ntrial
        c_start = pos(idx_trial);
        c_end = c_start + dur(idx_trial) - 1;
        c_signal = signal(c_start:c_end,:);

        % remove the fp and eog channel
        c_signal(:,[1,2,19]) = 0;
        if any(any(c_signal > threshold))
            trial_with_high_voltage(idx_trial) = 1;
            figure();
            imagesc(c_signal')
            yticks(1:nchannels)
            yticklabels(header.Label)
            colorbar;
            title(['trial: ' num2str(idx_trial)])
        end
    end
end