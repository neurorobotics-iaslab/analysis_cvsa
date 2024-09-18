function showTopoplot_ERDERS(chanlog_path, ERD, sampleRate, channels_label_selected, band, minDur, cueTYP, classes, nchannels, all_title)
disp('[INFO] plotting')
%% Visualization ERD
load(chanlog_path);
chanlocs_label = {chanlocs.labels};
fixPeriod = [1/sampleRate 2]*sampleRate;
cuePeriod = [2 3]*sampleRate;
a = 3:minDur/sampleRate;
b = 4:minDur/sampleRate;
b = cat(2, b, minDur/sampleRate);
cfPeriod = cat(1, a, b)' *sampleRate;
figure();

% fix
dataFix_1 = mean(mean(ERD(fixPeriod(1):fixPeriod(2), :, cueTYP == classes(1)), 3), 1);
dataFix_2 = mean(mean(ERD(fixPeriod(1):fixPeriod(2), :, cueTYP == classes(2)), 3), 1);
dataFix = dataFix_2 - dataFix_1;
c_fix = zeros(64, 1);
for i=1:length(chanlocs_label)
    for j = 1:nchannels
        if strcmpi(chanlocs_label{i}, channels_label_selected{j})
            if ~isnan(dataFix(j))
                c_fix(i) = dataFix(j);
            else
                c_fix(i) = 0;
            end

        end
    end
end

subplot(2, ceil((size(cfPeriod,1) + 2)/2), 1);
topoplot(squeeze(c_fix), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(c_fix)) max(abs(c_fix))]);
axis image;
title(['ERD/ERS fix (band [' num2str(band(1)) '-' num2str(band(2)) ']) - bottom right - bottom left']);
colorbar;

% cue
dataCue_1 = mean(mean(ERD(cuePeriod(1):cuePeriod(2), :, cueTYP == classes(1)), 3), 1);
dataCue_2 = mean(mean(ERD(cuePeriod(1):cuePeriod(2), :, cueTYP == classes(2)), 3), 1);
dataCue = dataCue_2 - dataCue_1;
c_cue = zeros(64, 1);
for i=1:length(chanlocs_label)
    for j = 1:nchannels
        if strcmpi(chanlocs_label{i}, channels_label_selected{j})
            if ~isnan(dataCue(j))
                c_cue(i) = dataCue(j);
            else
                c_cue(i) = 0;
            end

        end
    end
end

subplot(2, ceil((size(cfPeriod,1) + 2)/2), 2);
topoplot(squeeze(c_cue), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(c_cue)) max(abs(c_cue))]);
axis image;
title(['ERD/ERS cue (band [' num2str(band(1)) '-' num2str(band(2)) ']) - bottom right - bottom left']);
colorbar;

% cf in intervals
for idx_cf = 1:size(cfPeriod,1)
    c_cfPeriod = cfPeriod(idx_cf,:);
    dataCf_1 = mean(mean(ERD(c_cfPeriod(1):c_cfPeriod(2), :, cueTYP == classes(1)), 3), 1);
    dataCf_2 = mean(mean(ERD(c_cfPeriod(1):c_cfPeriod(2), :, cueTYP == classes(2)), 3), 1);
    dataCf = dataCf_2 - dataCf_1;
    c_cf = zeros(64,1);
    for i=1:length(chanlocs_label)
        for j = 1:nchannels
            if strcmpi(chanlocs_label{i}, channels_label_selected{j})
                if ~isnan(dataCf(j))
                    c_cf(i) = dataCf(j);
                else
                    c_cf(i) = 0;
                end
            end
        end
    end


    subplot(2, ceil((size(cfPeriod,1) + 2)/2), idx_cf + 2);
    topoplot(squeeze(c_cf), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(c_cf)) max(abs(c_cf))]);
    axis image;
    title(['ERD/ERS (band [' num2str(band(1)) '-' num2str(band(2)) ']) -- br - bl -- cf from ' num2str(ceil((cfPeriod(idx_cf,1) - cfPeriod(1,1))/sampleRate))...
        's to ' num2str(ceil((cfPeriod(idx_cf,2) - cfPeriod(1,1))/sampleRate)) 's']);
    colorbar;
end
sgtitle(all_title)

% cf all
c_cfPeriod = [3*sampleRate minDur];
dataCf_1 = mean(mean(ERD(c_cfPeriod(1):c_cfPeriod(2), :, cueTYP == classes(1)), 3), 1);
dataCf_2 = mean(mean(ERD(c_cfPeriod(1):c_cfPeriod(2), :, cueTYP == classes(2)), 3), 1);
dataCf = dataCf_2 - dataCf_1;
c_cf = zeros(64,1);
for i=1:length(chanlocs_label)
    for j = 1:nchannels
        if strcmpi(chanlocs_label{i}, channels_label_selected{j})
            if ~isnan(dataCf(j))
                c_cf(i) = dataCf(j);
            else
                c_cf(i) = 0;
            end
        end
    end
end


figure();
topoplot(squeeze(c_cf), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(c_cf)) max(abs(c_cf))], 'electrodes', 'labelpoint');
axis image;
title(['ERD/ERS (band [' num2str(band(1)) '-' num2str(band(2)) ']) -- br - bl -- cf from 0' ...
    's to ' num2str(ceil((c_cfPeriod(2) - c_cfPeriod(1))/sampleRate)) 's']);
colorbar;

% sgtitle(all_title)
end