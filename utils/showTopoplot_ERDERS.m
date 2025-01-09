%% plot in time the topoplot and the topoplot for only the cf. Therefore keep into consideration the period asked
function showTopoplot_ERDERS(chanlog_path, ERD, sampleRate, divisionSampleRate, channels_label_selected, band, minDurTrial, cueTYP, classes, nchannels, all_title)
disp('   [INFO] plotting')
c_cfPeriod = [3*sampleRate minDurTrial]; % for the last plot in order to have only cf
disp(c_cfPeriod);
start_cue = 2;
start_cf = 3;

%% Visualization ERD
load(chanlog_path);
chanlocs_label = {chanlocs.labels};

period1 = 1:ceil(sampleRate/divisionSampleRate):minDurTrial;
period2 = period1(2):ceil(sampleRate/divisionSampleRate):minDurTrial;
period2 = cat(2, period2, minDurTrial);
period = cat(1, period1, period2);

% show for each period the topoplot
figure();
for idx_period=1:size(period, 2)
    data_1 = mean(mean(ERD(period(1, idx_period):period(2, idx_period), :, cueTYP == classes(1)), 3), 1);
    data_2 = mean(mean(ERD(period(1, idx_period):period(2, idx_period), :, cueTYP == classes(2)), 3), 1);
    data = data_2 - data_1;
    chanlocs_data = zeros(size(chanlocs_label,2), 1);
    for i=1:length(chanlocs_label)
        for j = 1:nchannels
            if strcmpi(chanlocs_label{i}, channels_label_selected{j})
                if ~isnan(data(j))
                    chanlocs_data(i) = data(j);
                else
                    chanlocs_data(i) = 0;
                end

            end
        end
    end
    subplot(2, ceil((size(period, 2) + 2)/2), idx_period);
    topoplot(squeeze(chanlocs_data), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(chanlocs_data)) max(abs(chanlocs_data))]);
    axis image;
    colorbar;
    if period(2, idx_period) - 1 <= sampleRate*start_cue
        title(['fixation: ' num2str(period(1, idx_period)/sampleRate - 1/sampleRate) '-' num2str(period(2, idx_period)/sampleRate - 1/sampleRate) 's']);
    elseif period(2, idx_period) - 1 <= sampleRate*start_cf
        title(['cue: ' num2str(period(1, idx_period)/sampleRate - 1/sampleRate) '-' num2str(period(2, idx_period)/sampleRate - 1/sampleRate) 's']);
    else
        title(['cf: ' num2str(period(1, idx_period)/sampleRate - 1/sampleRate) '-' num2str(period(2, idx_period)/sampleRate - 1/sampleRate) 's']);
    end
end

% show the topoplot for all the cf (class 2 - class 1)
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

subplot(2, ceil((size(period, 2) + 2)/2), idx_period + 1);
topoplot(squeeze(c_cf), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(c_cf)) max(abs(c_cf))]);
axis image;
title(['all the cf: ' num2str(c_cfPeriod(1)/sampleRate - 1/sampleRate) '-' num2str(c_cfPeriod(2)/sampleRate - 1/sampleRate) 's']);
colorbar;

% show the topoplot for all cue and cf (class 2 - class 1)
dataCf_1 = mean(mean(ERD(start_cue*sampleRate:c_cfPeriod(2), :, cueTYP == classes(1)), 3), 1);
dataCf_2 = mean(mean(ERD(start_cue*sampleRate:c_cfPeriod(2), :, cueTYP == classes(2)), 3), 1);
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

subplot(2, ceil((size(period, 2) + 2)/2), idx_period + 2);
topoplot(squeeze(c_cf), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(c_cf)) max(abs(c_cf))], 'electrodes', 'labelpoint');
axis image;
title(['all cue and cf: ' num2str(start_cue) '-' num2str(c_cfPeriod(2)/sampleRate - 1/sampleRate) 's']);
colorbar;

all_title = [all_title '| br-bl | band: [' num2str(band(1)) ',' num2str(band(2)) ']'];
sgtitle(all_title)
end