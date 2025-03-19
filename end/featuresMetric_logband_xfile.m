%% compute the fischer to understand features using the log band
clc; clearvars;
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')

%% Initialization
a = 6:2:16;
b = a+2;
c = [a; b];
bands = [];
for i=1:length(a)
    bands{i} = [a(i), b(i)];
end
% bands = {[6 8], [8 10], [10 12], [12 14], [14 16]};
bands_str = cellfun(@(x) sprintf('%d-%d', x(1), x(2)), bands, 'UniformOutput', false);
nbands = length(bands);
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
signals = cell(1, nbands);
headers = cell(1, nbands);
for idx_band = 1:nbands
    headers{idx_band}.TYP = [];
    headers{idx_band}.DUR = [];
    headers{idx_band}.POS = [];
    signals{idx_band} = [];
end
classes = [730 731];
cf_event = 781;
fix_event = 786;
nchannels = 39;
nclasses = length(classes);
filterOrder = 4;
avg = 1;

%% Load file
[filenames, pathname] = uigetfile('*.gdf', 'Select GDF Files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end
subject = filenames{1}(1:2);
%% concatenate the files
nFiles = length(filenames);
for idx_file= 1: nFiles
    fullpath_file = fullfile(pathname, filenames{idx_file});
    disp(['file (' num2str(idx_file) '/' num2str(nFiles)  '): ', filenames{idx_file}]);
    [c_signal,header] = sload(fullpath_file);
    c_signal = c_signal(:,1:nchannels);
    channels_label = header.Label;

    for idx_band = 1:nbands
        band = bands{idx_band};

        signal_processed = proc_512hz(c_signal, header.SampleRate, band, filterOrder, avg);

        c_header = headers{idx_band};
        c_header.sampleRate = header.SampleRate;
        c_header.channels_labels = header.Label;
        c_header.TYP = cat(1, c_header.TYP, header.EVENT.TYP);
        c_header.DUR = cat(1, c_header.DUR, header.EVENT.DUR);
        c_header.POS = cat(1, c_header.POS, header.EVENT.POS);

        signals{idx_band} = signal_processed(:,:);
        headers{idx_band} = c_header;
    end



    %% Labelling data
    [~, channelsSelected] = ismember(channels_select, channels_label);
    nchannelsSelected = size(channelsSelected, 2);
    events = headers{1};
    sampleRate = events.sampleRate;
    cuePOS = events.POS(ismember(events.TYP, classes));
    cueDUR = events.DUR(ismember(events.TYP, classes));
    cueTYP = events.TYP(ismember(events.TYP, classes));

    fixPOS = events.POS(events.TYP == 786);
    fixDUR = events.DUR(events.TYP == 786);

    cfPOS = events.POS(events.TYP == 781);
    cfDUR = events.DUR(events.TYP == 781);

    %% Trial extraction
    n_trial = size(fixPOS, 1);
    trial_start = nan(n_trial, 1);
    trial_end = nan(n_trial, 1);
    trial_typ = nan(n_trial, 1);
    fix_start = nan(n_trial, 1);
    fix_end = nan(n_trial, 1);
    for idx_trial = 1:n_trial
        trial_start(idx_trial) = cuePOS(idx_trial);
        trial_typ(idx_trial) = cueTYP(idx_trial);
        trial_end(idx_trial) = cfPOS(idx_trial) + cfDUR(idx_trial) - 1;
        fix_start(idx_trial) = fixPOS(idx_trial);
        fix_end(idx_trial) = fixPOS(idx_trial) + fixDUR(idx_trial) - 1;
    end

    min_trial_data = min(trial_end - trial_start);
    min_fix_data = min(fix_end - fix_start);
    trial_data = nan(min_trial_data, nbands, nchannels, n_trial);
    basline_data = nan(min_fix_data, nbands, nchannels, n_trial);
    for idx_band = 1:nbands
        c_signal = signals{idx_band};
        for trial = 1:n_trial
            c_start = trial_start(trial);
            c_end = trial_start(trial) + min_trial_data - 1;
            trial_data(:,idx_band,:,trial) = c_signal(c_start:c_end,:);
            c_fix_start = fix_start(trial);
            c_fix_end = fix_start(trial) + min_fix_data - 1;
            basline_data(:,idx_band,:,trial) = c_signal(c_fix_start:c_fix_end,:);
        end
    end


    %% Fisher score in all the trial (cue + cf) and for each window
    nrows = 4;
    logband_hz = sampleRate;
    wlength = 0.10; %s
    wlength = ceil(wlength * logband_hz);
    overlap = 0.05; %s
    overlap = ceil(overlap * logband_hz);
    nwindows = round((min_trial_data - wlength)/overlap) + 1;
    windows_center = zeros(1, nwindows);
    plot_downsampling = 1;
    handles = [];
    cl = -inf;
    figure();
    fischers = zeros(nchannelsSelected, nbands, nwindows);
    typ = nan(nwindows, 1);
    for idx_w = 1:nwindows
        mean_c1 = nan(nchannelsSelected, nbands);
        mean_c2 = nan(nchannelsSelected, nbands);
        std_c1 = nan(nchannelsSelected, nbands);
        std_c2 = nan(nchannelsSelected, nbands);
        start_w = (idx_w-1)*overlap + 1;
        end_w = start_w + wlength - 1;
        if end_w > size(trial_data, 1)
            end_w = size(trial_data, 1);
        end

        for idx_ch=1:nchannelsSelected
            idx_chSel = channelsSelected(idx_ch);
            for idx_band=1:nbands
                mean_c1(idx_ch, idx_band) = mean(trial_data(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(1)), 'all');
                mean_c2(idx_ch, idx_band) = mean(trial_data(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(2)), 'all');
                std_c1(idx_ch, idx_band) = std(trial_data(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(1)), 0, 'all');
                std_c2(idx_ch, idx_band) = std(trial_data(start_w:end_w,idx_band,idx_chSel,trial_typ == classes(2)), 0, 'all');
            end
        end

        windows_center(idx_w) = (start_w + end_w) / 2;
        fischers(:,:,idx_w) = (mean_c1 - mean_c2) ./ sqrt(std_c1.^2 + std_c2.^2);
        if mod(idx_w, plot_downsampling) == 0
            % calculate fisher
            subplot(nrows, ceil(nwindows/nrows/plot_downsampling), idx_w/plot_downsampling)
            c_fisher = fischers(:,:,idx_w);

            % show the fisher features
            imagesc(c_fisher);
            ylabel('channels');
            ynew = 1:nchannelsSelected;
            set(gca, 'YTick',ynew, 'YTickLabel', channels_select);
            xlabel('frequencies');
            xnew = 1:nbands;
            x_label = bands_str;
            set(gca, 'XTick',xnew, 'XTickLabel',x_label);
            handles = [handles gca];
            cl = max(cl, max(abs(c_fisher), [], 'all'));
            if start_w > min(cueDUR)
                title(['cf | ' num2str(idx_w)])
            else
                title(['cue | ' num2str(idx_w)])
            end
        end

        if start_w > min(cueDUR)
            typ(idx_w) = cf_event;
        end
    end
    set(handles, 'clim', [-cl cl])
    sgtitle(['fisher score for some windows | ' subject ' | file: ' num2str(idx_file)]);

    % fisher for all cue - cf trial
    for idx_ch=1:nchannelsSelected
        idx_chSel = channelsSelected(idx_ch);
        for idx_band=1:nbands
            mean_c1(idx_ch, idx_band) = mean(trial_data(:,idx_band,idx_chSel,trial_typ == classes(1)), 'all');
            mean_c2(idx_ch, idx_band) = mean(trial_data(:,idx_band,idx_chSel,trial_typ == classes(2)), 'all');
            std_c1(idx_ch, idx_band) = std(trial_data(:,idx_band,idx_chSel,trial_typ == classes(1)), 0, 'all');
            std_c2(idx_ch, idx_band) = std(trial_data(:,idx_band,idx_chSel,trial_typ == classes(2)), 0, 'all');
        end
    end
    figure();
    fisher_all = (mean_c1 - mean_c2) ./ sqrt(std_c1.^2 + std_c2.^2);
    imagesc(fisher_all);
    ylabel('channels');
    ynew = 1:nchannelsSelected;
    set(gca, 'YTick',ynew, 'YTickLabel', channels_select);
    xlabel('frequencies');
    xnew = 1:nbands;
    x_label = bands_str;
    set(gca, 'XTick',xnew, 'XTickLabel',x_label);
    set(handles, 'clim', [-cl cl])
    title(['fisher score for cue + cf | ' subject ' | file: ' num2str(idx_file)]);

    %% metrics for the features selection
    n_features = 6;
    rho_fischer = nan(nchannelsSelected*nbands, nwindows);
    % take the best n features according to the overall trial
    [tot_sortedValues, tot_linearIndices] = maxk(fisher_all(:), n_features);
    [tot_rows, tot_cols] = ind2sub(size(fisher_all), tot_linearIndices);
    disp('ALL')
    disp(['   channels: ', strjoin(channels_select(tot_rows), ' ')]);
    disp(['   bands: ', strjoin(bands_str(tot_cols), ' ')])
    tot_select = zeros(size(fisher_all(:), 1), 1);
    tot_select(tot_linearIndices) = 1;
    for idx_w = 1:nwindows
        tmp = fischers(:,:,idx_w);
        rho_fischer(:,idx_w) = tmp(:);
    end

    overlap_ratio = zeros(1, nwindows);
    jaccard_similarity = zeros(1, nwindows);
    rho = zeros(1, nwindows);
    cos_similarity = zeros(1, nwindows);
    [~, fisher_Ranks] = sort(rho_fischer, 1, 'descend');

    tot_cos_similarity = zeros(1, nwindows);
    tot_jaccard_similarity = zeros(1, nwindows);
    tot_overlap_ratio = zeros(1, nwindows);
    disp('COMPARED')
    c_fischer = fischers(:,:,1); % take only interested values
    [c_sortedValues, c_linearIndices] = maxk(c_fischer(:), n_features);
    [c_rows, c_cols] = ind2sub(size(c_fischer), c_linearIndices);
    disp(['   window ' num2str(1) '/' num2str(nwindows)])
    disp(['      channels: ', strjoin(channels_select(c_rows), ' ')]);
    disp(['      bands: ', strjoin(bands_str(c_cols), ' ')])
    for idx_w=2:nwindows
        % take current fischer
        c_fischer = fischers(:,:,idx_w); % take only interested values
        [c_sortedValues, c_linearIndices] = maxk(c_fischer(:), n_features);
        [c_rows, c_cols] = ind2sub(size(c_fischer), c_linearIndices);
        disp(['   window ' num2str(idx_w) '/' num2str(nwindows)])
        disp(['      channels: ', strjoin(channels_select(c_rows), ' ')]);
        disp(['      bands: ', strjoin(bands_str(c_cols), ' ')])

        % take previous fischer
        p_fischer = fischers(:, :, idx_w - 1);
        [p_sortedValues, p_linearIndices] = maxk(p_fischer(:), n_features);

        % compute the difference in the consecutive windows
        n_intersection = length(intersect(c_linearIndices, p_linearIndices));
        n_union = length(union(c_linearIndices, p_linearIndices));
        c_select = zeros(size(c_fischer(:), 1), 1);
        c_select(c_linearIndices) = 1;
        n_select = zeros(size(p_fischer(:), 1), 1);
        n_select(p_linearIndices) = 1;

        % compute metrics for consecutive windows
        overlap_ratio(idx_w) = n_intersection / n_features;
        jaccard_similarity(idx_w) = n_intersection / n_union;
        rho(idx_w) = corr(fisher_Ranks(:,idx_w), fisher_Ranks(:,idx_w-1), 'Type', 'Spearman');
        cos_similarity(idx_w) = dot(c_select, n_select) / (norm(c_select) * norm(n_select));

        % compute metrics wrt the general good points
        tot_n_intersection = length(intersect(c_linearIndices, tot_linearIndices));
        tot_n_union = length(union(c_linearIndices, tot_linearIndices));
        tot_overlap_ratio(idx_w) = tot_n_intersection / n_features;
        tot_jaccard_similarity(idx_w) = tot_n_intersection / tot_n_union;
        tot_cos_similarity(idx_w) = dot(c_select, tot_select) / (norm(c_select) * norm(tot_select));
    end

    time_diff = zeros(1, nwindows);
    time_diff(2:end) = ((windows_center(1:end-1) + windows_center(2:end)) / 2) / sampleRate;

    figure();
    subplot(1,2,1);
    cf_start = find(typ == cf_event, 1);
    hold on;
    plot(overlap_ratio);
    plot(jaccard_similarity);
    plot(rho);
    plot(cos_similarity);
    plot([cf_start, cf_start], ylim);
    hold off
    xticks(1:5:size(time_diff,2))
    xticklabels(time_diff(1:5:end))
    xlabel('time [s]')
    legend({'overlap ration', 'jaccard similarity', ...
        'rho', 'cosine similarity', 'end cue'}, 'Location', 'best')
    title('difference consecutive windows')

    subplot(1,2,2)
    hold on
    plot(tot_overlap_ratio);
    plot(tot_jaccard_similarity);
    plot(tot_cos_similarity);
    plot([cf_start, cf_start], ylim);
    hold off
    xticks(1:5:size(time_diff,2))
    xticklabels(time_diff(1:5:end))
    xlabel('time [s]')
    legend({'overlap ration', 'jaccard similarity', ...
        'cosine similarity', 'end cue'}, 'Location', 'best')
    title('difference current windows with overall trial')

    sgtitle([subject ' | file: ' num2str(idx_file) ' | n features: ' num2str(n_features)])

    %% metrics for the features selection
    max_nfeatures = 10;
    nrows = 2;
    figure();
    for c_nf = 1:max_nfeatures
        wxw_img = zeros(nwindows);
        for idx_w1=1:nwindows
            c1_fischer = fischers(:,:,idx_w1);
            [c1_sortedValues, c1_linearIndices] = maxk(c1_fischer(:), c_nf);
            for idx_w2=1:nwindows
                c2_fischer = fischers(:,:,idx_w2);
                [c2_sortedValues, c2_linearIndices] = maxk(c2_fischer(:), c_nf);

                n_wxw_intersection = length(intersect(c1_linearIndices, c2_linearIndices));
                n_wxw_union = length(union(c1_linearIndices, c2_linearIndices));

                wxw_img(idx_w1,idx_w2) = n_wxw_intersection / c_nf;
            end
        end
        max_wxw = max(abs(wxw_img), [], 'all');
        wxw_img = wxw_img / max_wxw;
        max_wxw = max(abs(wxw_img), [], 'all');

        subplot(nrows, ceil(max_nfeatures/nrows), c_nf)
        hold on;
        imagesc(wxw_img);
        plot([cf_start, cf_start], ylim, 'k');
        hold off
        set(gca, 'clim', [0 max_wxw])
        xticks(1:7:size(windows_center,2))
        xticklabels(windows_center(1:7:end)/sampleRate)
        yticks(1:7:size(windows_center,2))
        yticklabels(windows_center(1:7:end)/sampleRate)
        title(['n features: ' num2str(c_nf)])
    end
    sgtitle([subject ' | overlap ratio between windows | file: ' num2str(idx_file)])
end