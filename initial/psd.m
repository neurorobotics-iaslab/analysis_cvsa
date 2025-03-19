%% compute the psd and show the fischer plot
clc; clearvars;

%%
signal = [];
TYP = [];
POS = [];
DUR = [];
nchannels = 39;

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

    events = header.EVENT;
    sampleRate = header.SampleRate;

    POS = [POS; events.POS + size(signal, 1)];
    TYP = [TYP; events.TYP];
    DUR = [DUR; events.DUR];
    signal = [signal; c_signal(1:events.POS(end)+events.DUR(end),:)];
end

%% Data processing
% apply PSD
disp('  [PROC] applaying pwelch');
% pwelch parameters 
wshift = 0.0625; % we have new signal after 0.0625 sec
pshift = 0.25; % shifting inside the windows
mlength = 1; % moving average
wlength = 0.5; % in order to have good estimation
wconv = 'backward';
f = 4:2:96; % in order to understand interesting features in a smaller
            %    number of frequencies
[PSD, freqs] = proc_spectrogram(signal, wlength, wshift,pshift,sampleRate, mlength);
[freq, idfreq] = intersect(f,freqs);
PSD = PSD(:,idfreq, :);

%% Update events from samples to windows
disp('  [PROC] updating events from samples to windows');
new_event.TYP = TYP;
new_event.POS = proc_pos2win(POS, wshift*sampleRate, wconv, mlength*sampleRate); 
new_event.conv = wconv;
new_event.SampleRate = sampleRate;
new_event.DUR = round(DUR/(wshift*sampleRate)) + 1; 

% other important informations
info.wshift = wshift;
info.pshift = pshift;
info.mlength = mlength;
info.wconv = wconv;
info.wlength = wlength;
info.freqs = freq;

%% General information
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
% channels_select = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
class_label = {'bl', 'br'};
classes = [730, 731];
n_classes = size(classes, 2);

n_trial = 20*nFiles; 
[~, channelsSelected] = ismember(channels_select, channels_label);

%% Labelling data 
events = new_event;
sampleRate = events.SampleRate;
cuePOS = events.POS(ismember(events.TYP, classes));
cueDUR = events.DUR(ismember(events.TYP, classes));
cueTYP = events.TYP(ismember(events.TYP, classes));

fixPOS = events.POS(events.TYP == 786);
fixDUR = events.DUR(events.TYP == 786);

cfPOS = events.POS(events.TYP == 781);
cfDUR = events.DUR(events.TYP == 781);

%% Trial extraction
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
n_freq = size(PSD, 2);
trial_data = nan(min_trial_data, n_freq, nchannels, n_trial);
basline_data = nan(min_fix_data, n_freq, nchannels, n_trial);
for trial = 1:n_trial
    c_start = trial_start(trial);
    c_end = trial_start(trial) + min_trial_data - 1;
    trial_data(:,:,:,trial) = PSD(c_start:c_end,:,:);
    c_fix_start = fix_start(trial);
    c_fix_end = fix_start(trial) + min_fix_data - 1;
    basline_data(:,:,:,trial) = PSD(c_fix_start:c_fix_end,:,:);
end


%% Calculate ERD/ERS 
ERD = nan(size(trial_data));
for trial=1:n_trial 
    for ch=1:nchannels
        for freq=1:n_freq
            baseline = mean(basline_data(:,freq,ch,trial));
            ERD(:,freq,ch,trial) = log(trial_data(:,freq,ch,trial)/baseline);
        end
    end
end

% %% Visualization
% figure();
% t = linspace(0, min_trial_data*info.wshift, min_trial_data);
% selclassLb = {'Bl', 'Br'};
% freqs = info.freqs;
% chandles = [];
% cl = -inf;
% for cId = 1:n_classes
%     climits = nan(2, length(channelsSelected));
%     for chId = 1:length(channelsSelected)
%         subplot(2, length(channelsSelected), (cId - 1)*length(channelsSelected) + chId);
%         cdata = mean(ERD(:, :, channelsSelected(chId), trial_typ == classes(cId)), 4);
%         imagesc(t, freqs, cdata');
%         set(gca,'YDir','normal');
%         climits(:, chId) = get(gca, 'CLim');
%         cl = max(cl, max(abs(cdata), [], 'all'));
%         chandles = cat(1, chandles, gca);
%         colormap(hot);
%         title([channels_label{channelsSelected(chId)} ' | ' selclassLb{cId}]);
%         xlabel('Time [s]');
%         ylabel('Frequency [Hz]');
%         line([1 1],get(gca,'YLim'),'Color',[0 0 0])
%         drawnow
%     end
% end
% colorbar;
% set(chandles, 'clim', [-cl cl])

%% Fisher score
% fisher score is channels x frequency
% calculate the mean for these session of the two classes
psd_hz = 16;
wlength = 0.5; %s
wlength = wlength * psd_hz;
overlap = 0.25; %s
overlap = overlap * psd_hz;
nwindows = round((min_trial_data - wlength)/overlap) + 1;
nrows = 2;
handles = [];
cl = -inf;

figure();
for idx_w = 1:nwindows
    start_w = (idx_w-1)*overlap + 1;
    end_w = start_w + wlength - 1;

    mean_c1 = nan(nchannels, n_freq);
    mean_c2 = nan(nchannels, n_freq);
    std_c1 = nan(nchannels, n_freq);
    std_c2 = nan(nchannels, n_freq);
    for ch=1:nchannels
        for freq=1:n_freq
            mean_c1(ch, freq) = mean(ERD(start_w:end_w,freq,ch,trial_typ == classes(1)), 'all');
            mean_c2(ch, freq) = mean(ERD(start_w:end_w,freq,ch,trial_typ == classes(2)), 'all');
            std_c1(ch, freq) = std(ERD(start_w:end_w,freq,ch,trial_typ == classes(1)), 0, 'all');
            std_c2(ch, freq) = std(ERD(start_w:end_w,freq,ch,trial_typ == classes(2)), 0, 'all');
        end
    end

    % calculate fisher
    fisher = (mean_c1 - mean_c2) ./ sqrt(std_c1.^2 + std_c2.^2);

    % show the fisher features
    start_freqs = find(freqs == 6);
    end_freqs = find(freqs == 16);
    subplot(nrows,ceil(nwindows/nrows),idx_w);
    imagesc(fisher(channelsSelected,start_freqs:end_freqs));
    ylabel('channels');
    ynew = 1:length(channelsSelected);
    set(gca, 'YTick',ynew, 'YTickLabel',channels_select);
    xlabel('frequencies');
    xnew = 1:n_freq;
    x_label = freqs;
    set(gca, 'XTick',xnew, 'XTickLabel',x_label);
    if start_w > min(cueDUR)
        title('cf')
    else
        title('cue')
    end
    cl = max(cl, max(abs(fisher(channelsSelected,start_freqs:end_freqs)), [], 'all'));
    handles = [handles gca];
end
set(handles, 'clim', [-cl cl])
sgtitle(['fisher score | ' subject]);