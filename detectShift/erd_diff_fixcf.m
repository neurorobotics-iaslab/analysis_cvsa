%% Look difference from erd between mean and std of fix with respect to cue and cf
% close all; clear all; clc;

%% Define variables
% file info
c_subject = 'f2'; % input
prompt = 'Enter "calibration" or "evaluation": ';
test_typ = input(prompt, 's');
day = ''; % input

path = ['/home/paolo/cvsa_ws/record/' c_subject day];
path_gdf = [path '/gdf/' test_typ];
chanlocs_path = '/home/paolo/new_chanlocs64.mat';

if strcmp(test_typ, "evaluation")
    feature_file = [path '/dataset/fischer_scores_ev.mat'];
    logband_file = [path '/dataset/logband_power_ev.mat'];
else
    feature_file = [path '/dataset/fischer_scores.mat'];
    logband_file = [path '/dataset/logband_power.mat'];
end

classes = [730,731];
nclasses = length(classes);

load(chanlocs_path);
files = dir(fullfile(path_gdf, '*.gdf'));  

band = {[8 14]};
nbands = length(band);

s=[]; events = struct('TYP',[],'POS',[],'SampleRate',512,'DUR',[]); Rk=[];

channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
       '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
 
% channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};


%% Concatenate files
for i=1:length(files)
    file = fullfile(path_gdf, files(i).name);

    %load(file);
    [signal,header] = sload(file);
    
    curr_s = signal(:,1:39);
    curr_h = header.EVENT;
    
    if strcmp(test_typ, "calibration") || (strcmp(test_typ, "evaluation"))
        start = find(curr_h.TYP == 1,1,'first');
        curr_h.TYP = curr_h.TYP(start:end);
        curr_h.POS = curr_h.POS(start:end);
        curr_h.DUR = curr_h.DUR(start:end);
    end
    % Create Rk vector (run)
    cRk = i*ones(size(curr_s,1),1);
    Rk = cat(1,Rk,cRk);
    % Concatenate events
    events.TYP = cat(1, events.TYP, curr_h.TYP);
    events.DUR = cat(1, events.DUR, curr_h.DUR);
    events.POS = cat(1, events.POS, curr_h.POS + size(s, 1));
    s = cat(1, s, curr_s);
end

%% Processing
% Create Vector labels
[nsamples,nchannels] = size(s);
[feedb_pos, feedb_dur, fix_dur, fix_pos, cue_dur, cue_pos, ntrials] = extract_info_label(events, 781, 786, [730 731]);

% Extract trial data
[TrialStart, TrialStop, FixStart, FixStop, Ck, Tk] = extract_trial_info(s, events, fix_pos, fix_dur, feedb_pos, feedb_dur, cue_pos, ntrials);

% Applay the filtering
s_processed = NaN(nsamples,nchannels,nbands);
for f_idx=1:nbands
    sel_band = band{f_idx}; %Hz
    t_window = 1; %[s]
    windowSize = events.SampleRate*t_window;
    filtOrder = 4;
    s_movavg = data_processing(s, nchannels, events.SampleRate, sel_band, filtOrder, t_window);
    s_processed(:,:,f_idx) = s_movavg;
end

% Trial extraction
trial_dur = min(TrialStop-TrialStart);
dataforTrial = NaN(trial_dur,nchannels,nbands,ntrials);
tCk = zeros(ntrials,1);
for trId=1:ntrials
    cstart = TrialStart(trId);
    cstop = cstart + trial_dur - 1;

    dataforTrial(:,:,:,trId) = s_processed(cstart:cstop,:,:);
    tCk(trId) = unique(nonzeros(Ck(cstart:cstop)));
end

% % Baseline extraction
% minFix_dur = min(FixStop - FixStart);
% reference = NaN(minFix_dur, nchannels, nbands, ntrials);
% for trId=1:ntrials
%     cstart = FixStart(trId); %=TrialStart(trId) o fix_pos(trId)
%     cstop = cstart+ minFix_dur - 1;
%     reference(:,:,:,trId) = s_processed(cstart:cstop,:,:);
% end
% baseline = repmat(mean(reference),[size(dataforTrial,1) 1 1]);

%% Compute LogBandPower [samples x channels x bands] 
% ERDfortrial = log(dataforTrial./baseline); % use the 4th dim for trials
ERDfortrial = log(dataforTrial);

%% Check difference from mean
covert_channels = find(~cellfun('isempty', channels_label));
cue_typ = events.TYP(events.TYP == 730 | events.TYP == 731);
for idx_b = 1:nbands
    fixation_period = mean(ERDfortrial(1:fix_dur, covert_channels, idx_b, :), 4);
    fix_mean = mean(fixation_period);
    fix_std = std(fixation_period);

    mean_1 = mean(ERDfortrial(:, covert_channels,idx_b,cue_typ==730),4);
    mean_2 = mean(ERDfortrial(:, covert_channels,idx_b,cue_typ==731),4);

    figure();
    
    
    for idx_ch = 1:size(covert_channels,2)
        subplot(5,4,idx_ch);
        l = {};
        hold on;
        plot(mean_1(:,idx_ch));
        l(end+1) = {channels_label{covert_channels(idx_ch)}};
        plot(xlim, [fix_mean(idx_ch), fix_mean(idx_ch)]);
        l(end+1) = {'mean on fix'};
        plot(xlim, [fix_mean(idx_ch) - 3 * fix_std(idx_ch), fix_mean(idx_ch)- 3 * fix_std(idx_ch)], 'k', 'LineWidth', 2);
        l(end+1) = {'std on fix'};
        plot(xlim, [fix_mean(idx_ch) + 3 * fix_std(idx_ch), fix_mean(idx_ch)+ 3 * fix_std(idx_ch)], 'k', 'LineWidth', 2);
        l(end+1) = {'std on fix'};
        plot([fix_dur(1), fix_dur(1)], ylim, 'k', 'LineWidth', 2);
        l(end+1) = {'cue'};
        plot([fix_dur(1)+cue_dur(1), fix_dur(1)+cue_dur(1)], ylim, 'k', 'LineWidth', 2);
        l(end+1) = {'cf'};
        plot([fix_dur(1)+cue_dur(1) + 0.625*512, fix_dur(1)+cue_dur(1) + 0.625*512], ylim, 'k--', 'LineWidth', 2);
        l(end+1) = {'cvsa'};
        legend(l)
        title('class 730')
        hold off;
    end

    figure();
    for idx_ch = 1:size(covert_channels,2)
        subplot(5,4,idx_ch);
        l = {};
        hold on;
        plot(mean_2(:,idx_ch));
        l(end+1) = {channels_label{covert_channels(idx_ch)}};
        plot(xlim, [fix_mean(idx_ch), fix_mean(idx_ch)]);
        l(end+1) = {'mean on fix'};
        plot(xlim, [fix_mean(idx_ch) - 3 * fix_std(idx_ch), fix_mean(idx_ch)- 3 * fix_std(idx_ch)], 'k', 'LineWidth', 2);
        l(end+1) = {'std on fix'};
        plot(xlim, [fix_mean(idx_ch) + 3 * fix_std(idx_ch), fix_mean(idx_ch)+ 3 * fix_std(idx_ch)], 'k', 'LineWidth', 2);
        l(end+1) = {'std on fix'};
        plot([fix_dur(1), fix_dur(1)], ylim, 'k', 'LineWidth', 2);
        l(end+1) = {'cue'};
        plot([fix_dur(1)+cue_dur(1), fix_dur(1)+cue_dur(1)], ylim, 'k', 'LineWidth', 2);
        l(end+1) = {'cf'};
        plot([fix_dur(1)+cue_dur(1) + 0.625*512, fix_dur(1)+cue_dur(1) + 0.625*512], ylim, 'k--', 'LineWidth', 2);
        l(end+1) = {'cvsa'};
        legend(l)
        title('class 731')
        hold off;
    end

    sgtitle(['band: ' num2str(band{idx_b}(1)) '-' num2str(band{idx_b}(2))]);
end