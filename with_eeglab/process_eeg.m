%% eeg_lab processing
clc; clearvars;
addpath(genpath('/home/paolo/Local/Matlab/EEG_preprocessing_toolbox'));

%% start
subject = 'c7';
day = '20250217';

DATA_PATH = ['/home/paolo/cvsa_ws/record/' subject '/' day '/calibration/gdf'];
files = dir([DATA_PATH '/*.gdf']);

load('/home/paolo/chanlocs39.mat');



settings.eeg.fs = 512; % [Hz]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%____________________________LOOK AT ME_______________________________ %
% settings.eeg.notch_cutoff = [50 100 150,200];
settings.eeg.notch_cutoff = [];
settings.eeg.notch_bw = 35;
settings.eeg.cutoff_h = 4;
settings.eeg.cutoff_l= 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%____________________________LOOK AT ME_______________________________ %
% settings.eeg.cutoff_h = 1;
% settings.eeg.cutoff_l= 200;     
settings.eeg.filter_order = 2; 
settings.eeg.ftype = 'cheby';
settings.eeg.drop_chans = [];

settings.eeg.asr.th = 20; % [std]
settings.eeg.asr.block_size = 10;
settings.eeg.asr.window_length = 2; % [sec]
settings.eeg.asr.window_overlap = 0.5; % [%]
settings.eeg.asr.max_bad_channels = 0.1; %0.2
settings.eeg.asr.auto_calib = true;
settings.eeg.asr.iterative = true;
settings.eeg.asr.iter_th = 0.05; % [%]
settings.eeg.asr.max_iter = 20; % [%]
settings.step_epoch = [-0 0.75];
settings.eeg.ica.method = 'ica'; 

merged_outfile = [DATA_PATH, '/processed_merged.mat'];
asr_outfile = [DATA_PATH, '/processed_asr.mat'];
ica_outfile = [DATA_PATH, '/processed_', settings.eeg.ica.method ,'.mat'];
ica_clean_outfile = [DATA_PATH, '/processed_', settings.eeg.ica.method ,'_clean.mat'];

filenames = find(cellfun(@(x) contains(x,'.'),{files.name}'));

%%
% Load Data and store it in a EEGLAB Struct
all_eegs = {};
for i=1:length(filenames)
    filename = files(filenames(i)).name;
    disp(['Loading ', filename, '...']);
    [s,h] = sload([DATA_PATH,'/',filename]);
    s = s(:,1:39);

    EEG = make_EEGLab_struct(s',chanlocs',settings);
    EEG.settings = settings;
    EEG.event_gdf = h.EVENT;
    
    all_eegs{i} = EEG;
end

%% MERGE
disp('Preprocessing and Merging files..')

eeg_len = sum(arrayfun(@(x) size(x{1}.data,2),all_eegs));
eeg_nch = size(all_eegs{1}.data,1);
settings = all_eegs{1}.settings;
chlocs = all_eegs{1}.chanlocs;
eeg = zeros(eeg_nch,eeg_len);

all_events.TYP = [];
all_events.POS = [];
all_events.DUR = [];
last_eeg_idx = 1;
for i=1:length(all_eegs)
    el = preprocess(all_eegs{i});
    eeg(:,last_eeg_idx:last_eeg_idx+size(el.data,2)-1) = el.data;
    all_events.TYP = [all_events.TYP; el.event_gdf.TYP];
    all_events.POS = [all_events.POS; (el.event_gdf.POS + last_eeg_idx-1)];
    all_events.DUR = [all_events.DUR; el.event_gdf.DUR];
    last_eeg_idx = last_eeg_idx + size(el.data,2);
    disp(' ');
end

EEG = make_EEGLab_struct(eeg,el.chanlocs,settings);
EEG.settings = settings;

EEG.event_gdf = all_events;

epoch = 1;
for i = 1:length(all_events.POS)
    events(i).latency = all_events.POS(i);       % Set latency
    events(i).duration = all_events.DUR(i);     % Set duration
    events(i).type = all_events.TYP(i);% Set type
    events(i).epoch = 1;
end
EEG.event_eeglab = events;

triggers.type  = arrayfun(@(x) nameEvents(x),EEG.event_gdf.TYP,'UniformOutput',false);
for i=1:length(triggers.type)
    if strcmp(triggers.type{i}, 'cf')
        if strcmp(triggers.type{i-1}, 'cue_br')
            triggers.type{i} = 'cf_br';
        else
            triggers.type{i} = 'cf_bl';
        end
    end

end
triggers.num = EEG.event_gdf.TYP;
triggers.duration = EEG.event_gdf.DUR;
triggers.time  = EEG.event_gdf.POS;
EEG.triggers = triggers;
EEG.event = EEG.event_eeglab;

disp(['Saving stage to ',merged_outfile]);
save(merged_outfile, 'EEG', '-v7.3');

%% Compute ASR
disp('Computing ASR..')
EEG = compute_asr(EEG, settings);

disp(['Saving stage to ',asr_outfile]);
save(asr_outfile, 'EEG', '-v7.3');

%% ICA
% load(asr_outfile);
EEG1 = EEG;
drk = getrank(EEG1.data);

if strcmp(settings.eeg.ica.method , 'ica')
    [EEG1.icaweights,EEG1.icasphere]=runica(EEG1.data(:,:),'sphering','on','lrate',1e-4,'maxsteps',100);
    EEG.icaweights = EEG1.icaweights;
    EEG.icasphere = EEG1.icasphere;
    EEG.icawinv=inv(EEG1.icaweights*EEG1.icasphere);
elseif strcmp(settings.eeg.ica.method , 'amica')
    [W,S,mods] = runamica15(double(EEG1.data(:,:)),...
                        'outdir',	'amicaout',...
                    	'num_chans',	EEG1.nbchan',...
                        'pcakeep', drk, ...
                        'max_threads',	4,...
                        'pdftype',	2 ...
                        );
    EEG.icaweights = W;
    EEG.icasphere = S(1:size(W,1),:);
    EEG.icawinv = mods.A(:,:,1);
    EEG.mods = mods;    
end

EEG = eeg_checkset(EEG, 'ica');

disp(['Saving stage to ',ica_outfile]);
save(ica_outfile, 'EEG', '-v7.3');

%% IC selection - part I

% close all;clc;
load(ica_outfile);

epoch = settings.step_epoch;
epoch = [-1.5 5];
labels = {'cf_bl', 'cf_br'};
data_epoch = extractepochs(EEG.data',EEG.triggers,labels,epoch,EEG.srate);

EEG1 = EEG;
EEG1.data = permute(data_epoch, [2 1 3]);
EEG1.trials = size(data_epoch,3);
EEG1.xmin = epoch(1);
EEG1.xmax = epoch(2)-1/EEG1.srate;
EEG1.pnts = size(EEG1.data,2);
for i= 1:length(EEG.event)
EEG1.event(i).epoch = i; 
end
EEG1.icaact = [];
EEG1.icaact = geticaact(EEG1);
EEG1 = iclabel(EEG1);
pop_viewprops(EEG1, 0);

%% IC selection - part II
prompt = 'Enter the components to remove:\n';
rm_comps = (str2num(input(prompt, "s")))';

EEG = pop_subcomp(EEG, rm_comps);
EEG.etc.removed_comps = rm_comps;

disp(['Saving stage to ',ica_clean_outfile]);
save(ica_clean_outfile, 'EEG', '-v7.3');

%% epoching
EEG730 = pop_epoch(EEG, {730}, [-1.5 5.5]);
EEG731 = pop_epoch(EEG, {731}, [-1.5 5.5]);
EEG730 = pop_rmbase(EEG730, [-1500 0]);
EEG731 = pop_rmbase(EEG731, [-1500 0]);

%% show erp
% pop_erpimage(EEG1, 1)
pop_plottopo(EEG730)

%% show erd
channels_all = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
        'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
channels_select = {'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
[~, idx_channels_select] = ismember(channels_select, channels_all);

for i = idx_channels_select
    figure()
    ch = i;
    [ersp, itc, times, freqs, ~] = newtimef(EEG730.data(ch, :, :), EEG730.pnts, ...
        [EEG730.xmin EEG730.xmax] * 1000, EEG730.srate, 0, 'baseline', [-1000 0], ...
        'freqs', [4 20], 'plotersp', 'on');
    sgtitle(['channel: ' channels_all{i} ' | 730'])

    figure()
    ch = i;
    [ersp, itc, times, freqs, ~] = newtimef(EEG731.data(ch, :, :), EEG731.pnts, ...
        [EEG731.xmin EEG731.xmax] * 1000, EEG731.srate, 0, 'baseline', [-1000 0], ...
        'freqs', [4 20], 'plotersp', 'on');
    sgtitle(['channel: ' channels_all{i} ' | 731'])
end

%% funciton
function [name] = nameEvents(type)
    if type == 1, name = 'trial_start'; 
    elseif type == 786, name = 'fixation'; 
    elseif type == 781, name = 'cf';
    elseif type == 730, name = 'cue_bl';
    elseif type == 731, name = 'cue_br';
    elseif type == 897, name = 'hit';
    elseif type == 898, name = 'miss';
    elseif type == 899, name = 'timeout';
    end
end