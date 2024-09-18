function fisher_score(signals, SelFreqs, SelChannels, psdSettings)

    nfiles = length(signals);
    %modalityId = [0 1];

    channelLb = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
        'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

    if nargin<4
        psdSettings.psdWlength = 0.5;
        psdSettings.wshift = 0.0625;
        psdSettings.mlength =  1;
        psdSettings.wconv =  'backward';
        psdSettings.pshift = 0.25;
        psdSettings.SampleRate = 512;
    end

    if nargin<3 || isempty(SelChannels)
        pickChannels = true(size(channelLb));
    else
        pickChannels = false(size(channelLb));
        for ch = SelChannels
             pickChannels = or(pickChannels, strcmpi(channelLb, ch));
        end
    end
    
    % Loading and concatenate PSD files and events
    psd = [];
    TYP = [];
    POS = [];
    DUR = [];
    Rk = []; 

    n_channels = length(channelLb);
    load('/home/paolo/lap_CVSA_39ch.mat');
    laplacianMask = laplacian;
    
    for fId = 1:nfiles
        
        data = signals{fId}; 
        signal_lap = data.s(:,1:n_channels)*laplacianMask;

        signal_lap = signal_lap(:,pickChannels);

        %% Evaluate PSD       
        
        [t_psd, freqs]    =  proc_spectrogram(signal_lap, psdSettings.psdWlength, psdSettings.wshift, psdSettings.pshift, psdSettings.SampleRate, psdSettings.mlength);
        t_psd         =  log(t_psd);
                
        events = events_raw2win(data.h.EVENT, psdSettings.wshift, psdSettings.SampleRate, psdSettings.wconv, psdSettings.mlength);
        
        % Extract the current events
        TYP = cat(1, TYP, events.TYP);
        DUR = cat(1, DUR, events.DUR);
        POS = cat(1, POS, events.POS + size(psd, 1));
        
        % Create Rk vector (run)
        t_Rk = fId*ones(size(t_psd, 1), 1);
        Rk = cat(1, Rk, t_Rk);

        psd = cat(1, psd, t_psd);

    end
    
    %% Extracting information from data
    events.TYP = TYP;
    events.DUR = DUR;
    events.POS = POS;
    
    nwindows  = size(psd, 1);

    
    %% Creating vector labels
    
    CFeedbackPOS = POS(TYP == 781);
    CFeedbackDUR = DUR(TYP == 781);
    
    CuePOS = POS(events.TYP == 730 | events.TYP == 731);
    CueDUR = DUR(events.TYP == 730 | events.TYP == 731);
    CueTYP = TYP(events.TYP == 730 | events.TYP == 731);
    
    FixPOS = POS(TYP == 786);
    FixDUR = DUR(TYP == 786);
    FixTYP = TYP(TYP == 786);
    
    NumTrials = length(CFeedbackPOS);
    
    % We consider the intersting period from Cue apperance to end of continuous feedback
    Ck = zeros(nwindows, 1);
    Tk = zeros(nwindows, 1);
    TrialStart = nan(NumTrials, 1);
    TrialStop  = nan(NumTrials, 1);
    for trId = 1:NumTrials
        cstart = CuePOS(trId);
        cstop  = CFeedbackPOS(trId) + CFeedbackDUR(trId) - 1;
        Ck(cstart:cstop) = CueTYP(trId);
        Tk(cstart:cstop) = trId;
        
        TrialStart(trId) = cstart;
        TrialStop(trId)  = cstop;
    end
    
    %% Apply log to the data (already done it)

    if nargin<2 || isempty(SelFreqs),   SelFreqs= freqs;    end
    
    [freqs, idfreqs] = intersect(freqs, SelFreqs);
    
    U = psd(:, idfreqs, :);
    
    NumWins  = size(U, 1);
    NumFreqs = size(U, 2);
    NumChans = size(U, 3);
    
    Runs = unique(Rk);
    NumRuns = length(Runs);
    
    %% Computing fisher score (for each run)
    disp('[proc] + Computing fisher score');
    Classes = [730 731];
    NumClasses = length(Classes);
    
    FisherScore = nan(NumFreqs, NumChans, NumRuns);
    FS2 = nan(NumFreqs*NumChans, NumRuns);
    for rId = 1:NumRuns
        rindex = Rk == Runs(rId); 
        
        cmu    = nan(NumFreqs, NumChans, 2);
        csigma = nan(NumFreqs, NumChans, 2);
        
        for cId = 1:NumClasses
            cindex = rindex & Ck == Classes(cId);
            cmu(:, :, cId) = squeeze(mean(U(cindex, :, :)));
            csigma(:, :, cId) = squeeze(std(U(cindex, :, :)));
        end
        
        FisherScore(:, :, rId) = abs(cmu(:, :, 2) - cmu(:, :, 1)) ./ sqrt( ( csigma(:, :, 1).^2 + csigma(:, :, 2).^2 ) );
    end
    
    %% Visualization Fisher score
    disp('[proc] |- Visualizing fisher score for offline runs');
    OfflineRuns = unique(Rk);
    NumCols = length(OfflineRuns);
    climits = [];
    handles = nan(length(OfflineRuns), 1);
    fig1 = figure;
    colormap('jet');
    for rId = 1:length(OfflineRuns)
        subplot(1, NumCols, rId);
        imagesc(FisherScore(:, :, OfflineRuns(rId))');
        axis square;
        colorbar;
        set(gca, 'XTick', 1:NumFreqs);
        set(gca, 'XTickLabel', freqs);
        set(gca, 'YTick', 1:NumChans);
        set(gca, 'YTickLabel', channelLb(pickChannels));
        xtickangle(90);
        xlabel('Hz');
        ylabel('channel');
        
        title(['Calibration run ' num2str(OfflineRuns(rId))]);
        
        climits = cat(2, climits, get(gca, 'CLim'));
        handles(OfflineRuns(rId)) = gca;
    end
    
    
    set(handles, 'clim', [0 max(max(climits))]);
    sgtitle('Fisher score');