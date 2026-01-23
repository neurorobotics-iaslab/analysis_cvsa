function [m] = computeMetrics(integratorCfg, artifact, gmm_prob, qda_prob, integrated_prob, event, event_start, gmm_classes, task_classes)
CODE_SX = 769; CODE_DX = 770; CODE_REST = 783;
CODE_HIT = 897; CODE_MISS = 898; CODE_TIMEOUT = 899;

cfPOS = event.POS(event.TYP == event_start);
cfDUR = event.DUR(event.TYP == event_start);
ntrials = length(cfPOS);

% Filtra solo gli eventi rilevanti per sincronizzare con i trial
all_cues = event.TYP(ismember(event.TYP, task_classes));
all_results = event.TYP(ismember(event.TYP, [CODE_HIT, CODE_MISS, CODE_TIMEOUT]));

% Safety Check
if length(all_cues) ~= ntrials || length(all_results) ~= ntrials
    ntrials = min([length(cfPOS), length(all_cues), length(all_results)]);
end

ic_index = find(gmm_classes == integratorCfg.ic_class_label);

% --- ACCUMULATORI ---
% Trial Counters
cnt_act_hit = 0; cnt_act_miss = 0; cnt_act_timeout = 0; cnt_act_total = 0;
cnt_rest_fp = 0; cnt_rest_ok = 0; cnt_rest_total = 0;

% Sample Counters
samp_corr_all = 0; samp_tot_all = 0;       % Tutti i trial attivi
samp_corr_noto = 0; samp_tot_noto = 0;     % Solo trial NO-TIMEOUT

% Time Lists (in samples)
list_time_hit = []; list_time_miss = [];
list_time_rest_fp = []; list_time_rest_ok = [];

% Stability Lists
list_wobble = []; list_wdr = []; list_rest_dev = [];

%% 2. CICLO DI ANALISI
for i = 1:ntrials
    idx_start = cfPOS(i);
    idx_end = cfPOS(i) + cfDUR(i) - 1;
    dur = idx_end - idx_start + 1;
    
    cue = all_cues(i);
    res = all_results(i);
    
    % Slice segnali
    s_art = artifact(idx_start:idx_end);
    s_qda = qda_prob(idx_start:idx_end, :);
    s_int = integrated_prob(idx_start:idx_end, 1); 
    
    % --- ACTIVE TASK (769/770) ---
    if cue == CODE_SX || cue == CODE_DX
        cnt_act_total = cnt_act_total + 1;
        
        % A. Trial Result & Time
        if res == CODE_HIT
            cnt_act_hit = cnt_act_hit + 1;
            list_time_hit = [list_time_hit; dur];
        elseif res == CODE_MISS
            cnt_act_miss = cnt_act_miss + 1;
            list_time_miss = [list_time_miss; dur];
        elseif res == CODE_TIMEOUT
            cnt_act_timeout = cnt_act_timeout + 1;
        end
        
        % B. Stability (Wobble & WDR)
        diffs = diff(s_int);
        list_wobble = [list_wobble; sum(abs(diffs))];
        
        is_target_sx = (cue == CODE_SX); % 1 se SX
        if is_target_sx, n_wrong = sum(diffs < -1e-6); else, n_wrong = sum(diffs > 1e-6); end
        list_wdr = [list_wdr; n_wrong / (length(diffs)+eps)];

        % C. Sample Accuracy (All vs No-Timeout)
        qda_idx = find(task_classes == cue);
        if qda_idx <= size(s_qda, 2)
            % Estrai samples validi (no artifact)
            valid_mask = (s_art == 0);
            if any(valid_mask)
                preds = s_qda(valid_mask, qda_idx) >= 0.5;
                n_corr = sum(preds);
                n_tot = length(preds);
                
                % 1. Accumula per "ALL"
                samp_corr_all = samp_corr_all + n_corr;
                samp_tot_all = samp_tot_all + n_tot;
                
                % 2. Accumula per "NO-TIMEOUT" solo se non è timeout
                if res ~= CODE_TIMEOUT
                    samp_corr_noto = samp_corr_noto + n_corr;
                    samp_tot_noto = samp_tot_noto + n_tot;
                end
            end
        end
        
    % --- REST TASK (783) ---
    elseif cue == CODE_REST
        cnt_rest_total = cnt_rest_total + 1;
        
        % A. Trial Result
        if res == CODE_TIMEOUT
            cnt_rest_ok = cnt_rest_ok + 1; % Successo (rimasto fermo)
            list_time_rest_ok = [list_time_rest_ok; dur];
        else
            cnt_rest_fp = cnt_rest_fp + 1; % Fallimento (attivato per sbaglio)
            list_time_rest_fp = [list_time_rest_fp; dur];
        end
        
        % B. Stability (Max Deviation)
        list_rest_dev = [list_rest_dev; max(abs(s_int - 0.5))];
    end
end

%% 3. CREAZIONE STRUTTURA OUTPUT

% --- Active Metrics ---
m.accuracy.trial.active.raw = cnt_act_hit / (cnt_act_total + eps);
m.accuracy.trial.active.no_timeout = cnt_act_hit / (cnt_act_hit + cnt_act_miss + eps);

m.accuracy.sample.active.all = samp_corr_all / (samp_tot_all + eps);
m.accuracy.sample.active.no_timeout = samp_corr_noto / (samp_tot_noto + eps);

m.time.active.hit_avg = mean(list_time_hit);
m.time.active.miss_avg = mean(list_time_miss);

m.stability.active.wobble_avg = mean(list_wobble);
m.stability.active.wdr_avg = mean(list_wdr);

% --- Rest Metrics ---
m.accuracy.trial.rest.acc = cnt_rest_ok / (cnt_rest_total + eps); % % Successo
m.accuracy.trial.rest.fpr = cnt_rest_fp / (cnt_rest_total + eps); % False Positive Rate

m.time.rest.fp_avg = mean(list_time_rest_fp); % Tempo medio all'errore
m.stability.rest.max_dev_avg = mean(list_rest_dev);

end