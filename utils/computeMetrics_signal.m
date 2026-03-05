function [m] = computeMetrics_signal(gmm_prob, qda_gmm_prob, qda_classic_prob, artifact, event, event_start, task_classes, gmm_classes, integratorCfg)
% COMPUTEMETRICS_COMPARISON - Confronta QDA Classico e QDA Focus-Trained su diverse soglie GMM.

    % 1. Parametri e Inizializzazione
    ic_index = find(gmm_classes == integratorCfg.ic_class_label);
    ths = 0:0.1:0.9; % Soglie da "All" a ">0.9"
    m.ths_labels = [{'All'}, arrayfun(@(x) ['>' num2str(x)], ths, 'UniformOutput', false)];
    n_ths = length(m.ths_labels);
    
    % Accumulatori: [Soglia]
    m.classic_correct = zeros(1, n_ths);
    m.gmm_inf_correct = zeros(1, n_ths);
    m.total_samples   = zeros(1, n_ths);

    % Estrazione Trial
    cfPOS = event.POS(event.TYP == event_start);
    cfDUR = event.DUR(event.TYP == event_start);
    all_cues = event.TYP(ismember(event.TYP, task_classes)); % Solo classi attive (SX/DX)

    % 2. Loop sui Trial (Solo Active Task per accuratezza QDA)
    for i = 1:length(cfPOS)
        idx_range = cfPOS(i) : min(cfPOS(i) + cfDUR(i) - 1, size(gmm_prob, 1));
        cue = all_cues(i);
        if ismember(cue, task_classes(1:2)) % active classes
            is_sx = (cue == task_classes(1));

            s_gmm = gmm_prob(idx_range, ic_index);
            s_qda_gmm = qda_gmm_prob(idx_range, 1);     % Assumiamo colonna 1 = SX
            s_qda_cls = qda_classic_prob(idx_range, 1);
            s_art = artifact(idx_range);

            for t = 1:length(s_gmm)
                if s_art(t) == 1 || isnan(s_gmm(t)), continue; end

                conf = s_gmm(t);
                % Correttezza dei due modelli
                corr_gmm_inf = (s_qda_gmm(t) >= 0.5) == is_sx;
                corr_classic = (s_qda_cls(t) >= 0.5) == is_sx;

                % Aggiornamento per ogni soglia
                for s = 1:n_ths
                    % La prima soglia (s=1) è "All", le altre seguono il vettore ths
                    if s == 1 || conf > ths(s-1)
                        m.total_samples(s) = m.total_samples(s) + 1;
                        m.gmm_inf_correct(s) = m.gmm_inf_correct(s) + corr_gmm_inf;
                        m.classic_correct(s) = m.classic_correct(s) + corr_classic;
                    end
                end
            end
        end
    end

    % 3. Calcolo Accuratezze Finali
    m.acc_classic = (m.classic_correct ./ m.total_samples) * 100;
    m.acc_gmm_inf = (m.gmm_inf_correct ./ m.total_samples) * 100;

    pe = 0.5; % Probabilità casuale per 2 classi bilanciate
    m.kappa_classic = (m.acc_classic/100 - pe) ./ (1 - pe);
    m.kappa_gmm_inf = (m.acc_gmm_inf/100 - pe) ./ (1 - pe);
    m.kappa_classic = max(0, m.kappa_classic);
    m.kappa_gmm_inf = max(0, m.kappa_gmm_inf);
    
    % Errore Standard (per le barre d'errore nel plot)
    m.std_classic = 100 * sqrt((m.acc_classic/100 .* (1 - m.acc_classic/100)) ./ m.total_samples);
    m.std_gmm_inf = 100 * sqrt((m.acc_gmm_inf/100 .* (1 - m.acc_gmm_inf/100)) ./ m.total_samples);
end