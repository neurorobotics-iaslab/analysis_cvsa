function [r2_values] = calc_r2_from_data(eeg_data, labels, varargin)
% CALC_R2_FROM_DATA Calcola Signed R-squared tra due classi.
%
% Usage:
%   r2 = calc_r2_from_data(X, y)
%   r2 = calc_r2_from_data(X, y, 'Plot', true, 'ChanLabels', {'C3','C4',...})
%
% Inputs:
%   eeg_data - Matrice [Samples x Channels]. I dati (potenza logaritmica).
%   labels   - Vettore [Samples x 1]. Le etichette di classe per ogni campione.
%              Deve contenere esattamente 2 classi (es. 730, 731).
%
% Optional Parameters (Name-Value pairs):
%   'Plot'       - true/false (default: false). Se true, genera il grafico.
%   'ChanLabels' - Cell array di stringhe. Nomi dei canali per l'asse X.
%   'title_data' - Stringa. Testo extra da aggiungere al titolo del grafico.
%
% Output:
%   r2_values - Vettore [1 x Channels]. Signed r^2 per ogni canale.
%               Positivo = Classe 1 > Classe 2.
%               Negativo = Classe 2 > Classe 1.

    %% 1. Parsing Input
    p = inputParser;
    addRequired(p, 'eeg_data', @isnumeric);
    addRequired(p, 'labels', @isnumeric);
    addParameter(p, 'Plot', false, @islogical);
    addParameter(p, 'ChanLabels', {}, @iscell);
    addParameter(p, 'title_data', '', @(x) ischar(x) || isstring(x));
    
    parse(p, eeg_data, labels, varargin{:});
    
    do_plot = p.Results.Plot;
    chan_labels = p.Results.ChanLabels;
    title_suffix = char(p.Results.title_data);

    % Controlli di base
    [n_samples, n_ch] = size(eeg_data);
    if length(labels) ~= n_samples
        error('Lunghezza labels (%d) non corrisponde alle righe di eeg_data (%d).', length(labels), n_samples);
    end

    % Identifica le due classi (rimuove 0 o NaN)
    classes = unique(labels);
    classes(classes == 0 | isnan(classes)) = [];
    
    if length(classes) ~= 2
        error('Input labels deve contenere esattamente 2 classi. Trovate: %d', length(classes));
    end
    
    c1 = classes(1); % es. 730
    c2 = classes(2); % es. 731
    
    % Dividi i dati
    X1 = eeg_data(labels == c1, :);
    X2 = eeg_data(labels == c2, :);
    
    n1 = size(X1, 1);
    n2 = size(X2, 1);
    
    if n1 < 5 || n2 < 5
        warning('Pochi campioni per classe (C1: %d, C2: %d). R^2 instabile.', n1, n2);
    end

    %% 2. Calcolo Matematico (Vettorizzato)
    % Point Biserial Correlation Coefficient squared (Signed)
    
    mu1 = mean(X1, 1);
    mu2 = mean(X2, 1);
    
    var1 = var(X1, 0, 1); % Varianza (normalizzata per N-1)
    var2 = var(X2, 0, 1);
    
    % Pooled Variance
    % var_pooled = [ (n1-1)*s1^2 + (n2-1)*s2^2 ] / (n1+n2-2)
    num = (n1 - 1) * var1 + (n2 - 1) * var2;
    den = n1 + n2 - 2;
    var_pooled = num / den;
    
    % Evita divisione per zero
    var_pooled(var_pooled < 1e-12) = 1e-12;
    
    % Point Biserial Correlation
    % r = (mu1 - mu2) / s_pooled * sqrt( (n1*n2) / (n1+n2)^2 )
    term_N = sqrt( (n1 * n2) / (n1 + n2)^2 );
    r = (mu1 - mu2) ./ sqrt(var_pooled) * term_N;
    
    % Signed R^2
    r2_values = sign(r) .* (r.^2);
    
    %% 3. Plotting (Opzionale)
    if do_plot
        figure('Color', 'w', 'Name', 'R2 Analysis Func');
        b = bar(r2_values);
        grid on; box on;
        ylabel(['Signed r^2 (Pos=' num2str(c1) ', Neg=' num2str(c2) ')']);
        title(['Class separability | ' title_suffix]);
        
        % Colora le barre: Rosso (Positivo/C1), Blu (Negativo/C2)
        b.FaceColor = 'flat';
        b.CData(r2_values > 0, :) = repmat([0.8 0.2 0.2], sum(r2_values > 0), 1);
        b.CData(r2_values < 0, :) = repmat([0.2 0.4 0.8], sum(r2_values < 0), 1);
        
        % Etichette canali
        if ~isempty(chan_labels)
            if length(chan_labels) == n_ch
                xticks(1:n_ch);
                xticklabels(strrep(chan_labels, 'EEG ', '')); % Pulisci label
                xtickangle(90);
                xlim([0, n_ch + 1]);
            end
        end
        
        % Linee soglia
        yline(0.05, 'k--', 'threshold 0.05');
        yline(-0.05, 'k--');
        xticks(1:length(chan_labels))
        xticklabels(chan_labels)
        
        % Stampa Top 5 in console
        [~, idx] = sort(abs(r2_values), 'descend');
        fprintf('\n--- TOP 5 CANALI ---\n');
        for i = 1:min(5, n_ch)
            ch = idx(i);
            val = r2_values(ch);
            if ~isempty(chan_labels), name = chan_labels{ch}; else, name = num2str(ch); end
            dir = ''; 
            if val > 0, dir = ['-> ' num2str(c1)]; else, dir = ['-> ' num2str(c2)]; end
            fprintf('%d. %s : r^2 = %+.4f %s\n', i, name, val, dir);
        end
    end
end