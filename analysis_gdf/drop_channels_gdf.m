clear all; clc;

% 1. Seleziona i file tramite finestra di dialogo
% 'MultiSelect', 'on' permette di selezionare più file tenendo premuto Ctrl o Shift
[files, folder] = uigetfile('*.gdf', 'Seleziona i file GDF da elaborare (taglio ultimi 3 canali)', 'MultiSelect', 'on');

% Se l'utente preme "Annulla", files sarà 0
if isequal(files, 0)
    disp('Selezione annullata.');
    return;
end

% Se viene selezionato un solo file, lo trasformiamo in cell per il ciclo for
if ischar(files)
    files = {files};
end

% Disabilita i warning fastidiosi per tutto il processo
warning('off', 'all');

% 2. Ciclo su ogni file selezionato
for i = 1:length(files)
    current_name = files{i};
    full_path = fullfile(folder, current_name);
    
    fprintf('Elaborazione di: %s...\n', current_name);
    
    % Caricamento
    [s, h] = sload(full_path);
    
    % Numero iniziale di canali
    ns_orig = h.NS;
    
    % Verifica che il file abbia almeno 3 canali
    if ns_orig <= 3
        fprintf(' -> ERRORE: il file ha solo %d canali, impossibile rimuoverne 3!\n', ns_orig);
        continue;
    end
    
    % 1. Taglia le ultime 3 colonne dal segnale
    s_new = s(:, 1:end-3);
    ns_new = ns_orig - 3;
    
    % 2. Inizializza un header PULITO per la scrittura (hdr_w)
    % Questo evita che Biosig si confonda con i parametri di lettura interni (come h.FID, h.FILE, ecc.)
    hdr_w = struct();
    
    % Assegna il nome del file di output fisso come richiesto dall'utente
    new_filename = fullfile(folder, 'prova_gdf.gdf');
    hdr_w.FileName = new_filename;
    hdr_w.TYPE = 'GDF';
    hdr_w.NS = ns_new;
    hdr_w.SampleRate = h.SampleRate;

    % Campi richiesti da sopen(...,'w') per il GDF:
    %   SPR  = samples per record (BIOSIG calcola HDR.AS.SampleRate da SPR; se manca → crash)
    %   NRec = -1 lascia che swrite lo determini scrivendo
    %   Dur  = durata di un record in secondi (SPR / SampleRate)
    %   T0   = data/ora di inizio (zeri se assente)
    if isfield(h, 'SPR') && ~isempty(h.SPR)
        hdr_w.SPR = h.SPR;
    else
        hdr_w.SPR = h.SampleRate;   % fallback: 1 record/sec
    end
    hdr_w.NRec = -1;
    if isfield(h, 'Dur') && ~isempty(h.Dur)
        hdr_w.Dur = h.Dur;
    else
        hdr_w.Dur = double(hdr_w.SPR) / double(h.SampleRate);
    end
    if isfield(h, 'T0') && ~isempty(h.T0)
        hdr_w.T0 = h.T0;
    else
        hdr_w.T0 = zeros(1, 6);
    end
    % GDFTYP: codice numerico (per canale) del tipo di campione (es. 16=int16, 17=float64)
    if isfield(h, 'GDFTYP') && length(h.GDFTYP) >= ns_new
        hdr_w.GDFTYP = h.GDFTYP(1:ns_new);
    else
        hdr_w.GDFTYP = repmat(16, 1, ns_new);   % default: int16
    end

    % Taglia e assegna le label dei canali
    hdr_w.Label = h.Label(1:ns_new);
    
    % Copia la tabella degli eventi
    if isfield(h, 'EVENT')
        hdr_w.EVENT = h.EVENT;
    end
    
    % Copia e taglia i parametri di calibrazione fisica e digitale
    fields_to_copy = {'PhysMax', 'PhysMin', 'DigMax', 'DigMin', 'Cal', 'Off', 'PhysDimCode', 'PhysDim'};
    for f_idx = 1:length(fields_to_copy)
        f_name = fields_to_copy{f_idx};
        if isfield(h, f_name)
            val = h.(f_name);
            if length(val) >= ns_new
                hdr_w.(f_name) = val(1:ns_new);
            end
        end
    end
    
    % Forza la definizione di PhysDimCode e PhysDim se non presenti o incompleti
    % (evita l'errore "HDR.PhysDimCode of the following channel(s) is(are) not defined")
    if ~isfield(hdr_w, 'PhysDimCode') || isempty(hdr_w.PhysDimCode) || length(hdr_w.PhysDimCode) < ns_new
        hdr_w.PhysDimCode = ones(1, ns_new) * 4275; % 4275 è il codice standard Biosig per microvolt (uV)
    end
    if ~isfield(hdr_w, 'PhysDim') || isempty(hdr_w.PhysDim) || length(hdr_w.PhysDim) < ns_new
        hdr_w.PhysDim = repmat({'uV'}, 1, ns_new);
    end
    
    % Copia e taglia i filtri
    if isfield(h, 'Filter')
        hdr_w.Filter = struct();
        filter_fields = {'HighPass', 'LowPass', 'Notch'};
        for f_idx = 1:length(filter_fields)
            f_name = filter_fields{f_idx};
            if isfield(h.Filter, f_name)
                val = h.Filter.(f_name);
                if length(val) >= ns_new
                    hdr_w.Filter.(f_name) = val(1:ns_new);
                end
            end
        end
    end
    
    % Apertura e Scrittura usando l'header pulito
    hdr_w = sopen(hdr_w, 'w'); 
    try
        swrite(hdr_w, s_new);
        hdr_w = sclose(hdr_w);
        fprintf(' -> Salvato con successo (%d -> %d canali) come: prova_gdf.gdf\n', ns_orig, ns_new);
    catch ME
        fprintf(' -> ERRORE durante il salvataggio: %s\n', ME.message);
        if isfield(hdr_w, 'FID') && hdr_w.FID > 0, sclose(hdr_w); end
    end
end

warning('on', 'all');
disp('--- Processo completato ---');
