clear all; clc;

% 1. Seleziona i file tramite finestra di dialogo
% 'MultiSelect', 'on' permette di selezionare più file tenendo premuto Ctrl o Shift
[files, folder] = uigetfile('*.gdf', 'Seleziona i file GDF da convertire', 'MultiSelect', 'on');

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
    
    % Modifica i trigger
    h.EVENT.TYP(h.EVENT.TYP == 769) = 771;
    h.EVENT.TYP(h.EVENT.TYP == 770) = 773;
    
    % Prepara nuovo nome con _copy
    [~, name, ext] = fileparts(current_name);
    new_filename = fullfile(folder, [name, '_copy', ext]);
    
    % Preparazione Header
    h.FileName = new_filename;
    h.TYPE = 'GDF';
    
    % Apertura e Scrittura
    h = sopen(h, 'w'); 
    try
        swrite(h, s);
        h = sclose(h);
        fprintf(' -> Salvato come: %s_copy%s\n', name, ext);
    catch ME
        fprintf(' -> ERRORE su %s: %s\n', current_name, ME.message);
        if isfield(h, 'FID') && h.FID > 0, sclose(h); end
    end
end

warning('on', 'all');
disp('--- Processo completato ---');