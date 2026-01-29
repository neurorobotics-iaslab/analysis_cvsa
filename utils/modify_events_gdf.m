% 1. Definisci il nome del file originale
clear all; clc;
filename = '/home/paolo/cvsa/ic_cvsa_ws/record_mi/c7/old/20260116/calibration/c7.20260116.145201.calibration.mi_bhbf.gdf';

[s, h] = sload(filename);

% 2. Modifica i trigger
h.EVENT.TYP(h.EVENT.TYP == 769) = 771;
h.EVENT.TYP(h.EVENT.TYP == 770) = 773;

% 3. Prepara il nuovo nome file
[path, name, ext] = fileparts(filename);
new_filename = fullfile(path, [name, '_copy', ext]);

% 4. Pulizia e Preparazione Header per la scrittura
% Alcuni campi vecchi possono mandare in confusione sopen in scrittura
h.FileName = new_filename;
h.TYPE = 'GDF';

% 5. Apertura, Scrittura e Chiusura
% Usiamo la sintassi HDR = SOPEN(HDR, 'w') che è la più compatibile
h = sopen(h, 'w'); 

try
    swrite(h, s);
    h = sclose(h);
    fprintf('File salvato correttamente: %s\n', new_filename);
catch ME
    if isfield(h, 'FID') && h.FID > 0, sclose(h); end
    rethrow(ME);
end
