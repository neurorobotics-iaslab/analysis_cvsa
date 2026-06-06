function [signal, header, basename] = load_gdf(gdf_path)
% LOAD_GDF  Read a GDF recording using BIOSIG's sload().
%   Returns the raw signal [N_samples x N_channels], the BIOSIG header
%   (with .SampleRate, .Label, .EVENT.POS/.TYP/.DUR), and the file
%   basename (no extension) — used to find the sibling YAML.
    if ~exist(gdf_path, 'file')
        error('load_gdf:notfound', 'GDF file not found:\n  %s', gdf_path);
    end
    log_step('load_gdf: reading "%s"', gdf_path);

    [signal, header] = sload(gdf_path);

    [~, basename, ~] = fileparts(gdf_path);
    log_step('load_gdf: %d samples x %d channels @ %.1f Hz, %d events, basename="%s"', ...
             size(signal, 1), size(signal, 2), header.SampleRate, ...
             numel(header.EVENT.POS), basename);
end
