function fromGdf2Mat(gdf_path, mat_path, all_merged, name_merged)
    files = dir(fullfile(gdf_path, '*.gdf'));

    if all_merged
        signal = [];
        header.EVENT.POS = [];
        header.EVENT.DUR = [];
        header.EVENT.TYP = [];
    end

    for i=1:length(files)
        disp(['[INFO] Loading file' files(i).name]);
        file = fullfile(gdf_path, files(i).name);
        [s,h] = sload(file);

        [~, filename, ~] = fileparts(file);
        filename = [filename '.mat'];
        sfile = fullfile(mat_path,filename);
        disp(['   [INFO] saving the file: ', sfile]);
        save(sfile, 's', 'h');

        if all_merged
            header.EVENT.DUR = cat(1, header.EVENT.DUR, h.EVENT.DUR);
            header.EVENT.TYP = cat(1, header.EVENT.TYP, h.EVENT.TYP);
            header.EVENT.POS = cat(1, header.EVENT.POS, h.EVENT.POS + size(signal, 1));
            signal = cat(1, signal, s);
        end
    end

    if all_merged
        merge_file_name = fullfile(mat_path, name_merged);
        disp(['[INFO] saving the file: ', merge_file_name]);
        save(merge_file_name, 'signal', 'header');
    end
end