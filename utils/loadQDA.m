function qda = loadQDA(path_file)
    c_qda = yaml.ReadYaml(path_file);

    % Accessing the parameters
    qda.nclasses    = c_qda.QdaCfg.params.nclasses;
    qda.nfeatures   = c_qda.QdaCfg.params.nfeatures;
    qda.subject     = c_qda.QdaCfg.params.subject;
    qda.samplerate  = c_qda.QdaCfg.params.sampleRate;
    qda.filterOrder = c_qda.QdaCfg.params.filterOrder;

    priors = nan(1, qda.nclasses);
    classes = nan(1, qda.nclasses);
    for i = 1:qda.nclasses
        priors(i) = c_qda.QdaCfg.params.priors{i};
        classes(i) = c_qda.QdaCfg.params.classlbs{i};
    end
    qda.priors = priors;
    qda.classes = classes;

    bands = zeros(qda.nfeatures, 2);
    idchans = zeros(1, qda.nfeatures);
    lbchans = cell(1, qda.nfeatures);
    for i = 1:qda.nfeatures
        bands(i, :) = [c_qda.QdaCfg.params.band{i, 1}, c_qda.QdaCfg.params.band{i, 2}];
        idchans(i) = c_qda.QdaCfg.params.idchans{i};
        lbchans{i} = c_qda.QdaCfg.params.chans{i};
    end
    qda.bands = bands;
    qda.idchans = idchans;
    qda.chans = lbchans;

    c_cov = zeros(qda.nfeatures);
    c_rotations = zeros(qda.nfeatures);
    qda.covs = cell(1, qda.nclasses);
    qda.rotations = cell(1, qda.nclasses);
    for id_class = 1:size(c_qda.QdaCfg.params.covs, 2)
        for i = 1:qda.nfeatures
            for j = 1:qda.nfeatures
                c_cov(i,j) = cell2mat(c_qda.QdaCfg.params.covs{id_class}(i,j));
                c_rotations(i,j) = cell2mat(c_qda.QdaCfg.params.rotations{id_class}(i,j));
            end
        end
        qda.covs{id_class} = c_cov;
        qda.rotations{id_class} = c_rotations;
    end

    c_scalings = zeros(1, qda.nfeatures);
    c_mean     = zeros(1, qda.nfeatures);
    qda.scalings = cell(1, qda.nclasses);
    qda.means    = cell(1, qda.nclasses);
    for id_class = 1:size(c_qda.QdaCfg.params.covs, 2)
        for i = 1:qda.nfeatures
                c_scalings(i) = cell2mat(c_qda.QdaCfg.params.scalings(id_class, i));
                c_mean(i)     = cell2mat(c_qda.QdaCfg.params.means(id_class, i));
        end
        qda.scalings{id_class} = c_scalings;
        qda.means{id_class}    = c_mean;
    end

end