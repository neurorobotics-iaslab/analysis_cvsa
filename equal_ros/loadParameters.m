function [ringBuffer, artifact, processing, gmm, qda, integrator] = loadParameters(path_file)
    param = ReadYaml(path_file);

    % ring buffer params
    ringBuffer.name = param.RingBufferCfg.name;
    ringBuffer.size = param.RingBufferCfg.params.size;

    % artifact params
    artifact = param.ArtifactCfg.params;

    % processing params
    processing.bands = str2num(param.processing_cvsa_node.filters_band);
    processing.filterOrder = param.processing_cvsa_node.filter_order;
    processing.chunkSize = param.processing_cvsa_node.chunkSize;

    % gmm params
    tmp = param.gmm_node.path_gmm_model;
    [~, name, ext] = fileparts(tmp);
    gmm.file_name = [name, ext];
    gmm.name = 'GMM';

    % qda params
    tmp = param.qda_node.path_qda_model;
    [~, name, ext] = fileparts(tmp);
    qda.file_name = [name, ext];
    qda.name = 'QDA';

    % integrator
    integrator.type = param.integrator.plugin(23:end);
    integrator.ic_threshold = param.integrator.ic_threshold;
    integrator.init_val = param.integrator.init_val;
    integrator.ic_class_label = param.integrator.ic_class_label;
    integrator.feedbackThs = cell2mat(param.trainingCVSA_node.thresholds);
    if all(integrator.type == 'Buffer')
        integrator.bufferSize = param.integrator.buffer_size;
    elseif all(integrator.type == 'Exponential')
        integrator.alpha = param.integrator.alpha;
    end
end