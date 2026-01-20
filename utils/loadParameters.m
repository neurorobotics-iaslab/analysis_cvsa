function [ringBuffer, artifact, processing, gmm, qda, integrator] = loadParameters(path_file)
    param = ReadYaml(path_file);

    % ring buffer params
    ringBuffer.name = param.RingBufferCfg.name;
    ringBuffer.size = param.RingBufferCfg.params.size;

    % artifact params
    artifact = param.ArtifactCfg.params;

    % processing params
    processing.bands = str2num(param.processing_bci_node.filters_band);
    processing.filterOrder = param.processing_bci_node.filter_order;
    processing.chunkSize = param.processing_bci_node.chunkSize;

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
    integrator.ic_class_label = param.integrator.ic_class_label;
    if all(param.protocol.task(1:2) == 'mi')
        integrator.feedbackThs = cell2mat(param.trainingwheel.thresholds);
    elseif all(param.protocol.task(1:4) == 'cvsa')
        integrator.feedbackThs = cell2mat(param.trainingCVSA_node.thresholds); 
    end
    
    integrator.k_gain = param.integrator.k_gain;
    integrator.increment_type = param.integrator.increment;
    
    if all(integrator.type == 'Buffer')
        integrator.init_val = param.integrator.init_val;
        integrator.bufferSize = param.integrator.buffer_size;
    elseif all(integrator.type == 'Exponential')
        integrator.init_val = param.integrator.init_percentual;
        integrator.alpha = param.integrator.alpha;
    end
end