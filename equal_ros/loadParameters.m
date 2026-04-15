function [ringBuffer, artifact, processing, qda_mi, qda_cvsa, integrator] = loadParameters(path_file)
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
    processing.do_hann = param.processing_bci_node.do_hann;

    % qda params
    if isfield(param, 'qda_node_mi')
        tmp = param.qda_node_mi.path_qda_model;
        [~, name, ext] = fileparts(tmp);
        qda_mi.file_name = [name, ext];
        qda_mi.name = 'QDA_MI';
    else
        warning('No qda_node_mi founded. Default values are used.');
        qda_mi.file_name = '';
        qda_mi.name = 'UNKNOWN';
    end
    if isfield(param, 'qda_node_cvsa')
        tmp = param.qda_node_cvsa.path_qda_model;
        [~, name, ext] = fileparts(tmp);
        qda_cvsa.file_name = [name, ext];
        qda_cvsa.name = 'QDA_CVSA';
    else
        warning('No qda_node_mi founded. Default values are used.');
        qda_cvsa.file_name = '';
        qda_cvsa.name = 'UNKNOWN';
    end 

    % integrator
    integrator.type = param.integrator.plugin(23:end);
    integrator.feedbackThs = cell2mat(param.integrator.thresholds);
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