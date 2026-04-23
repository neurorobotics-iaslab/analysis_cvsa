% you need yamlmatlab
function [ringBuffer, artifact, processing, qda_mi, qda_cvsa, integrator, paradigm] = loadParameters(path_file)
    param = yaml.ReadYaml(path_file);

    % ring buffer params
    if isfield(param, 'RingBufferCfg')
        ringBuffer.name = param.RingBufferCfg.name;
        ringBuffer.size = param.RingBufferCfg.params.size;
    else
        ringBuffer.name = 'RingBuffer';
        ringBuffer.size = 512;
    end

    % artifact params
    if isfield(param, 'ArtifactCfg')
        artifact = param.ArtifactCfg.params;
    else
        artifact = struct(); % Not available in calibration typically
    end

    % processing params
    if isfield(param, 'processing_power_node')
        processing.bands = str2num(param.processing_power_node.filters_band);
        processing.filterOrder = param.processing_power_node.filter_order;
        processing.chunkSize = param.processing_power_node.chunkSize;
        processing.do_hann = param.processing_power_node.do_hann;
    else
        % Default values for calibration
        processing.bands = [8 12];
        processing.filterOrder = 4;
        processing.chunkSize = 32;
        processing.do_hann = true;
    end

    % qda params
    paradigm = 'unknown';
    if isfield(param, 'qda_node_mi')
        tmp = param.qda_node_mi.path_qda_model;
        [~, name, ext] = fileparts(tmp);
        qda_mi.file_name = [name, ext];
        qda_mi.name = 'QDA_MI';
        qda_mi.path_to_model = tmp;
        paradigm = 'mi_lhrh';
    else
        warning('No qda_node_mi founded. Default values are used.');
        qda_mi.file_name = '';
        qda_mi.path_to_model = '';
        qda_mi.name = 'UNKNOWN';
    end
    if isfield(param, 'qda_node_cvsa')
        tmp = param.qda_node_cvsa.path_qda_model;
        [~, name, ext] = fileparts(tmp);
        qda_cvsa.file_name = [name, ext];
        qda_cvsa.name = 'QDA_CVSA';
        qda_cvsa.path_to_model = tmp;
        paradigm = 'cvsa_blbr';
    else
        warning('No qda_node_cvsa founded. Default values are used.');
        qda_cvsa.file_name = '';
        qda_cvsa.name = 'UNKNOWN';
    end

    if isfield(param, 'qda_node_mi') && isfield(param, 'qda_node_cvsa')
        paradigm = 'hybrid';
    end
    
    % If paradigm is still not defined via QDA, fallback to protocol (calibration case)
    if isfield(param, 'protocol') && isfield(param.protocol, 'task')
        if strcmp(param.protocol.task, 'cvsa_blbr')
            paradigm = 'cvsa_blbr';
        elseif strcmp(param.protocol.task, 'mi_lhrh')
            paradigm = 'mi_lhrh';
        elseif strcmp(param.protocol.task, 'hybrid_lr')
            paradigm = 'hybrid';
        end
    end

    % integrator
    if isfield(param, 'integrator')
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
    else
        integrator = struct();
    end
end