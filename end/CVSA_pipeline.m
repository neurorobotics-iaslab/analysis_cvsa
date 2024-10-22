%% Test the features extracted, the qda probs and the integrated probs only during CF
% works with only ne file
% clc; clear all; close all;

% loading data, gdf and classifier
subject = 'h8';
day = '/20240926'; % place / before the data
path_gdf = ['/home/paolo/cvsa_ws/record/' subject day '/gdf/evaluation'];
gdf_files = dir(fullfile(path_gdf, '*.gdf'));

yaml_QDA_path = ['/home/paolo/cvsa_ws/record/' subject day '/qda_' subject '.yaml'];
qda = loadQDA(yaml_QDA_path);

yaml_parameters_path = ['/home/paolo/cvsa_ws/record/' subject day '/bag'];

% variables for parameters given in the qda file + general parameters
filterOrder = qda.filterOrder;
sampleRate = qda.samplerate;
idx_selchs = qda.idchans;
bands = qda.bands;
frameSize = 32;
nchannels = 39;
bufferSize_data = 512;
sampleRate_ros = 16;
time_to_start_integrate = 0; % must be in sceonds
time_to_start_integrate = time_to_start_integrate * sampleRate;

% processing of gdf file
for idx_file=1:length(gdf_files)
    % load the gdf
    disp(['Loaded file: ' gdf_files(idx_file).name]);
    c_file = [path_gdf '/' gdf_files(idx_file).name];
    [s,header] = sload(c_file);
    s = s(:,1:nchannels);
    cueTYP = header.EVENT.TYP(ismember(header.EVENT.TYP, qda.classes));
    hit_miss_timeout = header.EVENT.TYP(ismember(header.EVENT.TYP, [897, 898, 899]));
    cfPOS_gdf = header.EVENT.POS(header.EVENT.TYP == 781);
    cfDUR_gdf = header.EVENT.DUR(header.EVENT.TYP == 781);
    ncf = size(cfDUR_gdf, 1);

    % load the parameters in the yaml file
    yaml_parameters_file = [yaml_parameters_path '/' gdf_files(idx_file).name(1:end-4) '.yaml'];
    yaml_parameters_file = yaml.ReadYaml(yaml_parameters_file);
    threshold = yaml_parameters_file.trainingCVSA_node.thresholds;
    threshold = [threshold{1}, threshold{2}];
    bufferSize_feedback = yaml_parameters_file.integrator.buffer_size;

    disp(['   Overall Accuracy: ' num2str(sum(header.EVENT.TYP == 897)) '/' num2str(sum(header.EVENT.TYP == 1))]);
    disp(['   Thresholds: (730) ' num2str(threshold(1)) ' - (731) ' num2str(threshold(2))])

    % define for the filtering
    zi_high = cell(1, length(bands));
    zi_low  = cell(1, length(bands));
    for i = 1:length(bands)
        zi_high{i} = [];
        zi_low{i}  = [];
    end

    % vairables for keeping the features and the pointers
    all_features = [];
    cfStart_features = []; % this is the same for the qda probs vector
    cfEnd_features = [];
    in_cf = false;
    idx_cf = 1;

    % define variable for simulate the buffer used in rosneuro
    buffer = nan(bufferSize_data, nchannels);
    
    % apply processing on all the data
    nchunk = floor(size(s,1)/frameSize);
    for idx_c=1:nchunk
        frame = s((idx_c-1)*frameSize+1:idx_c*frameSize,:);
        buffer(1:end-frameSize,:) = buffer(frameSize+1:end,:);
        buffer(end-frameSize+1:end, :) = frame;

        % check if the buffer is full
        if any(isnan(buffer))
            continue;
        end

        % iterate over bands
        features = nan(1, qda.nfeatures);
        for idx_band = 1:length(bands)
            c_band = bands(idx_band,:);
            [c_b_low, c_a_low] = butter(filterOrder, c_band(2)*(2/sampleRate),'low');
            [c_b_high, c_a_high] = butter(filterOrder, c_band(1)*(2/sampleRate),'high');
            c_zi_low = zi_low{idx_band};
            c_zi_high = zi_high{idx_band};

            % apply low and high pass filters
            [s_low, c_zi_low] = filter(c_b_low,c_a_low,buffer,c_zi_low);
            [tmp_data,c_zi_high] = filter(c_b_high,c_a_high,s_low,c_zi_high);
            zi_low{idx_band} = c_zi_low;
            zi_high{idx_band} = c_zi_high;

            % apply pow
            tmp_data = power(tmp_data, 2);

            % apply average
            tmp_data = mean(tmp_data, 1);

            % apply log
            tmp_data = log(tmp_data);

            idx_feature = idx_band; % we repeat the bands to have same index for band and features
            features(idx_feature) = tmp_data(idx_selchs(idx_feature));

        end
        % concatenate the features
        all_features = cat(1, all_features, features);

        if ~in_cf && idx_cf <= ncf
            if (idx_c-1)*frameSize+1 >= cfPOS_gdf(idx_cf) - frameSize + time_to_start_integrate
                in_cf = true;
                cfStart_features = cat(1, cfStart_features, size(all_features, 1) - 1); % the -1 is to anlign with rosneuro
            end
        else
            if idx_cf <= ncf
                if idx_c*frameSize >= cfPOS_gdf(idx_cf) + cfDUR_gdf(idx_cf) 
                    in_cf = false;
                    cfEnd_features = cat(1, cfEnd_features, size(all_features, 1) - 1);
                    idx_cf = idx_cf + 1;
                end
            end
        end
    end

    % variables to keep the qda prbs and the integrated probs
    matlab_integrated_onlyCF = [];
    matlab_qda_probs_onlyCF = [];
    matlab_new_pointer_startProbsCF = [];
    matlab_new_pointer_endProbsCF = [];

    % iterate over the continuous feedback
    for idx_cf=1:ncf
        % variables for the integrator
        idx_integrator = 1;
        integrator_buffer  = repmat([qda.classes(1), qda.classes(2)], 1, bufferSize_feedback/qda.nclasses);
        prob_integrated = [sum(integrator_buffer == qda.classes(1)), sum(integrator_buffer == qda.classes(2))] ./ bufferSize_feedback;
        matlab_integrated_onlyCF = cat(1, matlab_integrated_onlyCF, prob_integrated);

        % iterate over the features
        start_cf = cfStart_features(idx_cf);
        end_cf = cfEnd_features(idx_cf);
        matlab_new_pointer_startProbsCF = cat(1, matlab_new_pointer_startProbsCF, size(matlab_integrated_onlyCF, 1));
        for idx_feature=start_cf:end_cf
            feature = all_features(idx_feature,:);

            % apply the QDA
            prob = apply_qda(qda, feature);

            % update the integrator buffer
            if prob(1) >= prob(2)
                integrator_buffer(idx_integrator) = qda.classes(1);
            elseif prob(2) > prob(1)
                integrator_buffer(idx_integrator) = qda.classes(2);
            else
                disp('Found a nan value')
            end

            % update the index of the integrator
            if idx_integrator == bufferSize_feedback
                idx_integrator = 1;
            else
                idx_integrator = idx_integrator + 1;
            end

            % compute the integrate probability
            prob_integrated = [sum(integrator_buffer == qda.classes(1)), sum(integrator_buffer == qda.classes(2))] ./ bufferSize_feedback;

            % update the onlyCF vectors
            matlab_qda_probs_onlyCF = cat(1, matlab_qda_probs_onlyCF, prob);
            matlab_integrated_onlyCF = cat(1, matlab_integrated_onlyCF, prob_integrated);
        end
        matlab_new_pointer_endProbsCF = cat(1, matlab_new_pointer_endProbsCF, size(matlab_integrated_onlyCF, 1));
    end

    % plot the results
    x_trial = [matlab_new_pointer_startProbsCF matlab_new_pointer_startProbsCF]';
    y_trial = repmat([0, 1], size(matlab_new_pointer_startProbsCF, 1),1)';

    x_ths = repmat([0,matlab_new_pointer_endProbsCF(end)], 2, 1)';
    y_ths = [threshold; threshold];
    y_ths(:,2) = 1 - y_ths(:,2);

    y_qda_prob = matlab_qda_probs_onlyCF(:,1);
    y_qda_prob(y_qda_prob >= 0.5) = 1;
    y_qda_prob(y_qda_prob < 0.5)  = 0;
    x_qda_prob = linspace(1, matlab_new_pointer_endProbsCF(end), size(y_qda_prob, 1));

    figure()  
    hold on
    plot(matlab_integrated_onlyCF(:,1), 'Color','r')
    line(x_ths(:,1), y_ths(:,1), 'Color', 'b')
    line(x_ths(:,2), y_ths(:,2), 'Color', 'm')
    scatter(x_qda_prob, y_qda_prob, 20, [0, 0, 0], 'filled')
    line(x_trial,y_trial, 'Color', 'k')
    legend('Integrated probability', ['threshold ' num2str(qda.classes(1))], ['threshold ' num2str(qda.classes(2))],'QDA Probs', 'trials')
    for i=1:length(cueTYP)
        if hit_miss_timeout(i) == 897
            text_plot = [num2str(cueTYP(i)) ' | HIT'];
        elseif hit_miss_timeout(i) == 898
            text_plot = [num2str(cueTYP(i)) ' | MISS'];
        else
            text_plot = [num2str(cueTYP(i)) ' | TIMEOUT'];
        end
        if cueTYP(i) == 730
            text(matlab_new_pointer_startProbsCF(i), 0.95, text_plot)
        else
            text(matlab_new_pointer_startProbsCF(i), 0.05, text_plot)
        end
    end
%     For the time in the plot
%     dt = (matlab_new_pointer_endProbsCF - matlab_new_pointer_startProbsCF + 1);
%     time_trials = dt ./ sampleRate_ros;
    title('Integrated probs for each trial')

end







