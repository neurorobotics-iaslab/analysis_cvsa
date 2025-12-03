function [gmm, mu_data, sigma_data] = loadGMM(path)
try
    modelData = ReadYaml(path);
    params = modelData.GmmModelCfg.params;
    model_params = modelData.GmmModelCfg.model_params;
catch ME
    disp('Error in the loading of the YAML file. Is the file path correct? Do you have YAMLMatlab installed? Is the path correct?.');
    disp(ME.message);
    return;
end
mu_data = cellfun(@(x) x(1), params.mu);
sigma_data = cellfun(@(x) x(1), params.sigma);
nfeatures = model_params.nfeatures;

% model params
K = model_params.K;
weights = cell2mat(model_params.weights);
means_cell = model_params.means;
means = cell2mat(means_cell);
cov_cell = model_params.covariances; 
covariances = zeros(nfeatures, nfeatures, K);
for k = 1:K
    % Extract the k-th matrix
    matrix_cell = cov_cell{k}; 
    
    matrix_k = cell2mat(matrix_cell);
    covariances(:, :, k) = matrix_k;
end
gmm = gmdistribution(means, covariances, weights);
end