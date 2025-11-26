function prob = apply_qda_matrix(qda, data_in)
    % Get number of samples
    n_samples = size(data_in, 1);
    
    % Class 1 computations
    invCov1 = inv(qda.covs{1,1});  % Precompute inverse covariance
    detCov1 = det(qda.covs{1,1});  % Precompute determinant
    meanClass1 = qda.means{1,1};   % Mean vector
    
    meanDiff1 = data_in - meanClass1; % Subtract mean (broadcasting)
    exponent1 = sum((meanDiff1 * invCov1) .* meanDiff1, 2); % Mahalanobis distance
    likelihoodClass1 = (1 / ((2 * pi)^(size(data_in, 2) / 2) * sqrt(detCov1))) ...
                        * exp(-0.5 * exponent1);
    
    % Class 2 computations
    invCov2 = inv(qda.covs{1,2});
    detCov2 = det(qda.covs{1,2});
    meanClass2 = qda.means{1,2};
    
    meanDiff2 = data_in - meanClass2;
    exponent2 = sum((meanDiff2 * invCov2) .* meanDiff2, 2);
    likelihoodClass2 = (1 / ((2 * pi)^(size(data_in, 2) / 2) * sqrt(detCov2))) ...
                        * exp(-0.5 * exponent2);
    
    % Compute unnormalized posterior probabilities
    unnormalizedProbClass1 = likelihoodClass1 * qda.priors(1);
    unnormalizedProbClass2 = likelihoodClass2 * qda.priors(2);

    % Normalize to get posterior probabilities
    totalProb = unnormalizedProbClass1 + unnormalizedProbClass2;
    probClass1 = unnormalizedProbClass1 ./ totalProb;
    probClass2 = unnormalizedProbClass2 ./ totalProb;

    % Handle NaN values
    probClass1(isnan(probClass1)) = 1;
    probClass2(isnan(probClass2)) = 0;
    probClass1(isnan(probClass2)) = 1 - probClass2(isnan(probClass2));
    probClass2(isnan(probClass1)) = 1 - probClass1(isnan(probClass1));

    % Return probability matrix
    prob = [probClass1, probClass2];
end
