function prob = apply_qda(qda, data_in)
    % Class 1
    invCov1 = inv(qda.covs{1,1});
    detCov1 = det(qda.covs{1,1});
    meanClass1 = qda.means{1,1};
    meanDiff1 = data_in - meanClass1;
    likelihoodClass1 = (1 / ((2 * pi)^(length(data_in) / 2) * sqrt(detCov1))) * exp(-0.5 * (meanDiff1 * invCov1 * meanDiff1'));

    % Class 2 calculations
    invCov2 = inv(qda.covs{1,2});
    detCov2 = det(qda.covs{1,2});
    meanClass2 = qda.means{1,2};
    meanDiff2 = data_in - meanClass2;
    likelihoodClass2 = (1 / ((2 * pi)^(length(data_in) / 2) * sqrt(detCov2))) * exp(-0.5 * (meanDiff2 * invCov2 * meanDiff2'));

    % Compute the unnormalized posterior probabilities
    unnormalizedProbClass1 = likelihoodClass1 * qda.priors(1);
    unnormalizedProbClass2 = likelihoodClass2 * qda.priors(2);

    % Normalize to get posterior probabilities
    totalProb = unnormalizedProbClass1 + unnormalizedProbClass2;
    probClass1 = unnormalizedProbClass1 / totalProb;
    probClass2 = unnormalizedProbClass2 / totalProb;

    prob = [probClass1, probClass2];
end