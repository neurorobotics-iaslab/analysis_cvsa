%% file which test if the QDA results form python and matlab are equal
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/512hz/utils')
addpath('/home/paolo/cvsa_ws/src/analysis_cvsa/utils')

subject = 'g2';
day = '20250221';
type_qda = 'sustained';
chosen_qda = 3;
% Define path
path_data = ['/home/paolo/cvsa_ws/record/', subject, '/', day, '/calibration/two_classifier/dataset/', type_qda, '/', ...
    type_qda, '_data', num2str(chosen_qda), '.mat'];

path_qda = '/home/paolo/cvsa_ws/record/g2/20250221/calibration/two_classifier/classifier/sustained/cl_best_sustained.yaml';

% Load .mat file
data = load(path_data);

% Assume fix_mat is a function that extracts struct fields correctly
info = data.info;
typ = info.typ;
ntrial_test = info.nTest;

% Get training labels
typ_train = typ(1:end - ntrial_test);

% Trial start indices
start_trials = info.startTrial; 
X = data.X; % full data
y = data.y; % full labels
start_trials = [start_trials; size(X, 1)]; % append end index

% Training indices
ntrial_train = length(typ) - ntrial_test + 1;

% Initialize empty training sets
X_test = [];
y_test = [];

test_trial_idx = ntrial_train:length(typ);

for i = 1:ntrial_test
    idx_trial = test_trial_idx(i);
    c_start = start_trials(idx_trial); 
    c_end = start_trials(idx_trial+1) - 1;
    
    % Extract and append data
    X_test = [X_test; X(c_start:c_end, :)];
    y_test = [y_test; y(c_start:c_end, :)];
end

qda = loadQDA(path_qda);
y_pred = apply_qda_matrix(qda, X_test);
y_test(y_test == 730) = 0;
y_test(y_test == 731) = 1;

[xROC, yROC, ~, AUC] = perfcurve(y_test, y_pred(:,2), 1);
