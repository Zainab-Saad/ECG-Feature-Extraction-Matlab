display('Reading samples ECG signal from MIT-BIH Arrhythmia Database');
[ecg, Fs, tm] = rdsamp('mitdbase/106', 1);

display(['Reading and plotting annotations (human labels) of QRS complexes performed on the signals']);
[ann, type, subtype, chan, num, comments] = rdann('mitdbase/106', 'atr', 1);

%% PREPROCESSING
[R_peaks_val, R_peaks_ind, Q_peaks_ind, Q_peaks_val, ...
 S_peaks_ind, S_peaks_val, T_peaks_ind, T_peaks_val, delay] = pan_tompkin(ecg, Fs);

% RR Interval Calculation
count = 1;
RR_int = [];
for i = R_peaks_ind
    count = count + 1;
    if count > 2
        RR_int = [RR_int, (i - temp_i) / Fs];
    end
    temp_i = i;
end

% QS Interval Calculation
i = 1;
QS_int = [];
while i < length(S_peaks_ind)
    QS_int = [QS_int, (S_peaks_ind(i) - Q_peaks_ind(i)) / Fs];
    i = i + 1;
end

% Adjust Comments for Rhythm Identification
count = 1;
my_classes = {'N', 'VFL', 'VT'}; % Now includes VT
rhythm = comments(count);
while count < length(ann)
    if (type(count) == '+')
        rhythm = comments(count);
    end
    comments(count) = rhythm;
    count = count + 1;
end

%% Rhythm Segmentation and Feature Extraction
arrhythmiaData = struct();

count = 1;
while count <= length(comments)
    rhythm = cell2mat(comments(count));
    
    % Identify the rhythm type
    if length(rhythm) == 4 && all(rhythm == '(VFL')
        rhythmType = 'VFL';
    elseif length(rhythm) == 2 && all(rhythm == '(N')
        rhythmType = 'N';
    elseif length(rhythm) == 3 && all(rhythm == '(VT') % New class VT
        rhythmType = 'VT';
    else
        count = count + 1; % Skip unrecognized rhythms
        continue;
    end
    
    % Initialize fields if rhythmType doesn't exist
    if ~isfield(arrhythmiaData, rhythmType)
        arrhythmiaData.(rhythmType).R_peak_ind = [];
        arrhythmiaData.(rhythmType).Q_peak_ind = [];
        arrhythmiaData.(rhythmType).T_peak_ind = [];
        arrhythmiaData.(rhythmType).S_peak_ind = [];
    end
    
    % Find start and end of the rhythm section
    start_count = ann(count);
    while (count <= length(comments)) && ...
          (length(cell2mat(comments(count))) == length(rhythm)) && ...
          all(cell2mat(comments(count)) == rhythm)
        count = count + 1;
    end
    end_count = ann(count);
    
    % Update peak indices for the current rhythm
    arrhythmiaData.(rhythmType).R_peak_ind = [arrhythmiaData.(rhythmType).R_peak_ind, ...
        R_peaks_ind((R_peaks_ind > start_count) & (R_peaks_ind < end_count))];
    arrhythmiaData.(rhythmType).Q_peak_ind = [arrhythmiaData.(rhythmType).Q_peak_ind, ...
        Q_peaks_ind((Q_peaks_ind > start_count) & (Q_peaks_ind < end_count))];
    arrhythmiaData.(rhythmType).T_peak_ind = [arrhythmiaData.(rhythmType).T_peak_ind, ...
        T_peaks_ind((T_peaks_ind > start_count) & (T_peaks_ind < end_count))];
    arrhythmiaData.(rhythmType).S_peak_ind = [arrhythmiaData.(rhythmType).S_peak_ind, ...
        S_peaks_ind((S_peaks_ind > start_count) & (S_peaks_ind < end_count))];
end

% Calculate RR and QS intervals for each rhythm
rhythms = fieldnames(arrhythmiaData);
for i = 1:length(rhythms)
    rhythmType = rhythms{i};
    arrhythmiaData.(rhythmType).RR_int = calc_rr(arrhythmiaData.(rhythmType).R_peak_ind, Fs);
    arrhythmiaData.(rhythmType).QS_int = calc_qs(arrhythmiaData.(rhythmType).Q_peak_ind, ...
                                                 arrhythmiaData.(rhythmType).S_peak_ind, Fs);
end

%% Prepare Features and Labels
rhythmTypes = fieldnames(arrhythmiaData);
featureNames = {'R_peak_ind', 'Q_peak_ind', 'S_peak_ind', 'T_peak_ind', 'RR_int', 'QS_int'};
X = [];
Y = [];

% Assign numerical labels to each rhythm type
labelMapping = containers.Map(rhythmTypes, 1:length(rhythmTypes));

for i = 1:length(rhythmTypes)
    rhythmType = rhythmTypes{i};
    label = labelMapping(rhythmType);
    
    % Extract features for the current rhythm
    rhythmData = arrhythmiaData.(rhythmType);
    numObservations = length(rhythmData.R_peak_ind); % Assuming R_peak_ind determines observation count
    
    % Initialize a temporary feature matrix
    tempX = zeros(numObservations, length(featureNames));
    
    for j = 1:length(featureNames)
        feature = rhythmData.(featureNames{j});
        if ~isempty(feature)
            tempX(:, j) = feature(:);
        end
    end
    
    % Append the features and labels
    X = [X; tempX]; % Append features
    Y = [Y; label * ones(numObservations, 1)]; % Append corresponding labels
end

%% Split Data into Training and Testing Sets
rand_num = randperm(size(X, 1));
X = X(rand_num, :);
Y = Y(rand_num);

% Split into training and testing sets (80/20 split)
trainIdx = 1:round(0.8 * length(Y));
testIdx = round(0.8 * length(Y)) + 1:length(Y);

X_train = X(trainIdx, :);
y_train = Y(trainIdx);

X_test = X(testIdx, :);
y_test = Y(testIdx);


%% One-vs-All Multiclass Classification

% Unique classes in the labels
unique_classes = unique(Y);
num_classes = length(unique_classes);

% Initialize models and predictions
models = cell(num_classes, 1);
binary_predictions = zeros(length(Y), num_classes);

% Train a binary SVM for each class
for i = 1:num_classes
    fprintf('Training binary SVM for class %d vs all...\n', unique_classes(i));
    
    % Create binary labels: 1 for current class, -1 for all others
    binary_labels = -ones(size(Y));
    binary_labels(Y == unique_classes(i)) = 1;
    
    % Train SVM for the current class
    models{i} = fitcsvm(X_train, binary_labels(trainIdx), ...
        'KernelFunction', 'linear', ...
        'ClassNames', [-1, 1], ...
        'BoxConstraint', 1, ...
        'Standardize', true);
end

%% Testing and One-vs-All Prediction
binary_predictions = zeros(size(X_test, 1), num_classes); % Preallocate

for i = 1:num_classes
    % Predict scores for each class
    [~, scores] = predict(models{i}, X_test);

    % Assign only the positive class score
    binary_predictions(:, i) = scores(:, 2);
end

% Assign the class with the highest score as the prediction
[~, predicted_labels] = max(binary_predictions, [], 2);

%% Evaluate Performance
test_accuracy = sum(predicted_labels == Y(testIdx)) / length(Y(testIdx)) * 100;

fprintf('Test accuracy (One-vs-All): %.2f%%\n', test_accuracy);
