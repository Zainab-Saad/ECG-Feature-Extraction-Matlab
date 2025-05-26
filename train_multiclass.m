folder = "databases/";
fileList = dir(fullfile(folder, '*.hea'));
fileList = {fileList.name};

% Initialize variables to collect all data
all_X = [];
all_Y = [];
arrhythmiaData = struct();

% Iterate over each file
for i = 1:length(fileList)
    recordname = str2mat(fullfile(folder, fileList{i}(1:end-4))); % Remove file extension
    
    % Display file being processed
    display(['Reading ECG signal from file: ', recordname]);
    
    % Read ECG signal and annotations
    [ecg, Fs, tm] = rdsamp(recordname, 1);
    [ann, type, subtype, chan, num, comments] = rdann(recordname, 'atr', 1);
    
    %% Preprocessing (extract peaks and features)
    [R_peaks_val, R_peaks_ind, Q_peaks_ind, Q_peaks_val, ...
     S_peaks_ind, S_peaks_val, T_peaks_ind, T_peaks_val, delay] = pan_tompkin(ecg, Fs,0);

    % Calculate RR and QS intervals
    RR_int = []; QS_int = [];
    count = 1; temp_i = R_peaks_ind(1);
    for i = R_peaks_ind(2:end)
        RR_int = [RR_int, (i - temp_i) / Fs];
        temp_i = i;
    end
    for i = 1:length(S_peaks_ind)
        QS_int = [QS_int, (S_peaks_ind(i) - Q_peaks_ind(i)) / Fs];
    end

    % Rhythm Identification
    rhythm = comments(1);
    count = 1;
    my_classes = {'N', 'AFIB'};
    while count < length(ann)
        if (type(count) == '+')
            rhythm = comments(count);
        end
        comments(count) = rhythm;
        count = count + 1;
    end

    %% Feature Extraction and Labeling
    count = 1;
    while count <= length(comments)
        rhythm = cell2mat(comments(count));
        
        % Assign rhythm type
        if  length(rhythm) == 4 && all(rhythm == '(VFL')
            rhythmType = 'VFL';
        elseif length(rhythm) == 2 && all(rhythm == '(N')
            rhythmType = 'N';
        elseif length(rhythm) == 3 && all(rhythm == '(VT') % New class VT
            rhythmType = 'VT';
        elseif length(rhythm) == 5 && all(rhythm == '(AFIB') % New class VT
            rhythmType = 'AFIB';
        elseif length(rhythm) == 4 && all(rhythm == '(BII') % New class VT
            rhythmType = 'BII';
        else
            count = count + 1; % Skip unrecognized rhythms
            continue;
        end

    % Initialize fields if rhythmType doesn't exist
    if ~isfield(arrhythmiaData, rhythmType)
        arrhythmiaData.(rhythmType).R_peak_vals = [];
        arrhythmiaData.(rhythmType).R_peak_ind = [];
        arrhythmiaData.(rhythmType).Q_peak_vals = [];
        arrhythmiaData.(rhythmType).Q_peak_ind = [];
        arrhythmiaData.(rhythmType).T_peak_vals = [];
        arrhythmiaData.(rhythmType).S_peak_vals = [];
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
    arrhythmiaData.(rhythmType).R_peak_vals = [arrhythmiaData.(rhythmType).R_peak_vals, ...
        ecg(R_peaks_ind((R_peaks_ind > start_count) & (R_peaks_ind < end_count)))'];
    
    arrhythmiaData.(rhythmType).R_peak_ind = [arrhythmiaData.(rhythmType).R_peak_ind, ...
        R_peaks_ind((R_peaks_ind > start_count) & (R_peaks_ind < end_count))];
                
    arrhythmiaData.(rhythmType).Q_peak_vals = [arrhythmiaData.(rhythmType).Q_peak_vals, ...
        ecg(Q_peaks_ind((Q_peaks_ind > start_count) & (Q_peaks_ind < end_count)))'];

    arrhythmiaData.(rhythmType).Q_peak_ind = [arrhythmiaData.(rhythmType).Q_peak_ind, ...
        Q_peaks_ind((Q_peaks_ind > start_count) & (Q_peaks_ind < end_count))];
        
    arrhythmiaData.(rhythmType).T_peak_vals = [arrhythmiaData.(rhythmType).T_peak_vals, ...
        ecg(T_peaks_ind((T_peaks_ind > start_count) & (T_peaks_ind < end_count)))'];

    arrhythmiaData.(rhythmType).S_peak_vals = [arrhythmiaData.(rhythmType).S_peak_vals, ...
        ecg(S_peaks_ind((S_peaks_ind > start_count) & (S_peaks_ind < end_count)))'];

    arrhythmiaData.(rhythmType).S_peak_ind = [arrhythmiaData.(rhythmType).S_peak_ind, ...
        S_peaks_ind((S_peaks_ind > start_count) & (S_peaks_ind < end_count))];
end

% Calculate RR and QS intervals for each rhythm
rhythms = fieldnames(arrhythmiaData);
for i = 1:length(rhythms)
    rhythmType = rhythms{i};
    obs = min([length(arrhythmiaData.(rhythmType).R_peak_vals),length(arrhythmiaData.(rhythmType).Q_peak_vals),length(arrhythmiaData.(rhythmType).S_peak_vals),length(arrhythmiaData.(rhythmType).T_peak_vals)]);
        arrhythmiaData.(rhythmType).R_peak_vals = arrhythmiaData.(rhythmType).R_peak_vals(1:obs);
        arrhythmiaData.(rhythmType).R_peak_ind = arrhythmiaData.(rhythmType).R_peak_ind(1:obs);
        arrhythmiaData.(rhythmType).Q_peak_vals = arrhythmiaData.(rhythmType).Q_peak_vals(1:obs);
        arrhythmiaData.(rhythmType).S_peak_vals = arrhythmiaData.(rhythmType).S_peak_vals(1:obs);
        arrhythmiaData.(rhythmType).Q_peak_ind = arrhythmiaData.(rhythmType).Q_peak_ind(1:obs);
        arrhythmiaData.(rhythmType).S_peak_ind = arrhythmiaData.(rhythmType).S_peak_ind(1:obs);
        arrhythmiaData.(rhythmType).T_peak_vals = arrhythmiaData.(rhythmType).T_peak_vals(1:obs);
        arrhythmiaData.(rhythmType).RR_int = calc_rr(arrhythmiaData.(rhythmType).R_peak_ind, Fs);
        arrhythmiaData.(rhythmType).QS_int = calc_qs(arrhythmiaData.(rhythmType).Q_peak_ind, ...
                                                 arrhythmiaData.(rhythmType).S_peak_ind, Fs);
end

%% Prepare Features and Labels
rhythmTypes = fieldnames(arrhythmiaData);
featureNames = {'R_peak_vals', 'Q_peak_vals', 'S_peak_vals', 'T_peak_vals', 'RR_int', 'QS_int'};
X = [];
Y = [];

    labelMapping = containers.Map(rhythmTypes, 1:length(rhythmTypes));

    for i = 1:length(rhythmTypes)
        rhythmType = rhythmTypes{i};
        label = labelMapping(rhythmType);

        % Extract features
        rhythmData = arrhythmiaData.(rhythmType);
        numObservations = length(rhythmData.R_peak_ind);

        % Initialize temporary feature matrix
        tempX = zeros(numObservations, length(featureNames));

        for j = 1:length(featureNames)
            feature = rhythmData.(featureNames{j});
            if ~isempty(feature)
                tempX(:, j) = feature(:);
            end
        end

        % Append the features and labels
        X = [X; tempX]; % Append features
        Y = [Y; label * ones(numObservations, 1)]; % Append labels
    end

    %% Collect all features and labels from all files
    all_X = [all_X; X]; 
    all_Y = [all_Y; Y]; 
end


%% Split Data into Training and Testing Sets

X = all_X;
Y = all_Y;

% Get unique classes and find class 1 (assuming class 1 is the target class)
unique_classes = unique(Y);
target_class = 3; % We'll use class 1 as our target class
if ~ismember(target_class, unique_classes)
    % If the target class (AFIB) is not found in the *entire* dataset after
    % filtering databases, you might have removed the source of AFIB data.
    error('Target class %d (AFIB) not found in the combined data.', target_class);
end

% Create binary labels for AFIB vs all for the *entire* dataset
binary_labels = -ones(size(Y));
binary_labels(Y == target_class) = 1; % Assign 1 to the target class (AFIB)

%% Use Stratified Random Split to ensure class distribution is similar in train and test
fprintf('Performing stratified random split...\n');
% Create a partition object. 80% for training, stratified by the binary labels.
% Use the binary labels to ensure AFIB and Not AFIB are represented in both sets.
cvp = cvpartition(binary_labels, 'Holdout', 0.2, 'Stratify', true); % 20% for test

% Get the training and testing indices from the partition
trainIdx = training(cvp); % Logical indices for training
testIdx = test(cvp);     % Logical indices for testing

% Create the training and testing datasets using the logical indices
X_train = X(trainIdx, :);
y_train = Y(trainIdx); % Keep original multiclass labels for train if needed elsewhere
binary_labels_train = binary_labels(trainIdx); % Binary labels for training

X_test = X(testIdx, :);
y_test = Y(testIdx); % Keep original multiclass labels for test if needed elsewhere
binary_labels_test = binary_labels(testIdx); % Binary labels for testing

% Assign the binary labels to the variables used later
y_train_binary = binary_labels_train;
y_test_binary = binary_labels_test;

fprintf('Data split completed. Training samples: %d, Testing samples: %d\n', sum(trainIdx), sum(testIdx));
fprintf('Training set distribution: AFIB (%d): %d samples, Not AFIB (%d): %d samples\n', ...
        target_class, sum(y_train_binary == 1), -1, sum(y_train_binary == -1));
fprintf('Testing set distribution: AFIB (%d): %d samples, Not AFIB (%d): %d samples\n', ...
        target_class, sum(y_test_binary == 1), -1, sum(y_test_binary == -1));


%% Implement Undersampling with a Ratio
fprintf('Implementing undersampling with a ratio on the training data...\n');

% Find indices of minority class (AFIB, label 1) in the training data
minority_class_label = 1; % AFIB is labeled 1 in binary_labels
minority_indices = find(y_train_binary == minority_class_label);
num_minority_samples = length(minority_indices);

% Find indices of majority class (Not AFIB, label -1) in the training data
majority_class_label = -1; % Not AFIB is labeled -1 in binary_labels
majority_indices = find(y_train_binary == majority_class_label);
num_majority_samples = length(majority_indices);

fprintf('Training data before undersampling: Minority (%d): %d samples, Majority (%d): %d samples\n', ...
        minority_class_label, num_minority_samples, majority_class_label, num_majority_samples);

% Define the ratio of Majority to Minority samples in the undersampled training set
% E.g., 1 means 1:1 ratio (Majority count = Minority count)
% E.g., 2 means 1:2 ratio (Majority count = 2 * Minority count)
% Adjust this value based on desired training set size and performance
undersampling_ratio_majority_to_minority = 2; % <<< Start with a ratio (e.g., 2) and adjust

num_majority_to_select = round(num_minority_samples * undersampling_ratio_majority_to_minority);

% Ensure we don't select more majority samples than available
if num_majority_to_select > num_majority_samples
    warning('Cannot select %d majority samples for ratio %d, only %d available. Selecting all majority samples.', num_majority_to_select, undersampling_ratio_majority_to_minority, num_majority_samples);
    num_majority_to_select = num_majority_samples;
end

% Randomly select majority samples WITHOUT replacement
% Use randperm to get unique random indices
selected_majority_indices = majority_indices(randperm(num_majority_samples, num_majority_to_select));

% Combine the indices of the minority samples and the selected majority samples
balanced_train_indices = [minority_indices; selected_majority_indices];

% Shuffle the balanced indices to mix the classes
balanced_train_indices = balanced_train_indices(randperm(length(balanced_train_indices)));

% Create the undersampled training dataset
X_train_undersampled = X_train(balanced_train_indices, :);
y_train_binary_undersampled = y_train_binary(balanced_train_indices); % Use binary labels


fprintf('Training data after undersampling (ratio 1:%d): Minority (%d): %d samples, Majority (%d): %d samples\n', ...
        undersampling_ratio_majority_to_minority, ...
        minority_class_label, sum(y_train_binary_undersampled == minority_class_label), ...
        majority_class_label, sum(y_train_binary_undersampled == majority_class_label));


% Move undersampled training data and original test data to GPU
% Ensure data is single precision before moving to GPU
X_train_gpu = gpuArray(single(X_train_undersampled)); % Use undersampled training data
y_train_binary_gpu = gpuArray(single(y_train_binary_undersampled)); % Use undersampled binary labels
X_test_gpu = gpuArray(single(X_test)); % Use original test data

% Calculate class costs for cost matrix based on ORIGINAL TRAINING data IMBALANCE
% Note: target_class = 3 corresponds to binary label 1
% We will use the counts from the original y_train_binary before undersampling
num_afib_train_original = sum(y_train_binary == 1); % Count samples in original training set
num_not_afib_train_original = sum(y_train_binary == -1); % Count samples in original training set

% Calculate costs: Penalize misclassifying minority (AFIB, 1) as majority (-1)
% cost_fn is cost of predicting -1 when true is 1 (False Negative)
% Use counts from the ORIGINAL imbalanced data ratio for a strong penalty
cost_fn = num_not_afib_train_original / num_afib_train_original; % This will be a high value (~8.17 based on earlier numbers)
% cost_fp is cost of predicting 1 when true is -1 (False Positive)
cost_fp = 1; % A common starting point

% Define the Cost matrix for fitcsvm: Cost(i,j) is classifying into i when true is j
% ClassNames are [-1, 1]
% Cost matrix structure:
%           True -1   True 1
% Pred -1:    0       cost_fn
% Pred 1: cost_fp       0
Cost_Matrix = [0, cost_fn; cost_fp, 0];


% Train binary SVM for AFIB vs all
fprintf('Training binary SVM for AFIB vs all using undersampled data and cost-sensitive learning (original ratio cost)...\n');
afib_model = fitcsvm(X_train_gpu, y_train_binary_gpu, ... % Train on undersampled data
    'KernelFunction', 'linear', ...
    'ClassNames', [-1, 1], ...
    'BoxConstraint', 1, ...
    'Standardize', true, ...
    'Cost', Cost_Matrix); % Apply Cost matrix based on original imbalance
