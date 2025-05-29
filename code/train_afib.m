folder = "databases/";
fileList = dir(fullfile(folder, '*.hea'));
fileList = {fileList.name};

% Initialize variables to collect all data
all_X = [];
all_Y = [];
arrhythmiaData = struct();
total_len=0;
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

% Only keep AFIB and NORMAL classes
validRhythms = {'AFIB', 'N'};
X = [];
Y = [];

for i = 1:length(validRhythms)
    rhythmType = validRhythms{i};
    if isfield(arrhythmiaData, rhythmType)
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
        % Use 1 for AFIB, 0 for NORMAL
        if strcmp(rhythmType, 'AFIB')
            Y = [Y; ones(numObservations, 1)];
        else
            Y = [Y; zeros(numObservations, 1)];
        end
    end
end

all_X = [all_X; X]; 
all_Y = [all_Y; Y]; 
end
%%
X = all_X;
Y = all_Y;

% Print dataset statistics
fprintf('\nDataset Statistics:\n');
fprintf('Total samples: %d\n', length(Y));
fprintf('AFIB samples: %d (%.2f%%)\n', sum(Y==1), 100*sum(Y==1)/length(Y));
fprintf('NORMAL samples: %d (%.2f%%)\n', sum(Y==0), 100*sum(Y==0)/length(Y));

% Check for feature correlations
fprintf('\nFeature Correlations:\n');
correlation_matrix = corrcoef(X);
feature_names = {'R_peak_vals', 'Q_peak_vals', 'S_peak_vals', 'T_peak_vals', 'RR_int', 'QS_int'};
for i = 1:length(feature_names)
    for j = i+1:length(feature_names)
        fprintf('Correlation between %s and %s: %.4f\n', ...
            feature_names{i}, feature_names{j}, correlation_matrix(i,j));
    end
end

% Initialize performance metrics storage
num_folds = 5;
fold_accuracies = zeros(num_folds, 1);
fold_precisions = zeros(num_folds, 1);
fold_recalls = zeros(num_folds, 1);
fold_f1_scores = zeros(num_folds, 1);
fold_aucs = zeros(num_folds, 1);

% Create k-fold cross-validation partition with stratification
cv = cvpartition(Y, 'KFold', num_folds, 'Stratify', true);

fprintf('\nStarting %d-fold cross-validation...\n', num_folds);

% Perform k-fold cross-validation
for fold = 1:num_folds
    fprintf('\nFold %d/%d:\n', fold, num_folds);
    
    % Get training and test indices for this fold
    train_idx = training(cv, fold);
    test_idx = test(cv, fold);
    
    % Split data
    X_train = X(train_idx, :);
    y_train = Y(train_idx);
    X_test = X(test_idx, :);
    y_test = Y(test_idx);
    
    % Print fold statistics
    fprintf('Fold %d - Training set: %d samples (%d AFIB, %d NORMAL)\n', ...
        fold, length(y_train), sum(y_train==1), sum(y_train==0));
    fprintf('Fold %d - Test set: %d samples (%d AFIB, %d NORMAL)\n', ...
        fold, length(y_test), sum(y_test==1), sum(y_test==0));
    
    % Move training data to GPU
    X_train_gpu = gpuArray(single(X_train));
    y_train_gpu = gpuArray(single(y_train));
    
    % Calculate class weights to handle imbalance
    num_afib = sum(y_train == 1);
    num_normal = sum(y_train == 0);
    class_weights = zeros(size(y_train));
    class_weights(y_train == 1) = (num_normal/num_afib) * 1.5;  % Increased weight for AFIB class
    class_weights(y_train == 0) = 1.2;                         % Slightly increased weight for Normal class
    class_weights_gpu = gpuArray(single(class_weights));
    
    % Train SVM on GPU
    svm_model = fitcsvm(X_train_gpu, y_train_gpu, ...
        'KernelFunction', 'linear', ...
        'ClassNames', [0, 1], ...
        'BoxConstraint', 5, ...  % Increased from 1 to 5 for harder margin
        'Standardize', true, ...
        'Weights', class_weights_gpu);
    
    % Process test data in batches
    batch_size = 1000;  % Adjust this based on your GPU memory
    num_test_samples = size(X_test, 1);
    num_batches = ceil(num_test_samples / batch_size);
    
    % Initialize arrays for predictions and scores
    y_pred = zeros(num_test_samples, 1);
    scores = zeros(num_test_samples, 2);
    
    fprintf('Processing test data in %d batches...\n', num_batches);
    
    % Process each batch
    for batch = 1:num_batches
        % Calculate batch indices
        start_idx = (batch-1) * batch_size + 1;
        end_idx = min(batch * batch_size, num_test_samples);
        
        % Get current batch
        X_test_batch = X_test(start_idx:end_idx, :);
        
        % Move batch to GPU
        X_test_gpu = gpuArray(single(X_test_batch));
        
        % Get predictions for batch
        [y_pred_batch_gpu, scores_batch_gpu] = predict(svm_model, X_test_gpu);
        
        % Move predictions back to CPU
        y_pred(start_idx:end_idx) = gather(y_pred_batch_gpu);
        scores(start_idx:end_idx, :) = gather(scores_batch_gpu);
        
        % Clear GPU memory
        clear X_test_gpu y_pred_batch_gpu scores_batch_gpu;
        
        % Print progress
        fprintf('Processed batch %d/%d\n', batch, num_batches);
    end
    
    % Calculate metrics
    TP = sum(y_pred == 1 & y_test == 1);
    TN = sum(y_pred == 0 & y_test == 0);
    FP = sum(y_pred == 1 & y_test == 0);
    FN = sum(y_pred == 0 & y_test == 1);
    
    accuracy = (TP + TN) / (TP + TN + FP + FN);
    precision = TP / (TP + FP + eps);
    recall = TP / (TP + FN + eps);
    f1_score = 2 * (precision * recall) / (precision + recall + eps);
    
    % Calculate AUC
    [~,~,~,AUC] = perfcurve(y_test, scores(:,2), 1);
    
    % Store metrics
    fold_accuracies(fold) = accuracy;
    fold_precisions(fold) = precision;
    fold_recalls(fold) = recall;
    fold_f1_scores(fold) = f1_score;
    fold_aucs(fold) = AUC;
    
    % Print fold results
    fprintf('Fold %d Results:\n', fold);
    fprintf('Accuracy: %.2f%%\n', accuracy * 100);
    fprintf('Precision: %.4f\n', precision);
    fprintf('Recall: %.4f\n', recall);
    fprintf('F1 Score: %.4f\n', f1_score);
    fprintf('AUC: %.4f\n', AUC);
    
    % Plot confusion matrix for this fold
    figure;
    confusionchart(y_test, y_pred);
    title(sprintf('Confusion Matrix - Fold %d', fold));
    
    % Plot ROC curve for this fold
    figure;
    [X_roc,Y_roc,~,~] = perfcurve(y_test, scores(:,2), 1);
    plot(X_roc,Y_roc);
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    title(sprintf('ROC Curve - Fold %d (AUC = %.3f)', fold, AUC));
    grid on;
    
    % Clear GPU memory after each fold
    clear X_train_gpu y_train_gpu class_weights_gpu;
end

% Print overall results
fprintf('\nOverall Cross-Validation Results:\n');
fprintf('Mean Accuracy: %.2f%% ± %.2f%%\n', mean(fold_accuracies)*100, std(fold_accuracies)*100);
fprintf('Mean Precision: %.4f ± %.4f\n', mean(fold_precisions), std(fold_precisions));
fprintf('Mean Recall: %.4f ± %.4f\n', mean(fold_recalls), std(fold_recalls));
fprintf('Mean F1 Score: %.4f ± %.4f\n', mean(fold_f1_scores), std(fold_f1_scores));
fprintf('Mean AUC: %.4f ± %.4f\n', mean(fold_aucs), std(fold_aucs));

% Plot boxplot of performance metrics
figure;
boxplot([fold_accuracies, fold_precisions, fold_recalls, fold_f1_scores, fold_aucs], ...
    'Labels', {'Accuracy', 'Precision', 'Recall', 'F1 Score', 'AUC'});
title('Distribution of Performance Metrics Across Folds');
ylabel('Score');
grid on;

% Save the final trained model and essential variables
fprintf('\nSaving model and essential variables...\n');
save('afib_model.mat', 'svm_model', 'featureNames', 'Fs', 'class_weights');
fprintf('Model saved successfully as afib_model.mat\n');