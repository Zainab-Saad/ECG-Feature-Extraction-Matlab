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
    my_classes = {'N', 'VFL', 'VT'};
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

    % Print available rhythm types
    fprintf('\nAvailable rhythm types in data:\n');
    for i = 1:length(rhythmTypes)
        fprintf('%s\n', rhythmTypes{i});
    end

    % Only keep AFIB and NORMAL classes
    validRhythms = {'AFIB', 'N'};
    X = [];
    Y = [];

    fprintf('\nCollecting features for each rhythm type:\n');
    for i = 1:length(validRhythms)
        rhythmType = validRhythms{i};
        if isfield(arrhythmiaData, rhythmType)
            % Extract features
            rhythmData = arrhythmiaData.(rhythmType);
            numObservations = length(rhythmData.R_peak_ind);
            
            fprintf('Found %d samples for %s\n', numObservations, rhythmType);

            % Initialize temporary feature matrix
            tempX = zeros(numObservations, length(featureNames));

            for j = 1:length(featureNames)
                feature = rhythmData.(featureNames{j});
                if ~isempty(feature)
                    tempX(:, j) = feature(:);
                else
                    fprintf('Warning: Empty feature %s for %s\n', featureNames{j}, rhythmType);
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
        else
            fprintf('Warning: %s not found in data\n', rhythmType);
        end
    end

    % Print collected data statistics
    fprintf('\nCollected data statistics:\n');
    fprintf('Total samples collected: %d\n', length(Y));
    fprintf('AFIB samples: %d\n', sum(Y==1));
    fprintf('NORMAL samples: %d\n', sum(Y==0));

    % Validate data before proceeding
    if isempty(X) || isempty(Y)
        error('No data was collected. Please check the input data and rhythm types.');
    end

    if sum(Y==1) == 0 || sum(Y==0) == 0
        error('Missing one or both classes. AFIB samples: %d, NORMAL samples: %d', sum(Y==1), sum(Y==0));
    end

    %% Prepare test data
    X_test = X;  % Use collected data directly
    y_test = Y;

    % Print test data statistics
    fprintf('\nTest Data Statistics:\n');
    fprintf('Total test samples: %d\n', length(y_test));
    fprintf('AFIB samples: %d (%.2f%%)\n', sum(y_test==1), 100*sum(y_test==1)/length(y_test));
    fprintf('NORMAL samples: %d (%.2f%%)\n', sum(y_test==0), 100*sum(y_test==0)/length(y_test));

    % Validate test data
    if length(y_test) == 0
        error('Test data is empty. Please check data collection process.');
    end

    if sum(y_test==1) == 0 || sum(y_test==0) == 0
        error('Test data is missing one or both classes. AFIB: %d, NORMAL: %d', sum(y_test==1), sum(y_test==0));
    end

    %% Testing
    % Load the trained model and test data
    % Note: Make sure to run train_multiclass.m first to get the model and data

    % Process test data in batches
    batch_size = 1000;  % Same batch size as training
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

    % Print results
    fprintf('\nTest Set Performance:\n');
    fprintf('Accuracy: %.2f%%\n', accuracy * 100);
    fprintf('Precision: %.4f\n', precision);
    fprintf('Recall: %.4f\n', recall);
    fprintf('F1 Score: %.4f\n', f1_score);
    fprintf('AUC: %.4f\n', AUC);

    % Plot confusion matrix
    figure;
    confusionchart(y_test, y_pred);
    title('AFIB vs NORMAL Classification');

    % Plot ROC curve
    figure;
    [X_roc,Y_roc,~,~] = perfcurve(y_test, scores(:,2), 1);
    plot(X_roc,Y_roc);
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    title(sprintf('ROC Curve (AUC = %.3f)', AUC));
    grid on;

    % Evaluate different decision thresholds
    thresholds = [-2, -1, -0.5, 0, 0.5, 1];
    fprintf('\nEvaluating different decision thresholds:\n');
    fprintf('Threshold\tAccuracy\tPrecision\tRecall\tF1 Score\n');
    fprintf('--------------------------------------------------------\n');

    for threshold = thresholds
        % Apply threshold
        y_pred_threshold = scores(:,2) > threshold;
        
        % Calculate metrics
        TP = sum(y_pred_threshold == 1 & y_test == 1);
        TN = sum(y_pred_threshold == 0 & y_test == 0);
        FP = sum(y_pred_threshold == 1 & y_test == 0);
        FN = sum(y_pred_threshold == 0 & y_test == 1);
        
        acc = (TP + TN) / (TP + TN + FP + FN);
        prec = TP / (TP + FP + eps);
        rec = TP / (TP + FN + eps);
        f1 = 2 * (prec * rec) / (prec + rec + eps);
        
        fprintf('%.2f\t\t%.2f%%\t\t%.4f\t\t%.4f\t%.4f\n', ...
            threshold, acc*100, prec, rec, f1);
    end
end

