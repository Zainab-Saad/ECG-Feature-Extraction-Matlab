
%% Testing
% Get predictions and scores

%X_test_gpu_2 = gpuArray(single(X_test));        % <-- THIS LINE IS CRUCIAL FOR PREDICTION
%[~, afib_scores] = predict(afib_model, X_test_gpu_2);
num_test_samples = size(X_test_gpu, 1);
afib_scores_gpu = gpuArray(zeros(num_test_samples, 2, 'single')); % Preallocate scores on GPU (use single if X_test_gpu is single)

% Define batch size (e.g., process 1000 samples at a time, adjust based on trial and error)
batch_size = 1000; % You might need to adjust this number

% Process test data in batches
for start_idx = 1:batch_size:num_test_samples
    end_idx = min(start_idx + batch_size - 1, num_test_samples);
    
    % Extract the current batch from the GPU
    X_test_batch_gpu = X_test_gpu(start_idx:end_idx, :);
    
    % Predict on the current batch
    [~, batch_scores] = predict(afib_model, X_test_batch_gpu);
    
    % Store the scores for the current batch
    afib_scores_gpu(start_idx:end_idx, :) = batch_scores;
    
    fprintf('Processed batch %d to %d\n', start_idx, end_idx);
end


%%

% Move all gathered scores back to CPU at the end
afib_scores = gather(afib_scores_gpu);

% Define a range of decision thresholds to evaluate
% SVM scores can range, often from negative to positive. Start experimenting
% with thresholds around 0 and maybe negative values.
thresholds_to_evaluate = [-2, -1, -0.5, 0, 0.5, 1]; % <<< Adjust these thresholds

fprintf('\nEvaluating performance for different decision thresholds:\n');

for current_threshold = thresholds_to_evaluate
    fprintf('--- Threshold: %.2f ---\n', current_threshold);

    % Classify based on the current threshold
    afib_predictions = afib_scores(:, 2) > current_threshold; % Apply the new threshold

    % Convert boolean predictions to double (0/1) for calculations and confusionchart
    afib_predictions_double = double(afib_predictions);
    % y_test_binary_chart is already 0/1 from previous calculation

    % Calculate Confusion Matrix components (TP, TN, FP, FN) for class 1 (AFIB)
    % True Positives (TP): Actual 1, Predicted 1
    TP = sum(afib_predictions_double == 1 & y_test_binary_chart == 1);
    % True Negatives (TN): Actual 0, Predicted 0
    TN = sum(afib_predictions_double == 0 & y_test_binary_chart == 0);
    % False Positives (FP): Actual 0, Predicted 1
    FP = sum(afib_predictions_double == 1 & y_test_binary_chart == 0);
    % False Negatives (FN): Actual 1, Predicted 0
    FN = sum(afib_predictions_double == 0 & y_test_binary_chart == 1);

    % Calculate Precision for class 1 (AFIB)
    % Handle the case where TP + FP is zero to avoid division by zero
    Precision_AFIB = 0; % Default to 0 if no positive predictions are made
    if (TP + FP) > 0
        Precision_AFIB = TP / (TP + FP);
    end

    % Calculate Recall (Sensitivity) for class 1 (AFIB)
    % Handle the case where TP + FN is zero (no actual positive samples)
    Recall_AFIB = 0; % Default to 0 if no actual positive samples exist
    if (TP + FN) > 0
        Recall_AFIB = TP / (TP + FN);
    end

     % Calculate Overall Accuracy for this threshold
    Total_Samples = length(y_test_binary_chart);
    Overall_Accuracy = (TP + TN) / Total_Samples;


    % Display Metrics for the current threshold
    fprintf('Overall Accuracy: %.2f%%\n', Overall_Accuracy * 100);
    fprintf('Precision (AFIB): %.4f\n', Precision_AFIB);
    fprintf('Recall (AFIB): %.4f\n', Recall_AFIB);
    fprintf('----------------------\n');

    % Optional: Display confusion matrix for each threshold (can generate many figures)
    % figure;
    % confusionchart(y_test_binary_chart, afib_predictions_double);
    % title(['Confusion Matrix (Threshold: ', num2str(current_threshold, '%.2f'), ')']);

end % End of threshold loop

% You can keep the display of the confusion matrix for the *last* threshold evaluated
figure;
confusionchart(y_test_binary_chart, afib_predictions_double); % Note: this uses the results from the LAST threshold in the loop
title('AFIB vs All Classification Confusion Matrix (Last Threshold Evaluated)');


% The final reported accuracy/precision/recall will be from the last threshold
% If you want the values for a specific threshold printed at the end,
% you might need to store them in the loop and print after the loop.
% For now, they are printed inside the loop for each threshold.
