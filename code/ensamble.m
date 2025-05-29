% Load both AFIB and VT models
load('models/afib_model.mat', 'afib_model');
load('models/vt_model.mat', 'vt_model');

% Load test data
% Assuming X_test is your feature matrix from inference.m
% If you need to load it from a file, uncomment and modify the following line:
% load('test_data.mat', 'X_test');

% Get predictions from both models
[afib_pred, afib_scores] = predict(afib_model, X_test);
[vt_pred, vt_scores] = predict(vt_model, X_test);

% Initialize final prediction array
final_prediction = zeros(size(X_test, 1), 1);

% Apply priority logic: VT > AFIB > N
for i = 1:length(final_prediction)
    if vt_pred(i) == 1  % If VT is predicted
        final_prediction(i) = 2;  % VT class
    elseif afib_pred(i) == 1  % If AFIB is predicted
        final_prediction(i) = 1;  % AFIB class
    else
        final_prediction(i) = 0;  % Normal class
    end
end

% Calculate confidence scores
confidence_scores = zeros(size(X_test, 1), 1);
for i = 1:length(confidence_scores)
    if vt_pred(i) == 1
        confidence_scores(i) = vt_scores(i, 2);  % VT confidence
    elseif afib_pred(i) == 1
        confidence_scores(i) = afib_scores(i, 2);  % AFIB confidence
    else
        confidence_scores(i) = max(afib_scores(i, 1), vt_scores(i, 1));  % Normal confidence
    end
end

% Combine predictions and confidence scores
ensemble_results = [final_prediction, confidence_scores];

% Display results summary
fprintf('\nEnsemble Model Results:\n');
fprintf('Total samples: %d\n', length(final_prediction));
fprintf('VT predictions: %d (%.2f%%)\n', sum(final_prediction == 2), 100*sum(final_prediction == 2)/length(final_prediction));
fprintf('AFIB predictions: %d (%.2f%%)\n', sum(final_prediction == 1), 100*sum(final_prediction == 1)/length(final_prediction));
fprintf('Normal predictions: %d (%.2f%%)\n', sum(final_prediction == 0), 100*sum(final_prediction == 0)/length(final_prediction));

% Save results if needed
% save('ensemble_results.mat', 'ensemble_results', 'final_prediction', 'confidence_scores');
