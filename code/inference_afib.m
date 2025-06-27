%% Loop through each file and run the trained AFIB model, plot peaks, and report misclassifications
load('models/afib_model.mat');

folder = "databases/";
fileList = dir(fullfile(folder, '*.hea'));
fileList = {fileList.name};
afib_threshold = 0; % Set your custom threshold here

% Initialize arrays to store stats for all files
all_accuracy = [];
all_precision = [];
all_recall = [];
all_f1 = [];
all_auc = [];

for i = 1:length(fileList)
    recordname = str2mat(fullfile(folder, fileList{i}(1:end-4))); % Remove file extension
    display(['Processing file: ', recordname]);
    [ecg, Fs, tm] = rdsamp(recordname, 1);
    [ann, type, subtype, chan, num, comments] = rdann(recordname, 'atr', 1);
    % Rhythm Identification (propagate using + marker)
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
    % Feature extraction and labeling (match train_afib.m)
    [R_peaks_val, R_peaks_ind, Q_peaks_ind, Q_peaks_val, S_peaks_ind, S_peaks_val, T_peaks_ind, T_peaks_val, delay] = pan_tompkin(ecg, Fs, 0);
    % Build arrhythmiaData struct as in training
    arrhythmiaData = struct();
    count = 1;
    found_any = false;
    while count <= length(comments)
        rhythm = cell2mat(comments(count));
        % Assign rhythm type
        if  length(rhythm) == 4 && all(rhythm == '(VFL')
            rhythmType = 'VFL';
        elseif  length(rhythm) == 4 && all(rhythm == '(AFL')
            rhythmType = 'AFL';
        elseif length(rhythm) == 2 && all(rhythm == '(N')
            rhythmType = 'N';
        elseif length(rhythm) == 2 && all(rhythm == '(P')
            rhythmType = 'P';
        elseif length(rhythm) == 3 && all(rhythm == '(VT')
            rhythmType = 'VT';
        elseif length(rhythm) == 5 && all(rhythm == '(AFIB')
            rhythmType = 'AFIB';
        elseif length(rhythm) == 4 && all(rhythm == '(BII')
            rhythmType = 'BII';
        else
            count = count + 1;
            continue;
        end
        % Only AFIB and N are valid for this script
        if strcmp(rhythmType, 'AFIB') || strcmp(rhythmType, 'N')
            found_any = true;
        end
        if ~isfield(arrhythmiaData, rhythmType)
            arrhythmiaData.(rhythmType).R_peak_vals = [];
            arrhythmiaData.(rhythmType).R_peak_ind = [];
            arrhythmiaData.(rhythmType).Q_peak_vals = [];
            arrhythmiaData.(rhythmType).Q_peak_ind = [];
            arrhythmiaData.(rhythmType).T_peak_vals = [];
            arrhythmiaData.(rhythmType).S_peak_vals = [];
            arrhythmiaData.(rhythmType).S_peak_ind = [];
        end
        start_count = ann(count);
        while (count <= length(comments)) && ...
              (length(cell2mat(comments(count))) == length(rhythm)) && ...
              all(cell2mat(comments(count)) == rhythm)
            count = count + 1;
        end
        if count > length(ann)
            end_count = length(ecg);
        else
            end_count = ann(count);
        end
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
    % Skip file if no (AFIB) or (N) segments found
    if ~found_any
        fprintf('No (AFIB) or (N) segments found in file: %s\n', fileList{i});
        continue;
    end
    % Calculate RR and QS intervals for each rhythm
    rhythms = fieldnames(arrhythmiaData);
    for j = 1:length(rhythms)
        rhythmType = rhythms{j};
        obs = min([length(arrhythmiaData.(rhythmType).R_peak_vals),length(arrhythmiaData.(rhythmType).Q_peak_vals),length(arrhythmiaData.(rhythmType).S_peak_vals),length(arrhythmiaData.(rhythmType).T_peak_vals)]);
        arrhythmiaData.(rhythmType).R_peak_vals = arrhythmiaData.(rhythmType).R_peak_vals(1:obs);
        arrhythmiaData.(rhythmType).R_peak_ind = arrhythmiaData.(rhythmType).R_peak_ind(1:obs);
        arrhythmiaData.(rhythmType).Q_peak_vals = arrhythmiaData.(rhythmType).Q_peak_vals(1:obs);
        arrhythmiaData.(rhythmType).S_peak_vals = arrhythmiaData.(rhythmType).S_peak_vals(1:obs);
        arrhythmiaData.(rhythmType).Q_peak_ind = arrhythmiaData.(rhythmType).Q_peak_ind(1:obs);
        arrhythmiaData.(rhythmType).S_peak_ind = arrhythmiaData.(rhythmType).S_peak_ind(1:obs);
        arrhythmiaData.(rhythmType).T_peak_vals = arrhythmiaData.(rhythmType).T_peak_vals(1:obs);
        arrhythmiaData.(rhythmType).RR_int = calc_rr(arrhythmiaData.(rhythmType).R_peak_ind, Fs);
        arrhythmiaData.(rhythmType).QS_int = calc_qs(arrhythmiaData.(rhythmType).Q_peak_ind, arrhythmiaData.(rhythmType).S_peak_ind, Fs);
    end
    % Prepare Features and Labels (only AFIB and N)
    validRhythms = {'AFIB', 'N'};
    featureNames = {'R_peak_vals', 'Q_peak_vals', 'S_peak_vals', 'T_peak_vals', 'RR_int', 'QS_int'};
    X = [];
    true_labels = [];
    for j = 1:length(validRhythms)
        rhythmType = validRhythms{j};
        if isfield(arrhythmiaData, rhythmType)
            rhythmData = arrhythmiaData.(rhythmType);
            numObservations = length(rhythmData.R_peak_ind);
            tempX = zeros(numObservations, length(featureNames));
            for k = 1:length(featureNames)
                feature = rhythmData.(featureNames{k});
                if ~isempty(feature)
                    tempX(:, k) = feature(:);
                end
            end
            X = [X; tempX];
            if strcmp(rhythmType, 'AFIB')
                true_labels = [true_labels; ones(numObservations, 1)];
            else
                true_labels = [true_labels; zeros(numObservations, 1)];
            end
        end
    end
    % Predict using trained model
    [~, scores] = predict(afib_svm_model, X);
    y_pred = scores(:,2) > afib_threshold;
    % Find misclassifications
    mis_idx = find(y_pred ~= true_labels);
    % Compute metrics
    TP = sum(y_pred == 1 & true_labels == 1);
    TN = sum(y_pred == 0 & true_labels == 0);
    FP = sum(y_pred == 1 & true_labels == 0);
    FN = sum(y_pred == 0 & true_labels == 1);
    accuracy = (TP + TN) / (TP + TN + FP + FN);
    precision = TP / (TP + FP + eps);
    recall = TP / (TP + FN + eps);
    f1_score = 2 * (precision * recall) / (precision + recall + eps);
    % Compute AUC if both classes are present
    if any(true_labels==1) && any(true_labels==0)
        [~,~,~,AUC] = perfcurve(true_labels, scores(:,2), 1);
    else
        AUC = NaN;
    end
    % Plot ECG and peaks, highlight misclassified and correctly classified beats
    figure;
    plot(ecg, 'k'); hold on;
    scatter(R_peaks_ind, ecg(R_peaks_ind), 'rv', 'filled', 'DisplayName', 'R peaks');
    scatter(Q_peaks_ind, ecg(Q_peaks_ind), 'bo', 'filled', 'DisplayName', 'Q peaks');
    scatter(S_peaks_ind, ecg(S_peaks_ind), 'gs', 'filled', 'DisplayName', 'S peaks');
    scatter(T_peaks_ind, ecg(T_peaks_ind), 'md', 'filled', 'DisplayName', 'T peaks');
    % Robustly map indices for plotting
    n_N = 0; n_AFIB = 0;
    if isfield(arrhythmiaData, 'N')
        n_N = length(arrhythmiaData.N.R_peak_ind);
    end
    if isfield(arrhythmiaData, 'AFIB')
        n_AFIB = length(arrhythmiaData.AFIB.R_peak_ind);
    end
    % Plot correctly classified peaks
    correct_idx = find(y_pred == true_labels);
    for m = 1:length(correct_idx)
        idx_m = correct_idx(m);
        if idx_m <= n_N && n_N > 0
            r_idx = arrhythmiaData.N.R_peak_ind(idx_m);
        elseif idx_m > n_N && n_AFIB > 0 && idx_m <= (n_N + n_AFIB)
            r_idx = arrhythmiaData.AFIB.R_peak_ind(idx_m - n_N);
        else
            r_idx = NaN;
        end
        if ~isnan(r_idx)
            y = ecg(r_idx);
            if true_labels(idx_m) == 1
                scatter(r_idx, y, 60, 'g', 'o', 'filled', 'DisplayName', 'Correct AFIB');
            else
                scatter(r_idx, y, 60, 'b', 'o', 'filled', 'DisplayName', 'Correct Normal');
            end
        end
    end
    % Plot misclassified peaks
    for m = 1:length(mis_idx)
        idx_m = mis_idx(m);
        if idx_m <= n_N && n_N > 0
            r_idx = arrhythmiaData.N.R_peak_ind(idx_m);
        elseif idx_m > n_N && n_AFIB > 0 && idx_m <= (n_N + n_AFIB)
            r_idx = arrhythmiaData.AFIB.R_peak_ind(idx_m - n_N);
        else
            r_idx = NaN;
        end
        if ~isnan(r_idx)
            y = ecg(r_idx);
            if true_labels(idx_m) == 1
                true_color = 'r';
            else
                true_color = 'b';
            end
            if y_pred(idx_m) == 1
                pred_marker = 'x';
            else
                pred_marker = 'o';
            end
            scatter(r_idx, y, 100, true_color, pred_marker, 'LineWidth', 2, 'DisplayName', sprintf('Misclass: True %d, Pred %d', true_labels(idx_m), y_pred(idx_m)));
            text(r_idx, y, sprintf('T:%d/P:%d', true_labels(idx_m), y_pred(idx_m)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8, 'Color', true_color);
        end
    end
    title(sprintf('ECG and Peaks for %s (Misclassifications: %d)', fileList{i}, length(mis_idx)));
    xlabel('Sample Index'); ylabel('ECG Amplitude');
    legend('show'); grid on; hold off;
    % Print misclassification summary and metrics
    fprintf('File: %s | Total beats: %d | Misclassified: %d\n', fileList{i}, length(true_labels), length(mis_idx));
    fprintf('Accuracy: %.4f | Precision: %.4f | Recall: %.4f | F1: %.4f | AUC: %.4f\n', accuracy, precision, recall, f1_score, AUC);
    num_true_afib = sum(true_labels == 1);
    fprintf('Number of true AFIB beats: %d\n', num_true_afib);
    % Store stats
    all_accuracy = [all_accuracy; accuracy];
    all_precision = [all_precision; precision];
    all_recall = [all_recall; recall];
    all_f1 = [all_f1; f1_score];
    all_auc = [all_auc; AUC];
end

% Print mean stats at the end
fprintf('\n==== Mean Inference Statistics Across All Files ====');
fprintf('\nMean Accuracy: %.4f', mean(all_accuracy, 'omitnan'));
fprintf('\nMean Precision: %.4f', mean(all_precision, 'omitnan'));
fprintf('\nMean Recall: %.4f', mean(all_recall, 'omitnan'));
fprintf('\nMean F1 Score: %.4f', mean(all_f1, 'omitnan'));
fprintf('\nMean AUC: %.4f', mean(all_auc, 'omitnan'));

% Optionally, print std as well
fprintf('\nStd Accuracy: %.4f', std(all_accuracy, 'omitnan'));
fprintf('\nStd Precision: %.4f', std(all_precision, 'omitnan'));
fprintf('\nStd Recall: %.4f', std(all_recall, 'omitnan'));
fprintf('\nStd F1 Score: %.4f', std(all_f1, 'omitnan'));
fprintf('\nStd AUC: %.4f\n', std(all_auc, 'omitnan'));
