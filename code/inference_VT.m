%% Loop through each file and run the trained VT model, plot peaks, and report misclassifications
load('vt_model.mat');

folder = "databases/";
fileList = dir(fullfile(folder, '*.hea'));
fileList = {fileList.name};

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
    [R_peaks_val, R_peaks_ind, Q_peaks_ind, Q_peaks_val, S_peaks_ind, S_peaks_val, T_peaks_ind, T_peaks_val, ~] = pan_tompkin(ecg, Fs, 0);
    % Calculate RR and QS intervals
    RR_int = [];
    QS_int = [];
    temp_i = R_peaks_ind(1);
    for idx = R_peaks_ind(2:end)
        RR_int = [RR_int, (idx - temp_i) / Fs];
        temp_i = idx;
    end
    for idx = 1:length(S_peaks_ind)
        QS_int = [QS_int, (S_peaks_ind(idx) - Q_peaks_ind(idx)) / Fs];
    end
    % Rhythm Identification and Propagation
    [~, filename, ~] = fileparts(recordname);
    if startsWith(filename, 'cu')
        vt_start = [];
        vt_end = [];
        for idx = 1:length(type)
            if type(idx) == '['
                vt_start = [vt_start, ann(idx)];
            elseif type(idx) == ']'
                vt_end = [vt_end, ann(idx)];
            end
        end
        if length(vt_start) > length(vt_end)
            remaining_vt = length(vt_start) - length(vt_end);
            vt_end = [vt_end, repmat(length(ecg), 1, remaining_vt)];
        end
        for idx = 1:length(ann)
            comments{idx} = '(N';
        end
        for idx = 1:length(vt_start)
            start_idx = find(ann >= vt_start(idx), 1, 'first');
            end_idx = find(ann <= vt_end(idx), 1, 'last');
            if ~isempty(start_idx) && ~isempty(end_idx)
                for j = start_idx:end_idx
                    comments{j} = '(VT';
                end
            end
        end
    else
        rhythm = comments(1);
        count = 1;
        while count < length(ann)
            if (type(count) == '+')
                rhythm = comments(count);
            end
            comments(count) = rhythm;
            count = count + 1;
        end
    end
    % Feature Extraction and Labeling
    arrhythmiaData = struct();
    count = 1;
    found_any = false;
    while count <= length(comments)
        rhythm = cell2mat(comments(count));
        if  length(rhythm) == 4 && all(rhythm == '(VFL')
            rhythmType = 'VFL';
        elseif length(rhythm) == 2 && all(rhythm == '(N')
            rhythmType = 'N';
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
        % Only VT and N are valid for this script
        if strcmp(rhythmType, 'VT') || strcmp(rhythmType, 'N')
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
        if count <= length(ann)
            end_count = ann(count);
        else
            end_count = length(ecg);
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
    % Skip file if no (VT) or (N) segments found
    if ~found_any
        fprintf('No (VT) or (N) segments found in file: %s\n', fileList{i});
        continue;
    end
    rhythms = fieldnames(arrhythmiaData);
    for idx = 1:length(rhythms)
        rhythmType = rhythms{idx};
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
    featureNames = {'R_peak_vals', 'Q_peak_vals', 'S_peak_vals', 'T_peak_vals', 'RR_int', 'QS_int'};
    X = [];
    true_labels = [];
    validRhythms = {'VT', 'N'};
    for idx = 1:length(validRhythms)
        rhythmType = validRhythms{idx};
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
            if strcmp(rhythmType, 'VT')
                true_labels = [true_labels; ones(numObservations, 1)];
            else
                true_labels = [true_labels; zeros(numObservations, 1)];
            end
        end
    end
    obs = size(X,1);
    % Predict using trained model
    [y_pred, scores] = predict(vt_svm_model, X);
    % Find misclassifications
    mis_idx = find(y_pred ~= true_labels);
    % Find correctly classified indices
    correct_idx = find(y_pred == true_labels);
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
    % Plot ECG and peaks, highlight correctly and misclassified beats
    figure;
    plot(ecg, 'k'); hold on;
    scatter(R_peaks_ind, ecg(R_peaks_ind), 'rv', 'filled', 'DisplayName', 'R peaks');
    scatter(Q_peaks_ind, ecg(Q_peaks_ind), 'bo', 'filled', 'DisplayName', 'Q peaks');
    scatter(S_peaks_ind, ecg(S_peaks_ind), 'gs', 'filled', 'DisplayName', 'S peaks');
    scatter(T_peaks_ind, ecg(T_peaks_ind), 'md', 'filled', 'DisplayName', 'T peaks');
    % Plot correctly classified peaks (only first instance gets DisplayName for legend clarity)
    plotted_vt = false; plotted_n = false;
    for m = 1:length(correct_idx)
        idx_m = correct_idx(m);
        x = R_peaks_ind(idx_m);
        y = ecg(x);
        if true_labels(idx_m) == 1
            if ~plotted_vt
                scatter(x, y, 60, 'g', 'o', 'filled', 'DisplayName', 'Correct VT'); plotted_vt = true;
            else
                scatter(x, y, 60, 'g', 'o', 'filled', 'DisplayName', '');
            end
        else
            if ~plotted_n
                scatter(x, y, 60, 'b', 'o', 'filled', 'DisplayName', 'Correct Normal'); plotted_n = true;
            else
                scatter(x, y, 60, 'b', 'o', 'filled', 'DisplayName', '');
            end
        end
    end
    % Plot misclassified peaks (only first instance gets DisplayName for legend clarity)
    plotted_mis = false;
    for m = 1:length(mis_idx)
        idx_m = mis_idx(m);
        x = R_peaks_ind(idx_m);
        y = ecg(x);
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
        if ~plotted_mis
            scatter(x, y, 100, true_color, pred_marker, 'LineWidth', 2, 'DisplayName', 'Misclassified'); plotted_mis = true;
        else
            scatter(x, y, 100, true_color, pred_marker, 'LineWidth', 2, 'DisplayName', '');
        end
        text(x, y, sprintf('T:%d/P:%d', true_labels(idx_m), y_pred(idx_m)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8, 'Color', true_color);
    end
    title(sprintf('ECG and Peaks for %s (Misclassifications: %d)', fileList{i}, length(mis_idx)));
    xlabel('Sample Index'); ylabel('ECG Amplitude');
    legend('show'); grid on; hold off;
    % Print misclassification summary and metrics
    fprintf('File: %s | Total beats: %d | Misclassified: %d\n', fileList{i}, obs, length(mis_idx));
    fprintf('Accuracy: %.4f | Precision: %.4f | Recall: %.4f | F1: %.4f | AUC: %.4f\n', accuracy, precision, recall, f1_score, AUC);
    num_true_vt = sum(true_labels == 1);
    fprintf('Number of true VT beats: %d\n', num_true_vt);
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

fprintf('\nStd Accuracy: %.4f', std(all_accuracy, 'omitnan'));
fprintf('\nStd Precision: %.4f', std(all_precision, 'omitnan'));
fprintf('\nStd Recall: %.4f', std(all_recall, 'omitnan'));
fprintf('\nStd F1 Score: %.4f', std(all_f1, 'omitnan'));
fprintf('\nStd AUC: %.4f\n', std(all_auc, 'omitnan'));