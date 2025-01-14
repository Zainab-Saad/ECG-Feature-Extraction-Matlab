%%

display('Reading samples ECG signal from MIT-BIH Arrhythmia Database')
[ecg,Fs,tm]=rdsamp('mitdbase/207',1);

display('Reading and plotting annotations (human labels) of QRS complexes performend on the signals')
[ann,type,subtype,chan,num]=rdann('mitdbase/207','atr',1);

%%

[R_peaks_val, R_peaks_ind,Q_peaks_ind, Q_peaks_val,S_peaks_ind, S_peaks_val,T_peaks_ind, T_peaks_val, delay] = pan_tompkin(ecg, Fs);


%% Calculating RR Interval
count = 1;
RR_int = [];
for i = R_peaks_ind
    count = count +1;
    if count >2
        RR_int = [RR_int , (i-temp_i)/Fs];
    end
    temp_i = i;

end

%% Calculating QS interval
i =1 ;
QS_int = [];
while i < length(S_peaks_ind)
    QS_int = [QS_int (S_peaks_ind(i)-Q_peaks_ind(i))/Fs];
    i = i+1;
end
%%
figure
plot(tm,ecg);hold on;grid on
plot(tm(ann+1),ecg(ann+1),'ro');

%%
display('Ploting 3D version of signal and labels')
[RR,tms]=ann2rr('mitdbase/208','atr');

%%
my_classes = ['N' 'V'];
count = 1;
V_start = [];
V_stop = [];
count =1;
while count<length(type)
    if type(count) == '['
        V_start = [V_start count];
    end
    if type(count) == ']'
        V_stop = [V_stop count];
    end
    count = count +1;
end
%%
classes=unique(type);
ms=length(classes);
SVMModels=cell(ms,1);
for j = 1:numel(classes)
    indx=strcmp(Y,classes(j)); % Create binary classes for each classifier
    SVMModels{j}=fitcsvm(X,indx,'ClassNames',[false true],'Standardize',true,...
        'KernelFunction','polynomial');
end

%%
folder = "mitdbase/";
fileList = dir(fullfile(folder, '*.hea'));
fileList = {fileList.name};
for i = fileList
    recordname = cell2mat(folder + cell2mat(i));
    recordname = recordname(1:end-4);
    display('Reading samples ECG signal from MIT-BIH Arrhythmia Database Sample Number:')
    display(recordname)
    [ecg,Fs,tm]=rdsamp(recordname,1);
    [ann,type,subtype,chan,num]=rdann(recordname,'atr',1);

end
%%


folder = "mitdbase/";
fileList = dir(fullfile(folder, '*.hea'));
fileList = {fileList.name};

% Initialize empty arrays for features (X) and labels (Y)
X = [];
Y = [];

% Loop through each file in the directory and incrementally train the SVM
for i = 1:length(fileList)
    recordname = fullfile(folder, fileList{i});
    recordname = recordname(1:end-4); % Remove the file extension (e.g., '.hea')

    display('Reading samples ECG signal from MIT-BIH Arrhythmia Database Sample Number:');
    display(recordname);
    
    % Read ECG signal and annotations
    [ecg, Fs, tm] = rdsamp(recordname, 1);
    [ann, type, subtype, chan, num] = rdann(recordname, 'atr', 1);

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
    tempX = [];
    tempY = [];

    % Assign numerical labels to each rhythm type
    labelMapping = containers.Map(rhythmTypes, 1:length(rhythmTypes));

    for i = 1:length(rhythmTypes)
        rhythmType = rhythmTypes{i};
        label = labelMapping(rhythmType);

        % Extract features for the current rhythm
        rhythmData = arrhythmiaData.(rhythmType);
        numObservations = length(rhythmData.R_peak_ind); % Assuming R_peak_ind determines observation count

        % Initialize a temporary feature matrix
        currentX = zeros(numObservations, length(featureNames));

        for j = 1:length(featureNames)
            feature = rhythmData.(featureNames{j});
            if ~isempty(feature)
                currentX(:, j) = feature(:);
            end
        end

        % Append the features and labels
        tempX = [tempX; currentX]; % Append features
        tempY = [tempY; label * ones(numObservations, 1)]; % Append corresponding labels
    end

    % Combine current file's data with overall data for training
    X = [X; tempX]; % Append features from current file
    Y = [Y; tempY]; % Append labels from current file

    %% Train the SVM Model (using all the data up until now)
    display('Training SVM model with current data...');
    
        % Continue training the existing model with the new data
    model = fitcsvm(X, Y, 'KernelFunction', 'linear', 'Standardize', true);
    
    % Store the updated model after training
    models{i} = model;

    %% Evaluate Performance on the current file (for demonstration)
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

    % Testing the model
    [predicted_labels, score] = predict(model, X_test);

    test_accuracy = sum(predicted_labels == y_test) / length(y_test) * 100;
    fprintf('Test accuracy for %s: %.2f%%\n', fileList{i}, test_accuracy);
end
