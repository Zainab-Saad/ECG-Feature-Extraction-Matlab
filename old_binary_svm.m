display('Reading samples ECG signal from MIT-BIH Arrhythmia Database')
[ecg,Fs,tm]=rdsamp('mitdbase/207',1);

display(['Reading and plotting annotations (human labels) of QRS complexes performend on the signals'])
[ann,type,subtype,chan,num,comments] =rdann('mitdbase/207','atr',1);
%% PRE PROCESSING
[R_peaks_val, R_peaks_ind,Q_peaks_ind, Q_peaks_val,S_peaks_ind, S_peaks_val,T_peaks_ind, T_peaks_val, delay] = pan_tompkin(ecg, Fs);

% RR_int
count = 1;
RR_int = [];
for i = R_peaks_ind
    count = count +1;
    if count >2
        RR_int = [RR_int , (i-temp_i)/Fs];
    end
    temp_i = i;

end
% Calculating QS interval
i =1 ;
QS_int = [];
while i < length(S_peaks_ind)
    QS_int = [QS_int (S_peaks_ind(i)-Q_peaks_ind(i))/Fs];
    i = i+1;
end


% Using Comments
comments2 = comments;
count =1;
my_classes = ['N' , 'VFL'];
rhythm = comments(count);
while count<length(ann)
    if (type(count) == '+')
        rhythm = comments(count);
    end
    comments(count) = rhythm;

    count = count +1;
end

%%
start_count = 1;
end_count = 1;
count = 1;
VFL_R_peak_ind = []; % Initialize to store R peak indexes
VFL_Q_peak_ind = []; % Initialize to store Q peak indexes
VFL_T_peak_ind = []; % Initialize to store T peak indexes
VFL_S_peak_ind = [];

N_R_peak_ind = []; % Initialize to store R peak indexes
N_Q_peak_ind = []; % Initialize to store Q peak indexes
N_T_peak_ind = []; % Initialize to store T peak indexes
N_S_peak_ind = [];
while count <= length(comments)
    rhythm = cell2mat(comments(count));
    
    % Check if rhythm matches '(VFL'
    if (length(rhythm) == 4) && all(rhythm == '(VFL')
        start_count = ann(count);
        
        % Find the end of the VFL rhythm
        while (count <= length(comments)) && ...
              (length(cell2mat(comments(count))) == 4) && ...
              all(cell2mat(comments(count)) == '(VFL')
            count = count + 1;
        end
        
        
        end_count = ann(count); % Mark the end of VFL section
        
        % Get R peak indices within the VFL section
        VFL_R_peak_ind = [VFL_R_peak_ind, R_peaks_ind((R_peaks_ind > start_count) & (R_peaks_ind < end_count))];
        VFL_Q_peak_ind = [VFL_Q_peak_ind, Q_peaks_ind((Q_peaks_ind > start_count) & (Q_peaks_ind < end_count))];
        VFL_T_peak_ind = [VFL_T_peak_ind, T_peaks_ind((T_peaks_ind > start_count) & (T_peaks_ind < end_count))];
        VFL_S_peak_ind = [VFL_S_peak_ind, S_peaks_ind((S_peaks_ind > start_count) & (S_peaks_ind < end_count))];
    
    else 
        if (length(rhythm) == 2) && all(rhythm == '(N')
        start_count = ann(count);
        
        % Find the end of the VFL rhythm
        while (count <= length(comments)) && ...
              (length(cell2mat(comments(count))) == 2) && ...
              all(cell2mat(comments(count)) == '(N')
            count = count + 1;
        end
        
        end_count = ann(count); % Mark the end of VFL section
        
        % Get R peak indices within the VFL section
        N_R_peak_ind = [N_R_peak_ind, R_peaks_ind((R_peaks_ind > start_count) & (R_peaks_ind < end_count))];
        N_Q_peak_ind = [N_Q_peak_ind, Q_peaks_ind((Q_peaks_ind > start_count) & (Q_peaks_ind < end_count))];
        N_T_peak_ind = [N_T_peak_ind, T_peaks_ind((T_peaks_ind > start_count) & (T_peaks_ind < end_count))];
        N_S_peak_ind = [N_S_peak_ind, S_peaks_ind((S_peaks_ind > start_count) & (S_peaks_ind < end_count))];
        else count = count + 1; % Move to the next comment
    end
    
    end
   

end
N_RR_int= calc_rr (N_R_peak_ind, Fs);
VFL_RR_int = calc_rr (VFL_R_peak_ind,Fs);
N_QS_int= calc_qs (N_Q_peak_ind,N_S_peak_ind, Fs);
VFL_QS_int = calc_qs (VFL_Q_peak_ind,VFL_S_peak_ind,Fs);

%% Structured Approach
% Initialize the main data structure to hold arrhythmia features
arrhythmiaData = struct();

count = 1;
while count <= length(comments)
    rhythm = cell2mat(comments(count));
    
    % Identify the rhythm type
    if (length(rhythm) == 4) && all(rhythm == '(VFL')
        rhythmType = 'VFL';
    elseif (length(rhythm) == 2) && all(rhythm == '(N')
        rhythmType = 'N';
    else
        count = count + 1; % Skip unrecognized rhythms
        continue;
    end
    
    % If rhythmType doesn't exist in arrhythmiaData, initialize it
    if ~isfield(arrhythmiaData, rhythmType)
        arrhythmiaData.(rhythmType).R_peak_ind = [];
        arrhythmiaData.(rhythmType).Q_peak_ind = [];
        arrhythmiaData.(rhythmType).T_peak_ind = [];
        arrhythmiaData.(rhythmType).S_peak_ind = [];
    end
    
    % Find the start and end of the current rhythm section
    start_count = ann(count);
    while (count <= length(comments)) && ...
          (length(cell2mat(comments(count))) == length(rhythm)) && ...
          all(cell2mat(comments(count)) == rhythm)
        count = count + 1;
    end
    end_count = ann(count); % Mark the end of the section
    
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


%% Feature making
X = [[VFL_R_peak_ind N_R_peak_ind]' , [VFL_Q_peak_ind N_Q_peak_ind]'... 
    , [VFL_S_peak_ind N_S_peak_ind]' , [VFL_T_peak_ind N_T_peak_ind]'...
    [VFL_RR_int N_RR_int]', [VFL_QS_int N_QS_int]'];


Y  = (1:length(X));
Y (1: length(VFL_T_peak_ind)) = 1;
Y (length(VFL_T_peak_ind):end) = 2;


%% ========================= SVM =========================

rand_num = randperm(size(X,1));
X_train = X(rand_num(1:round(0.8*length(rand_num))),:);
y_train = Y(rand_num(1:round(0.8*length(rand_num))));

X_test = X(rand_num(round(0.8*length(rand_num))+1:end),:);
y_test = Y(rand_num(round(0.8*length(rand_num))+1:end));
%% CV partition
y_train = y_train';
y_test = y_test';
c = cvpartition(y_train,'k',5);
%% feature selection

opts = statset('display','iter');
classf = @(train_data, train_labels, test_data, test_labels)...
    sum(predict(fitcsvm(train_data, train_labels,'KernelFunction','rbf'), test_data) ~= test_labels);

[fs, history] = sequentialfs(classf, X_train, y_train, 'cv', c, 'options', opts,'nfeatures',4);
%% Best hyperparameter

X_train_w_best_feature = X_train(:,fs);

Md1 = fitcsvm(X_train_w_best_feature,y_train,'KernelFunction','linear','OptimizeHyperparameters','auto',...
      'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
      'expected-improvement-plus','ShowPlots',true)); % Bayes' Optimization 사용.


%% Final test with test set
X_test_w_best_feature = X_test(:,fs);
test_accuracy_for_iter = sum((predict(Md1,X_test_w_best_feature) == y_test))/length(y_test)*100

%% hyperplane 확인

figure;
hgscatter = gscatter(X_train_w_best_feature(:,3),X_train_w_best_feature(:,4),y_train);
hold on;
h_sv=plot(Md1.SupportVectors(:,1),Md1.SupportVectors(:,2),'ko','markersize',8);


% test set의 data를 하나 하나씩 넣어보자.

gscatter(X_test_w_best_feature(:,1),X_test_w_best_feature(:,2),y_test,'rb','xx')

% decision plane
XLIMs = get(gca,'xlim');
YLIMs = get(gca,'ylim');
[xi,yi] = meshgrid([XLIMs(1):1000000:XLIMs(2)],[YLIMs(1):1000000:YLIMs(2)]);
dd = [xi(:), yi(:)];
pred_mesh = predict(Md1, dd);
redcolor = [1, 0.8, 0.8];
bluecolor = [0.8, 0.8, 1];
pos = find(pred_mesh == 1);
h1 = plot(dd(pos,1), dd(pos,2),'s','color',redcolor,'Markersize',5,'MarkerEdgeColor',redcolor,'MarkerFaceColor',redcolor);
pos = find(pred_mesh == 2);
h2 = plot(dd(pos,1), dd(pos,2),'s','color',bluecolor,'Markersize',5,'MarkerEdgeColor',bluecolor,'MarkerFaceColor',bluecolor);
uistack(h1,'bottom');
uistack(h2,'bottom');
legend([hgscatter;h_sv],{'setosa','versicolor','support vectors'})
