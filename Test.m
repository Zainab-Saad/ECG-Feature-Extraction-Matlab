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
