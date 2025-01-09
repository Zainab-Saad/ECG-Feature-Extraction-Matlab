display('Reading samples ECG signal from MIT-BIH Arrhythmia Database')
[ecg,Fs,tm]=rdsamp('mitdbase/207',1);

display('Reading and plotting annotations (human labels) of QRS complexes performend on the signals')
[ann,type,subtype,chan,num]=rdann('mitdbase/207','atr',1);



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
    c