display('Reading samples ECG signal from MIT-BIH Arrhythmia Database')
[ecg,Fs,tm]=rdsamp('mitdb/102',1);

display('Reading and plotting annotations (human labels) of QRS complexes performend on the signals')
[ann,type,subtype,chan,num]=rdann('mitdb/102','atr',1);
%%

[qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(ecg(1:5000), Fs);

