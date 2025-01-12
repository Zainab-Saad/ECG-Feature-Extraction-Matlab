function [QS_int] = calc_qs(Q_peaks_ind, S_peaks_ind,Fs)
i =1 ;
QS_int = [];
while i <= length(S_peaks_ind)
    QS_int = [QS_int (S_peaks_ind(i)-Q_peaks_ind(i))/Fs];
    i = i+1;
end
QS_int (1) = mean(QS_int);
end