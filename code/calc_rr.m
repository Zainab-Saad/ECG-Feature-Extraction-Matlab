function [RR_int] = calc_rr(R_peaks_ind,Fs)
count = 1;
RR_int = [];
for i = R_peaks_ind
    count = count +1;
    if count >2
        RR_int = [RR_int , (i-temp_i)/Fs];
    end
    temp_i = i;
end
RR_int = [mean(RR_int) RR_int];

end