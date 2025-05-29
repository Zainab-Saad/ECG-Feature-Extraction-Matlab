ECG Arrhythmia Classification Project

This project implements an ECG signal processing and classification pipeline in MATLAB. It focuses on extracting key features from ECG waveforms and using a Linear Support Vector Machine (SVM) to classify different types of arrhythmias.

Project Structure

The core MATLAB code is primarily located within the code/ directory. The main scripts and functions are described below.

Code Files (code/)

This directory contains the essential MATLAB functions and the main script for the project.

pan_tompkin.m

Purpose: Implements the Pan-Tompkins algorithm for robust detection of QRS complexes (R peaks) in an ECG signal. It also includes functionality to locate the associated Q, S, and T peaks.

Inputs:
- ecg: The raw ECG signal vector.
- fs: The sampling frequency of the ECG signal (in Hz).
- gr: A flag (0 or 1) to enable or disable plotting of intermediate and final results.

Outputs:
- qrs_amp_raw: Amplitudes of the detected R peaks in the bandpass-filtered signal.
- qrs_i_raw: Indices of the detected R peaks in the bandpass-filtered signal.
- Q_peaks, Q_peaks_val: Indices and values of detected Q peaks.
- S_peaks, S_peaks_val: Indices and values of detected S peaks.
- T_peaks, T_peaks_val: Indices and values of detected T peaks.
- delay: Processing delay due to filtering.

Key Functionality:
- Applies various filters (bandpass, derivative, moving average) to the ECG signal.
- Uses adaptive thresholding based on signal and noise levels to identify QRS locations.
- Searches for Q, S, and T peaks relative to the detected R peaks.
- Includes options for visual verification through plotting.

calc_rr.m

Purpose: Calculates the time intervals between consecutive R peaks (RR intervals).

Inputs:
- R_peaks_ind: A vector of indices for the detected R peaks.
- Fs: The sampling frequency.

Output:
- RR_int: A vector of calculated RR intervals in seconds.

Key Functionality:
- Computes the difference in time between successive R peak locations.
- Note: The function prepends the mean of all calculated RR intervals to the output vector.

calc_qs.m

Purpose: Calculates the time intervals between detected Q and S peaks (QS intervals).

Inputs:
- Q_peaks_ind: A vector of indices for the detected Q peaks.
- S_peaks_ind: A vector of indices for the detected S peaks.
- Fs: The sampling frequency.

Output:
- QS_int: A vector of calculated QS intervals in seconds.

Key Functionality:
- Computes the time difference between each S peak and its corresponding Q peak.
- Note: The function sets the first element of the output vector QS_int to the mean of all calculated QS intervals.

multiclass_SVM.m

Purpose: This is a main script that orchestrates the entire process: loading ECG data, extracting features (peak values, RR/QS intervals), preparing data for machine learning, and training/testing a multiclass SVM classifier using a One-vs-All strategy.

Inputs: Reads ECG data and annotations from files in the mitdbase/ directory. Relies on pan_tompkin.m, calc_rr.m, and calc_qs.m.

Outputs: Prints the test accuracy of the trained SVM classifier. Generates various data variables in the MATLAB workspace.

Key Functionality:
- Iterates through ECG files in the specified database folder.
- Uses rdsamp and rdann to read signals and annotations.
- Calls pan_tompkin to find ECG peaks.
- Calls calc_rr and calc_qs to compute intervals.
- Processes annotations to assign rhythm labels ('N', 'VFL', 'VT', 'AFIB', 'BII', etc.) to segments of the signal.
- Extracts and organizes features (peak values, intervals) per rhythm type.
- Aggregates features and numerical labels from all files into large matrices (all_X, all_Y).
- Splits the data into training and testing sets (80/20 split, currently sequential).
- Implements a One-vs-All multiclass SVM using fitcsvm, training a binary classifier for each class against all others.
- Predicts labels on the test set using the One-vs-All approach and evaluates accuracy.
- Includes commented-out sections for alternative classification/tuning methods, including preparing data for GPU processing. 