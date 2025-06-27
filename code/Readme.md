# ECG Arrhythmia Classification Project

This project implements an ECG signal processing and classification pipeline in MATLAB. It focuses on extracting key features from ECG waveforms and using a Linear Support Vector Machine (SVM) to classify different types of arrhythmias.

## Project Structure

The core MATLAB code is primarily located within the `code/` directory. The main scripts and functions are described below.

## Code Files (`code/`)

This directory contains the essential MATLAB functions and the main script for the project.

### pan_tompkin.m

**Purpose:** Implements the Pan-Tompkins algorithm for robust detection of QRS complexes (R peaks) in an ECG signal. It also includes functionality to locate the associated Q, S, and T peaks.

**Inputs:**
- `ecg`: The raw ECG signal vector.
- `fs`: The sampling frequency of the ECG signal (in Hz).
- `gr`: A flag (0 or 1) to enable or disable plotting of intermediate and final results.

**Outputs:**
- `qrs_amp_raw`: Amplitudes of the detected R peaks in the bandpass-filtered signal.
- `qrs_i_raw`: Indices of the detected R peaks in the bandpass-filtered signal.
- `Q_peaks`, `Q_peaks_val`: Indices and values of detected Q peaks.
- `S_peaks`, `S_peaks_val`: Indices and values of detected S peaks.
- `T_peaks`, `T_peaks_val`: Indices and values of detected T peaks.
- `delay`: Processing delay due to filtering.

**Key Functionality:**
- Applies various filters (bandpass, derivative, moving average) to the ECG signal.
- Uses adaptive thresholding based on signal and noise levels to identify QRS locations.
- Searches for Q, S, and T peaks relative to the detected R peaks.
- Includes options for visual verification through plotting.

### calc_rr.m

**Purpose:** Calculates the time intervals between consecutive R peaks (RR intervals).

**Inputs:**
- `R_peaks_ind`: A vector of indices for the detected R peaks.
- `Fs`: The sampling frequency.

**Output:**
- `RR_int`: A vector of calculated RR intervals in seconds.

**Key Functionality:**
- Computes the difference in time between successive R peak locations.
- Note: The function prepends the mean of all calculated RR intervals to the output vector.

### calc_qs.m

**Purpose:** Calculates the time intervals between detected Q and S peaks (QS intervals).

**Inputs:**
- `Q_peaks_ind`: A vector of indices for the detected Q peaks.
- `S_peaks_ind`: A vector of indices for the detected S peaks.
- `Fs`: The sampling frequency.

**Output:**
- `QS_int`: A vector of calculated QS intervals in seconds.

**Key Functionality:**
- Computes the time difference between each S peak and its corresponding Q peak.
- Note: The function sets the first element of the output vector QS_int to the mean of all calculated QS intervals.

### train_VT.m

**Purpose:** This script trains a Support Vector Machine (SVM) classifier to detect Ventricular Tachycardia (VT) from ECG data. It reads ECG files, extracts features (R, Q, S, T peak values, RR and QS intervals), and assigns VT labels using annotation logic, with special handling for CUDB files (using '[' and ']' markers). The script performs k-fold cross-validation, prints performance metrics (accuracy, precision, recall, F1 score, AUC), and saves the trained VT model for later inference.

**Inputs:** Reads ECG data and annotations from files in the `databases/` directory. Relies on `pan_tompkin.m`, `calc_rr.m`, and `calc_qs.m` for feature extraction and peak detection.

**Outputs:** Prints cross-validation results and performance metrics for each fold and overall. Saves the trained VT SVM model and related variables to `vt_model.mat` for use in inference scripts.

**Key Functionality:**
- Iterates through ECG files in the specified database folder.
- Uses `rdsamp` and `rdann` to read signals and annotations.
- Calls `pan_tompkin` to detect R peaks and extract Q, S, T peaks.
- Calls `calc_rr` and `calc_qs` to compute RR and QS intervals.
- Assigns VT and normal rhythm labels using annotation comments and CUDB-specific logic.
- Aggregates features and labels from all files into matrices for training.
- Performs k-fold cross-validation with stratification, calculates and prints accuracy, precision, recall, F1, and AUC for each fold.
- Allows for custom decision threshold for VT detection.
- Plots confusion matrices and ROC curves for each fold.
- Saves the final trained model and feature information for later use.

### train_afib.m

**Purpose:** This script trains a Support Vector Machine (SVM) classifier to detect Atrial Fibrillation (AFIB) from ECG data. It reads ECG files, extracts features (R, Q, S, T peak values, RR and QS intervals), and assigns AFIB labels by propagating rhythm annotations using the '+' marker logic. The script performs k-fold cross-validation, prints performance metrics (accuracy, precision, recall, F1 score, AUC), and saves the trained AFIB model for later inference.

**Inputs:** Reads ECG data and annotations from files in the `databases/` directory. Relies on `pan_tompkin.m`, `calc_rr.m`, and `calc_qs.m` for feature extraction and peak detection.

**Outputs:** Prints cross-validation results and performance metrics for each fold and overall. Saves the trained AFIB SVM model and related variables to `afib_model.mat` for use in inference scripts.

**Key Functionality:**
- Iterates through ECG files in the specified database folder.
- Uses `rdsamp` and `rdann` to read signals and annotations.
- Calls `pan_tompkin` to detect R peaks and extract Q, S, T peaks.
- Calls `calc_rr` and `calc_qs` to compute RR and QS intervals.
- Propagates rhythm labels using annotation comments and the '+' marker to assign AFIB and normal labels.
- Aggregates features and labels from all files into matrices for training.
- Performs k-fold cross-validation with stratification, calculates and prints accuracy, precision, recall, F1, and AUC for each fold.
- Allows for custom decision threshold for AFIB detection.
- Plots confusion matrices and ROC curves for each fold.
- Saves the final trained model and feature information for later use.

### inference_afib.m

**Purpose:** This script performs per-file inference for Atrial Fibrillation (AFIB) detection using a pre-trained SVM model. It reads ECG files, extracts features (R, Q, S, T peak values, RR and QS intervals), loads the AFIB model, and predicts AFIB presence for each beat. The script visualizes predictions alongside the ECG signal and provides summary statistics.

**Inputs:** Reads ECG data and annotations from files in the `databases/` directory. Requires the trained AFIB SVM model (`afib_model.mat`) and uses `pan_tompkin.m`, `calc_rr.m`, and `calc_qs.m` for feature extraction and peak detection.

**Outputs:** For each processed ECG file, displays a plot of the ECG signal with predicted AFIB/Normal labels, and prints summary statistics (e.g., number of AFIB/Normal beats, accuracy if ground truth is available). Handles errors and missing data robustly.

**Key Functionality:**
- Iterates through specified ECG files for inference.
- Uses `rdsamp` and `rdann` to read signals and annotations.
- Calls `pan_tompkin` to detect R peaks and extract Q, S, T peaks.
- Calls `calc_rr` and `calc_qs` to compute RR and QS intervals.
- Loads the trained AFIB SVM model and feature information from `afib_model.mat`.
- Extracts features for each beat and applies the SVM model to predict AFIB or Normal rhythm.
- Visualizes predictions on the ECG plot, marking AFIB and Normal beats distinctly.
- Prints summary statistics for each file, including counts of predicted classes and (if available) accuracy versus ground truth.
- Includes error handling for missing data, annotation mismatches, or feature extraction issues.

### inference_VT.m

**Purpose:** This script performs per-file inference for Ventricular Tachycardia (VT) detection using a pre-trained SVM model. It reads ECG files, extracts features (R, Q, S, T peak values, RR and QS intervals), loads the VT model, and predicts VT presence for each beat. The script visualizes predictions alongside the ECG signal and provides summary statistics.

**Inputs:** Reads ECG data and annotations from files in the `databases/` directory. Requires the trained VT SVM model (`vt_model.mat`) and uses `pan_tompkin.m`, `calc_rr.m`, and `calc_qs.m` for feature extraction and peak detection.

**Outputs:** For each processed ECG file, displays a plot of the ECG signal with predicted VT/Normal labels, and prints summary statistics (e.g., number of VT/Normal beats, accuracy if ground truth is available). Handles errors and missing data robustly.

**Key Functionality:**
- Iterates through specified ECG files for inference.
- Uses `rdsamp` and `rdann` to read signals and annotations.
- Calls `pan_tompkin` to detect R peaks and extract Q, S, T peaks.
- Calls `calc_rr` and `calc_qs` to compute RR and QS intervals.
- Loads the trained VT SVM model and feature information from `vt_model.mat`.
- Extracts features for each beat and applies the SVM model to predict VT or Normal rhythm.
- Visualizes predictions on the ECG plot, marking VT and Normal beats distinctly.
- Prints summary statistics for each file, including counts of predicted classes and (if available) accuracy versus ground truth.
- Includes error handling for missing data, annotation mismatches, or feature extraction issues.
