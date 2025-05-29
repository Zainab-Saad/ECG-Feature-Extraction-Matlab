ECG Feature Extraction and Classification Project

This project processes ECG signals to detect arrhythmias using feature extraction and machine learning.

Project Structure:
1. databases/ - Contains ECG signal files (.hea and .dat files)
2. train_multiclass.m - Main training script that:
   - Reads ECG signals
   - Extracts features using Pan-Tompkins algorithm
   - Trains SVM model using 5-fold cross-validation
   - Saves trained model as 'afib_model.mat'

3. inference.m - Testing script that:
   - Loads the trained model
   - Processes test data
   - Evaluates model performance
   - Shows confusion matrix and metrics

How to Use:
1. Place your ECG signal files in the 'databases' folder
2. Run train_multiclass.m to train the model
3. Run inference.m to test the model

Required MATLAB Toolboxes:
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- WFDB Toolbox (for reading ECG files)

Note: Make sure all required toolboxes are installed before running the scripts.
