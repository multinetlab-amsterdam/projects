#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 14:45:09 2025

@author: ekoderman

status: in-progress

REVIEW
reviewed by: Sebastien Dam
review date: January 2026

This script is used for performing the Mondrian conformal prediction (MCP) to obtain model prediction uncertainty estimates.
Ref: https://mapie.readthedocs.io/en/latest/theoretical_description_mondrian.html

"""
from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict
from math import ceil
from utils import summary_mondrian_conformal
import joblib
from utils import plot_patient_prediction

# Paths
df_path = Path("/folder/data/processed/combined/combined_report_05112025.csv")
train_idx_path = Path("/folder/subj_ids/train_test_IDs/random_stratified/development_idx.npy")
test_idx_path = Path("/folder/subj_ids/train_test_IDs/random_stratified/external_test_set_idx.npy")

oof_results_path = Path("/folder/models/ordinal_reg_xgb/top3_class1w_iter150_27112025.pkl")
final_model_path = Path("/folder/models/ordinal_reg_xgb/final_model_for_mcp_01122025.pkl")
external_results_path = Path('/folder/models/ordinal_reg_xgb/final_model_held-out_results_01122025.pkl')

output_dir = Path('/folder/models/mcp')

# Load
combined_df = pd.read_csv(df_path)
combined_df.drop(columns='Unnamed: 0', inplace=True)
train_idx = np.load(train_idx_path, allow_pickle=True)
test_idx = np.load(test_idx_path, allow_pickle=True)

results_oof = joblib.load(oof_results_path)
final_model = joblib.load(final_model_path)
results_external = joblib.load(external_results_path)

## make 3 classed outcome
conditions = [combined_df['KPS_post-op_value'] == 0,
              combined_df['KPS_post-op_value'] < 70,
              combined_df['KPS_post-op_value'] >= 70]
values = [0,1,2]

combined_df['KPS_post_op_binned'] = np.select(conditions,values)

# Select development set and reset index
dev_df = combined_df.iloc[train_idx].reset_index(drop=True)
test_df = combined_df.iloc[test_idx].reset_index(drop=True)

#%%
print(f'DEVELOPMENT: {dev_df["KPS_post_op_binned"].value_counts(normalize=True)*100}')
print(f'TEST: {test_df["KPS_post_op_binned"].value_counts(normalize=True)*100}')

#%%

# in this first step we compute empirical distribution of model error to compute thresholds (quantiles) that guarantee
# that in the set of classes chosen the classes appeared also in the true distribution with a 90% (determined by alpha) certainty

# OOF = using only the folds that were used for testing
# Providing the empirical distribution of model error for each class (how wrong was the model)
# nonconformity score = how wrong the model is for this sample
# Using OOF scores ensures the error distribution is unbiased.
# The quantiles are computed per class (Mondrian) to handle imbalance.
# These thresholds can then be applied to external predictions to generate conformal prediction sets.

##  1) Build OOF nonconformity scores for mcp
n_classes = results_oof[0]['y_proba'].shape[1]
nonconf_by_class = defaultdict(list)
for fold in results_oof:
    y_true = np.asarray(fold['y_true'])
    y_proba = np.asarray(fold['y_proba'])   # shape (n_fold_samples, n_classes)
    for i, y in enumerate(y_true):
        #p_true = y_proba[i, y]
        #p_other_classes = np.delete(y_proba[i], y)
        #y_expected = np.sum(np.arange(n_classes) * y_proba[i])
        #s = abs(y_expected - y)
        #s = -np.log(np.clip(p_true, 1e-12, None))
        #s = max(p_other_classes) - p_true
        p_true = y_proba[i, y]
        #y_expected = np.sum(np.arange(n_classes) * y_proba[i])
        #s = 0.5*(1 - p_true) + 0.8*abs(y_expected - y)
        s = 1 - p_true
        nonconf_by_class[y].append(s)

# nonconf_by_clas: empirical distribution of model error
# the error distribution reflects “unseen” performance, which is what Mondrian CP needs to compute proper thresholds.

# below gives us a empirical distribution of errors from which quantiles (thresholds) can be calculated

# --- 2) Compute Mondrian q_c per class with finite-sample correction ---

def compute_q_from_scores(nonconf_by_class, alpha_dict=None):
    """
    Compute Mondrian thresholds q_c per class using class-specific alpha.
    
    Parameters
    ----------
    nonconf_by_class : dict
        {class: list of nonconformity scores}
    alpha_dict : dict, optional
        {class: alpha} for per-class miscoverage. If None, uses 0.1 for all.
        
    Returns
    -------
    q_dict : dict
        Threshold per class
    n_dict : dict
        Number of OOF calibration examples per class
    """
    q_dict = {} # the thresholds
    # Any predicted probability for class c that gives 1 - prob <= q_dict[c] is included in the prediction set
    # Intuition: it’s the (1 - alpha) quantile of the class-specific error distribution
    n_dict = {} ## number of OOF calibration examples available for class c (distribution)
    for c in range(n_classes):
        scores = np.array(nonconf_by_class.get(c, []))
        n_c = len(scores)
        n_dict[c] = n_c
        
        if n_c == 0:
            q = 1.0 # This happens if there were no OOF samples for this class, e.g., underrepresented class 1 in some folds.
            # q is the maximum allowed nonconformity for a class to be included in the prediction set
            # It’s “conservative” because you don’t exclude a class just because you lack calibration examples
        else:
            alpha_c = alpha_dict.get(c, 0.1) if alpha_dict else 0.1
            k = int(ceil((n_c + 1) * (1 - alpha_c))) - 1
            k = max(0, min(k, n_c - 1))
            q = np.sort(scores)[k]
        q_dict[c] = q # Threshold for class c (Mondrian quantile) to achieve coverage 1-alpha
    return q_dict, n_dict

# Compute per-class thresholds
# alpha = fraction of times the true label is allowed to fall outside the predicted set.
# 1-alpha: desired coverage

# various alpha parameters were tried out empirically, after which we decided to go with the same alpha value for all classes
# 0.15 seems a good enough cut-off for clinical practice as there are currently no standards but future studies could search for better ways to define this threshold

alpha_per_class = {0: 0.15, 1: 0.15, 2: 0.15} # desired miscoverage 
q_dict, n_per_class = compute_q_from_scores(nonconf_by_class, alpha_dict=alpha_per_class)
print("Calibration counts per class (OOF):", n_per_class) ## number of OOF calibration examples available for class c (distribution)
print("q per class (OOF):", q_dict) # Threshold for class c (Mondrian quantile) to achieve coverage 1-alpha

#%%

# --- 6) Apply Mondrian thresholds to external set ---
def mondrian_predict_sets_from_q(final_model, X, q_dict):
    probs = np.asarray(final_model.predict_proba(X))
    n, C = probs.shape
    pred_sets = []
    # optionally return per-class nonconformity s too
    for i in range(n):
        incl = []
        for c in range(C):
            s = 1.0 - probs[i, c]
            if s <= q_dict[c]:
                incl.append(c)
        pred_sets.append(incl)
    return pred_sets, probs

X_test = results_external['X_test']
y_test = results_external['y_true']

pred_sets_ext, probs_ext = mondrian_predict_sets_from_q(final_model, X_test, q_dict)

# --- 5) Compute empirical coverage per class on external set ---
cover_counts = defaultdict(int)
total_counts = defaultdict(int)

y_ext_arr = np.asarray(y_test)  # your external true labels
for i, ytrue in enumerate(y_ext_arr):
    total_counts[ytrue] += 1
    if ytrue in pred_sets_ext[i]:
        cover_counts[ytrue] += 1

for c in range(n_classes):
    n = total_counts.get(c, 0)
    covered = cover_counts.get(c, 0)
    print(f"class {c}: coverage = {covered}/{n} = {(covered/n if n>0 else float('nan')):.3f} (n={n})")

proba_list = [row for row in results_external['y_proba']]

#print(pred_sets_ext[['y_true', 'y_pred', 'prediction_set']].head(10))
df_pred_sets = pd.DataFrame({
    'y_true': y_test,
    'y_pred': results_external['y_pred'],
    'y_proba': proba_list,
    'prediction_set': pred_sets_ext
})

# Add nonconformity scores per sample to df_pred_sets
df_pred_sets["scores"] = df_pred_sets["y_proba"].apply(lambda p: 1 - np.array(p))
df_pred_sets["tau"] = [q_dict] * len(df_pred_sets)   # store thresholds
summary_mondrian_conformal(df_pred_sets)

#%%

def extract_patient_features(results, patient_index):
    """
    Extract the feature names and their values for a specific patient from the results.

    Parameters:
    -----------
    results : dict
        Dictionary containing the keys 'X_test' and 'feature_set_name', which hold
        the patient feature data and the feature names.
    patient_index : int
        The index of the patient in the test set.

    Returns:
    --------
    feature_names : list
        List of feature names.
    feature_values : list
        List of feature values for the specified patient.
    """
    X_test = results['X_test']  # Assuming X_test is a DataFrame or ndarray
    feature_set_name = results.get('feature_set_name', [])

    if isinstance(X_test, pd.DataFrame):
        # If X_test is a DataFrame, columns should be the feature names
        feature_names = X_test.columns.tolist()
        feature_values = X_test.iloc[patient_index].tolist()
    elif isinstance(X_test, np.ndarray):
        # If X_test is a numpy array, we may need to access columns by index
        feature_names = feature_set_name
        feature_values = X_test[patient_index].tolist()
    else:
        raise ValueError("X_test is neither a DataFrame nor a numpy array")

    return feature_names, feature_values

# just to play around:
# IM0379 = idx 36
idx = 36
feature_names_idx, feature_values_idx = extract_patient_features(results_external, idx)

patient_row = df_pred_sets.iloc[idx]
plot_patient_prediction(patient_row, class_labels=[0,1,2], feature_names=feature_names_idx, feature_values=feature_values_idx)

#%% # store

class_mapping = {
    0: "KPS=0",
    1: "KPS<70",
    2: "KPS>=70"
}

mondrian_cp = {
    "model": final_model,
    "q_dict": q_dict,
    "alpha_per_class": alpha_per_class,
    "class_mapping": class_mapping,
    "calibration_counts": n_per_class,
    "method": "Mondrian CP, OOF, score=1-p_true",
}

joblib.dump(mondrian_cp, output_dir / 'mondrian_cp.pkl')