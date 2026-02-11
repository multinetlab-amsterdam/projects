#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 11:38:05 2025

@author: ekoderman

status: final

REVIEW
reviewed by: Sebastien Dam
review date: January 2026

This script is used for testing the chosen model on the 30% held-out set.

"""
from pathlib import Path
import pandas as pd
import numpy as np
from utils import bootstrap_auc_multiclass, plot_multiclass_roc_cv, compute_foldwise_metrics, XGBOrdinalWrapper, evaluate_across_folds
import joblib
from collections import Counter
import shap
import matplotlib.pyplot as plt
import json

# Paths
df_path = Path("/folder/data/processed/combined/combined_report_05112025.csv")
train_idx_path = Path("/folder/subj_ids/train_test_IDs/random_stratified/development_idx.npy")
test_idx_path = Path("/folder/subj_ids/train_test_IDs/random_stratified/external_test_set_idx.npy")
results_path = Path("/folder/models/ordinal_reg_xgb/top3_class1w_iter150_27112025.pkl")

# Load
combined_df = pd.read_csv(df_path)
combined_df.drop(columns='Unnamed: 0', inplace=True)
train_idx = np.load(train_idx_path, allow_pickle=True)
test_idx = np.load(test_idx_path, allow_pickle=True)
results_list = joblib.load(results_path)

## make 3 classed outcome
conditions = [combined_df['KPS_post-op_value'] == 0,
              combined_df['KPS_post-op_value'] < 70,
              combined_df['KPS_post-op_value'] >= 70]
values = [0,1,2]

combined_df['KPS_post_op_binned'] = np.select(conditions,values)

# Select development set and reset index
dev_df = combined_df.iloc[train_idx].reset_index(drop=True)
test_df = combined_df.iloc[test_idx].reset_index(drop=True)

print(f'DEVELOPMENT: {dev_df["KPS_post_op_binned"].value_counts(normalize=True)*100}')
print(f'TEST: {test_df["KPS_post_op_binned"].value_counts(normalize=True)*100}')

#%%
# final check of the development set performance
eval_results = evaluate_across_folds(results_list, class_labels=[0,1,2])
summary, roc_scores, pr_scores = bootstrap_auc_multiclass(results_list, n_iterations=100, plot_ci=False)

#%% # 1. pick the best params from all the folds

# Compute weighted ROC AUC for each fold

# Collect all hyperparams
param_lists = [fold['best_params'] for fold in results_list]

# Suppose you have only numeric hyperparameters for simplicity
final_params = {}
for key in param_lists[0].keys():
    values = [p[key] for p in param_lists]
    if isinstance(values[0], (int, float)):
        val = np.median(values)
        # cast to int if this param should be integer
        if key in ['n_estimators', 'max_depth', 'min_child_weight']:
            val = int(val)
        final_params[key] = val
    else:
        # For categorical params, pick the most common
        final_params[key] = Counter(values).most_common(1)[0][0]

print(final_params)

#%% # 2. train the final model on the full development set

X_dev = dev_df[['age_at_surgery', 'KPS_pre-op_value', 'enhancing_part_volume']] # based on the top3 features that the pre-trained model was trained on
y_dev = dev_df['KPS_post_op_binned']

#%%
# compute sample weights as done in training
class_counts = Counter(y_dev)
total = sum(class_counts.values())
weight_for_class = {cls: total / class_counts[cls] for cls in class_counts}
sample_weight = np.array([weight_for_class[v] for v in y_dev])
sample_weight = sample_weight / np.mean(sample_weight)

# train final XGBOrdinal model
final_model_wrapper = XGBOrdinalWrapper(random_state=42, aggregation='weighted', norm=True, **final_params)
final_model_wrapper.fit(X_dev, y_dev, sample_weight=sample_weight)
final_model = final_model_wrapper.model

#%% 3. make predictions on the held-out set

X_test = test_df[['age_at_surgery', 'KPS_pre-op_value', 'enhancing_part_volume']]
y_test = test_df['KPS_post_op_binned']

y_pred_test = final_model.predict(X_test)
y_proba_test = final_model.predict_proba(X_test)

# 4. compute SHAP for the external set
shap_per_class_test = []
for i, clf in final_model.clfs.items():
    explainer = shap.TreeExplainer(clf)
    shap_values = explainer.shap_values(X_test)
    shap_per_class_test.append(shap_values)

# 5. collect results
external_results = {
    'y_true': y_test.values,
    'y_pred': y_pred_test,
    'y_proba': y_proba_test,
    'shap_per_class': shap_per_class_test,
    'X_test': X_test
}

#%%

summary, roc_scores, pr_scores = bootstrap_auc_multiclass([external_results], n_iterations=500, plot_ci=True)

#%%

eval_results = evaluate_across_folds([external_results], class_labels=[0,1,2])

#%%

plot_multiclass_roc_cv([external_results], class_labels=[0,1,2])

#%%

metrics_df = compute_foldwise_metrics([external_results], class_labels=[0,1,2])
print(metrics_df)

#%%

# Extract number of classifiers (one per ordinal threshold)
n_classifiers = len(external_results['shap_per_class'])
feature_names = external_results['X_test'].columns

# Loop over classifiers and plot
for i, class_shap in enumerate(external_results['shap_per_class']):
    plt.figure(figsize=(10, 6))
    shap.summary_plot(class_shap, external_results['X_test'], feature_names=feature_names, show=False)
    plt.title(f'SHAP summary for classifier threshold {i} (class > {i})')
    plt.tight_layout()
    plt.show()

#%%
## store
output_results = '/folder/models/ordinal_reg_xgb/final_model_held-out_results_01122025.pkl'
joblib.dump(external_results, output_results)

output_model = '/folder/models/ordinal_reg_xgb/final_model_for_mcp_01122025.pkl'
joblib.dump(final_model, output_model)

output_final_params = '/folder/models/ordinal_reg_xgb/final_model_best_params_01122025.json'
with open(output_final_params, 'w') as f:
    json.dump(final_params, f, indent=4)