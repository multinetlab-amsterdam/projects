#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 12:20:49 2025

@author: ekoderman

status: final

REVIEW
reviewed by: Sebastien Dam
review date: January 2026

This script is used for training models (XGB ordinal regression) on the 70% training set.
It starts with the function used to do that: train_and_evaluate_ordinal_shap
It is meant to train several models by re-running with different feature sets.
Additionally, it investigate the SHAP values per classifier and global.

A custom wrapper (XGBOrdinalWrapper) was made for exposing XGBOrdinal to BayesSearchCV.

"""
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from skopt import BayesSearchCV
import shap
from utils import evaluate_across_folds
import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter
import joblib
from utils import XGBOrdinalWrapper, compute_foldwise_metrics, plot_multiclass_roc_cv
from datetime import datetime
import json

def train_and_evaluate_ordinal_shap(X, y, cv_outer, param_space, n_iter, feature_set_name):
    """
    Trains XGBOrdinal using outer CV with Bayesian search for hyperparams.
    Collects per-class and global SHAP values.
    Returns:
        - results_list: dict with fold results
        - per_class_shap_all_folds: list of concatenated SHAP values per classifier
        - feature_names
        - classes
    """
    
    results_list = []
    per_class_shap_all_folds = []
    feature_names = X.columns.tolist()
    classes = np.unique(y)

    for fold_idx, (train_idx, test_idx) in enumerate(cv_outer.split(X, y)):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
                
       # Higher weight → the model tries harder to fit that sample, because misclassifying it increases the loss more.
        #Lower weight → the model “cares less” about getting that sample right.
        
        class_counts = Counter(y_train)
        total = sum(class_counts.values())
        weight_for_class = {cls: total / class_counts[cls] for cls in class_counts}
        sample_weight = np.array([weight_for_class[v] for v in y_train])
        sample_weight = sample_weight / np.mean(sample_weight)
        
        inner_cv = StratifiedKFold(n_splits=3) 
        
        # with current setup = 8 cores needed
        # ---- Bayesian search for this fold ----
        xgb_wrapped = XGBOrdinalWrapper(random_state=42, n_jobs=2, aggregation='weighted', norm=True, tree_method='hist') #hist = better CPU performance
        optimizer = BayesSearchCV(
            xgb_wrapped,
            param_space,
            n_iter=n_iter,
            cv=inner_cv,
            scoring="roc_auc_ovr_weighted",
            n_jobs=4,
            refit=True,
            verbose=2,
            random_state=42
        )
        
        optimizer.fit(X_train, y_train, sample_weight=sample_weight)
        best_model = optimizer.best_estimator_.model  # the trained XGBOrdinal

        # ---- Predictions ----
        y_pred = best_model.predict(X_test)
        y_proba = best_model.predict_proba(X_test)

        # ---- SHAP per-class ----
        shap_per_class = []
        for i, clf in best_model.clfs.items():
            explainer = shap.TreeExplainer(clf)
            shap_values = explainer.shap_values(X_test)
            shap_per_class.append(shap_values)

        # Store per-class SHAP values across folds
        if not per_class_shap_all_folds:
            # Initialize list of lists
            per_class_shap_all_folds = [[] for _ in range(len(shap_per_class))]
        for i, class_shap in enumerate(shap_per_class):
            per_class_shap_all_folds[i].append(class_shap)

        # ---- Save results ----
        results_list.append({
            'feature_set_name': feature_set_name,
            'fold': fold_idx,
            'best_params': optimizer.best_params_,
            'best_model': best_model,
            'y_true': y_test.values,
            'y_pred': y_pred,
            'y_proba': y_proba,
            'shap_per_class': shap_per_class,
            'X_test': X_test
        })

        print(f"Fold {fold_idx} done. Best params: {optimizer.best_params_}")

    # Concatenate per-class SHAP across folds
    per_class_shap_all_folds = [np.vstack(shap_list) for shap_list in per_class_shap_all_folds]

    return results_list, per_class_shap_all_folds, feature_names, classes

#%% ## Prepare data

# Paths
df_path = Path("/folder/data/processed/combined/combined_df_strata_groups_22122025.csv")
train_idx_path = Path("/folder/KPS_prediction_code_review/subj_ids/train_test_IDs/random_stratified/development_idx.npy")

# Load
combined_df = pd.read_csv(df_path)
combined_df.drop(columns='Unnamed: 0', inplace=True)
train_idx = np.load(train_idx_path, allow_pickle=True)

# Select development set and reset index
dev_df = combined_df.iloc[train_idx].reset_index(drop=True)

## make 3 classed outcome
conditions = [dev_df['KPS_post-op_value'] == 0,
              dev_df['KPS_post-op_value'] < 70,
              dev_df['KPS_post-op_value'] >= 70]
values = [0,1,2]

dev_df['KPS_post_op_binned'] = np.select(conditions,values)
print(f'{dev_df["KPS_post_op_binned"].value_counts(normalize=True)*100}')

#### only keep relevant columns
dev_df_og = dev_df.copy()
dev_df.columns

subj_ids_kps = dev_df[['PRECOG_ID', 'KPS_post-op_value']]
drop_cols = ['PRECOG_ID', 'strata', 'source', 'KPS_post-op_value', 'survival_categories']

dev_df.drop(columns=drop_cols, inplace=True)
dev_df.columns

# %%
# Drop columns that are not features
combined_df = dev_df.copy()
cols_to_drop = ['event', 'time_to_event', 'KPS_post_op_binned']
combined_df = combined_df.drop(columns=cols_to_drop, errors='ignore')
# our main outcome is not present in the combined_df (KPS_post-op_binned)

# Define feature sets
feature_sets = {
    "clinical": ['Sex','Age at resection','Preoperative seizures','Preoperative KPS'],
    "combined": list(combined_df.columns),  # convert Index -> list
    'top3': ['Age at resection', 'Enhancing component volume', 'Preoperative KPS'],
    'top4': ['Age at resection', 'Enhancing component volume', 'Preoperative KPS', 'Cognitive networks overlap']
}

# Use the feature set
X = dev_df[feature_sets["top4"]]
y = dev_df['KPS_post_op_binned']
cv_outer = StratifiedKFold(n_splits=10)

#%%
# TRAIN THE MODEL

#ranges
param_space = {
    "max_depth": (3, 7),
    "learning_rate": (0.01, 0.2, "log-uniform"),
    "n_estimators": (100, 600),
    "reg_alpha": (0, 10),
    "reg_lambda": (0, 10),
    "min_child_weight": (1, 20),
    'class1_scale': [0.95, 1, 1.2, 1.5], 
    #"tree_method": ["hist", "approx", "auto"] # avoid this - not very CPU efficient - instead set to 'hist' when initializing
}

results_top4, shap_values, feature_names, classes = train_and_evaluate_ordinal_shap(
    X, y, cv_outer, param_space, n_iter=150, feature_set_name='top4'
)

#%%

metrics_df = compute_foldwise_metrics(results_top4, class_labels=classes)
print(metrics_df)

#%%
### the function below provides metrics aggregated across folds - only means
eval_results = evaluate_across_folds(results_top4, class_labels=[0,1,2])

#%%

plot_multiclass_roc_cv(results_top4, class_labels=classes)

#%% ### Investigate SHAPs per class

n_classifiers = len(results_top4[0]['shap_per_class'])
feature_names = X.columns

# Collect SHAP values and corresponding X_test per classifier
shap_per_classifier = [[] for _ in range(n_classifiers)]
X_per_classifier = [[] for _ in range(n_classifiers)]

for fold in results_top4:
    fold_shap = fold['shap_per_class']
    X_test_fold = fold['X_test']  # now directly stored
    
    for i, class_shap in enumerate(fold_shap):
        shap_per_classifier[i].append(class_shap)
        X_per_classifier[i].append(X_test_fold)

# Concatenate SHAP values and matching features per classifier
shap_per_classifier_concat = [np.vstack(shap_list) for shap_list in shap_per_classifier]
X_per_classifier_concat = [pd.concat(X_list, axis=0) for X_list in X_per_classifier]

# Plot SHAP summary per classifier
for i, (class_shap_all_folds, X_all_folds) in enumerate(zip(shap_per_classifier_concat, X_per_classifier_concat)):
    plt.figure(figsize=(10, 6))
    shap.summary_plot(class_shap_all_folds, X_all_folds, feature_names=feature_names, show=False)
    plt.title(f'SHAP summary for classifier threshold {i} (class > {i})')

#%% # GLOBAL SHAP values
# Aggregation across classifiers! and then select top 5 from these

mean_abs_shap_per_classifier = []

for shap_vals in shap_per_classifier_concat:
    # mean |SHAP| over samples
    mean_abs = np.mean(np.abs(shap_vals), axis=0)
    mean_abs_shap_per_classifier.append(mean_abs)

mean_abs_shap_per_classifier = np.vstack(mean_abs_shap_per_classifier)
aggregated_mean_abs_shap = mean_abs_shap_per_classifier.mean(axis=0)

shap_importance_df = (
    pd.DataFrame({
        "feature": feature_names,
        "mean_abs_shap": aggregated_mean_abs_shap
    })
    .sort_values("mean_abs_shap", ascending=False)
    .reset_index(drop=True)
)

shap_importance_df.head(5)

#%% ## investigate some hyperparameters

#print(results_list[0]['best_model'])
best_params_df = pd.DataFrame([r["best_params"] for r in results_top4])
print(best_params_df["class1_scale"])

#%% ## store results & training info ##

output_dir = '/folder/models/ordinal_reg_xgb/top4_23122025_same_split.pkl'
joblib.dump(results_top4, output_dir)

metadata = {}

metadata["param_space"] = param_space                         # already JSON serializable
metadata["feature_set_name"] = "top4"
metadata["features"] = list(X.columns)                        # convert Index → list
metadata["n_outer_folds"] = cv_outer.get_n_splits()
metadata["n_iter"] = 150
metadata["timestamp"] = datetime.utcnow().isoformat() + "Z"
metadata['scoring'] = 'roc_auc_ovr_weighted'

with open("/folder/models/ordinal_reg_xgb/top4_23122025_same_split.json", "w") as f:
    json.dump(metadata, f, indent=2)