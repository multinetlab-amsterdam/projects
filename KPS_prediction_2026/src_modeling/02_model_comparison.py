#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 11:16:33 2025

@author: ekoderman

status: final

REVIEW
reviewed by: Sebastien Dam
review date: January 2026

This script is used to compare results across models obtained in the 01_model_train_xgbordreg.py.
The best model is chosen based on visual inspection of these results and on consideration of clinical usability of the model.
THe best model is then tested on external test set in the next script (03_model-test_held-out_set.py)

"""
# Imports
import numpy as np
import pandas as pd
import joblib
from sklearn.metrics import f1_score, mean_absolute_error, roc_auc_score, cohen_kappa_score, matthews_corrcoef, balanced_accuracy_score
from sklearn.preprocessing import label_binarize
from collections import defaultdict
import matplotlib.pyplot as plt
from utils import compute_foldwise_metrics
import os

# Paths
clinical_mri = joblib.load('/folder/models/ordinal_reg_xgb/combined_23122025_same_split.pkl')
clinical = joblib.load('/folder/models/ordinal_reg_xgb/clinical_23122025_same_split.pkl')
top3 = joblib.load('/folder/models/ordinal_reg_xgb/top3_23122025_same_split.pkl')
top4 = joblib.load('/folder/models/ordinal_reg_xgb/top4_23122025_same_split.pkl')

#%%

#compare results across models
all_results = [clinical_mri, clinical, top3, top4]

# f1_scores[model_name][class_idx] -> list of F1 scores (across folds)
f1_scores = defaultdict(lambda: defaultdict(list))

for model in all_results:
    for res in model:
        model_name = res['feature_set_name']
        y_true = res['y_true']
        y_pred = res['y_pred']
    
        f1_per_class = f1_score(y_true, y_pred, average=None)
    
        for c, f1 in enumerate(f1_per_class):
            f1_scores[model_name][c].append(f1)

# aggregate across folds (mean F1 per class per model)
models = list(f1_scores.keys())
classes = sorted(next(iter(f1_scores.values())).keys())


# aggregate
f1_mean = {m: [] for m in models}
f1_std  = {m: [] for m in models}

for m in models:
    for c in classes:
        f1_mean[m].append(np.mean(f1_scores[m][c]))
        f1_std[m].append(np.std(f1_scores[m][c], ddof=1))  # sample SD

# ---- plotting ----
x = np.arange(len(models))
width = 0.25

colors = ['tab:blue', 'tab:green', 'tab:orange']
alpha_bar = 0.6
alpha_err = 0.8

fig, ax = plt.subplots(figsize=(8, 5))

class_labels = {
    0: 'KPS = 0',
    1: 'KPS < 70',
    2: 'KPS ≥ 70'
}

for i, c in enumerate(classes):
    heights = [f1_mean[m][i] for m in models]
    errors  = [f1_std[m][i] for m in models]

    bars = ax.bar(
        x + i * width,
        heights,
        width,
        color=colors[i],
        alpha=alpha_bar,
        edgecolor='black',
        linewidth=0.5,
        label=class_labels[c]
    )

    ax.errorbar(
        x + i * width,
        heights,
        yerr=errors,
        fmt='none',
        ecolor=colors[i],
        elinewidth=1.5,
        capsize=6,
        alpha=alpha_err
    )

ax.set_xticks(x + width)
ax.set_xticklabels(models, rotation=30, ha='right')
ax.set_ylabel('F1 score')
ax.set_xlabel('Model')
ax.set_ylim(0, 1)
ax.legend(frameon=False)
plt.tight_layout()
plt.show()


#%%

classes = [0,1,2]
metrics_df = compute_foldwise_metrics(top4, class_labels=classes)
print(metrics_df)

#%%

# Map each results list to a model name
model_mapping = {
    'clinical_mri': clinical_mri,
    'clinical': clinical,
    'top3': top3,
    'top4': top4
}
all_results_flat = []

for model_name, model_res in model_mapping.items():
    for r in model_res:  # each r is a dict
        r_copy = r.copy()
        r_copy['model'] = model_name
        all_results_flat.append(r_copy)

all_results = all_results_flat

# ---------------------------------------------------------
# Convert to a DataFrame of per-fold scores
# ---------------------------------------------------------
rows = []

# Determine all unique classes across datasets
all_classes = np.unique(np.concatenate([r['y_true'] for r in all_results]))

# Binarize y for ROC-AUC per class
for r in all_results:
    y_true = np.asarray(r['y_true'])
    y_pred = np.asarray(r['y_pred'])
    
    # One-hot encode true labels for ROC-AUC
    y_true_bin = label_binarize(y_true, classes=all_classes)
    
    # Compute per-class ROC-AUC if probabilities available
    if 'y_proba' in r:
        y_proba = np.asarray(r['y_proba'])
        per_class_rocauc = {}
        for i, cls in enumerate(all_classes):
            try:
                per_class_rocauc[f'roc_auc_class_{cls}'] = roc_auc_score(y_true_bin[:, i], y_proba[:, i])
            except ValueError:
                per_class_rocauc[f'roc_auc_class_{cls}'] = np.nan
    else:
        per_class_rocauc = {f'roc_auc_class_{cls}': np.nan for cls in all_classes}
    
    rows.append({
        'model': r['model'],
        'fold': r.get('fold', r.get('fold_idx', None)),
        'f1_weighted': f1_score(y_true, y_pred, average='weighted'),
        'bal_accuracy': balanced_accuracy_score(y_true, y_pred),
        'mae': mean_absolute_error(y_true, y_pred),
        'qwk': cohen_kappa_score(y_true, y_pred, weights='quadratic'),
        'mcc': matthews_corrcoef(y_true, y_pred),
        **per_class_rocauc
    })
df_scores = pd.DataFrame(rows)

#%%
# Compute mean of each metric across folds for each model
df_summary = df_scores.groupby('model').mean(numeric_only=True).reset_index()

# Optional: reorder columns if needed
cols = ['model'] + [c for c in df_summary.columns if c != 'model']
df_summary = df_summary[cols]

#%%

output_dir = '/folder/models/ordinal_reg_xgb/'
# Make sure the directory exists
os.makedirs(output_dir, exist_ok=True)

# Full path to the CSV file
csv_path = os.path.join(output_dir, 'model_selection_scores_across_folds.csv')
df_scores.to_csv(csv_path, index=False)

# Save the summary
summary_csv_path = os.path.join(output_dir, 'model_selection_summary.csv')
df_summary.to_csv(summary_csv_path, index=False)