#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 11:12:27 2026

@author: ekoderman
"""
from pathlib import Path
import joblib
from sklearn.metrics import f1_score, mean_absolute_error, roc_auc_score, cohen_kappa_score, matthews_corrcoef, balanced_accuracy_score
from utils import bootstrap_auc_multiclass, evaluate_across_folds
import shap
from utils import bootstrap_metrics_global

external_results_path = Path('/folder/models/ordinal_reg_xgb/final_model_held-out_results_09012026.pkl')
results_external = joblib.load(external_results_path)

#%%

qwk =  cohen_kappa_score(results_external['y_true'], results_external['y_pred'], weights='quadratic')
mcc = matthews_corrcoef(results_external['y_true'], results_external['y_pred'])

#%%

summary, roc_scores, pr_scores = bootstrap_auc_multiclass([results_external], n_iterations=500, plot_ci=True)

#%%
eval_results = evaluate_across_folds([results_external], class_labels=[0,1,2])

#%%

y_true_all = results_external['y_true'] 
y_pred_all = results_external['y_pred']

class_labels = [0,1,2]

ci_results = bootstrap_metrics_global(y_true_all, y_pred_all, class_labels=class_labels, n_iterations=500)
print(ci_results)
