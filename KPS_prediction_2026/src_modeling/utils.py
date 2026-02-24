#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 17:32:43 2025

@author: ekoderman
"""
from sklearn.base import BaseEstimator, ClassifierMixin
from xgbordinal import XGBOrdinal
import numpy as np
from collections import Counter
from sklearn.metrics import recall_score

# --- Wrapper to make XGBOrdinal compatible with BayesSearchCV ---
class XGBOrdinalWrapper(BaseEstimator, ClassifierMixin):
    
    def __init__(self, aggregation='weighted', norm=True, random_state=None, class1_scale=1.0, **kwargs):
        self.aggregation = aggregation
        self.norm = norm
        self.random_state = random_state
        self.class1_scale = class1_scale
        self.kwargs = kwargs
        self.model = None
        self.classes_ = [0,1,2]

    def fit(self, X, y, sample_weight=None):
        from collections import Counter
        
        self.model = XGBOrdinal(
            aggregation=self.aggregation,
            norm=self.norm,
            random_state=self.random_state,
            **self.kwargs
        )
        
        # Compute sample_weight inside the wrapper if not passed
        if sample_weight is None:
            class_counts = Counter(y)
            total = sum(class_counts.values())
            sample_weight = np.array([total / class_counts[v] for v in y])
        
        # Apply class1_scale to minority class
        sample_weight[y == 1] *= self.class1_scale
        
        # Normalize so mean weight = 1
        sample_weight = sample_weight / np.mean(sample_weight)
        
        # Initialize and train the XGBOrdinal model
        self.model = XGBOrdinal(
            aggregation=self.aggregation,
            norm=self.norm,
            random_state=self.random_state,
            **self.kwargs
        )
        
        # Fit the model
        self.model.fit(X, y, sample_weight=sample_weight)

        return self

    def predict(self, X):
        return self.model.predict(X)

    def predict_proba(self, X):
        return self.model.predict_proba(X)

    def get_params(self, deep=True):
        params = dict(self.kwargs)
        params.update({
            'aggregation': self.aggregation,
            'norm': self.norm,
            'random_state': self.random_state,
            'class1_scale': self.class1_scale
        })
        return params

    def set_params(self, **params):
        for k, v in params.items():
            if k in ['aggregation', 'norm', 'random_state', 'class1_scale']:
                setattr(self, k, v)
            else:
                self.kwargs[k] = v
        return self

def weighted_mae_sklearn(y_true, y_pred):
    counts = Counter(y_true)
    weights = np.array([1 / counts[c] for c in y_true])
    weights = weights / np.mean(weights)  # normalize to keep scale similar
    return np.mean(weights * np.abs(y_true - y_pred))

# Option 1: class 1 recall
def recall_class_1(y_true, y_pred):
    return recall_score(y_true, y_pred, labels=[1], average='macro')  # or 'binary' if only class 1

#recall_class_1_scorer = make_scorer(recall_class_1, greater_is_better=True)


def compute_foldwise_metrics(fold_results, class_labels=None):
    """
    Compute metrics per fold and aggregate mean ± std across folds.

    Parameters
    ----------
    fold_results : list of dict
        Each dict must contain keys: 'y_true', 'y_pred'.
        Optionally, 'y_proba' can be included for ROC-AUC (not used here).
    class_labels : list-like, optional
        Class names. Defaults to sorted unique values across folds.

    Returns
    -------
    metrics_summary : pd.DataFrame
        Rows are metrics (macro or per class), columns: 'mean', 'std'.
    """
    import numpy as np
    import pandas as pd
    from sklearn.metrics import f1_score, precision_score, recall_score, balanced_accuracy_score, confusion_matrix, classification_report
    from sklearn.metrics import matthews_corrcoef, cohen_kappa_score
    
    # Determine class labels
    if class_labels is None:
        class_labels = np.unique(np.concatenate([fold["y_true"] for fold in fold_results]))
    
    # Storage
    f1_macro_list = []
    f1_per_class_list = []
    precision_per_class_list = []
    recall_per_class_list = []
    bal_acc_list = []
    mcc_list = []
    qwk_list = []
    
    # Compute metrics per fold
    for fold in fold_results:
        y_true = fold["y_true"]
        y_pred = fold["y_pred"]
        
        f1_macro_list.append(f1_score(y_true, y_pred, average='macro'))
        f1_per_class_list.append(f1_score(y_true, y_pred, average=None))
        precision_per_class_list.append(precision_score(y_true, y_pred, average=None))
        recall_per_class_list.append(recall_score(y_true, y_pred, average=None))
        bal_acc_list.append(balanced_accuracy_score(y_true, y_pred))
        mcc_list.append(matthews_corrcoef(y_true, y_pred))
        qwk_list.append(cohen_kappa_score(y_true, y_pred, weights='quadratic'))
    
    # Aggregate mean ± std
    metrics_dict = {
        "f1_macro": (np.mean(f1_macro_list), np.std(f1_macro_list)),
        "balanced_accuracy": (np.mean(bal_acc_list), np.std(bal_acc_list)),
        "mcc": (np.mean(mcc_list), np.std(mcc_list)),
        "qwk": (np.mean(qwk_list), np.std(qwk_list))
    }
    
    # Per-class metrics
    for i, cls in enumerate(class_labels):
        metrics_dict[f"f1_class_{cls}"] = (np.mean([f[i] for f in f1_per_class_list]), 
                                           np.std([f[i] for f in f1_per_class_list]))
        metrics_dict[f"precision_class_{cls}"] = (np.mean([p[i] for p in precision_per_class_list]),
                                                  np.std([p[i] for p in precision_per_class_list]))
        metrics_dict[f"recall_class_{cls}"] = (np.mean([r[i] for r in recall_per_class_list]),
                                               np.std([r[i] for r in recall_per_class_list]))
    
    # Convert to DataFrame
    metrics_summary = pd.DataFrame(metrics_dict, index=["mean","std"]).T
    metrics_summary = metrics_summary[["mean","std"]]  # ensure column order
    
    # Confusion Matrix with concatenated y_pred and y_true across folds
    # Concatenate all fold results
    y_true_all = np.concatenate([fold["y_true"] for fold in fold_results])
    y_pred_all = np.concatenate([fold["y_pred"] for fold in fold_results])
    
    cm = confusion_matrix(y_true_all, y_pred_all)
    print("Confusion Matrix:\n", cm)
    print("\nClassification Report:\n", classification_report(y_true_all, y_pred_all))
    
    return metrics_summary

def evaluate_across_folds(fold_results, class_labels=None, plot=True):
    """
    Aggregate y_true, y_pred, y_pred_proba across folds and evaluate.
    
    Parameters
    ----------
    fold_results : list of dict
        Each dict should contain keys: 'y_true', 'y_pred', 'y_pred_proba'.
    class_labels : list-like, optional
        Class names/labels.
    plot : bool
        Whether to plot metrics.
        
    Returns
    -------
    results : dict
        Metrics returned by evaluate_multiclass.
    """
    import numpy as np
    # Concatenate all fold results
    y_true_all = np.concatenate([fold["y_true"] for fold in fold_results])
    y_pred_all = np.concatenate([fold["y_pred"] for fold in fold_results])
    #y_pred_proba_all = np.vstack([fold["y_pred_proba"] for fold in fold_results])
    y_pred_proba_all = np.vstack([fold["y_proba"] for fold in fold_results])
    
    # Call the evaluation function
    return evaluate_multiclass(
        y_test=y_true_all, 
        y_pred=y_pred_all, 
        y_pred_proba=y_pred_proba_all,
        class_labels=class_labels,
        plot=plot
    )

def evaluate_multiclass(y_test, y_pred, y_pred_proba, class_labels=None, plot=True):
    """
    Evaluate multiclass predictions with standard and ordinal metrics.
    
    Parameters
    ----------
    y_test : array-like, shape (n_samples,)
        True labels.
    y_pred : array-like, shape (n_samples,)
        Predicted labels.
    y_pred_proba : array-like, shape (n_samples, n_classes)
        Predicted probabilities.
    class_labels : list-like, optional
        Class names/labels for plotting. Defaults to sorted unique y_test.
    plot : bool, default=True
        Whether to plot confusion matrix, per-class metrics, and ROC curves.
    
    Returns
    -------
    results : dict
        Dictionary with metrics:
        - 'confusion_matrix'
        - 'classification_report'
        - 'f1_per_class', 'precision_per_class', 'recall_per_class'
        - 'quadratic_weighted_kappa'
        - 'mean_absolute_error'
        - 'balanced_accuracy'
        - 'mcc'
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.metrics import (
        confusion_matrix, classification_report, f1_score, precision_score,
        recall_score, roc_curve, auc, cohen_kappa_score, mean_absolute_error,
        balanced_accuracy_score, matthews_corrcoef
    )
    from itertools import combinations
    
    if class_labels is None:
        class_labels = np.unique(y_test)
    
    # Confusion Matrix
    cm = confusion_matrix(y_test, y_pred)
    print("Confusion Matrix:\n", cm)
    print("\nClassification Report:\n", classification_report(y_test, y_pred))
    print("\nPredicted probabilities (first 5 samples):\n", y_pred_proba[:5])
    
    # Per-class metrics
    f1_per_class = f1_score(y_test, y_pred, average=None)
    precision_per_class = precision_score(y_test, y_pred, average=None)
    recall_per_class = recall_score(y_test, y_pred, average=None)
    
    # Ordinal metrics
    qwk = cohen_kappa_score(y_test, y_pred, weights='quadratic')
    mae = mean_absolute_error(y_test, y_pred)
    
    # Global metrics
    bal_acc = balanced_accuracy_score(y_test, y_pred)
    mcc = matthews_corrcoef(y_test, y_pred)
    
    print("Quadratic Weighted Kappa (QWK):", round(qwk,3))
    print("Mean Absolute Error (Ordinal Distance):", round(mae,3))
    print("Balanced Accuracy:", round(bal_acc,3))
    print("Matthews Correlation Coefficient (MCC):", round(mcc,3))
    
    if plot:
        # Confusion matrix heatmap
        plt.figure(figsize=(6,5))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=class_labels,
                    yticklabels=class_labels)
        plt.xlabel("Predicted")
        plt.ylabel("True")
        plt.title("Confusion Matrix")
        plt.show()
        
        # Per-class metrics barplot
        metrics_df = pd.DataFrame({
            'Class': class_labels,
            'F1': f1_per_class,
            'Precision': precision_per_class,
            'Recall': recall_per_class
        }).melt(id_vars='Class', var_name='Metric', value_name='Score')
        
        plt.figure(figsize=(8,5))
        sns.barplot(data=metrics_df, x='Class', y='Score', hue='Metric', palette='viridis')
        plt.ylim(0,1)
        plt.title('Per-Class Metrics')
        plt.show()
    
    results = {
        "confusion_matrix": cm,
        "classification_report": classification_report(y_test, y_pred, output_dict=True),
        "f1_per_class": f1_per_class,
        "precision_per_class": precision_per_class,
        "recall_per_class": recall_per_class,
        "quadratic_weighted_kappa": qwk,
        "mean_absolute_error": mae,
        "balanced_accuracy": bal_acc,
        "mcc": mcc
    }
    
    return results

def plot_multiclass_roc_cv(fold_results, class_labels):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc

    """
    Computes and plots ROC mean ± std across folds for each class (one-vs-rest).
    
    fold_results: list of dicts, each containing:
        - 'y_true'
        - 'y_proba'
    class_labels: list or array of class names
    """

    n_classes = len(class_labels)
    mean_fpr = np.linspace(0, 1, 200)  # common FPR grid for interpolation

    plt.figure(figsize=(8, 6))

    for class_idx in range(n_classes):
        tprs = []
        aucs = []

        for fold_res in fold_results:
            y_true = fold_res["y_true"]
            y_proba = fold_res["y_proba"][:, class_idx]

            # binarize for one-vs-rest
            y_binary = (y_true == class_idx).astype(int)

            fpr, tpr, _ = roc_curve(y_binary, y_proba)
            roc_auc = auc(fpr, tpr)

            # interpolate TPR
            tpr_interp = np.interp(mean_fpr, fpr, tpr)
            tpr_interp[0] = 0.0

            tprs.append(tpr_interp)
            aucs.append(roc_auc)

        # aggregate
        mean_tpr = np.mean(tprs, axis=0)
        std_tpr = np.std(tprs, axis=0)
        mean_auc = auc(mean_fpr, mean_tpr)

        # plotting
        label = f"Class {class_labels[class_idx]} (AUC={mean_auc:.2f} ± {np.std(aucs):.2f})"
        plt.plot(mean_fpr, mean_tpr, label=label)

        # SD shading
        tpr_upper = np.minimum(mean_tpr + std_tpr, 1)
        tpr_lower = np.maximum(mean_tpr - std_tpr, 0)

        plt.fill_between(mean_fpr, tpr_lower, tpr_upper, alpha=0.2)

    plt.plot([0, 1], [0, 1], 'k--', lw=1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Mean ROC Curves Across Folds (One-vs-Rest)")
    plt.legend()
    plt.show()


def evaluate_post_calibration_mcp(fold_results, class_labels=None):
    """
    Visualize model performance after Platt/Isotonic calibration and MCP.
    
    Parameters
    ----------
    fold_results : list of dict
        Each dict should contain keys:
        'y_true', 'y_pred', 'y_proba' (calibrated), 'prediction_set' (MCP)
    class_labels : list-like, optional
        Names of classes. Defaults to sorted unique labels.
        
    Returns
    -------
    results : dict
        Coverage, accuracy, per-class inclusion rates, and global metrics.
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from collections import Counter
    
    # Aggregate across folds
    y_true_all = np.concatenate([f['y_true'] for f in fold_results])
    y_pred_all = np.concatenate([f['y_pred'] for f in fold_results])
    y_proba_all = np.vstack([f['y_proba'] for f in fold_results])
    pred_sets_all = np.concatenate([f['prediction_set'] for f in fold_results])
    
    if class_labels is None:
        class_labels = np.unique(y_true_all)
    
    n_classes = len(class_labels)
    n_samples = len(y_true_all)
    
    # Compute coverage: fraction of samples where true class is in MCP prediction set
    coverage = np.mean([y_true_all[i] in pred_sets_all[i] for i in range(n_samples)])
    print(f"Global MCP coverage: {coverage:.3f} (should be ~ 1-alpha)")
    
    # Per-class coverage
    per_class_coverage = {}
    for cls in class_labels:
        idx = np.where(y_true_all == cls)[0]
        per_class_coverage[cls] = np.mean([y_true_all[i] in pred_sets_all[i] for i in idx])
    print("Per-class MCP coverage:", per_class_coverage)
    
    # Plot predicted probability distributions before and after calibration
    plt.figure(figsize=(10,6))
    for cls in range(n_classes):
        sns.histplot(y_proba_all[:, cls], bins=20, label=f'Class {cls}', alpha=0.5)
    plt.title("Calibrated predicted probabilities per class")
    plt.xlabel("Probability")
    plt.ylabel("Count")
    plt.legend()
    plt.show()
    
    # Plot MCP prediction set size distribution
    set_sizes = [len(s) for s in pred_sets_all]
    plt.figure(figsize=(6,4))
    sns.histplot(set_sizes, bins=np.arange(1, n_classes+2)-0.5, discrete=True)
    plt.xlabel("Size of MCP prediction set")
    plt.ylabel("Number of samples")
    plt.title("Distribution of MCP prediction set sizes")
    plt.show()
    
    # Confusion matrix for point predictions
    from sklearn.metrics import confusion_matrix
    cm = confusion_matrix(y_true_all, y_pred_all, labels=class_labels)
    plt.figure(figsize=(6,5))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=class_labels, yticklabels=class_labels)
    plt.xlabel("Predicted (argmax)")
    plt.ylabel("True class")
    plt.title("Confusion Matrix after Calibration")
    plt.show()
    
    results = {
        "global_coverage": coverage,
        "per_class_coverage": per_class_coverage,
        "set_sizes": set_sizes,
        "y_true": y_true_all,
        "y_pred": y_pred_all,
        "y_proba": y_proba_all
    }
    
    return results



def plot_prediction_set_sizes(prediction_sets_df):
    """
    Plot distribution of prediction set sizes (number of classes in each set).
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    
    set_sizes = prediction_sets_df['prediction_set'].apply(len)
    plt.figure(figsize=(6,4))
    sns.histplot(set_sizes, bins=np.arange(set_sizes.min(), set_sizes.max()+2)-0.5, 
                 discrete=True, kde=False, color='skyblue')
    plt.xlabel("Prediction Set Size")
    plt.ylabel("Number of Patients")
    plt.title("Distribution of Mondrian Conformal Prediction Set Sizes")
    plt.xticks(range(set_sizes.min(), set_sizes.max()+1))
    plt.show()
    print(f'percentage distribution of set sizes: {set_sizes.value_counts(normalize=True)*100}')


def plot_coverage_per_class(prediction_sets_df):
    """
    Compute and plot coverage per true class: fraction of samples where true class is in prediction set.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    coverage_dict = {}
    classes = sorted(prediction_sets_df['y_true'].unique())
    
    for cls in classes:
        cls_df = prediction_sets_df[prediction_sets_df['y_true'] == cls]
        coverage = cls_df.apply(lambda row: row['y_true'] in row['prediction_set'], axis=1).mean()
        coverage_dict[cls] = coverage
    
    plt.figure(figsize=(6,4))
    sns.barplot(x=list(coverage_dict.keys()), y=list(coverage_dict.values()), palette="Set2")
    plt.ylim(0,1)
    plt.xlabel("True Class")
    plt.ylabel("Coverage (fraction of samples where true class in set)")
    plt.title("Mondrian Conformal Prediction Coverage by Class")
    plt.show()
    

    
def plot_patient_prediction(patient_row, class_labels=None, feature_names=None, feature_values=None):
    """
    Visualize predictions and conformal prediction set for a single patient.
    
    Parameters
    ----------
    patient_row : pd.Series
        Row from prediction_sets_df containing:
        - 'y_true': true class
        - 'y_proba': array of predicted probabilities
        - 'prediction_set': list or set of classes included in conformal prediction
    class_labels : list, optional
        List of class labels corresponding to probabilities. If None, uses range(len(y_proba))
    feature_names : list, optional
        List of feature names
    feature_values : list, optional
        List of feature values corresponding to the patient
    """
    import matplotlib.pyplot as plt
    
    y_true = patient_row['y_true']
    y_proba = patient_row['y_proba']
    prediction_set = patient_row['prediction_set']
    
    if class_labels is None:
        class_labels = list(range(len(y_proba)))
    
    plt.figure(figsize=(8,5))
    bars = plt.bar(class_labels, y_proba, color='lightgray', alpha=0.7)
    
    # Highlight the classes in the conformal prediction set
    for i, label in enumerate(class_labels):
        if label in prediction_set:
            bars[i].set_color('dodgerblue')
            bars[i].set_alpha(0.9)
    
    # Mark the true class
    if y_true in class_labels:
        plt.scatter([y_true], [y_proba[y_true]], color='red', s=120, marker='*', label='True class')
    
    # Annotate feature names and their values
    if feature_names and feature_values:
        for i, (feature, value) in enumerate(zip(feature_names, feature_values)):
            # Adding the feature name and its value on the plot
            plt.text(i, 1.00, f'{feature}: {value}', ha='center', va='bottom', fontsize=10)#, rotation=45)
    
    plt.xticks(class_labels)
    plt.ylabel('Predicted probability')
    plt.title('Per-patient prediction with conformal set')
    plt.ylim(0, 1.05)
    plt.legend()
    plt.show()

def summary_mondrian_conformal(prediction_sets_df):
    """
    Generate a summary report for Mondrian Conformal Predictions:
    - Average prediction set size
    - Coverage overall
    - Coverage per class
    """
    set_sizes = prediction_sets_df['prediction_set'].apply(len)
    overall_coverage = prediction_sets_df.apply(
        lambda row: row['y_true'] in row['prediction_set'], axis=1).mean()
    
    classes = sorted(prediction_sets_df['y_true'].unique())
    coverage_per_class = {}
    for cls in classes:
        cls_df = prediction_sets_df[prediction_sets_df['y_true'] == cls]
        coverage_per_class[cls] = cls_df.apply(lambda row: row['y_true'] in row['prediction_set'], axis=1).mean()
    
    print(f"Overall coverage: {overall_coverage:.3f}")
    print(f"Average prediction set size: {set_sizes.mean():.2f}")
    print("Coverage per class:")
    for cls, cov in coverage_per_class.items():
        print(f"  Class {cls}: {cov:.3f}")
    
    # Optional plots
    plot_prediction_set_sizes(prediction_sets_df)
    plot_coverage_per_class(prediction_sets_df)
    
    


def mondrian_conformal_prediction(results_list, alpha=0.1):
    """
    Implements Mondrian conformal prediction for ordinal classification using your trained results_list.
    
    Parameters
    ----------
    results_list : list of dicts
        Output from train_and_evaluate_ordinal_shap (per-fold results). Must contain 'y_true', 'y_proba', 'X_test'.
    alpha : float
        Miscoverage rate. 0.1 -> 90% confidence prediction sets.
    
    Returns
    -------
    all_prediction_sets : pd.DataFrame
        DataFrame with index matching concatenated test data, columns:
        - 'y_true': true labels
        - 'prediction_set': list of classes included in prediction set
        - 'y_pred': point prediction (argmax of probability)
        - 'y_proba': predicted probabilities
    """
    import numpy as np
    import pandas as pd
    from collections import defaultdict
    
    # Step 1: Concatenate all folds
    all_y_true = np.concatenate([fold['y_true'] for fold in results_list])
    all_y_proba = np.vstack([fold['y_proba'] for fold in results_list])
    all_X_test = pd.concat([fold['X_test'] for fold in results_list], axis=0).reset_index(drop=True)
    
    n_samples, n_classes = all_y_proba.shape
    classes = np.arange(n_classes)
    
    # Step 2: Compute Mondrian nonconformity scores
    # Using 1 - probability of true class as nonconformity
    nonconformity_scores = 1 - all_y_proba[np.arange(n_samples), all_y_true]
    
    # Step 3: Stratify by class (Mondrian)
    class_indices = defaultdict(list)
    for idx, y in enumerate(all_y_true):
        class_indices[y].append(idx)
    
    class_quantiles = {}
    for cls, idxs in class_indices.items():
        scores_cls = nonconformity_scores[idxs]
        class_quantiles[cls] = np.quantile(scores_cls, 1 - alpha)  # (1-alpha) quantile per class
    
    # Step 4: Compute prediction sets
    prediction_sets = []
    point_preds = []
    
    for i in range(n_samples):
        y_proba_i = all_y_proba[i]
        y_true_i = all_y_true[i]
        
        # Include class k in prediction set if 1 - p_k <= threshold for that class
        included_classes = [k for k in classes if 1 - y_proba_i[k] <= class_quantiles[k]]
        prediction_sets.append(included_classes)
        
        # Optional: point prediction as argmax
        point_preds.append(np.argmax(y_proba_i))
    
    # Step 5: Combine into DataFrame
    all_prediction_sets = all_X_test.copy()
    all_prediction_sets['y_true'] = all_y_true
    all_prediction_sets['prediction_set'] = prediction_sets
    all_prediction_sets['y_pred'] = point_preds
    all_prediction_sets['y_proba'] = list(all_y_proba)  # store full probability vector
    for k in classes:
        all_prediction_sets[f'prob_class_{k}'] = all_y_proba[:, k]
    
    return all_prediction_sets

######## THE FUNCTION BELOW WORKS WHEN THERE IS NO STORING/RELOADING OF THE MODELS:

def calculate_mean_shap_signal(results, feature_set: str):
    """
    Calculate and return the mean SHAP signal for a given set of XGBoost folds.
    Compatible with both array and shap.Explanation outputs.
    """
    import numpy as np
    import pandas as pd
    import shap

    xgb_folds = results[feature_set]
    all_shap_values = []
    all_X = []

    for fold in xgb_folds:
        model = fold["best_model"]
        X = fold["X_test"]

        explainer = shap.TreeExplainer(model)
        shap_vals = explainer.shap_values(X)

        # Ensure output is a shap.Explanation
        if isinstance(shap_vals, np.ndarray):
            shap_vals = shap.Explanation(
                values=shap_vals,
                base_values=explainer.expected_value,
                data=X.values,
                feature_names=X.columns
            )

        all_shap_values.append(shap_vals)
        all_X.append(X)

    # Combine all folds
    combined_values = np.vstack([sv.values for sv in all_shap_values])
    combined_data = np.vstack([sv.data for sv in all_shap_values])
    combined_X = pd.concat(all_X, axis=0).reset_index(drop=True)

    combined_shap = shap.Explanation(
        values=combined_values,
        base_values=np.mean([sv.base_values for sv in all_shap_values], axis=0),
        data=combined_data,
        feature_names=all_shap_values[0].feature_names
    )

    mean_abs_shap = np.abs(combined_shap.values).mean(axis=0)
    mean_signal_shap = np.mean(mean_abs_shap)
    feature_importance = dict(zip(combined_X.columns, mean_abs_shap))

    shap.summary_plot(combined_shap, combined_X, show=False)

    print(f"Mean SHAP signal across all features: {mean_signal_shap:.4f}")
    return mean_signal_shap, feature_importance, combined_shap, combined_X    

##### THIS FUNCITON BELOW WORKS WHEN THE MODELS GET STORED/RELOADED

# def calculate_mean_shap_signal(results, feature_set: str):
#     """
#     Compute mean SHAP signal across XGBoost LOOCV/MCP folds.
#     """
#     import numpy as np
#     import pandas as pd
#     import shap

#     xgb_folds = results[feature_set]
#     all_values = []
#     all_base_values = []
#     all_data = []
#     all_X = []

#     for fold in xgb_folds:
#         model = fold['best_model']
#         X = fold['X_test']

#         # Generic Explainer using predict_proba
#         #explainer = shap.Explainer(model.predict_proba, X, output_names=['0','1'])
#         #shap_values = explainer(X)
#         #model = fold['best_model'].get_booster()  # gets xgb.Booster object

#         # if base_score string problem persists
#         booster = model.get_booster()  

#         explainer = shap.TreeExplainer(booster)
#         shap_values = explainer.shap_values(X) 

#         # Binary classification: take class 1
#         all_values.append(shap_values.values[:, :, 1])
#         # base_values shape: (1, n_classes) or scalar
#         if isinstance(shap_values.base_values, np.ndarray) and shap_values.base_values.ndim > 0:
#             all_base_values.append(shap_values.base_values[:, 1])
#         else:
#             all_base_values.append(np.array([shap_values.base_values]*X.shape[0]))
#         all_data.append(shap_values.data)
#         all_X.append(X)

#     # Combine all folds
#     combined_values = np.vstack(all_values)
#     combined_base_values = np.hstack(all_base_values)
#     combined_data = np.vstack(all_data)
#     combined_X = pd.concat(all_X, axis=0).reset_index(drop=True)

#     # Create combined Explanation object
#     combined_shap = shap.Explanation(
#         values=combined_values,
#         base_values=combined_base_values,
#         data=combined_data,
#         feature_names=all_X[0].columns
#     )

#     # Compute mean absolute SHAP signal
#     mean_abs_shap = np.abs(combined_shap.values).mean(axis=0)
#     mean_signal_shap = np.mean(mean_abs_shap)

#     feature_importance = dict(zip(combined_X.columns, mean_abs_shap))

#     # Summary plot
#     shap.summary_plot(combined_shap, combined_X, show=False)

#     print(f"Mean SHAP signal across all features: {mean_signal_shap:.4f}")
#     return mean_signal_shap, feature_importance, combined_shap, combined_X



def train_and_evaluate_all_algorithms(X, y, cv_outer, model, feature_set: str):
    import time
    import numpy as np
    from sklearn.base import clone
    from sklearn.preprocessing import LabelEncoder
    from sklearn.impute import SimpleImputer

    print(f"Training and evaluating {feature_set} ({model.__class__.__name__})...")
    results_list = []

    # Encode y if needed (some classifiers need integers)
    if y.dtype == 'O' or y.dtype.name == 'category':
        le = LabelEncoder()
        y = le.fit_transform(y)

    # Identify numeric columns
    numeric_cols = X.select_dtypes(include=[np.number]).columns
    binary_cols = [c for c in numeric_cols if X[c].dropna().nunique() == 2]
    continuous_cols = [c for c in numeric_cols if c not in binary_cols]

    # Imputers
    imputer_binary = SimpleImputer(strategy="most_frequent")
    imputer_cont = SimpleImputer(strategy="median")

    for fold, (train_idx, test_idx) in enumerate(cv_outer.split(X, y), 1):
        X_train, X_test = X.iloc[train_idx].copy(), X.iloc[test_idx].copy()
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # Impute
        if binary_cols:
            X_train[binary_cols] = imputer_binary.fit_transform(X_train[binary_cols])
            X_test[binary_cols] = imputer_binary.transform(X_test[binary_cols])
        if continuous_cols:
            X_train[continuous_cols] = imputer_cont.fit_transform(X_train[continuous_cols])
            X_test[continuous_cols] = imputer_cont.transform(X_test[continuous_cols])

        start = time.time()
        clf = clone(model)
        clf.fit(X_train, y_train)
        elapsed = time.time() - start

        # Prediction handling
        if hasattr(clf, "predict_proba"):
            y_pred_proba = clf.predict_proba(X_test)[:, 1]
        elif hasattr(clf, "decision_function"):
            decision = clf.decision_function(X_test)
            y_pred_proba = (decision - decision.min()) / (decision.max() - decision.min())
        else:
            y_pred_proba = np.zeros_like(y_test, dtype=float)

        y_pred = clf.predict(X_test)

        results_list.append({
            "best_model": clf,
            "best_score": np.nan,  # no inner CV score
            "y_true": y_test,
            "y_pred": y_pred,
            "y_pred_proba": y_pred_proba,
            "fit_time": elapsed,
            "X_test": X_test
        })

        print(f"Fold {fold}: done in {elapsed:.2f}s")

    return results_list

def plot_cv_results_all(models_results: dict):
    """
    Plot cross-validated ROC curves and performance metrics for multiple models.
    
    models_results: dict of {model_name: folds_list}, where folds_list is output of train_and_evaluate
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.metrics import roc_curve, auc, matthews_corrcoef, f1_score, balanced_accuracy_score
    import pandas as pd

    # === Metrics ===
    all_metrics = []

    for model_name, folds in models_results.items():
        for fold in folds:
            y_true = fold["y_true"]
            y_pred = fold["y_pred"]
            all_metrics.append({
                "Model": model_name,
                "MCC": matthews_corrcoef(y_true, y_pred),
                "F1": f1_score(y_true, y_pred),
                "Balanced Accuracy": balanced_accuracy_score(y_true, y_pred)
            })

    df_metrics = pd.DataFrame(all_metrics)
    df_long = df_metrics.melt(id_vars="Model", var_name="Metric", value_name="Score")

    plt.figure(figsize=(10,6))
    sns.barplot(data=df_long, x="Metric", y="Score", hue="Model", palette="tab10", ci="sd")
    plt.ylim(0,1)
    plt.title("Cross-Validated Performance Metrics Across Models")
    plt.legend(title="Model")
    plt.tight_layout()
    plt.show()

    # === ROC Curves ===
    plt.figure(figsize=(8,6))
    for model_name, folds in models_results.items():
        tprs = []
        aucs = []
        mean_fpr = np.linspace(0,1,100)

        for fold in folds:
            y_true = fold["y_true"]
            y_proba = fold["y_pred_proba"]
            fpr, tpr, _ = roc_curve(y_true, y_proba)
            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(auc(fpr,tpr))
        
        mean_tpr = np.mean(tprs, axis=0)
        mean_auc = np.mean(aucs)
        std_auc = np.std(aucs)

        plt.plot(mean_fpr, mean_tpr, lw=2, label=f"{model_name} (AUC = {mean_auc:.2f} ± {std_auc:.2f})")

    plt.plot([0,1],[0,1], linestyle='--', color='gray')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curves Across Models")
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_cv_results_across_tiers(folds_dict, mcp=False):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc, f1_score, matthews_corrcoef, balanced_accuracy_score
    import pandas as pd
    import seaborn as sns

    model_colors = sns.color_palette("tab10", n_colors=len(folds_dict))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax_roc, ax_bar = axes

    all_metrics = []

    for color, (model_name, folds) in zip(model_colors, folds_dict.items()):
        tprs, aucs = [], []
        mean_fpr = np.linspace(0, 1, 100)

        for fold in folds:
            y_true = fold["y_true"]
            y_proba_full = fold["y_pred_proba"]
            
            # Always take true model probabilities for ROC
            y_proba = y_proba_full[:, 1]

            # Compute ROC and interpolate for averaging
            fpr, tpr, _ = roc_curve(y_true, y_proba)
            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(auc(fpr, tpr))

            # Use MCP or normal predictions for metrics
            y_pred = fold["y_pred"] if mcp else fold.get("y_pred_original", fold["y_pred"])

            all_metrics.append({
                "Model": model_name,
                "MCC": matthews_corrcoef(y_true, y_pred),
                "F1": f1_score(y_true, y_pred),
                "Balanced Accuracy": balanced_accuracy_score(y_true, y_pred)
            })

        mean_tpr = np.mean(tprs, axis=0)
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)

        ax_roc.plot(mean_fpr, mean_tpr, color=color, lw=2,
                    label=f"{model_name} (AUC = {mean_auc:.2f} ± {std_auc:.2f})")
        ax_roc.fill_between(mean_fpr,
                            np.maximum(0, mean_tpr - np.std(tprs, axis=0)),
                            np.minimum(1, mean_tpr + np.std(tprs, axis=0)),
                            color=color, alpha=0.15)

    # ---- ROC formatting ----
    ax_roc.plot([0, 1], [0, 1], linestyle="--", color="gray")
    ax_roc.set_xlabel("False Positive Rate")
    ax_roc.set_ylabel("True Positive Rate")
    ax_roc.set_title("ROC Curves Across Models (Model Probabilities)")
    ax_roc.legend(loc="lower right")
    ax_roc.grid(True)

    # ---- Metrics barplot ----
    metrics_df = pd.DataFrame(all_metrics)
    metrics_long = metrics_df.melt(id_vars="Model", var_name="Metric", value_name="Score")
    sns.barplot(data=metrics_long, x="Metric", y="Score", hue="Model",
                palette=model_colors, errorbar="sd", ax=ax_bar)
    ax_bar.set_ylim(0, 1)
    ax_bar.set_title("Cross-Validated Metrics (MCP vs Base if applicable)")
    ax_bar.grid(True)

    fig.tight_layout()
    plt.show()


    
def plot_loocv_results(results_list, model_name="Model", mcp=False):
    """
    Plot metrics, confusion matrix, and ROC curve for LOOCV results.
    
    Parameters
    ----------
    results_list : list of dicts
        Each dict contains 'y_true', 'y_pred', 'y_pred_proba', etc.
    model_name : str
        Name of the model (used in plot titles)
    mcp : bool
        If True, treats 'y_pred_proba' as Mondrian Conformal Prediction output (2D)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.metrics import (f1_score, balanced_accuracy_score,
                                 matthews_corrcoef, confusion_matrix,
                                 roc_curve, auc)
    
    # Aggregate predictions and labels
    y_true_all = np.concatenate([r["y_true"] for r in results_list])
    y_pred_all = np.concatenate([r["y_pred"] for r in results_list])
    
    # Probabilities for ROC curve
    if mcp:
        # MCP: y_pred_proba is (n_samples, n_classes)
        y_proba_all = np.concatenate([np.array(r["y_pred_proba"])[:, 1] for r in results_list])
    else:
        # Regular thresholded probabilities
        y_proba_all = np.array([r["y_pred_proba"][1] if len(r["y_pred_proba"]) > 1 
                                else r["y_pred_proba"][0] for r in results_list])
    
    # Metrics
    bal_acc = balanced_accuracy_score(y_true_all, y_pred_all)
    f1 = f1_score(y_true_all, y_pred_all, average='macro')
    mcc = matthews_corrcoef(y_true_all, y_pred_all)
    
    print(f"{model_name} LOOCV Metrics:")
    print("Balanced Accuracy:", bal_acc)
    print("F1 Score (macro):", f1)
    print("MCC:", mcc)
    
    # Confusion matrix
    cm = confusion_matrix(y_true_all, y_pred_all)
    plt.figure(figsize=(6,5))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=np.unique(y_true_all),
                yticklabels=np.unique(y_true_all))
    plt.xlabel("Predicted")
    plt.ylabel("True")
    plt.title(f"{model_name} LOOCV Confusion Matrix")
    plt.show()
    
    # ROC curve (binary)
    fpr, tpr, thresholds = roc_curve(y_true_all, y_proba_all)
    roc_auc = auc(fpr, tpr)
    
    plt.figure()
    plt.plot(fpr, tpr, color='blue', lw=2, label=f"ROC (AUC = {roc_auc:.2f})")
    plt.plot([0,1],[0,1], color='gray', linestyle='--')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"{model_name} LOOCV ROC Curve")
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.show()

def plot_cv_results(folds, model_name):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import (
        roc_curve, auc, f1_score, matthews_corrcoef,
        balanced_accuracy_score, confusion_matrix, ConfusionMatrixDisplay
    )
    import pandas as pd
    import seaborn as sns

    # ROC curve
    tprs, aucs = [], []
    mean_fpr = np.linspace(0, 1, 100)

    for fold in folds:
        y_true = fold["y_true"]
        y_proba = fold["y_pred_proba"]
        fpr, tpr, _ = roc_curve(y_true, y_proba)
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(auc(fpr, tpr))

    mean_tpr = np.mean(tprs, axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)

    plt.figure()
    plt.plot(mean_fpr, mean_tpr, color='blue', lw=2,
             label=f"Mean ROC (AUC = {mean_auc:.2f} ± {std_auc:.2f})")
    plt.fill_between(mean_fpr,
                     np.maximum(0, mean_tpr - np.std(tprs, axis=0)),
                     np.minimum(1, mean_tpr + np.std(tprs, axis=0)),
                     color='blue', alpha=0.2)
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve Across Folds: {model_name}")
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.tight_layout()

    # Metrics histogram
    fold_metrics = []
    all_conf_mats = []

    for i, fold in enumerate(folds):
        y_true_fold = fold["y_true"]
        y_pred_fold = fold["y_pred"]

        fold_metrics.append({
            "MCC": matthews_corrcoef(y_true_fold, y_pred_fold),
            "F1": f1_score(y_true_fold, y_pred_fold),
            "Balanced Accuracy": balanced_accuracy_score(y_true_fold, y_pred_fold)
        })

        #confusion matrix for each fold
        cm = confusion_matrix(y_true_fold, y_pred_fold)
        all_conf_mats.append(cm)

        #disp = ConfusionMatrixDisplay(confusion_matrix=cm)
        #disp.plot(cmap='Blues', values_format='d')
        #plt.title(f"Confusion Matrix - Fold {i+1}: {model_name}")
        #plt.tight_layout()

    # Average confusion matrix
    mean_cm = np.mean(all_conf_mats, axis=0)
    disp = ConfusionMatrixDisplay(confusion_matrix=mean_cm)
    disp.plot(cmap='Blues')
    plt.title(f"Mean Confusion Matrix Across Folds: {model_name}")
    plt.tight_layout()

    # Aggregate metrics bar plot
    metrics_df = pd.DataFrame(fold_metrics)
    metrics_long = metrics_df.melt(var_name="Metric", value_name="Score")

    plt.figure(figsize=(8, 5))
    sns.barplot(data=metrics_long, x="Metric", y="Score", palette="coolwarm", errorbar="sd")
    plt.ylim(0, 1)
    plt.title(f"Cross-Validated Performance Metrics (mean ± SD): {model_name}")
    plt.grid()
    plt.tight_layout()

def plot_cv_results_mcp(folds, model_name):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import (
        roc_curve, auc, f1_score, matthews_corrcoef,
        balanced_accuracy_score, confusion_matrix, ConfusionMatrixDisplay
    )
    import pandas as pd
    import seaborn as sns

    # ROC curve
    tprs, aucs = [], []
    mean_fpr = np.linspace(0, 1, 100)

    for fold in folds:
        y_true = fold["y_true"]
        y_proba = fold["y_pred_proba"]

        # Always take true model probabilities for ROC
        y_proba_pos = y_proba[:, 1]

        # Handle MCP output: y_proba could be (n_samples, n_classes)
        # if y_proba.ndim == 2 and y_proba.shape[1] > 1:
        #     # Take probability of positive class
        #     y_proba_pos = y_proba[:, 1]
        # else:
        #     y_proba_pos = np.ravel(y_proba)  # flatten if needed

        fpr, tpr, _ = roc_curve(y_true, y_proba_pos)
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(auc(fpr, tpr))

    mean_tpr = np.mean(tprs, axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)

    plt.figure()
    plt.plot(mean_fpr, mean_tpr, color='blue', lw=2,
             label=f"Mean ROC (AUC = {mean_auc:.2f} ± {std_auc:.2f})")
    plt.fill_between(mean_fpr,
                     np.maximum(0, mean_tpr - np.std(tprs, axis=0)),
                     np.minimum(1, mean_tpr + np.std(tprs, axis=0)),
                     color='blue', alpha=0.2)
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve Across Folds: {model_name}")
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.tight_layout()

    # Metrics histogram and confusion matrices
    fold_metrics = []
    all_conf_mats = []

    for i, fold in enumerate(folds):
        y_true_fold = fold["y_true"]
        y_pred_fold = fold["y_pred"]
        #y_pred_fold = fold["y_pred_original"]

        fold_metrics.append({
            "MCC": matthews_corrcoef(y_true_fold, y_pred_fold),
            "F1": f1_score(y_true_fold, y_pred_fold),
            "Balanced Accuracy": balanced_accuracy_score(y_true_fold, y_pred_fold)
        })

        # Confusion matrix for each fold
        cm = confusion_matrix(y_true_fold, y_pred_fold)
        all_conf_mats.append(cm)

        #disp = ConfusionMatrixDisplay(confusion_matrix=cm)
        #disp.plot(cmap='Blues', values_format='d')
        #plt.title(f"Confusion Matrix - Fold {i+1}: {model_name}")
        #plt.tight_layout()

    # Average confusion matrix across folds
    mean_cm = np.mean(all_conf_mats, axis=0)
    disp = ConfusionMatrixDisplay(confusion_matrix=mean_cm)
    disp.plot(cmap='Blues')
    plt.title(f"Mean Confusion Matrix Across Folds: {model_name}")
    plt.tight_layout()

    # Aggregate metrics bar plot
    metrics_df = pd.DataFrame(fold_metrics)
    metrics_long = metrics_df.melt(var_name="Metric", value_name="Score")

    plt.figure(figsize=(8, 5))
    sns.barplot(data=metrics_long, x="Metric", y="Score", palette="coolwarm", errorbar="sd")
    plt.ylim(0, 1)
    plt.title(f"Cross-Validated Performance Metrics (mean ± SD): {model_name}")
    plt.grid()
    plt.tight_layout()
    
    

def compute_nonconformity_scores(y_true, y_proba, inverse_class_weight):
    """
    Compute nonconformity scores for ordinal MCP.

    Parameters
    ----------
    y_true : pd.Series or np.array
        True class labels
    y_proba : np.array
        Predicted probabilities, shape (n_samples, n_classes)
    class_penalty : dict or None
        Optional penalty per class, e.g., {0: 0.5, 1: 0.1, 2: 0.5}

    Returns
    -------
    scores : np.array
        Nonconformity scores, shape (n_samples, n_classes)
    """
    import numpy as np

    n_samples, n_classes = y_proba.shape
    scores = np.zeros_like(y_proba)

    for i in range(n_samples):
        for c in range(n_classes):
            # Ordinal distance + probability term + class penalty
            # scores[i, c] = (abs(c - y_true.iloc[i]) / (n_classes - 1)
            #                 + (1 - y_proba[i, c])
            #                 + class_penalty.get(c, 0.0))
            scores[i, c] = (1 - y_proba[i, c]) * inverse_class_weight[c]
    return scores


def mcp_predict(y_proba_test, class_nonconf_scores, class_penalty=None, alpha_per_class=None):
    """
    Compute MCP predictions for test set using precomputed class-conditional scores.
    Class penalties are consistently applied.

    Parameters
    ----------
    y_proba_test : np.array
        Test predicted probabilities, shape (n_samples, n_classes)
    class_nonconf_scores : dict
        Precomputed training nonconformity scores per class {c: np.array(scores_c)}
    class_penalty : dict or None
        Penalty per class
    alpha_per_class : dict
        Per-class alpha thresholds for coverage (lower alpha = stricter)

    Returns
    -------
    y_pred_mcp : np.array
        MCP-predicted classes
    p_values_list : list of dicts
        Per-class p-values for each test sample
    """
    import numpy as np

    n_samples, n_classes = y_proba_test.shape
    y_pred_mcp = []
    p_values_list = []

    if class_penalty is None:
        class_penalty = {c: 0.0 for c in range(n_classes)}
    if alpha_per_class is None:
        alpha_per_class = {c: 0.1 for c in range(n_classes)}

    for i in range(n_samples):
        # Compute test nonconformity score with class penalties
        distance = np.abs(np.arange(n_classes) - np.argmax(y_proba_test[i])) / (n_classes - 1)
        prob_term = 1 - y_proba_test[i]
        score_i = distance + prob_term + np.array([class_penalty.get(c, 0.0) for c in range(n_classes)])

        # Compute p-values per class
        p_i = {c: np.mean(class_nonconf_scores[c] >= score_i[c]) for c in range(n_classes)}
        p_values_list.append(p_i)

        # Determine eligible classes based on alpha thresholds
        eligible = [c for c, p in p_i.items() if p > alpha_per_class.get(c, 0.1)]
        if not eligible:
            # No class passes threshold: pick class with highest p-value
            y_pred_mcp.append(max(p_i, key=p_i.get))
        elif len(eligible) == 1:
            y_pred_mcp.append(eligible[0])
        else:
            # Multiple eligible: pick the one with highest predicted probability
            eligible_probs = y_proba_test[i, eligible]
            y_pred_mcp.append(eligible[np.argmax(eligible_probs)])

    return np.array(y_pred_mcp), p_values_list    

def compute_sample_weights(y_train):
    from collections import Counter
    import numpy as np
    
    # inverse-frequency component
    classes = np.unique(y_train)
    counts = Counter(y_train)
    total = len(y_train)
    class_weight = {c: total / (len(classes) * counts[c]) for c in classes}

    # combine with ordinal distance penalty
    median_class = np.median(classes)
    num_classes = len(classes)

    sample_weights = []
    for y in y_train:
        freq_term = class_weight[y]
        distance_term = 1 + abs(y - median_class) / (num_classes - 1)
        weight = freq_term * distance_term
        sample_weights.append(weight)

    return np.array(sample_weights)

def train_and_evaluate_ordinal_mcp(X, y, cv_outer, optimizer, feature_set:str, alpha=0.1):
    """
    Train and evaluate XGBOrdinal with Mondrian Conformal Prediction (class-conditional thresholds)
    across outer CV folds.

    Parameters
    ----------
    X : pd.DataFrame
        Feature matrix.
    y : pd.Series
        Target labels (ordinal encoded as integers).
    cv_outer : cross-validation splitter
        Outer CV, e.g., StratifiedKFold or LeaveOneOut().
    optimizer : RandomizedSearchCV
        Inner CV hyperparameter optimizer (fitted on XGBClassifier parameters).
    feature_set : str
        Name of the feature set.
    alpha : float
        Desired miscoverage rate (1 - confidence). Default 0.1 for 90% confidence.

    Returns
    -------
    results_list : list of dicts
        Each dict contains predictions, probabilities, best_model, MCP info, etc.
    """
    import time
    import numpy as np
    from xgbordinal import XGBOrdinal
    from collections import Counter
    from xgboost import XGBClassifier


    print(f"Training {feature_set} with ordinal MCP (alpha={alpha})...")
    results_list = []
    classes = np.unique(y)

    for train_idx, test_idx in cv_outer.split(X, y):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # # ----- Sample weights per class (more aggressive for rare classes) -----
        # classes = sorted(set(y_train))
        # counts = Counter(y_train)
        # total = len(y_train)
        # num_classes = len(classes)
        
        # # Base inverse-frequency weights (smoothed)
        # #base_weights = {c: min(2.0, total / (num_classes * counts[c])) for c in classes}
        # base_weights = ({c: total / (num_classes * counts[c]) for c in classes})
        
        # # Ordinal penalty: weight increases with distance from median class
        # median_class = np.median(classes)
        
        # ordinal_weights = {}
        # for c in classes:
        #     distance_from_median = abs(c - median_class) / (num_classes - 1)  # normalized distance
        #     #Setting alpha < 1 reduces the bias toward extremes, letting rare middle classes retain higher weight.
        #     #ordinal_weights[c] = base_weights[c] * (0.5 + alpha * distance_from_median)
        #     ordinal_weights[c] = base_weights[c] * (1 + distance_from_median)
        
        # # Assign sample weights
        # sample_weights = np.array([ordinal_weights[c] for c in y_train])
        
        #Compute inverse-frequency weights
        classes = sorted(set(y_train))
        counts = Counter(y_train)
        total = len(y_train)
        
        # Weight each class inversely proportional to its frequency
        class_weights = {c: total / (len(classes) * counts[c]) for c in classes}
        
        # # Assign weights to each sample
        # # Optional: amplify weight of class 0
        weight_multiplier = {0: 0.8, 1: 2.0, 2: 1.5}  # double weight for class 0
        sample_weights = np.array([class_weights[c] * weight_multiplier[c] for c in y_train])
        #sample_weights = np.array([class_weights[c] for c in y_train])
        
        #sample_weights = compute_sample_weights(y_train)

        # ----- Fit optimizer (inner CV) using XGBClassifier template -----
        start = time.time()
        optimizer.fit(X_train, y_train, sample_weight=sample_weights)
        elapsed = time.time() - start

        best_params = optimizer.best_params_

        # ----- Train XGBOrdinal with best params -----
        ordinal_model = XGBOrdinal(**best_params)
        ordinal_model.fit(X_train, y_train, sample_weight=sample_weights)
        
        # 1. Compute training nonconformity scores
        y_proba_train = ordinal_model.predict_proba(X_train)
        scores_train = compute_nonconformity_scores(y_train, y_proba_train, sample_weights)
        
        # 2. Split by class for Mondrian calibration
        classes = np.unique(y_train)
        class_nonconf_scores = {c: scores_train[y_train == c, c] for c in classes}
        
        # 3. Compute MCP predictions on test set
        y_proba_test = ordinal_model.predict_proba(X_test)
        y_pred_og = ordinal_model.predict(X_test)
        
        alpha_per_class = {0: 0.05, 1: 0.2, 2: 0.2}  # lower alpha = stricter coverage (higher alpha = class more likely to be chosen)
        class_penalty = {0: 0, 1: 0.8, 2: 0.0}  # optional
        y_pred_mcp, p_values_list = mcp_predict(y_proba_test, class_nonconf_scores, class_penalty, alpha_per_class)

        results_list.append({
            "best_model": ordinal_model,
            "best_params": best_params,
            "y_true": y_test.values,
            "y_pred_original": y_pred_og,
            "y_pred": np.array(y_pred_mcp),
            "y_pred_proba": y_proba_test,
            "p_values": p_values_list,
            "fit_time": elapsed
        })

    return results_list

def train_and_evaluate_ordinal_mcp_v2(X, y, cv_outer, optimizer, feature_set: str, alpha=0.1):
    """
    Outer CV training + MCP prediction for XGBOrdinal.
    Works with your RandomizedSearchCV setup where parameters are NOT nested.
    """

    import time
    import numpy as np
    from collections import Counter
    from xgbordinal import XGBOrdinal

    print(f"Training {feature_set} with ordinal MCP (alpha={alpha})...")
    results_list = []
    classes = np.unique(y)

    def unflatten_params_to_xgbordinal(params):
        """Convert flat parameter dict into {'xgb_params': {...}} for XGBOrdinal."""
        xgb_params = {}
        for k, v in params.items():
            # these are all real XGBoost parameters in your search space
            xgb_params[k] = v

        return {"xgb_params": xgb_params}

    for train_idx, test_idx in cv_outer.split(X, y):

        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # ----- Sample weights (your version, unchanged) -----
        classes = sorted(set(y_train))
        counts = Counter(y_train)
        total = len(y_train)

        class_weights = {c: total / (len(classes) * counts[c]) for c in classes}
        weight_multiplier = {0: 0.8, 1: 2.0, 2: 1.5}
        sample_weights = np.array([class_weights[c] * weight_multiplier[c] for c in y_train])

        # ----- Inner CV tuning -----
        start = time.time()
        optimizer.fit(X_train, y_train, sample_weight=sample_weights)
        elapsed = time.time() - start

        # Convert optimizer params → XGBOrdinal-compatible params
        best_params = unflatten_params_to_xgbordinal(optimizer.best_params_)

        # ----- Train final ordinal model -----
        ordinal_model = XGBOrdinal(
            aggregation="weighted",
            norm=True,
            random_state=42,
            **best_params              # <-- this now works
        )
        ordinal_model.fit(X_train, y_train, sample_weight=sample_weights)

        # ----- Nonconformity scores -----
        y_proba_train = ordinal_model.predict_proba(X_train)
        scores_train = compute_nonconformity_scores(y_train, y_proba_train, sample_weights)

        classes_train = np.unique(y_train)
        class_nonconf_scores = {
            c: scores_train[y_train == c, c] for c in classes_train
        }

        # ----- MCP predictions -----
        y_proba_test = ordinal_model.predict_proba(X_test)
        y_pred_original = ordinal_model.predict(X_test)

        alpha_per_class = {0: 0.05, 1: 0.2, 2: 0.2}
        class_penalty = {0: 0, 1: 0.8, 2: 0.0}

        y_pred_mcp, p_values_list = mcp_predict(
            y_proba_test,
            class_nonconf_scores,
            class_penalty,
            alpha_per_class
        )

        # ----- Store results -----
        results_list.append({
            "best_model": ordinal_model,
            "best_params": optimizer.best_params_,
            "y_true": y_test.values,
            "y_pred_original": y_pred_original,
            "y_pred": np.array(y_pred_mcp),
            "y_pred_proba": y_proba_test,
            "p_values": p_values_list,
            "fit_time": elapsed,
        })

    return results_list


def plot_perm_importance_across_folds(results_list, n_repeats=10, scoring='accuracy', random_state=42):
    """
    Compute and plot permutation importance across multiple CV folds.

    Parameters
    ----------
    results_list : list of dicts
        Each dict must contain:
          - 'X_test': pd.DataFrame
          - 'y_true': array-like
          - 'best_model': fitted model
    n_repeats : int
        Number of repeats for permutation importance.
    scoring : str
        Scoring metric for permutation importance.
    random_state : int
        Random seed for reproducibility.
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    from sklearn.inspection import permutation_importance
    import numpy as np
    
    feature_names = results_list[0]["X_test"].columns
    all_importances = {f: [] for f in feature_names}

    # Compute permutation importance for each fold
    for fold_idx, res in enumerate(results_list):
        X_test = res["X_test"]
        y_test = res["y_true"]
        model = res["best_model"]

        perm_imp = permutation_importance(
            model, X_test, y_test, n_repeats=n_repeats, random_state=random_state, scoring=scoring
        )

        for i, f in enumerate(feature_names):
            all_importances[f].extend(perm_imp.importances[i])

    # Aggregate across folds
    importance_df = pd.DataFrame({
        "feature": feature_names,
        "mean_importance": [np.mean(all_importances[f]) for f in feature_names],
        "std_importance": [np.std(all_importances[f]) for f in feature_names]
    }).sort_values(by="mean_importance", ascending=True)

    # Plot
    plt.figure(figsize=(8,6))
    plt.barh(importance_df["feature"], importance_df["mean_importance"], xerr=importance_df["std_importance"], color="skyblue")
    plt.xlabel("Permutation Importance (mean ± std across folds)")
    plt.title("Feature Importance via Permutation (All Folds)")
    plt.tight_layout()
    plt.show()

    return importance_df

def plot_shap_across_folds(results_list, X_columns):
    import shap
    import matplotlib.pyplot as plt
    from xgboost import XGBClassifier

    for fold_idx, res in enumerate(results_list):
        model = res["best_model"]
        X_test = res["X_test"]

        print(f"\n--- Fold {fold_idx + 1} ---")

        # Check if model.clfs exists
        if hasattr(model, "clfs"):
            for idx, clf in enumerate(model.clfs):
                if isinstance(clf, XGBClassifier):
                    print(f"\nFold {fold_idx+1}, Classifier {idx} SHAP summary:")
                    explainer = shap.TreeExplainer(clf)
                    shap_values = explainer.shap_values(X_test)
                    shap.summary_plot(shap_values, X_test, feature_names=X_columns, show=True)
                else:
                    print(f"Skipping internal clf at index {idx}: not a fitted XGBClassifier")
        else:
            print(f"Fold {fold_idx+1}: No internal classifiers found in XGBOrdinal")

# def plot_shap_across_folds(results_list, X_columns):
#     """
#     Plot SHAP values for each fold in results_list (from train_and_evaluate_ordinal_mcp).
#     """
#     import shap
#     import matplotlib.pyplot as plt
#     from xgboost import XGBClassifier

#     for fold_idx, res in enumerate(results_list):
#         model = res["best_model"]
#         X_test = res["X_test"]

#         print(f"\n--- Fold {fold_idx + 1} ---")

#         # Use the main internal fitted classifier
#         if hasattr(model, "clf") and isinstance(model.clf, XGBClassifier):
#             clf = model.clf
#         else:
#             raise ValueError(f"No proper fitted XGBClassifier found in fold {fold_idx + 1}")

#         # Prepare classifier for SHAP
#         explainer = shap.TreeExplainer(clf)
#         shap_values = explainer.shap_values(X_test)

#         # If multi-class, shap_values is a list
#         if isinstance(shap_values, list):
#             for c, sv in enumerate(shap_values):
#                 print(f"Class {c} SHAP summary:")
#                 shap.summary_plot(sv, X_test, feature_names=X_columns, show=True)
#         else:
#             print("Scalar output SHAP summary:")
#             shap.summary_plot(shap_values, X_test, feature_names=X_columns, show=True)

def get_shap_ready_xgb(clf):
    """
    Returns a standard XGBClassifier with a properly formatted booster for SHAP.
    clf must be a fitted XGBClassifier instance.
    """
    from xgboost import XGBClassifier

    if not hasattr(clf, "get_booster"):
        raise ValueError("clf must be a fitted XGBClassifier instance, got type: {}".format(type(clf)))

    # Extract booster
    booster = clf.get_booster()  # <-- must be instance, not class

    # Create a new XGBClassifier with the same parameters
    params = clf.get_xgb_params()
    plain_xgb = XGBClassifier(**params)
    plain_xgb._Booster = booster
    if hasattr(clf, "_le"):
        plain_xgb._le = clf._le

    return plain_xgb


def train_and_evaluate_ordinal_mcp_v1(X, y, cv_outer, optimizer, feature_set:str, alpha=0.1):
    """
    Train and evaluate XGBOrdinal with Mondrian Conformal Prediction (class-conditional thresholds)
    across outer CV folds.

    Parameters
    ----------
    X : pd.DataFrame
        Feature matrix.
    y : pd.Series
        Target labels (ordinal encoded as integers).
    cv_outer : cross-validation splitter
        Outer CV, e.g., StratifiedKFold or LeaveOneOut().
    optimizer : RandomizedSearchCV
        Inner CV hyperparameter optimizer (fitted on XGBClassifier parameters).
    feature_set : str
        Name of the feature set.
    alpha : float
        Desired miscoverage rate (1 - confidence). Default 0.1 for 90% confidence.

    Returns
    -------
    results_list : list of dicts
        Each dict contains predictions, probabilities, best_model, MCP info, etc.
    """
    import time
    import numpy as np
    from xgbordinal import XGBOrdinal
    from collections import Counter

    print(f"Training {feature_set} with ordinal MCP (alpha={alpha})...")
    results_list = []
    classes = np.unique(y)

    for train_idx, test_idx in cv_outer.split(X, y):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # ----- Sample weights per class (more aggressive for rare classes) -----
        classes = sorted(set(y_train))
        counts = Counter(y_train)
        total = len(y_train)
        num_classes = len(classes)
        
        # Base inverse-frequency weights (smoothed)
        base_weights = {c: min(2.0, total / (num_classes * counts[c])) for c in classes}
        
        # Ordinal penalty: weight increases with distance from median class
        median_class = np.median(classes)
        
        ordinal_weights = {}
        for c in classes:
            distance_from_median = abs(c - median_class) / (num_classes - 1)  # normalized distance
            #Setting alpha < 1 reduces the bias toward extremes, letting rare middle classes retain higher weight.
            ordinal_weights[c] = base_weights[c] * (0.5 + alpha * distance_from_median)
            #ordinal_weights[c] = base_weights[c] * (1 + distance_from_median)
        
        # Assign sample weights
        sample_weights = np.array([ordinal_weights[c] for c in y_train])

        # ----- Fit optimizer (inner CV) using XGBClassifier template -----
        start = time.time()
        optimizer.fit(X_train, y_train, sample_weight=sample_weights)
        elapsed = time.time() - start

        best_params = optimizer.best_params_

        # ----- Train XGBOrdinal with best params -----
        ordinal_model = XGBOrdinal(**best_params)
        ordinal_model.fit(X_train, y_train, sample_weight=sample_weights)
        
        # 1. Compute training nonconformity scores
        class_penalty = {0: 0.8, 1: 0.1, 2: 0.5}  # optional
        y_proba_train = ordinal_model.predict_proba(X_train)
        scores_train = compute_nonconformity_scores(y_train, y_proba_train, class_penalty)
        
        # 2. Split by class for Mondrian calibration
        classes = np.unique(y_train)
        class_nonconf_scores = {c: scores_train[y_train == c, c] for c in classes}
        
        # 3. Compute MCP predictions on test set
        y_proba_test = ordinal_model.predict_proba(X_test)
        y_pred_og = ordinal_model.predict(X_test)
        
        alpha_per_class = {0: 0.1, 1: 0.2, 2: 0.1}  # lower alpha = stricter coverage
        y_pred_mcp, p_values_list = mcp_predict(y_proba_test, class_nonconf_scores, alpha_per_class)

        results_list.append({
            "best_model": ordinal_model,
            "best_params": best_params,
            "y_true": y_test.values,
            "y_pred_original": y_pred_og,
            "y_pred": np.array(y_pred_mcp),
            "y_pred_proba": y_proba_test,
            "p_values": p_values_list,
            "fit_time": elapsed,
            "X_test": X_test
        })

    return results_list

def train_and_evaluate_3class_sample_weights(X, y, continuous_features, binary_features,
                              cv_outer, cv_inner, param_distributions,
                              feature_set:str, n_iter=30, random_state=42):
    """
    Train and evaluate a 3-class ordinal classifier with preprocessing and sample weights.
    No SMOTE is used; class imbalance is handled via XGBOrdinal scale_pos_weight.

    Parameters
    ----------
    X : pd.DataFrame
        Features.
    y : pd.Series
        Target (3-class ordinal).
    continuous_features : list[str]
        Names of continuous columns.
    binary_features : list[str]
        Names of binary/categorical columns.
    cv_outer : sklearn.model_selection.StratifiedKFold
        Outer CV splitter.
    cv_inner : sklearn.model_selection.BaseCrossValidator
        Inner CV splitter for RandomizedSearchCV.
    param_distributions : dict
        Parameter distributions for RandomizedSearchCV.
    feature_set : str
        Name of the feature tier (for logging).
    n_iter : int
        Number of iterations for RandomizedSearchCV.
    random_state : int
        Random state.

    Returns
    -------
    results_list : list[dict]
        Best model, best score, predictions, fit times, and test data per fold.
    """
    import time
    import numpy as np
    from sklearn.pipeline import Pipeline
    from sklearn.compose import ColumnTransformer
    from sklearn.impute import SimpleImputer
    from sklearn.model_selection import RandomizedSearchCV
    from xgbordinal import XGBOrdinal
    from sklearn.model_selection import StratifiedShuffleSplit

    print(f"Training and evaluating {feature_set} model...")
    results_list = []

    # Preprocessing
    preprocessor = ColumnTransformer([
        ('cont', SimpleImputer(strategy='median'), continuous_features),
        ('bin', SimpleImputer(strategy='most_frequent'), binary_features)
    ])

    for fold, (train_idx, test_idx) in enumerate(cv_outer.split(X, y), 1):
        print(f"\nOuter fold {fold}/{cv_outer.get_n_splits()}")
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # Class counts
        class_counts = np.bincount(y_train)
        count_0, count_1, count_2 = class_counts
        scale_pos_weight = {0: 1, 1: 0.9*(count_0 + count_2) / count_1 if count_1 > 0 else 1, 2: 1}
        #scale_pos_weight = {0:1, 1: 0.7*(count_0+count_2)/count_1, 2:1}
        sample_weight = np.array([scale_pos_weight[cls] for cls in y_train])

        print("Class distribution:", class_counts)
        print("Scale_pos_weight applied:", scale_pos_weight)

        # Define pipeline
        xgb = XGBOrdinal(aggregation='weighted', norm=True,
                         random_state=random_state, n_jobs=-1)
        pipeline = Pipeline([
            ('preprocessor', preprocessor),
            ('classifier', xgb)
        ])
        
        # Convert y_train to binary for stratifying on class 1
        y_train_class1 = np.where(y_train == 1, 1, 0)
        cv_inner_class1 = StratifiedShuffleSplit(n_splits=cv_inner.get_n_splits(), 
                                                 test_size=cv_inner.test_size, 
                                                 random_state=cv_inner.random_state)

        # RandomizedSearchCV
        optimizer = RandomizedSearchCV(
            estimator=pipeline,
            param_distributions=param_distributions,
            n_iter=n_iter,
            cv=cv_inner_class1.split(X_train, y_train_class1),  # <-- use this split
            #cv=cv_inner,
            n_jobs=2,
            verbose=1,
            scoring='f1_macro',
            random_state=random_state,
            refit=True
        )

        start = time.time()
        optimizer.fit(X_train, y_train, classifier__sample_weight=sample_weight)
        elapsed = time.time() - start

        best_pipeline = optimizer.best_estimator_
        best_score = optimizer.best_score_
        y_pred = best_pipeline.predict(X_test)
        y_pred_proba = best_pipeline.predict_proba(X_test)
        
        # Apply threshold tweak for class 1
        y_pred = np.argmax(y_pred_proba, axis=1)
        y_pred[(y_pred == 1) & (y_pred_proba[:,1] < 0.3)] = 0  # only predict class 1 if probability > 0.4

        results_list.append({
            "fold": fold,
            "best_model": best_pipeline,
            "best_score": best_score,
            "y_true": y_test.values,
            "y_pred": y_pred,
            "y_pred_proba": y_pred_proba,
            "fit_time": elapsed,
            "X_test": X_test
        })

        print(f"Fold {fold} done. Best inner CV score: {best_score:.4f}, fit time: {elapsed:.1f}s")

    return results_list


def train_and_evaluate_3class(X, y, continuous_features, binary_features,
                              cv_outer, cv_inner, param_distributions,
                              feature_set:str, n_iter=30, random_state=42):
    """
    Train and evaluate a 3-class ordinal classifier with preprocessing, SMOTE, and nested CV.

    Parameters:
    - X: pd.DataFrame, features
    - y: pd.Series, target (3-class ordinal)
    - continuous_features: list[str], continuous column names
    - binary_features: list[str], binary/categorical column names
    - cv_outer: outer CV splitter (StratifiedKFold)
    - cv_inner: inner CV splitter (for RandomizedSearchCV)
    - param_distributions: dict of hyperparameters for RandomizedSearchCV
    - feature_set: str, feature tier name for logging
    - n_iter: int, number of iterations for RandomizedSearchCV
    - random_state: int, random seed

    Returns:
    - results_list: list of dicts with best_model, best_score, y_true, y_pred, y_pred_proba, fit_time, X_test
    """
    import time
    import numpy as np
    from sklearn.impute import SimpleImputer
    from sklearn.compose import ColumnTransformer
    from imblearn.over_sampling import SMOTE
    from imblearn.pipeline import Pipeline as ImbPipeline
    from sklearn.model_selection import RandomizedSearchCV
    from xgbordinal import XGBOrdinal  # assuming this is your ordinal classifier

    results_list = []
    print(f"Training and evaluating {feature_set} model...")

    # Preprocessing step
    preprocessor = ColumnTransformer([
        ('cont', SimpleImputer(strategy='median'), continuous_features),
        ('bin', SimpleImputer(strategy='most_frequent'), binary_features)
    ])

    for fold, (train_idx, test_idx) in enumerate(cv_outer.split(X, y), 1):
        print(f"\nOuter fold {fold}/{cv_outer.get_n_splits()}")
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # Print class distribution before SMOTE
        print("Class distribution before SMOTE:", np.bincount(y_train))
        

        xgb = XGBOrdinal(aggregation='weighted', norm=True,
                                  random_state=random_state, n_jobs=-1)
        
        # Define pipeline with preprocessing + SMOTE + classifier
        pipeline = ImbPipeline([
            ('preprocessor', preprocessor),
            ('smote', SMOTE(sampling_strategy='not majority', random_state=random_state)),
            ('classifier', xgb)
        ])

        # RandomizedSearchCV wrapping the pipeline
        optimizer = RandomizedSearchCV(
            estimator=pipeline,
            param_distributions=param_distributions,
            n_iter=n_iter,
            cv=cv_inner,
            n_jobs=2,
            verbose=1,
            scoring='neg_log_loss',
            random_state=random_state,
            refit=True
        )

        start = time.time()
        optimizer.fit(X_train, y_train)
        elapsed = time.time() - start

        # Get best pipeline
        best_pipeline = optimizer.best_estimator_

        # Optional: print class distribution after SMOTE (using preprocessed X)
        X_train_imputed = best_pipeline.named_steps['preprocessor'].transform(X_train)
        X_resampled, y_resampled = best_pipeline.named_steps['smote'].fit_resample(X_train_imputed, y_train)
        print("Class distribution after SMOTE (preview):", np.bincount(y_resampled))

        # Predictions
        y_pred = best_pipeline.named_steps['classifier'].predict(
            best_pipeline.named_steps['preprocessor'].transform(X_test)
        )
        y_pred_proba = best_pipeline.named_steps['classifier'].predict_proba(
            best_pipeline.named_steps['preprocessor'].transform(X_test)
        )

        results_list.append({
            "fold": fold,
            "best_model": best_pipeline.named_steps['classifier'],
            "best_score": optimizer.best_score_,
            "y_true": y_test.values,
            "y_pred": y_pred,
            "y_pred_proba": y_pred_proba,
            "fit_time": elapsed,
            "X_test": X_test
        })

        print(f"Fold {fold} done. Best inner CV score: {optimizer.best_score_:.4f}, fit time: {elapsed:.1f}s")

    return results_list


def train_and_evaluate(X, y, cv_outer, optimizer, feature_set:str):
    import time
    
    print(f"Training and evaluating {feature_set} model...")
    results_list = []
    
    for train_idx, test_idx in cv_outer.split(X, y):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        start = time.time()
        optimizer.fit(X_train, y_train)
        elapsed = time.time() - start

        print(f"Best score (inner CV) for fold: {optimizer.best_score_:.4f}")
        best_model = optimizer.best_estimator_

        y_pred = best_model.predict(X_test)
        y_pred_proba = best_model.predict_proba(X_test)[:, 1]

        results_list.append({
            "best_model": best_model,
            "best_score": optimizer.best_score_,
            "y_true": y_test.values,
            "y_pred": y_pred,
            "y_pred_proba": y_pred_proba,
            "fit_time": elapsed,
            "X_test": X_test
        })
    
    return results_list

def train_and_evaluate_mcp(X, y, cv_outer, optimizer, feature_set:str, alpha):
    """
    Train and evaluate XGBoost with Mondrian Conformal Prediction (class-conditional thresholds)
    across outer CV folds (LOOCV or K-fold).

    Parameters
    ----------
    X : pd.DataFrame
        Feature matrix.
    y : pd.Series
        Target labels (binary: minority=1, majority=0).
    cv_outer : cross-validation splitter
        Outer CV, e.g., LeaveOneOut() or StratifiedKFold().
    optimizer : RandomizedSearchCV
        Inner CV hyperparameter optimizer.
    feature_set : str
        Name of the feature set.
    alpha : float
        Desired miscoverage rate (1 - confidence). Default 0.1 for 90% confidence.

    Returns
    -------
    results_list : list of dicts
        Each dict contains predictions, probabilities, best_model, MCP info, etc.
    """
    import time
    import numpy as np
    from sklearn.utils import check_array
    from sklearn.utils import class_weight
    
    
    print(f"Training and evaluating {feature_set} model with MCP (alpha={alpha})...")
    results_list = []

    for train_idx, test_idx in cv_outer.split(X, y):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # Dynamically compute scale_pos_weight for imbalance
        #pos_class_weight = (len(y_train) - np.sum(y_train)) / np.sum(y_train)
        #optimizer.estimator.set_params(scale_pos_weight=pos_class_weight)
        #sample_weights = class_weight.compute_sample_weight(class_weight='balanced', y=y_train)

        # Fit optimizer (inner CV)
        start = time.time()
        #optimizer.fit(X_train, y_train, sample_weights=sample_weights)
        optimizer.fit(X_train, y_train)
        elapsed = time.time() - start

        best_model = optimizer.best_estimator_

        # Predict probabilities on training set (needed for MCP nonconformity scores)
        train_probs = best_model.predict_proba(X_train)
        train_probs = check_array(train_probs)

        # Nonconformity score: 1 - probability of true class
        nonconf_scores = np.array([1 - train_probs[i, y_train.iloc[i]] for i in range(len(y_train))])

        # Class-conditional nonconformity scores (Mondrian)
        classes = np.unique(y_train)
        class_nonconf = {c: nonconf_scores[y_train == c] for c in classes}

        # Predict probabilities on test set
        test_probs = best_model.predict_proba(X_test)
        test_probs = check_array(test_probs)

        # ---- Original predictions (0.5 threshold) ----
        y_pred_original = best_model.predict(X_test)

        # Compute MCP p-values for each class
        p_values = []
        for i in range(len(X_test)):
            p_i = {}
            for c in classes:
                # Fraction of training scores >= 1 - predicted probability
                p_i[c] = np.mean(class_nonconf[c] >= (1 - test_probs[i, c]))
            p_values.append(p_i)

        # Assign predicted class: pick any class with p-value > alpha
        y_pred_mcp = []
        for i, pv in enumerate(p_values):
            eligible_classes = [c for c, p in pv.items() if p > alpha]
            if len(eligible_classes) == 0:
                # No class passes threshold: pick class with max p-value
                y_pred_mcp.append(max(pv, key=pv.get))
            elif len(eligible_classes) == 1:
                y_pred_mcp.append(eligible_classes[0])
            else:
                # Multiple classes eligible: pick class with highest probability
                y_pred_mcp.append(max(eligible_classes, key=lambda c: test_probs[i, c]))

        results_list.append({
            "best_model": best_model,
            "best_score": optimizer.best_score_,
            "y_true": y_test.values,
            "y_pred": np.array(y_pred_mcp),          # MCP-adjusted
            "y_pred_original": y_pred_original,      # Original 0.5-threshold
            "y_pred_proba": test_probs,
            "p_values": p_values,
            "fit_time": elapsed,
            "X_test": X_test
        })

    return results_list

def visualize_mcp_strength(results_list, alpha_per_class=None):
    """
    Visualize the strength of MCP outcomes (prediction set sizes and coverage).
    """

    import numpy as np
    import matplotlib.pyplot as plt

    n_classes = len(results_list[0]['y_pred_proba'][0])
    if alpha_per_class is None:
        alpha_per_class = {c: 0.1 for c in range(n_classes)}

    set_sizes = []
    correct_cover = []
    all_y_true = []

    for fold_res in results_list:
        y_true = fold_res['y_true']
        p_values = fold_res['p_values']
        for i, true_label in enumerate(y_true):
            included = [c for c in range(n_classes)
                        if p_values[i][c] > alpha_per_class.get(c, 0.1)]
            set_sizes.append(len(included))
            correct_cover.append(true_label in included)
            all_y_true.append(true_label)

    set_sizes = np.array(set_sizes)
    correct_cover = np.array(correct_cover)

    # --- Plot 1: distribution of prediction set sizes ---
    plt.figure(figsize=(6, 4))
    plt.hist(set_sizes, bins=np.arange(0.5, n_classes + 1.5, 1), align='mid',
             color='steelblue', rwidth=0.8)
    plt.xticks(range(1, n_classes + 1))
    plt.xlabel('Prediction set size (#classes included)')
    plt.ylabel('Frequency')
    plt.title('Distribution of MCP prediction set sizes')
    plt.show()

    # --- Plot 2: coverage by set size ---
    unique_sizes = np.unique(set_sizes)
    coverage_by_size = [correct_cover[set_sizes == s].mean() for s in unique_sizes]

    plt.figure(figsize=(6, 4))
    plt.plot(unique_sizes, coverage_by_size, marker='o', color='darkorange')
    plt.xlabel('Prediction set size')
    plt.ylabel('Empirical coverage')
    plt.title('Coverage vs. prediction set size (MCP strength)')
    plt.ylim(0, 1)
    plt.show()

    # --- Summary ---
    avg_size = set_sizes.mean()
    overall_coverage = correct_cover.mean()
    print(f"Average prediction set size: {avg_size:.2f}")
    print(f"Overall empirical coverage: {overall_coverage:.3f}")

def visualize_mcp_results(results_list, alpha_per_class=None):
    """
    Visualize effect of Mondrian Conformal Prediction (MCP) using original predictions
    stored in 'y_pred_original' and MCP-adjusted predictions in 'y_pred'.
    Includes class-wise probability distributions and MCP intervals.

    Parameters
    ----------
    results_list : list of dicts
        Output from `train_and_evaluate_ordinal_mcp`.
    alpha_per_class : dict, optional
        Class-specific MCP miscoverage rates, e.g., {0: 0.05, 1: 0.2, 2: 0.2}.
        If None, uses global alpha=0.1.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import confusion_matrix

    all_probs_before = []
    all_probs_after = []
    all_mcp_intervals = []

    n_classes = len(results_list[0]['y_pred_proba'][0])  # assume same across folds

    if alpha_per_class is None:
        alpha_per_class = {c: 0.1 for c in range(n_classes)}

    for fold_res in results_list:
        test_probs = fold_res['y_pred_proba']
        p_values = fold_res['p_values']

        for i in range(len(test_probs)):
            # Original probabilities
            all_probs_before.append(test_probs[i])

            # MCP-adjusted probability: pick class with max p-value
            y_pred_mcp = np.argmax([p_values[i][c] for c in range(n_classes)])
            all_probs_after.append(test_probs[i, y_pred_mcp])

            # MCP interval per sample using class-specific alpha
            # Interval: include classes where p_value > alpha[class]
            lower = min([test_probs[i,c] for c in range(n_classes) if p_values[i][c] > alpha_per_class.get(c, 0.1)] + [0])
            upper = max([test_probs[i,c] for c in range(n_classes) if p_values[i][c] > alpha_per_class.get(c, 0.1)] + [0])
            all_mcp_intervals.append((lower, upper))

    all_probs_before = np.array(all_probs_before)
    all_probs_after = np.array(all_probs_after)
    all_mcp_intervals = np.array(all_mcp_intervals)

    # 1. Class-wise histogram (original)
    plt.figure(figsize=(10,6))
    for c in range(n_classes):
        plt.hist(all_probs_before[:, c], bins=20, alpha=0.5, label=f'Class {c} probs')
    plt.xlabel('Predicted probability')
    plt.ylabel('Frequency')
    plt.title('Class-wise probability distribution (original)')
    plt.legend()
    plt.show()

    # 2. MCP-adjusted histogram
    plt.figure(figsize=(10,4))
    plt.hist(all_probs_after, bins=20, alpha=0.7, color='orange', label='MCP-adjusted (max p-value)')
    plt.xlabel('Probability')
    plt.ylabel('Frequency')
    plt.title('MCP-adjusted probabilities (max p-value class)')
    plt.legend()
    plt.show()

    # 3. MCP intervals
    plt.figure(figsize=(10,4))
    plt.fill_between(np.arange(len(all_mcp_intervals)),
                     all_mcp_intervals[:,0], all_mcp_intervals[:,1],
                     color='skyblue', alpha=0.5, label='MCP interval per sample')
    plt.scatter(np.arange(len(all_probs_after)), all_probs_after, color='red', s=10, label='MCP-selected probability')
    plt.xlabel('Sample index')
    plt.ylabel('Probability')
    plt.title('Mondrian Conformal Prediction intervals (class-specific alpha)')
    plt.legend()
    plt.show()
    
    # ==== MCP vs Original Comparison ====
    y_true_all = np.concatenate([r["y_true"] for r in results_list])
    y_pred_orig_all = np.concatenate([r["y_pred_original"] for r in results_list])
    y_pred_mcp_all = np.concatenate([r["y_pred"] for r in results_list])
    
    total_changed = np.sum(y_pred_mcp_all != y_pred_orig_all)
    corrected = np.sum((y_pred_mcp_all == y_true_all) & (y_pred_orig_all != y_true_all))
    worsened = np.sum((y_pred_mcp_all != y_true_all) & (y_pred_orig_all == y_true_all))
    
    print("\n=== MCP vs Original Threshold Comparison ===")
    print(f"Total predictions changed by MCP: {total_changed}/{len(y_true_all)} "
          f"({100 * total_changed / len(y_true_all):.1f}%)")
    print(f"Corrected by MCP (was wrong, now right): {corrected}")
    print(f"Worsened by MCP (was right, now wrong): {worsened}")
    
    print("\nConfusion matrices:")
    print("Original:\n", confusion_matrix(y_true_all, y_pred_orig_all))
    print("MCP-adjusted:\n", confusion_matrix(y_true_all, y_pred_mcp_all))
    

def train_and_evaluate_with_imputation(X, y, cv_outer, estimator, feature_set: str, use_optimizer=True):
    """
    Train and evaluate a model with imputation within folds.
    
    Parameters
    ----------
    X : pd.DataFrame
        Feature matrix.
    y : pd.Series
        Target variable.
    cv_outer : cross-validation splitter
    estimator : estimator or optimizer
        If use_optimizer=True, should be BayesSearchCV / GridSearchCV; else a standard estimator.
    feature_set : str
        Name of the feature set for printing.
    use_optimizer : bool
        Whether estimator is a hyperparameter optimizer or a plain estimator.
    
    Returns
    -------
    results_list : list of dicts
    """
    from sklearn.impute import SimpleImputer
    from sklearn.compose import ColumnTransformer
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import OneHotEncoder
    import numpy as np
    import pandas as pd
    import time
    
    print(f"Training and evaluating {feature_set} model with imputation...")
    results_list = []

    # Identify categorical vs numeric columns
    cat_cols = X.select_dtypes(include=['object', 'category']).columns.tolist()
    num_cols = X.select_dtypes(include=[np.number]).columns.tolist()

    for train_idx, test_idx in cv_outer.split(X, y):
        X_train, X_test = X.iloc[train_idx].copy(), X.iloc[test_idx].copy()
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # Imputation transformers
        num_imputer = SimpleImputer(strategy='median')
        cat_imputer = SimpleImputer(strategy='most_frequent')

        preprocessor = ColumnTransformer(
            transformers=[
                ('num', num_imputer, num_cols),
                ('cat', cat_imputer, cat_cols)
            ],
            remainder='passthrough'
        )

        # Wrap in pipeline
        pipeline = Pipeline([
            ('imputer', preprocessor),
            ('model', estimator)
        ])

        start = time.time()
        pipeline.fit(X_train, y_train)
        elapsed = time.time() - start

        # Extract best estimator
        if use_optimizer:
            best_model = pipeline.named_steps['model'].best_estimator_
            best_score = pipeline.named_steps['model'].best_score_
        else:
            best_model = pipeline.named_steps['model']
            best_score = None  # No cross-validated inner score

        # Transform test set
        X_test_imputed = pd.DataFrame(preprocessor.transform(X_test),
                                      columns=num_cols + cat_cols,
                                      index=X_test.index)

        y_pred = best_model.predict(X_test_imputed)
        y_pred_proba = best_model.predict_proba(X_test_imputed)[:, 1]

        results_list.append({
            "best_model": best_model,
            "best_score": best_score,
            "y_true": y_test.values,
            "y_pred": y_pred,
            "y_pred_proba": y_pred_proba,
            "fit_time": elapsed,
            "X_test": X_test_imputed
        })

    return results_list


def train_and_evaluate_smote(X, y, cat_cols, num_cols, cv_inner, cv_outer, preprocessor, tier_name="clinical_features"):
    from imblearn.over_sampling import SMOTENC, SMOTE
    from sklearn.impute import SimpleImputer
    from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
    from sklearn.pipeline import Pipeline
    from sklearn.compose import ColumnTransformer
    from sklearn.preprocessing import StandardScaler
    from xgboost import XGBClassifier
    from skopt import BayesSearchCV
    from skopt.space import Integer, Real
    from sklearn.model_selection import RandomizedSearchCV
    import numpy as np
    from scipy.stats import uniform, randint
    import time
    
    print(f"Training and evaluating {tier_name} model...")
    
    results_list = []
    
    # Apply SMOTENC before cross-validation splits
    smote = SMOTENC(categorical_features=cat_cols, random_state=42)
    
    # Base model (XGBoost classifier)
    pos_class_weight = (len(y) - np.sum(y)) / np.sum(y)
    xgb_clf = XGBClassifier(eval_metric='auc',
                            scale_pos_weight=pos_class_weight,
                            max_delta_step=1,
                            random_state=42,
                            n_jobs=2)

    # Parameter distributions for RandomizedSearchCV
    param_distributions = {
        'max_depth': randint(3, 5),
        'learning_rate': uniform(0.01, 0.09),
        'n_estimators': randint(100, 151),
        'reg_alpha': uniform(1e-5, 1.0),  # L1 regularization
        'reg_lambda': uniform(1e-5, 1.0), # L2 regularization
        'subsample': uniform(0.5, 0.5)     # Row subsampling
    }
    
    # Randomized search optimizer
    optimizer = RandomizedSearchCV(
        estimator=xgb_clf,
        param_distributions=param_distributions,
        n_iter=30,
        cv=cv_inner,
        n_jobs=2,
        verbose=1,
        scoring="roc_auc",
        random_state=42,
        refit=True
    )
    
    # Outer cross-validation loop
    for train_idx, test_idx in cv_outer.split(X, y):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
    
        # Preprocess data (impute missing values in train and test data)
        X_train_processed = preprocessor.fit_transform(X_train)
        X_test_processed = preprocessor.transform(X_test)

        # Apply SMOTENC to the training data only
        X_resampled, y_resampled = smote.fit_resample(X_train_processed, y_train)
    
        # Train the model on resampled, preprocessed data
        start = time.time()
        optimizer.fit(X_resampled, y_resampled)
        elapsed = time.time() - start
        
        print(f"Best score (inner CV) for fold: {optimizer.best_score_:.4f}")
        best_model = optimizer.best_estimator_
    
        # Evaluate the best model on the test set (without SMOTE)
        #X_test_transformed = preprocessor.transform(X_test)
        y_pred = best_model.predict(X_test_processed)
        y_pred_proba = best_model.predict_proba(X_test_processed)[:, 1]
    
        results_list.append({
            "best_model": best_model,
            "best_score": optimizer.best_score_,
            "y_true": y_test.values,
            "y_pred": y_pred,
            "y_pred_proba": y_pred_proba,
            "fit_time": elapsed,
            "X_test": X_test
        })
    
    return results_list

def bootstrap_auc_multiclass(results_list, n_iterations=100, random_state=42, plot_ci=True):
    #26112026
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_auc_score, roc_curve
    from sklearn.utils import resample
    from sklearn.preprocessing import label_binarize

    np.random.seed(random_state)

    # Get all classes
    all_classes = np.unique(np.concatenate([r['y_true'] for r in results_list]))

    # Initialize storage
    auc_roc_scores = {cls: [] for cls in all_classes}
    auc_pr_scores = {cls: [] for cls in all_classes}
    fpr_per_class = {cls: [] for cls in all_classes}
    tpr_per_class = {cls: [] for cls in all_classes}

    # Bootstrapping
    for result in results_list:
        y_true = np.array(result['y_true'])
        y_proba = np.array(result['y_proba'])  # shape (n_samples, n_classes)
        y_true_bin = label_binarize(y_true, classes=all_classes)

        for iteration in range(n_iterations):
            indices = resample(np.arange(len(y_true)), n_samples=len(y_true),
                               random_state=random_state + iteration)
            y_true_resampled = y_true_bin[indices]
            y_proba_resampled = y_proba[indices]

            for i, cls in enumerate(all_classes):
                if y_true_resampled[:, i].sum() in [0, len(y_true_resampled)]:
                    continue  # skip invalid resamples

                # ROC-AUC
                auc_roc = roc_auc_score(y_true_resampled[:, i], y_proba_resampled[:, i])
                auc_roc_scores[cls].append(auc_roc)
                fpr, tpr, _ = roc_curve(y_true_resampled[:, i], y_proba_resampled[:, i])
                fpr_per_class[cls].append(fpr)
                tpr_per_class[cls].append(tpr)

    # Plotting ROC curves
    plt.figure(figsize=(10, 8))
    for cls in all_classes:
        if not auc_roc_scores[cls]:
            continue
        mean_auc = np.mean(auc_roc_scores[cls])

        # Interpolate TPR at common FPR grid
        mean_fpr = np.linspace(0, 1, 100)
        interp_tprs = [np.interp(mean_fpr, fpr, tpr) 
                       for fpr, tpr in zip(fpr_per_class[cls], tpr_per_class[cls])]
        mean_tpr = np.mean(interp_tprs, axis=0)
        lower_tpr = np.percentile(interp_tprs, 2.5, axis=0)
        upper_tpr = np.percentile(interp_tprs, 97.5, axis=0)

        plt.plot(mean_fpr, mean_tpr, label=f'Class {cls} ROC (AUC={mean_auc:.3f})')
        
        if plot_ci:
            plt.fill_between(mean_fpr, lower_tpr, upper_tpr, alpha=0.2)

    plt.plot([0, 1], [0, 1], linestyle='--', color='gray')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Bootstrapped One-vs-Rest ROC Curves')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.show()
    
    # Summary statistics
    summary = {cls: {
        'roc_mean': np.mean(auc_roc_scores[cls]),
        'roc_2.5%': np.percentile(auc_roc_scores[cls], 2.5),
        'roc_97.5%': np.percentile(auc_roc_scores[cls], 97.5)
    } for cls in all_classes}
    
    for cls, stats in summary.items():
        print(f"Class {cls}: ROC-AUC={stats['roc_mean']:.3f} [{stats['roc_2.5%']:.3f}, {stats['roc_97.5%']:.3f}]")

    return summary, auc_roc_scores, auc_pr_scores
