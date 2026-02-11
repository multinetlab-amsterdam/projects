#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 11:50:57 2025

@author: ekoderman (e.koderman@amsterdamumc.nl)
status: final

REVIEW
reviewed by: Sebastien Dam
review date: January 2026

Split 30% of the data stratified on clinical variables and main outcomes
output:
    - train/test subject ID files
    - combined_dataframe containing all subjects and stratas used

EDA is performed for comparison between the development/external (train/test set)
In the following scripts, the train/test ID files can be used to index into the combined_dataframe.
"""

# Imports
from pathlib import Path
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit
from utils import EDA_visualization, validate_merged_df
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
import numpy as np

# Paths
combined_path = Path("/folder/data/processed/combined/combined_ID_input_target_features_22122025.csv")
output_dir = Path("/folder/subj_ids/train_test_IDs/random_stratified")
clinical_path = Path("/folder/data/processed/clinical_data/clinical_data_filtered.csv")

# Load
combined_df = pd.read_csv(combined_path)
clinical_data = pd.read_csv(clinical_path)

# stratify on the main outcomes:
    # - event (censored (0-alive) or uncensored (1-death))
    # - survival within three categories based on alex paper (short term (up to 12m), midd term (12m-15m), long term (15m -))
    # - KPS post-op 1-funct independent or 0-funct dependent

##### SURVIVAL
# Create the 'survival_categories' only for patients with event = 1
combined_df.loc[combined_df['event'] == 1, 'survival_categories'] = pd.cut(
    combined_df.loc[combined_df['event'] == 1, 'time_to_event'], 
    bins=[-1, 12, 15, float('inf')], 
    labels=["short-term", "mid-term", "long-term"]
)

# Add "censored" as an allowed category
combined_df["survival_categories"] = (
    combined_df["survival_categories"]
    .astype("category")
    .cat.add_categories(["censored"])
)

# Now safe to assign
combined_df.loc[combined_df['event'] == 0, 'survival_categories'] = "censored"

##### KPS BINNED
# Already in the combined_df but keep raw KPS post-op values for future as well

# Create composite stratification label combining survival interval and KPS
combined_df["strata"] = (
    combined_df["survival_categories"].astype(str) + "_" +
    combined_df["KPS_post_op_binned"].astype(str)
)

#%% # If some composite strata are very small, collapse or remove rare strata (otherwise StratifiedShuffleSplit will fail).

# ---- Drop rare strata ----
min_count = 2  # need at least a few per group for stratification
strata_counts = combined_df["strata"].value_counts()

valid_strata = strata_counts[strata_counts >= min_count].index
invalid_strata = strata_counts[strata_counts < min_count]

# Subset the dataframe to only valid strata
filtered_df = combined_df[combined_df["strata"].isin(valid_strata)].reset_index(drop=True)

# Optional: inspect or save invalid combinations
invalid_combinations = (
    combined_df[combined_df["strata"].isin(invalid_strata.index)]
    .assign(strata_count=lambda x: x["strata"].map(strata_counts))
)

# Example: view the rare strata combinations
print("\n Rare strata combinations (below threshold):")
print(invalid_combinations["strata"].value_counts())

# Good - no rare strata, we keep all the data!

#%%
# ---- Split stratified ----
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.3, random_state=42)
train_idx, test_idx = next(sss.split(filtered_df, filtered_df["strata"]))

dev_df = filtered_df.iloc[train_idx].reset_index(drop=True)
external_test_df = filtered_df.iloc[test_idx].reset_index(drop=True)

# Add a 'dataset' column to distinguish the two sets
filtered_df["source"] = np.where(
    filtered_df["PRECOG_ID"].isin(dev_df["PRECOG_ID"]),
    "Development",
    np.where(filtered_df["PRECOG_ID"].isin(external_test_df["PRECOG_ID"]), "External test", "Excluded")
)

print("Development set:", dev_df.shape)
print("External test set:", external_test_df.shape)
print("Number of strata used:", filtered_df['strata'].nunique())

# only X (number of strata) of those combinations actually exist (i.e., appear in your data with at least 3 patients each).
# More strata = finer control over balance,
# but too many small strata = risk of singleton failures
# 6-20: Reasonable diversity; good for stratified splitting

#%% ## EDA: check what were working with

for col in ['survival_categories','KPS_post_op_binned']:
    print('test')
    print(filtered_df.iloc[test_idx][col].value_counts(normalize=True))
    print('development')
    print(filtered_df.iloc[train_idx][col].value_counts(normalize=True))

# Train and external test sets should have similar distributions!

fig, axes = plt.subplots(1, 2, figsize=(15, 5))

sns.countplot(
    data=filtered_df,
    x='survival_categories', hue='source', ax=axes[0], legend=False
)
axes[0].set_title("Survival category distribution")

sns.countplot(
    data=filtered_df,
    x='KPS_post_op_binned', hue='source', ax=axes[1], legend=False
)
axes[1].set_title("KPS post-op distribution")

plt.tight_layout()
plt.show()

# Summary stats
print("\n--- TRAIN SUMMARY ---")
print(dev_df.describe(include='all').T)
validate_merged_df(dev_df, key_col='PRECOG_ID')
EDA_visualization(dev_df, 'development set')

print("\n--- TEST SUMMARY ---")
print(external_test_df.describe(include='all').T)
validate_merged_df(external_test_df, key_col='PRECOG_ID')
EDA_visualization(external_test_df, 'external test set')

#%% # EDA: Overview of surgery years in the held-out set

external_test_df = filtered_df[filtered_df['source'] == 'External test']
clinical_external = clinical_data[clinical_data['PRECOG_ID'].isin(external_test_df['PRECOG_ID'])]

# Extract year from datetime and drop NaNs just in case
years = clinical_external['surgery_date'].str[:4].astype(int)

# Plot histogram
sns.histplot(years, bins=range(int(years.min()), int(years.max()) + 1))
plt.xlabel("Surgery Year")
plt.ylabel("Count")
plt.title("Distribution of Surgery Years")

#%% # EDA: further outcome distribution checking between training/test set (development/external)

clinical_external = clinical_data[clinical_data['PRECOG_ID'].isin(external_test_df['PRECOG_ID'])]
clinical_development = clinical_data[clinical_data['PRECOG_ID'].isin(dev_df['PRECOG_ID'])]

survivors_above_12m_dev = clinical_development[clinical_development['months_resection_death'] >= 12]

kps_above_0_dev = clinical_development[clinical_development['KPS_post-op_value'] > 0]

conditions = [kps_above_0_dev['KPS_post-op_value'] < 70,
              kps_above_0_dev['KPS_post-op_value'] >= 70]
values = [1,2]

kps_above_0_dev['kps_post_op_binary'] = np.select(conditions,values)

#distribution_plot(kps_above_0_dev,'kps_post_op_binary',['KPS < 70', 'KPS ≥ 70'],'Dev set: KPS distributions of patients survival at least 12m')

survivors_above_12m_ext = clinical_external[clinical_external['months_resection_death'] >= 12]
kps_above_0_ext = clinical_external[clinical_external['KPS_post-op_value'] > 0]

conditions = [kps_above_0_ext['KPS_post-op_value'] < 70,
              kps_above_0_ext['KPS_post-op_value'] >= 70]
values = [1,2]

kps_above_0_ext['kps_post_op_binary'] = np.select(conditions,values)

#distribution_plot(kps_above_0_ext,'kps_post_op_binary',['KPS < 70', 'KPS ≥ 70'],'External test set: KPS distributions of patients survival at least 12m')

# --- Plot side by side ---
sns.set(style='whitegrid')
fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)

datasets = [
    (kps_above_0_dev, f"Dev set: N={len(kps_above_0_dev)}"),
    (kps_above_0_ext, f"External test set: N={len(kps_above_0_ext)}")
]
labels = ['KPS < 70', 'KPS ≥ 70']

for ax, (df, title) in zip(axes, datasets):
    total = len(df)
    ax = sns.histplot(df['kps_post_op_binary'], discrete=True, ax=ax, shrink=0.8)

    # Annotate
    for p in ax.patches:
        count = int(p.get_height())
        percent = (count / total) * 100 if total > 0 else 0
        x = p.get_x() + p.get_width() / 2
        y = p.get_height()
        ax.annotate(f'{count} ({percent:.1f}%)', (x, y), ha='center', va='bottom', fontsize=9)
    
    # Tick labels matching available categories
    unique_vals = sorted(df['kps_post_op_binary'].dropna().unique())
    ax.set_xticks(unique_vals)
    ax.set_xticklabels([labels[i-1] for i in unique_vals])
    
    ax.set_title(title)
    ax.set_xlabel('')
    ax.set_ylabel('Count')

fig.suptitle("KPS Post-Op Distributions (≥12-month survivors)", fontsize=13)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

#%% STORING DATA

### store the full combined df with the strata rule
filtered_df.to_csv('/folder/data/processed/combined/combined_df_strata_groups_22122025.csv')

### store the numpy indices
np.save(output_dir / 'development_idx_22122025.npy', train_idx)
np.save(output_dir / 'external_test_set_idx_22122025.npy', test_idx)

### store the pipeline used for splitting
joblib.dump(sss, output_dir / 'sss_22122025.joblib')
