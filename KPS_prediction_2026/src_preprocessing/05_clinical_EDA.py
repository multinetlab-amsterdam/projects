#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 15:30:27 2025

@author: ekoderman (e.koderman@amsterdamumc.nl)
status: final

REVIEW
reviewed by: Sebastien Dam
review date: 20260112

This script performs exploratory data analysis (EDA) on filtered clinical data of PRECOG subjects.
It generates HTML reports using ydata-profiling and creates a grouped KPS post-op distribution plot.

"""
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from ydata_profiling import ProfileReport

#Paths
clinical_path = Path("/folder/data/processed/clinical_data/clinical_data_filtered.csv")
overlap_subjs = Path("/folder/subj_ids/overlap/overlap_subjs_23092025.csv")
output_dir = Path("/folder/subj_ids/overlap")

#Load
merged_df = pd.read_csv(clinical_path)
subj_ids = pd.read_csv(overlap_subjs)

filtered_df = merged_df[merged_df['PRECOG_ID'].isin(subj_ids['PRECOG_ID'])]
filtered_df.drop(columns=['Unnamed: 0'], inplace=True)

#########  Generate an EDA overview report ####################
## 1) FUll dataset
profile = ProfileReport(
    filtered_df,
    title='PRECOG: Filtered subjects patient characteristics',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)
profile.to_file(output_dir / "patient_characteristics_overview.html")

## 2) survivors with 12months and above
survivors = filtered_df[
    (filtered_df['months_resection_death']>=12) |
    (filtered_df['months_resection_death'].isna())]

profile = ProfileReport(
    survivors,
    title='PRECOG: Filtered subjects survivors',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)
profile.to_file(output_dir / "patient_characteristics_survivors12m.html")

### 3) Patients who passed away
dead_before12m = filtered_df[filtered_df['months_resection_death'] < 12]
profile = ProfileReport(
    dead_before12m,
    title='PRECOG: Filtered subjects dead before 12 months',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)
profile.to_file(output_dir / "patient_characteristics_dead_before_12m.html")

### 4) Functionally dependent with KPS less than 70
func_dependence = survivors[(survivors['KPS_post-op_value'] < 70)]
profile = ProfileReport(
    func_dependence,
    title='PRECOG: Filtered subjects functional dependence',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)
profile.to_file(output_dir / "patient_characteristics_func_dependent.html")


### 5) FUnctionally independent with KPS 70 or more
func_independence = survivors[(survivors['KPS_post-op_value'] >= 70)]

profile = ProfileReport(
    func_independence,
    title='PRECOG: Filtered subjects functional independence',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)

profile.to_file(output_dir / "patient_characteristics_func_independent.html")

print(f'Made Profile reports (.html) and stored: {output_dir}')

######## Outcome distribution (KPS) plot with percentages ######
print(f'check nans in kps post-op: {filtered_df["KPS_post-op_value"].isna().sum()}')

conditions = [
    (filtered_df['KPS_post-op_value'] == 0),
    (filtered_df['KPS_post-op_value'].between(10,60)),
    (filtered_df['KPS_post-op_value'].between(70,100))
    ]

choices = [0,1,2]

filtered_df['KPS_post-op_grouped'] = np.select(conditions, choices, default=np.nan)
print(f'check nans in kps post-op: {filtered_df["KPS_post-op_grouped"].isna().sum()}')

# plot
sns.set(style='whitegrid')
ax = sns.histplot(filtered_df['KPS_post-op_grouped'], discrete=True)
total = len(filtered_df['KPS_post-op_grouped'])

# Annotate bars
for p in ax.patches:
    count = int(p.get_height())
    percent = (count / total) * 100
    x = p.get_x() + p.get_width() / 2
    y = p.get_height()
    ax.annotate(f'{count} ({percent:.1f}%)', (x, y), ha='center', va='bottom', fontsize=10)

# Custom tick labels
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['KPS = 0', 'KPS < 70', 'KPS ≥ 70'])

# Set labels and title
ax.set_xlabel('KPS post-op')
ax.set_ylabel('Count')
ax.set_title('KPS post-op grouped based on ranges')
plt.ylim([0,275]) #change to fit but keep it consistent with the plot below

plt.tight_layout()
#plt.show()

plt.savefig(output_dir / "review_KPS_postop_grouped_distribution.png", dpi=300, bbox_inches="tight")
print(f'Made KPS grouped distribution plot figure (.png) and stored: {output_dir}')
