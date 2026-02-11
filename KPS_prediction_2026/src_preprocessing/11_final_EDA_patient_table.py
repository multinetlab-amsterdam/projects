#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 16:34:41 2025

@author: ekoderman

status: final

REVIEW
reviewed by: Sebastien Dam
review date: January 2026

This script performs detailed exploratory data analysis (EDA) using the package ProfileReport.
This package outputs a .html object which can be investigated further in a browser and shows all summary statistics
with distribution plots for each variable.

4 reports are created:
    - KPS = 0: group of patients who died before 12 months
    - KPS = 10 - 60: group of functionally dependent patients
    - KPS = 70 - 100: group of functionally independent patients
    - Overall subjects: all patients

These reports are used to construct the Patient characteristics table.

"""
import pandas as pd
from pathlib import Path
import numpy as np
from ydata_profiling import ProfileReport

output_dir = Path("/folder/subj_ids/overlap")
df_path = Path("/folder/combined_df_strata_groups_22122025.csv")
combined_df = pd.read_csv(df_path)

drop_cols = ['Unnamed: 0', 'strata', 'survival_categories', 'source']
combined_df.drop(columns=drop_cols, inplace=True)

## make 3 classed outcome
conditions = [combined_df['KPS_post-op_value'] == 0,
              combined_df['KPS_post-op_value'] < 70,
              combined_df['KPS_post-op_value'] >= 70]
values = [0,1,2]

combined_df['KPS_post_op_binned'] = np.select(conditions,values)

kps_0 = combined_df[combined_df['KPS_post_op_binned'] == 0]
kps_1 = combined_df[combined_df['KPS_post_op_binned'] == 1]
kps_2 = combined_df[combined_df['KPS_post_op_binned'] == 2]

## 1) KPS = 0
profile = ProfileReport(
    kps_0,
    title='Patient characterists: KPS = 0',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)

profile.to_file(output_dir / "all_features_kps_0.html")

### 2) KPS = 10 - 60
profile = ProfileReport(
    kps_1,
    title='Patient characterists: KPS = 10-60',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)

profile.to_file(output_dir / "all_features_kps_10-60.html")

### 3) KPS = 70 - 100
profile = ProfileReport(
    kps_2,
    title='Patient characterists: KPS = 70-100',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)

profile.to_file(output_dir / "all_features_kps_70-100.html")

### 4) Overall (N = 552))
profile = ProfileReport(
    combined_df,
    title='Patient characterists: Overall (N = 552)',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)

profile.to_file(output_dir / "all_features_overall_patients.html")

print(f'Made Profile reports (.html) and stored: {output_dir}')
