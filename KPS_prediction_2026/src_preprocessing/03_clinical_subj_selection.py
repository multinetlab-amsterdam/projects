#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 16:58:36 2025

@author: ekoderman (e.koderman@amsterdamumc.nl)
status: final

REVIEW
reviewed by: Sebastien Dam
review date: 20260112

This script performs filtering for clinical data based on the merged clinical dataframe.
It excludes:
    1) subjects with missing vital status
    2) subjects with missing KPS post-op in the desired month range (9-15 months)
    3) subjects with missing KPS post-op because there wasn't enough info in the clinical records for KPS inference

It also checks for overlapping subjects between these 3 groups (one subject can be missing for many reasons).
    
It stores:
    1) json reports with the excluded subejcts and reasons for exclusion
    2) final clinical dataframe with the desired subjects
"""

### Imports
from pathlib import Path
import pandas as pd
from utils import add_excluded_subjects
import json

###  Paths 
data_path = Path("/folder/subj_ids/clinical/clinical_info_all_subjs.xlsx")
output_dir_subjs = Path("/folder/subj_ids/clinical")
output_dir_data = Path("/folder/data/processed/clinical_data")

###  Load 
merged_df = pd.read_excel(data_path, sheet_name='Sheet1')
merged_df_filtered = merged_df.copy()
excluded_subjs_dict = {}

### missing vital status
missing_vital_status = merged_df[merged_df['vital_status'].isna()]
merged_df_filtered = merged_df_filtered[~merged_df_filtered['PRECOG_ID'].isin(missing_vital_status['PRECOG_ID'])]
add_excluded_subjects(excluded_subjs_dict, 'missing vital status', missing_vital_status['PRECOG_ID'])

### KPS_resection_motnhs outside of the month range
outside_month_range = merged_df[(merged_df['months_resection_kps']<9) | (merged_df['months_resection_kps']>15)]
merged_df_filtered = merged_df_filtered[~merged_df_filtered['PRECOG_ID'].isin(outside_month_range['PRECOG_ID'])]
add_excluded_subjects(excluded_subjs_dict, 'kps post-op outside the month range', outside_month_range['PRECOG_ID'])

#### missing KPS post-op value because not enough info to infer it from
not_enough_info = merged_df[merged_df['KPS_post-op_value'].isna()]
merged_df_filtered = merged_df_filtered[~merged_df_filtered['PRECOG_ID'].isin(not_enough_info['PRECOG_ID'])]
add_excluded_subjects(excluded_subjs_dict, 'not enough info for KPS post-op inference', not_enough_info['PRECOG_ID'])

### check for overlapping subjs

# --- Extract sets of IDs from your dict ---
excluded_sets = {
    reason: set(info['subjs_id'])
    for reason, info in excluded_subjs_dict.items()
}

# 1) Total unique excluded subjects
all_excluded_ids = set().union(*excluded_sets.values())
print(f"Total unique excluded subjects: {len(all_excluded_ids)}")

# 2) Check overlaps between groups
print("\nOverlap between exclusion groups:")
reasons = list(excluded_sets.keys())
for i in range(len(reasons)):
    for j in range(i+1, len(reasons)):
        r1, r2 = reasons[i], reasons[j]
        overlap = excluded_sets[r1] & excluded_sets[r2]
        if overlap:
            print(f"  - {r1} ∩ {r2}: {len(overlap)} overlapping subjects")

add_excluded_subjects(excluded_subjs_dict, 'total unique excluded subjects', all_excluded_ids)

### save 
# Save the dictionary to a JSON file
with open(output_dir_subjs / 'clinical_excluded_subjs.json', 'w') as f:
    json.dump(excluded_subjs_dict, f, indent=4)
    
merged_df_filtered.to_csv(output_dir_data / 'review_clinical_data_final_subjs.csv')
