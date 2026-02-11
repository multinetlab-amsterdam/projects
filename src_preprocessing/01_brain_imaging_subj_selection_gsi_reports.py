#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 11:28:21 2025

@author: ekoderman (e.koderman@amsterdamumc.nl)
status: final

REVIEW
reviewed by: Sebastien Dam
review date: 20260112
    
Manual quality assessment of all scans and tumor masks was performed offline in an excel file. This scripts extracts subjects that passed
the QA. It also filters for subjects that don't have a cerebellar tumor. This information is provided in the gsi-rads reports so it loads
the reports across subjects.

This script consists of several steps:
    1) Filters for subjects that have passed the 1) QA check and are 2) enhancing
    2) Calls the script 00_combine_gsi_rads_reports.py to obtain gsi-rads for the desired subjects
    3) Excludes subjects with 3) cerebellar overlap as this info is provided in the gsi-rads reports

Outputs:
    1) 2 csv files with subject IDs of enhancing and non-enhancing subjects
    2) a dataframe with gsi-rads report for the enhancing subjects
    3) a json report on subject exclusion and reasons

"""

### Imports
import pandas as pd
from ydata_profiling import ProfileReport
from pathlib import Path
import importlib
combine_and_save = getattr(importlib.import_module("00_combine_gsi_rads_reports"), "combine_and_save")
from utils import add_excluded_subjects
import os
import json

###  Paths 
raw_data_path = Path("/folder/data/raw/brain_imaging_segmentation_qa.csv")
gsi_input_dir = Path("/path/to/tumor_masks")
gsi_output_dir = Path("/folder/data/processed/gsi_reports")
subjs_output_dir = Path("/folder/subj_ids/brain_imaging")

### Load data 
imago_seg_file = pd.read_csv(raw_data_path)

#########  Generate an EDA overview report ####################
profile = ProfileReport(
    imago_seg_file,
    title='PRECOG: Brain Imaging Overview before subj selection',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)
profile.to_file(subjs_output_dir / "brain_imaging_overview.html")

######### Filter for subjects that have passed the 1) QA check and are 2) enhancing ###########
subjs_enh = imago_seg_file[
    (imago_seg_file['SEGMENTATION'] == 'success') &
    (imago_seg_file['ENH Y/N'] == True) &
    (imago_seg_file['SEGMENTATION QA'] == 'Sufficient')
]

# we store excluded subjects to keep track
subjs_seg_failed = imago_seg_file[imago_seg_file['SEGMENTATION'] == 'failed']
subjs_seg_insufficient = imago_seg_file[imago_seg_file['SEGMENTATION QA'] == 'Insufficient']

excluded_subjs_dict = {}
add_excluded_subjects(excluded_subjs_dict, 'segmentation failed', subjs_seg_failed['IMAGO SUBJECT NR'])
add_excluded_subjects(excluded_subjs_dict, 'segmentation insufficient', subjs_seg_insufficient['IMAGO SUBJECT NR'])

# Non-enhancing subjects are not the focus of this study so the segmentation QA was not manually checked
subjs_non_enh = imago_seg_file[
    (imago_seg_file['SEGMENTATION'] == 'success') &
    (imago_seg_file['ENH Y/N'] == False)]

# Extract just the IDs 
subjs_enh_ID = subjs_enh['IMAGO SUBJECT NR']
subjs_non_enh_ID = subjs_non_enh['IMAGO SUBJECT NR']

######### 3) Filter for subjects that do not have a cerebellar overlap #########

# first we need to import gsi-rads report to get insight into the cerebellar overlap
combined_df, csv_path = combine_and_save(gsi_input_dir, gsi_output_dir)
print(f'Duplicates after combining gsi-rads reports: {combined_df.duplicated().sum()}')
enh_gsi = combined_df[combined_df['subject_id'].isin(subjs_enh_ID)]

# Manually add the subject IM1220 because its not present in the gsi_input_dir
gsi_IM1220 = pd.read_csv('/data/IM1220/report.csv')
gsi_IM1220.insert(0, 'subject_id', 'IM1220')
final_gsi = pd.concat([enh_gsi, gsi_IM1220], ignore_index=True)

# Combine left and right cerebellum overlap
cerebellum_overlap = final_gsi[['subject_id', 'MNI__cerebellum_left_main_overlap', 'MNI__cerebellum_right_main_overlap']]
cerebellum_overlap.loc[:, ['total_cerebellum_overlap']] = cerebellum_overlap['MNI__cerebellum_left_main_overlap'] + cerebellum_overlap['MNI__cerebellum_right_main_overlap']

nonzero_vals = cerebellum_overlap['total_cerebellum_overlap'][cerebellum_overlap['total_cerebellum_overlap'] > 0]
mean_val = nonzero_vals.mean()
std_val = nonzero_vals.std()
k = 2  # how many standard deviations

cerebellum_overlap.loc[:, ['is_outlier']] = (
    (cerebellum_overlap['total_cerebellum_overlap'] - mean_val).abs() > k * std_val
)

outliers = cerebellum_overlap[cerebellum_overlap['is_outlier']]
add_excluded_subjects(excluded_subjs_dict, 'cerebellum overlap', outliers['subject_id'])

filtered_final_gsi = final_gsi[~final_gsi['subject_id'].isin(outliers['subject_id'])].copy()

######### Saving values #########

# Save subject IDs
final_subjs = filtered_final_gsi['subject_id']

os.makedirs(subjs_output_dir, exist_ok=True)
final_subjs.to_csv(subjs_output_dir / 'IMAGO_enh_subjs_final.csv', index=False)
subjs_non_enh_ID.to_csv(subjs_output_dir / 'IMAGO_non_enh_subjs.csv', index=False)

# Save the dictionary to a JSON file
with open(subjs_output_dir / 'brain_imaging_excluded_subjs.json', 'w') as f:
    json.dump(excluded_subjs_dict, f, indent=4)

### Save gsi-rads reports based on subject selection
filtered_final_gsi.to_csv(gsi_output_dir / 'gsi_reports_final_subjs.csv')