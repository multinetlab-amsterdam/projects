#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 12:29:44 2025

@author: ekoderman (e.koderman@amsterdamumc.nl)
status: final

REVIEW
reviewed by: Sebastien Dam
review date: January 2026

This script creates a merged dataframe that contains clinical, volumetric, and radiomics features.
This is the dataframe that will be used as input for splitting the train/test split.

It performs the following dataframe modifications:
    - merges clinical & gsi-rads & volumetric outputs
    - drops irrelevant columns
    - performs exploratory data analysis (correlation matrix)
    - label mapping
    - creates time_to_event and event which are used in the next step (10_external_holdout_clinical_set.py)
    - renames column names to a more readable format

"""
# Imports
import pandas as pd
from pathlib import Path
from utils import validate_merged_df, convert_ids, calculate_month_difference, EDA_visualization
import seaborn as sns
from ydata_profiling import ProfileReport
import matplotlib.pyplot as plt
import numpy as np

# Paths
vol_components_path = Path("/folder/data/processed/volume_components/PRECOG_tumor_vol_components_final.csv")
clinical_path = Path("/folder/data/processed/clinical_data/clinical_data_filtered.csv")
gsi_path = Path("/folder/data/processed/gsi_reports/concise_gsi_report_overlap_sum_22122025.csv")
overlap_subjs_path = Path("/folder/subj_ids/overlap/overlap_subjs_23092025.csv")

output_dir = Path('/folder/data/processed/combined')

# Load
vol_df = pd.read_csv(vol_components_path)
overlap_ids = pd.read_csv(overlap_subjs_path)

clinical_df = pd.read_csv(clinical_path)
clinical_overlap = clinical_df[clinical_df["PRECOG_ID"].isin(overlap_ids['PRECOG_ID'])]

gsi_df = pd.read_csv(gsi_path)
gsi_df.drop(columns='Unnamed: 0', inplace=True)
gsi_df = convert_ids(gsi_df, 'subject_id','PRECOG_ID')
gsi_df.drop(columns = 'subject_id', inplace=True)

#%% CHeck and validate dataframes

print(set(clinical_overlap['PRECOG_ID']) == set(gsi_df['PRECOG_ID']))
print(set(clinical_overlap['PRECOG_ID']) == set(vol_df['PRECOG_ID']))
print(set(vol_df['PRECOG_ID']) == set(gsi_df['PRECOG_ID']))

validate_merged_df(clinical_overlap, key_col='PRECOG_ID')
validate_merged_df(vol_df, key_col='PRECOG_ID')
validate_merged_df(gsi_df, key_col='PRECOG_ID')

#%% merge dataframes into one and drop irrelevant columns

merged_df = (
    clinical_overlap
    .merge(vol_df, on='PRECOG_ID', how='inner')
    .merge(gsi_df, on='PRECOG_ID', how='inner')
)

drop_cols = ['Unnamed: 0_x', 'KPS_post-op:I/N', 'KPS_post-op_notes', 'KPS_post_op_date', 'KPS_pre-op:I/N',
             'KPS_pre-op_date', 'KPS_pre-op_notes', 'KPS_post-op_database', 'KPS_pre-op_database', 'rated_by',
             'Unnamed: 0_y', 'date_death']

merged_df.drop(columns = drop_cols, inplace=True)

#%% ### prepare vals to be ML ready and investigate input features

map_dict = {'M': 0, 'F': 1}
merged_df['Sex'].replace(map_dict, inplace=True)

cols_input_features = ['Sex', 'age_at_surgery', 'epilepsy_history', 'KPS_pre-op_value', 'T2 hyperintensity',
                       'Necrotic core', 'Enhancing component', 'Multifocality', 'subcortical_overlap_sum',
                       'MNI_total_parietal_overlap_sum', 'MNI_total_temporal_overlap_sum', 'MNI_total_frontal_overlap_sum',
                       'MNI_total_occipital_overlap_sum', 'Yeo7_cognitive_overlap_sum', 'Yeo7_sensory_limbic_overlap_sum',
                       'Laterality index', 'ExpectedResidualVolume (ml)', 'ResectionIndex']

input_df = merged_df[cols_input_features]

#### EDA - check distributions and missing vals
print(input_df.info())
print(input_df.describe(include='all').T)

EDA_visualization(input_df,'check')

#%%

sns.scatterplot(data=input_df,x='T2 hyperintensity', y='ExpectedResidualVolume (ml)')

#%% Laterality index & Resection Index inspection

# 1 is left, -1 right
sns.scatterplot(data=merged_df, x='Laterality index', y='ResectionIndex')
plt.grid()

top_5 = merged_df['Laterality index'].value_counts().nlargest(5).index
high_LI_df = merged_df[merged_df['Laterality index'].isin(top_5)]
high_LI_df = high_LI_df[['PRECOG_ID','Laterality index', 'ResectionIndex']]

all_LI_df = merged_df[['PRECOG_ID','Laterality index', 'ResectionIndex']]

convert_ids(high_LI_df, 'PRECOG_ID', 'IMAGO_ID')
convert_ids(all_LI_df, 'PRECOG_ID', 'IMAGO_ID')
convert_ids(gsi_df, 'PRECOG_ID', 'IMAGO_ID')

#%%

profile = ProfileReport(
    input_df,
    title='PRECOG: combined input features overview',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)
profile.to_file(output_dir / 'combined_input_features_22122025.html')

#%% Prepare target variables to be ML ready

# KPS post-op one column with 3 bins: 0 = dead before 12m, 1 = functionally dependent, 2 = functionally independent

merged_df['KPS_post-op_value'] = merged_df['KPS_post-op_value'].astype(int)

conditions = [merged_df['KPS_post-op_value'] == 0,
              (merged_df['KPS_post-op_value'] >= 10) & (merged_df['KPS_post-op_value'] < 70),
                merged_df['KPS_post-op_value'] >= 70]

values = [0,1,2]

merged_df['KPS_post_op_binned'] = np.select(conditions, values)
print('nans in KPS post-op binned target variable?')
print(merged_df['KPS_post_op_binned'].isna().sum())

#%% time_to_event & event

merged_df.rename(columns ={'vital_status':'event'},inplace=True)
merged_df['event'].replace({'passed away':1, 'alive':0}, inplace=True)
merged_df.rename(columns = {'months_resection_death':'time_to_event'}, inplace=True)

## add the date_vital_status to the time_to_event! - these are censored patients (still alive)
calculate_month_difference(merged_df, 'surgery_date', 'date_vital_status', 'months_resection_vital_status')
merged_df['time_to_event'].fillna(merged_df['months_resection_vital_status'], inplace=True)

#%% ## clean it up

#drop highly correlated features with (expected residual volume) and variables that are not longer needed (replaced by time_to_event)
merged_df.columns
drop_cols = ['date_vital_status', 'surgery_date', 'months_resection_kps', 'months_resection_vital_status', 'ExpectedResidualVolume (ml)', 'Volume in MNI (ml)']
merged_df.drop(columns=drop_cols, inplace=True)

#%% ## rename column names

merged_df = merged_df.rename(columns={
    'T2 hyperintensity': 'T2 hyperintensity volume',
    'KPS_pre-op_value': 'Preoperative KPS',
    'age_at_surgery': 'Age at resection',
    'epilepsy_history': 'Preoperative seizures',
    'Necrotic core': 'Necrotic core volume',
    'Enhancing component': 'Enhancing component volume',
    'subcortical_overlap_sum': 'Subcortical overlap',
    'MNI_total_parietal_overlap_sum': 'Parietal overlap',
    'MNI_total_temporal_overlap_sum': 'Temporal overlap',
    'MNI_total_frontal_overlap_sum': 'Frontal overlap',
    'MNI_total_occipital_overlap_sum': 'Occipital overlap',
    'Yeo7_cognitive_overlap_sum': 'Cognitive networks overlap',
    'Yeo7_sensory_limbic_overlap_sum': 'Sensory limbic networks overlap',
    'ResectionIndex': 'Resection index'})

print(merged_df.columns)

#%% ## store this intermediate full df = full clean csv that includes all ID+features+targets

merged_df.to_csv(output_dir / 'combined_ID_input_target_features_22122025.csv', index=False)

