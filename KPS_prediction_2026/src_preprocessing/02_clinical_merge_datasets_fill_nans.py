#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 12:23:06 2025

@author: ekoderman (e.koderman@amsterdamumc.nl)
status: final

REVIEW
reviewed by: Sebastien Dam
review date: 20260112

This notebook combines 2 excel files containing clinical information with the main outcomes (KPS post-op) that were created offline.
For additional clinical information it also loads an external database.

It performs the following steps:
    - merges first and second batch
    - drops unneeded columns
    - checks and normalizes date format
    - fills missing data from the external clinical file (dbtr)
    - normalizes categorical columns
    - cleans KPS columns (formatting + if KPS is a range it takes the max value)
    - calculates months different (resection-death and resection-KPS)
    - final validation: checks for inconsistencies and duplicates

It outputs a 1) profile report and 2) merged clinical dataframe that contains only relevant clinical information from all these databases.
"""

#imports
import pandas as pd
from ydata_profiling import ProfileReport
from pathlib import Path
from utils import validate_merged_df, check_date_column, fill_missing_from_dbtr, calculate_month_difference, clean_kps_column
import numpy as np

###  Paths 
first_batch_path = Path("/folder/data/raw/KPS_overview.xlsx")
second_batch_path = Path("/folder/data/raw/KPS_overview_IMAGO.xlsx")
dbtr = Path("folder/data/raw/clinical_file.csv")
output_dir = Path("/folder/subj_ids/clinical")

# Load data
first_batch = pd.read_excel(first_batch_path, sheet_name="MAIN")
second_batch = pd.read_excel(second_batch_path, sheet_name="Final")
dbtr_df = pd.read_csv(dbtr, encoding="latin-1")

# Merge together
merged_df = pd.merge(first_batch, second_batch, on="PRECOG_ID", how='outer')

#%%

print(merged_df.columns)

# drop columns that we don't need
drop_columns = ['ziekenhuis radiotherapie', 'ziekenhuis chemotherapie',
            'Pre_op_MRI', 'grade', 'enhancing', 'days_diff_surgery_KPS_post-op (274-457 = 9-15m)',
            'months_diff_surgery_death', 'months_diff_surgery_KPS_post-op', 'months_difference_surgery_death',
            'months_surgery_to_postop_KPS']

# we will re-calculate the differences in month to ensure these are all done consistently across the subjects
merged_df.drop(columns=drop_columns, inplace=True)
print(merged_df.columns)

#%% # Test for potentially spurious values
print('Test 1:')
validate_merged_df(merged_df, key_col="PRECOG_ID")

# Address the [FAIL] with conflicting column paris found (x,y)
for col in ['KPS_post-op:I/N', 'KPS_post-op_notes', 'KPS_post-op_value', 'KPS_post_op_date', 'KPS_pre-op:I/N',
            'KPS_pre-op_date', 'KPS_pre-op_notes', 'date_vital_status', 'epilepsy_history', 'surgery_date', 'vital_status',
            'KPS_post-op_database', 'KPS_pre-op_database','KPS_pre-op_value', 'date_death', 'rated_by']:
    merged_df[col] = merged_df[col+'_x'].combine_first(merged_df[col+'_y']) #This keeps _x values unless they’re null, in which case it takes _y
    merged_df.drop([col+'_x', col+'_y'], axis=1, inplace=True)

# Test again
print('Test 2:')
validate_merged_df(merged_df, key_col="PRECOG_ID")

#%% ##### age_at_resection

# ensure resection_date is consistent - first check what we are dealing with
print('surgery date:')
check_date_column(merged_df, 'surgery_date')

# Fill missing age
merged_df = fill_missing_from_dbtr(
    merged_df, dbtr_df,
    target_col='age_at_surgery',
    dbtr_col='leeftijd ten tijde van diagnose',
    sentinel_values = ['999', 999]
)

#%% ###### sex

print(merged_df['Sex'].value_counts(dropna=False))
# Fill missing sex
merged_df = fill_missing_from_dbtr(
    merged_df, dbtr_df,
    target_col='Sex',
    dbtr_col='geslacht',
    sentinel_values = ['999', 999]
)

# check if values match
print(merged_df['Sex'].value_counts()) #they don't - let's normalize them

sex_map = {
    'man': 'M',
    'vrouw': 'F',
    '1': 'M' #this one person was additionally checked in the Multinet database
    }

merged_df['Sex'] = merged_df['Sex'].replace(sex_map)

# check if values match
print(merged_df['Sex'].value_counts(dropna=False))

#%% #### seizure history

print('epilepsy history initally (nans):')
print(merged_df['epilepsy_history'].isna().sum())

# check if values match
print(merged_df['epilepsy_history'].value_counts())

merged_df['epilepsy_history'] = merged_df['epilepsy_history'].replace(['999',999], np.nan)

# check if values match
print(merged_df['epilepsy_history'].value_counts())

#%% #### vital status

# check vital status and whether or not the date_death is missing info or the patient is still alive
print('first check: vital status')
print(merged_df['vital_status'].value_counts())

vital_status_map = {
    'overleden': 'passed away',
    'in leven': 'alive',
    '999': np.nan,
    999: np.nan
    }

merged_df['vital_status'] = merged_df['vital_status'].replace(vital_status_map)

print('second check: vital status')
print(merged_df['vital_status'].value_counts())

#%% ### date_death

print('missing vital status:')
print(merged_df['vital_status'].isna().sum())

merged_df['date_death'] = merged_df['date_death'].replace([999, '999'], pd.NaT)
merged_df['date_death'] = pd.to_datetime(merged_df['date_death'], errors='coerce')
merged_df['date_death'].value_counts()

print('missing death_date:')
print(merged_df['date_death'].isna().sum())

# nans are OK - they mean that the patient is alive
#%% #### months: resection-death

merged_df = calculate_month_difference(
    merged_df,
    start_col='surgery_date',
    end_col='date_death',
    new_col='months_resection_death'
)

print('Months resection death nans:')
print(merged_df['months_resection_death'].isna().sum())

#%% #### KPS

merged_df.columns

cols_to_clean = ['KPS_pre-op_value', 'KPS_post-op_value']
for col in cols_to_clean:
    merged_df[col] = merged_df[col].replace([999, '999'], np.nan)

cols_to_date = ['KPS_post_op_date', 'KPS_pre-op_date']
for col in cols_to_date:
    merged_df[col] = merged_df[col].replace([999, '999'], pd.NaT)
    merged_df[col] = pd.to_datetime(merged_df[col], errors='coerce')

# take the highest value of the patients that have either above or below 70 KPS status
merged_df = clean_kps_column(merged_df, 'KPS_pre-op_value')
merged_df = clean_kps_column(merged_df, 'KPS_post-op_value')

#%% ### KPS months calculation

merged_df = calculate_month_difference(
    merged_df,
    start_col='surgery_date',
    end_col='KPS_post_op_date',
    new_col='months_resection_kps'
)

print(merged_df['months_resection_kps'].describe())

#%% #### put all 999 values to nans

affected_cols = [col for col in merged_df.columns if merged_df[col].isin([999, '999']).any()]
print("Columns containing 999:", affected_cols)

merged_df = merged_df.replace([999, '999'], np.nan)

#%% ### final validation test

validate_merged_df(merged_df, key_col="PRECOG_ID")
duplicates = merged_df[merged_df.duplicated(subset='PRECOG_ID', keep=False)]
merged_df = merged_df.drop_duplicates(subset='PRECOG_ID', keep='first')

print('Final check:')
validate_merged_df(merged_df, key_col="PRECOG_ID")

#%% #########  Generate an EDA overview report ####################
exploratory_cols = ['epilepsy_history', 'months_resection_death', 'months_resection_kps',
                    'KPS_post-op_value', 'KPS_pre-op_value', 'age_at_surgery', 'Sex', 'vital_status']

exploratory_df = merged_df[exploratory_cols]

profile = ProfileReport(
    exploratory_df,
    title='PRECOG: clinical data EDA before subj selection',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)

#%% ## Save

merged_df.to_excel(output_dir / 'clinical_info_all_subjs.xlsx', index=False)
profile.to_file(output_dir / 'clinical_info_overview.html')