#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 12:29:20 2025

@author: ekoderman (e.koderman@amsterdamumc.nl)
status: final

REVIEW
reviewed by: Sebastien Dam
review date: January 2026

As the gsi-rads report provides close to 200 features, this notebook aims to aggregate and extract the relevant ones based on expert opinion.

"""
#Imports
import pandas as pd
from pathlib import Path
from utils import validate_merged_df
import seaborn as sns
from ydata_profiling import ProfileReport

#Path
reports_path = Path('/folder/data/processed/gsi_reports/combined_gsi_report_20250916_150607.csv')
overlap_subjs_path = Path('/folder/subj_ids/overlap/overlap_subjs_23092025_IMAGO_ID.csv')
output_dir = Path('/folder/data/processed/gsi_reports')

#Load
gsi_reports = pd.read_csv(reports_path)
overlap = pd.read_csv(overlap_subjs_path)

#manually add XXX because of a re-ran gsi-rads after corrected resection date
gsi_IM1220 = pd.read_csv('/XXX/report.csv')
gsi_IM1220.insert(0, 'subject_id', 'XXX')
gsi_reports = pd.concat([gsi_reports, gsi_XXX], ignore_index=True)

#index
gsi_report_overlap = gsi_reports[gsi_reports['subject_id'].isin(overlap['IMAGO_ID'])]
gsi_report_overlap.drop(columns=['Unnamed: 0'], inplace=True)

#%%

sns.scatterplot(data = gsi_report_overlap, x='Volume in MNI (ml)', y = 'Volume original (ml)')
print(gsi_report_overlap['Volume in MNI (ml)'].mean())
print(gsi_report_overlap['Volume original (ml)'].mean())
# the per component volumes are in MNI so lets stick with that

#%% ### checks

validate_merged_df(gsi_report_overlap)

#%% ### laterality index merged

gsi_report_overlap['Laterality index'] = (gsi_report_overlap['Left laterality (%)'] - gsi_report_overlap['Right laterality (%)']) / 100
gsi_report_overlap['Laterality index'].describe()

# 1 is left, -1 right
#sns.histplot(gsi_report_overlap['Laterality index'],kde=True)

sns.histplot(
    gsi_report_overlap['Laterality index'],
    kde=True,        # add smooth density curve
    bins=30,         # adjust number of bins
    color="steelblue",
    edgecolor="white"
)

#%%  ### Yeo7 networks

# agreggate YEO7 cognitive networks & primary networks

# Cognitively Relevant Networks:
# Default Mode Network (DMN) → Memory, self-referential thinking, mind-wandering, and social cognition
# Frontoparietal Network (FPN) → Executive function, working memory, cognitive flexibility, and decision-making
# Salience/Ventral Attention Network (SVAN) → Attention switching, detecting important stimuli, and emotional processing
# Dorsal Attention Network (DAN) → Top-down attention control, visuospatial attention, and goal-directed actions

# Less Cognitive, More Sensory/Motor Networks:
# Somatomotor Network (SMN) → Motor control, sensory processing
# Limbic Network → Emotion regulation, autonomic function
# Visual Network → Visual perception and processing

schaefer7_columns = [col for col in gsi_report_overlap.columns if 'schaefer7' in col]
print(schaefer7_columns)

sensory_limbic = ['Schaefer7_schaefer7_visual_main_overlap', 'Schaefer7_schaefer7_somatomotor_main_overlap',
                  'Schaefer7_schaefer7_limbic_main_overlap']

gsi_report_overlap['Yeo7_sensory_limbic_overlap_sum'] = gsi_report_overlap[sensory_limbic].sum(axis=1)

cognitive = ['Schaefer7_schaefer7_salienceventralattention_main_overlap',
             'Schaefer7_schaefer7_dorsalattention_main_overlap',
             'Schaefer7_schaefer7_frontoparietalcontrol_main_overlap',
             'Schaefer7_schaefer7_default_main_overlap']

gsi_report_overlap['Yeo7_cognitive_overlap_sum'] = gsi_report_overlap[cognitive].sum(axis=1)

#%% #### MNI lobes aggregation

parietal = ['MNI__parietal_left_main_overlap', 'MNI__parietal_right_main_overlap']
temporal = ['MNI__temporal_left_main_overlap', 'MNI__temporal_right_main_overlap']
occipital = ['MNI__occipital_left_main_overlap', 'MNI__occipital_right_main_overlap']
frontal = ['MNI__frontal_left_main_overlap', 'MNI__frontal_right_main_overlap']

gsi_report_overlap['MNI_total_parietal_overlap_sum'] = gsi_report_overlap[parietal].sum(axis=1)
gsi_report_overlap['MNI_total_temporal_overlap_sum'] = gsi_report_overlap[temporal].sum(axis=1)
gsi_report_overlap['MNI_total_occipital_overlap_sum'] = gsi_report_overlap[occipital].sum(axis=1)
gsi_report_overlap['MNI_total_frontal_overlap_sum'] = gsi_report_overlap[frontal].sum(axis=1)

#%% ### subcortical white matter fiber tracts

#structures below are chosen by a neuro-surgeon (expert-based):

subcortical_cols = [
    'BCBCorpus_callosum_overlap',
    'BCBArcuate_Long_Segment_Left_overlap',
    'BCBArcuate_Long_Segment_Right_overlap',
    'BCBCortico_Spinal_Left_overlap',
    'BCBCortico_Spinal_Right_overlap',
    'BCBArcuate_Anterior_Segment_Right_overlap',
    'BCBArcuate_Anterior_Segment_Left_overlap',
    'BCBCingulum_Left_anterior_overlap',
    'BCBCingulum_Right_Anterior_overlap',
    'BCBSuperior_Londgitudinal_Fasciculus_III_Right_overlap',
    'BCBSuperior_Londgitudinal_Fasciculus_III_Left_overlap',
    'BCBSuperior_Londgitudinal_Fasciculus_II_Right_overlap',
    'BCBSuperior_Londgitudinal_Fasciculus_II_Left_overlap',
    'BCBSuperior_Londgitudinal_Fasciculus_I_Right_overlap',
    'BCBSuperior_Londgitudinal_Fasciculus_I_Left_overlap',
    'BCBArcuate_Posterior_Segment_Left_overlap',
    'BCBArcuate_Posterior_Segment_Right_overlap',
    'BCBFrontal_Aslant_tract_Right_overlap',
    'BCBFrontal_Aslant_Tract_Left_overlap',
    'BCBInferior_Fronto_Occipital_fasciculus_Left_overlap',
    'BCBInferior_Fronto_Occipital_fasciculus_Right_overlap',
    'BCBFronto_Striatal_Left_overlap',
    'BCBFronto_Striatal_Right_overlap',
    'BCBFrontal_Superior_Longitudinal_Left_overlap',
    'BCBFrontal_Superior_Longitudinal_Right_overlap',
    'BCBCingulum_Left_overlap',
    'BCBCingulum_Right_overlap',
    'BCBCingulum_Right_Posterior_overlap',
    'BCBCingulum_Left_posterior_overlap'
]

gsi_report_overlap['subcortical_overlap_sum'] = gsi_report_overlap[subcortical_cols].sum(axis=1)

#%% ## make a concise df

relevant_sum_columns = ['subject_id', 'Multifocality','subcortical_overlap_sum',
                    'MNI_total_parietal_overlap_sum','MNI_total_temporal_overlap_sum','MNI_total_frontal_overlap_sum',
                    'MNI_total_occipital_overlap_sum','Yeo7_cognitive_overlap_sum', 'Yeo7_sensory_limbic_overlap_sum',
                    'Laterality index', 'Volume in MNI (ml)', 'ExpectedResidualVolume (ml)', 'ResectionIndex']

concise_report_overlap_sum = gsi_report_overlap[relevant_sum_columns]

#%% ## make the profileReport

profile = ProfileReport(
    concise_report_overlap_sum,
    title='PRECOG: gsi reports concise overlap sum',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)

profile.to_file(output_dir / 'concise_gsi_report_overlap_sum_22122025.html')

#%% ## store csv

concise_report_overlap_sum.to_csv(output_dir / 'concise_gsi_report_overlap_sum_22122025.csv')
