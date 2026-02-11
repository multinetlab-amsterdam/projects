#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 17:57:06 2025

@author: ekoderman (e.koderman@amsterdamumc.nl)
status: final

REVIEW
reviewed by: Sebastien Dam
review date: 20260112

This script loads the final selection of subjects based on clinical data and a selection of subjects based on brain imaging data.
It checks for overlap between them and makes the final selection of subejcts for whom we have all the data available.
It also investigates the subjects that ended up outside the overlap and why.

It stores the subjects IDs from the overlap, meaning the clinical and the brain imaging data is available for these subjects.

"""
# Imports
import pandas as pd
from pathlib import Path
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import json
from utils import convert_exclusion_dict

# Paths
brain_imaging_path = Path("/folder/subj_ids/brain_imaging/PRECOG_enh_subjs_final.csv")
clinical_path = Path("/folder/data/processed/clinical_data/clinical_data_filtered.csv")
output_dir = Path("/folder/subj_ids/overlap")

excluded_brain_path = Path("/folder/subj_ids/brain_imaging/brain_imaging_excluded_subjs_PRECOG.json")
excluded_clinical_path = Path("/folder/subj_ids/clinical/clinical_excluded_subjs.json")
non_enh_path = Path("/folder/subj_ids/brain_imaging/PRECOG_non_enh_subjs.csv")

# Load
brain_imaging_data = pd.read_csv(brain_imaging_path)
clinical_data = pd.read_csv(clinical_path)
non_enh_subjs = pd.read_csv(non_enh_path)

# Load JSON into a Python dict
with open(excluded_clinical_path, "r") as f:
    excluded_clinical = json.load(f)

with open(excluded_brain_path, "r") as f:
    excluded_brain = json.load(f)

# Ensure PRECOG_ID columns are strings and stripped of whitespace
brain_imaging_data['PRECOG_ID'] = brain_imaging_data['PRECOG_ID'].astype(str).str.strip()
clinical_data['PRECOG_ID'] = clinical_data['PRECOG_ID'].astype(str).str.strip()

brain_imaging_subjs = set(brain_imaging_data['PRECOG_ID'])
clinical_subjs = set(clinical_data['PRECOG_ID'])

overlap = brain_imaging_subjs & clinical_subjs

# Plot the Venn diagram
plt.figure(figsize=(6, 6))
venn2([clinical_subjs, brain_imaging_subjs],
      set_labels=(f'Clinical\n N={len(clinical_subjs)}',f'Brain Imaging\n N={len(brain_imaging_subjs)}', ))

plt.title('Overlap of Subjects')
plt.show()

#%% ### investigate the ones outside the overlap for subjs unique to brain_imaging

brain_only = brain_imaging_subjs - clinical_subjs
# they were excluded from the clinical_subjs

unique_brain = convert_exclusion_dict(excluded_clinical, brain_only)
duplicates = unique_brain[unique_brain.duplicated(subset='subj_id', keep=False)]
unique_brain = unique_brain.drop_duplicates(subset='subj_id', keep='first')

# this explains 73 subjs, what about the other 3?

unique_set = set(unique_brain['subj_id'])
brain_outsiders = brain_only - unique_set

# they were excluded offline in excel sheets for reasons documented elsewhere (protected due to patient sensitive data).

#%% ### investigate the ones outside the overlap for subjs unique to clinical

clinical_only = clinical_subjs - brain_imaging_subjs # N = 135
# for clinical_only = are all these non-enhancing subjs? - double-check to confirm

overlap_clinical_non_enh = clinical_only & set(non_enh_subjs['PRECOG_ID']) # N = 120
unique_clinical_brain_reasons = convert_exclusion_dict(excluded_brain, clinical_only) # N = 14
unique_set = set(unique_clinical_brain_reasons['subj_id'])


#%% ## save
overlap_df = pd.DataFrame(overlap, columns=['PRECOG_ID'])
overlap_df.to_csv(output_dir / 'review_overlap_subjs_23092025.csv')
