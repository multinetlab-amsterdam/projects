#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 10:36:25 2025

@author: ekoderman (e.koderman@amsterdamumc.nl)
status: final

REVIEW
reviewed by: Sebastien Dam
review date: January 2026

This script is more like a notebook (cell-based used in Spyder). It takes as input:
    1) subject list with IDs (csv))
    2) dataframe with interim volume components in mL
    3) file with paths to patient specific scan paths (json)
    
Then it manually corrects 2 subjects because different paths need to be used than for the rest (they were manually corrected outside of this python workflow)
It also creates a Profile Report for the final volume dataframe and stores them in the output_dir provided.

"""
from pathlib import Path
import pandas as pd
import nibabel as nib
import numpy as np
from ydata_profiling import ProfileReport

#Paths
overlap_subjs_path = Path("/folder/subj_ids/overlap/overlap_subjs_23092025.csv")
volume_components_path = Path("/folder/data/processed/volume_components/PRECOG_overlap_vol_components_MNI_input_interim.csv")
json_dir_base = Path('/path/to/gsi_rads_jobs/preop')
output_dir = Path("/folder/data/processed/volume_components")

#Load
subj_ids = pd.read_csv(overlap_subjs_path)
vol_components = pd.read_csv(volume_components_path)
#%% #### add subject XXX
# not in the initial pool of subjects because of a wrong resection date

unique_subj = subj_ids[~subj_ids['PRECOG_ID'].isin(vol_components['PRECOG_ID'])]
# the extra is XXX  - we need to add this one manually because initially the wrong resection date was used

# make it nice and appendable to the main df
subject_id = unique_subj['PRECOG_ID']
subject_id = subject_id.iloc[0] if hasattr(subject_id, "iloc") else subject_id
sub = str(subject_id)

# retrieve the segmentation
subj_IM1220_path = Path('/XXX/registration/input_segmentation_to_MNI.nii.gz')
img = nib.load(subj_XXX_path)
data = img.get_fdata()
voxel_size = np.prod(img.header.get_zooms()[:3])  # voxel volume (mm³)

unique_labels = np.unique(data)
volume_dict = {"PRECOG_ID": sub}

#calculate components
for label in unique_labels:
    if label == 0:  # skip background
        continue
    num_voxels = np.sum(data == label)
    volume_mm3 = num_voxels * voxel_size
    volume_cm3 = round(volume_mm3 / 1000, 2)  # convert to cm³
    volume_dict[f"Label_{int(label)}"] = volume_cm3

# convert to DataFrame and append as the 552nd row
row_df = pd.DataFrame([volume_dict])
vol_components = pd.concat([vol_components, row_df], ignore_index=True)

#%% ### YYY - different segmentation mask
# the enhancing part is much smaller in the corrected version than in the initial volume mask

# correct the volume components of YYY to include the corrected segmentation file
print(f"initial calculation of YYY: {vol_components[vol_components['PRECOG_ID'] == 'YYY']}")
# MNI
corrected_YYY_path = Path('/YYY/corrected_tumor_mask/input_segmentation_to_MNI_YYY.nii.gz')

sub = 'YYY'
img = nib.load(corrected_P0267_ath)
data = img.get_fdata()
voxel_size = np.prod(img.header.get_zooms()[:3])  # voxel volume (mm³)

unique_labels = np.unique(data)
volume_dict = {"PRECOG_ID": sub}

#calculate components
for label in unique_labels:
    if label == 0:  # skip background
        continue
    num_voxels = np.sum(data == label)
    volume_mm3 = num_voxels * voxel_size
    volume_cm3 = round(volume_mm3 / 1000, 2)  # convert to cm³
    volume_dict[f"Label_{int(label)}"] = volume_cm3

# convert to DataFrame
row_df = pd.DataFrame([volume_dict])
# drop the old YYY row
vol_components = vol_components[vol_components["PRECOG_ID"] != sub]
# append the new one
vol_components = pd.concat([vol_components, row_df], ignore_index=True)
print(f"calculation of YYY after correction: {vol_components[vol_components['PRECOG_ID'] == 'YYY']}")

#%% Rename labels into actual components

rename_cols = {'Label_1': 'T2 hyperintensity', 'Label_2': 'Necrotic core', 'Label_3': 'Enhancing component'}
vol_components = vol_components.rename(columns=rename_cols)

profile = ProfileReport(
    vol_components,
    title='PRECOG: per component volume report',
    explorative=True,        # adds extra visuals
    minimal=True              # speeds up for large datasets
)
profile.to_file(output_dir / 'per_component_volume_report.html')

#%%
# save
vol_components.to_csv(output_dir / 'PRECOG_tumor_vol_components_final.csv')
