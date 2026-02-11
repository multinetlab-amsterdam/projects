#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 16:01:59 2025

@author: ekoderman (e.koderman@amsterdamumc.nl)
status: final

REVIEW
reviewed by: Sebastien Dam
review date: January 2026

This script takes as input the desired subject ID list (overlap_subjs) and the directory with the paths to these subjects' scans.
It extracts the tumor mask in MNI that was used as input to the gsi-rads algorithm.
Then it calculates the volume in mL for each of the tumor mask component (enhancing, T2 hyperintensity, and necrotic core).

NB: This script automatically outputs an interim dataframe of these volume components. This dataframe needs to 
be manually corrected in the next preprocessing step for devious subjects (07_tumor_volume_components_corrections.py).

"""
from pathlib import Path
import pandas as pd
import json
import glob
import os
from utils import process_subject_component_volumes
from joblib import Parallel, delayed

#Paths
overlap_subjs_path = Path("/folder/subj_ids/subjs_overlap/overlap_subjs_23092025_IMAGO_ID.csv")
json_dir_base = Path('/path/to/gsi_rads_jobs.json')
output_dir= Path('/folder/data/processed/volume_components')

#Load
subj_ids = pd.read_csv(overlap_subjs_path)
subj_ids_set = set(subj_ids['IMAGO_ID'].tolist())  # convert to set for fast lookup

# find all JSON files under that base directory
json_files = sorted(glob.glob(os.path.join(json_dir_base, "*.json")))

if not json_files:
    raise FileNotFoundError(f"No json files found in {json_dir_base}")
    
records = []
for in_json in json_files:
    with open(in_json, 'r') as file:
        data = json.load(file)

    subject_id = data.get("subject")  # e.g. IM1134

    # Skip if subject not in CSV
    if subject_id not in subj_ids_set:
        continue

    # build expected search path for this subject
    subject_dir = f"/path/to/gsi_rads_jobs/{subject_id}/preop"
    search_pattern = os.path.join(subject_dir, "wdir", "*", "gsi_rads", "*", "registration", "input_segmentation_to_MNI.nii.gz")

    seg_files = glob.glob(search_pattern)
    if not seg_files:
        print(f"No input_segmentation_to_MNI.nii.gz found for {subject_id}")
        seg_file = None
    else:
        seg_file = seg_files[0]  # take first match

    records.append({
        "subject": subject_id,
        "segmentation_mni": seg_file,
        "json_file": in_json
    })

# turn into dataframe for easy handling
df_paths = pd.DataFrame(records)

# ---- run in parallel ----
if __name__ == "__main__":
    results = Parallel(n_jobs=2)(delayed(process_subject_component_volumes)(row) for _, row in df_paths.iterrows())
    volume_df = pd.DataFrame(results)
    print(volume_df.head())
    volume_df.to_csv(output_dir / 'overlap_vol_components_MNI_input_interim.csv', index=False)






