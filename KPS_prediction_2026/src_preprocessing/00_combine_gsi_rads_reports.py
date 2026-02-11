#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 11:29:38 2025

@author: ekoderman (e.koderman@amsterdamumc.nl)
status: final

REVIEW
reviewed by: Sebastien Dam
review date: 20260109

The gsi_rads toolbox (for creating automatic reports on the tumor masks) outputs the per-subject report into the per-subject folder.
This script aims to aggregate the per-subject reports into a one dataframe that contains these reports across subjects.    

Run in a command line the following: python combine_gsi_reports.py 'base_directory_of_gsi_reports_across_subjects' 'desired_output_dir' 
outputs one dataframe with combined gsi_rads reports across subjects

This script currently only checks for preop folders. It extracts the most recent gsi-rads report created.
It can be used as a stand alone script for any aggregation of gsi-rads report or called as a function as it's
used in the next preprocessing step (01_brain_imaging_subj_selection_gsi_reports.py)

"""
import pandas as pd
import argparse
import os
from datetime import datetime
from pathlib import Path
from tqdm import tqdm # to be removed

def collect_report_files(base_path):
    report_files = []

    for root, _, files in tqdm(os.walk(base_path)):
        if "report.csv" not in files:
            continue

        # Check if this is a 'preop' folder
        parts = Path(root).parts
        if len(parts) < 5 or parts[-5] != "preop":
            continue  # skip anything not under preop

        parent_folder = os.path.dirname(root)
        subfolders = [f for f in os.listdir(parent_folder)
                      if os.path.isdir(os.path.join(parent_folder, f))]

        # Find the most recent timestamp folder
        most_recent_timestamp = None
        most_recent_folder = None
        for subfolder in subfolders:
            try:
                timestamp = datetime.strptime(subfolder, "%d%m%Y_%H%M%S")
                if most_recent_timestamp is None or timestamp > most_recent_timestamp:
                    most_recent_timestamp = timestamp
                    most_recent_folder = os.path.join(parent_folder, subfolder)
            except ValueError:
                continue

        # Append only **one report per subject**
        if most_recent_folder:
            report_path = os.path.join(most_recent_folder, "report.csv")
            report_files.append(report_path)

    return report_files

def combine_reports(base_path):
    """Combine all report.csv files into one DataFrame, adding a subject_id column."""
    report_files = collect_report_files(base_path)
    df_list = []

    for file in tqdm(report_files):
        # Extract subject_id from folder structure (pre-filtering version)
        # Example: subject_id = file.split(os.sep)[-7]
        subject_id = file.split(os.sep)[-7]
        df = pd.read_csv(file)
        df.insert(0, "subject_id", subject_id)
        df_list.append(df)

    # SD: perhaps add a try/raise block to check if df_list is empty or not (which would mean that the input directory is not the good one)
    combined_df = pd.concat(df_list, ignore_index=True)
    return combined_df

def combine_and_save(base_directory, output_dir):
    """Combine reports and save to output_dir."""
    combined_dataframe = combine_reports(base_directory)

    time_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f'combined_gsi_report_{time_stamp}.csv')
    combined_dataframe.to_csv(output_path, index=False)

    return combined_dataframe, output_path


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Combine GSI-RADS report.csv files from a directory.")
    parser.add_argument('base_directory', type=str,
                        help="Directory where the GSI-RADS subject folders are located.")
    parser.add_argument('output_directory', type=str,
                        help="Directory where the combined CSV will be saved.")

    args = parser.parse_args()
    combined_df, saved_path = combine_and_save(args.base_directory, args.output_directory)
    print(f"Combined report saved to: {saved_path}")