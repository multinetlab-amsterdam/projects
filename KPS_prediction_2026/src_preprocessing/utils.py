#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 14:29:55 2025

@author: ekoderman
"""
import pandas as pd
import numpy as np
import os
import nibabel as nib
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
from dateutil import parser

def EDA_visualization(df, title:str):
    # --- 2. Missing values heatmap ---
    plt.figure(figsize=(10,6))
    sns.heatmap(df.isnull(), cbar=False)
    plt.title(f"Missing Values Heatmap: {title}")

    # --- 3. Distribution of numerical features ---
    num_cols = df.select_dtypes(include=['int64','float64']).columns
    df[num_cols].hist(bins=30, figsize=(18, 13))
    plt.suptitle(f"Distributions of Numerical Features: {title}")
    plt.show()

    # --- 4. Correlation matrix ---
    corr = df[num_cols].corr()

    plt.figure(figsize=(12, 8))
    sns.heatmap(corr, annot=False, cmap="coolwarm", center=0)
    plt.title(f"Correlation Matrix of Numerical Features: {title}")

def plot_violins(df, column1:str, column2:str, plot_label:str):
    sns.violinplot(
        data=df,
        x=column1,
        y=column2,
        inner='quartile',   # shows median and quartiles inside the violin
        palette='Set2'      # optional: nicer colors
    )
    
    sns.stripplot(
        data=df,
        x=column1,
        y=column2,
        color='k',          # black points
        size=2,             # size of points
        jitter=True         # spread points horizontally
    )
    
    plt.xticks(rotation=50)
    plt.title(label = plot_label)
    plt.xlabel('')

def convert_ids(df, current_id_column, convert_to_id: str):
    """
    Convert subject identifiers between IMAGO_ID and PRECOG_ID using a lookup table.
    
    The function adds a new ID column to the input DataFrame by mapping values from
    the existing ID column to the requested identifier type. Rows are not dropped
    or duplicated; unmapped IDs result in NaN values.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing subject identifiers.
    current_id_column : str
        Name of the column with the current subject IDs.
    convert_to_id : str
        Target ID type ('PRECOG_ID' or 'IMAGO_ID').
    
    Returns
    -------
    df : pd.DataFrame
        DataFrame with an additional column containing the converted IDs.
    """
    id_file_path = Path('/data/anw/anw-work/MULTINET/ekoderman/KPS_prediction/subj_ids/imago_precog_id_key.csv')
    id_df = pd.read_csv(id_file_path)

    if convert_to_id == 'PRECOG_ID' and current_id_column != 'PRECOG_ID':
        # Create a mapping from IMAGO_ID -> PRECOG_ID
        mapping = dict(zip(id_df['IMAGO_ID'], id_df['PRECOG_ID']))
        df['PRECOG_ID'] = df[current_id_column].map(mapping)
        
    elif convert_to_id == 'IMAGO_ID' and current_id_column != 'IMAGO_ID':
        # Create a mapping from PRECOG_ID -> IMAGO_ID
        mapping = dict(zip(id_df['PRECOG_ID'], id_df['IMAGO_ID']))
        df['IMAGO_ID'] = df[current_id_column].map(mapping)
    
    return df

def process_subject_component_volumes(row):
    """
    Compute volumetric measurements for a subject's segmentation image.

    This function processes a single subject's segmentation file (NIfTI format)
    and calculates the volume of each labeled region (ignoring the background label 0)
    in cubic centimeters (cm³). It is designed to be used with a DataFrame row
    containing subject information and the path to the segmentation file.

    Parameters
    ----------
    row : pandas.Series
        A row from a DataFrame containing at least the following fields:
        - 'subject': str, the subject identifier.
        - 'segmentation_mni': str, file path to the NIfTI segmentation image
          in MNI space.

    Returns
    -------
    dict
        A dictionary containing:
        - 'Subject': str, the subject identifier.
        - 'Label_X': float, volume in cm³ for each label X in the segmentation
          (excluding background label 0).
        - 'Error': str, optional key present if an error occurred (e.g., file not found,
          corrupted NIfTI).

    Notes
    -----
    - Voxel volume is computed from the product of the voxel dimensions in the image header.
    - The function rounds volumes to 2 decimal places in cubic centimeters.
    - Designed to be compatible with parallel processing of multiple subjects.

    Example
    -------
    >>> row = {'subject': 'IM1134', 'segmentation_mni': '/path/to/seg.nii.gz'}
    >>> process_subject(row)
    {'Subject': 'IM1134', 'Label_1': 12.34, 'Label_2': 45.67}
    """

    sub = row['subject']
    #file_path = row['segmentation_native']
    file_path = row['segmentation_mni']

    if not os.path.exists(file_path):
        return {"Subject": sub, "Error": "File not found"}
    
    try:
        img = nib.load(file_path)
        data = img.get_fdata()
        voxel_size = np.prod(img.header.get_zooms()[:3])  # voxel volume (mm³)

        unique_labels = np.unique(data)
        volume_dict = {"Subject": sub}

        for label in unique_labels:
            if label == 0:  # skip background
                continue
            num_voxels = np.sum(data == label)
            volume_mm3 = num_voxels * voxel_size
            volume_cm3 = round(volume_mm3 / 1000, 2)  # convert to cm³
            volume_dict[f"Label_{int(label)}"] = volume_cm3

        return volume_dict

    except Exception as e:
        return {"Subject": sub, "Error": str(e)}


def convert_exclusion_dict(excluded_dict:dict, set_ids_to_compare:set):
    """
    Convert an exclusion dictionary into a pandas DataFrame of subjects and reasons.
    
    Parameters
    ----------
    excluded_dict : dict
        Dictionary where each key is an exclusion reason and each value is a dictionary
        containing at least a key 'subjs_id' mapping to a list of excluded subject IDs.
    set_ids_to_compare : set
        Set of subject IDs to check against the exclusion dictionary.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - 'subj_id': IDs from set_ids_to_compare that appear in the exclusion dictionary
        - 'reason': the corresponding reason for exclusion
    """
    matches = []
    for reason, info in excluded_dict.items():
        for subj in info['subjs_id']:
            if subj in set_ids_to_compare:
                matches.append({'subj_id': subj, 'reason': reason})
    return pd.DataFrame(matches)


def clean_kps_column(df, col_name):
    """
    Clean KPS column:
    - If value is a range like '10-60', take the maximum (60 in this case).
    - Convert all values to integers.
    """
    def extract_max(val):
        if pd.isna(val):
            return np.nan
        if isinstance(val, str) and '-' in val:
            # Split on '-', convert to numbers, take max
            try:
                numbers = [int(x.strip()) for x in val.split('-')]
                return max(numbers)
            except:
                return np.nan
        try:
            return int(val)
        except:
            return np.nan
    
    df[col_name] = df[col_name].apply(extract_max).astype('Int64')  # nullable integer
    
    return df


def calculate_month_difference(
    df, start_col, end_col, new_col, round_result=True
):
    """
    Calculate the difference in months between two date columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the date columns.
    start_col : str
        Column name for the start date (e.g., 'Surgery_date').
    end_col : str
        Column name for the end date (e.g., 'DOD').
    new_col : str
        Name of the new column to store the month difference.
    round_result : bool
        Whether to round the result to nearest integer. Default True.

    Returns
    -------
    df : pd.DataFrame
        DataFrame with the new column added.
    """

    # Convert to datetime (everything else)
    df[start_col] = df[start_col].apply(parse_date_mixed)
    df[end_col]   = df[end_col].apply(parse_date_mixed)

    #df[start_col] = pd.to_datetime(df[start_col], errors='coerce', dayfirst=True)
    #df[end_col] = pd.to_datetime(df[end_col], errors='coerce', dayfirst=True)

    # Calculate difference in months
    df[new_col] = ((df[end_col] - df[start_col]).dt.days / 30.4375) #average days in month

    # Round if requested
    if round_result:
        df[new_col] = df[new_col].round().astype('Int64')  # nullable integer

    return df


def parse_date_mixed(x):
    """Parse dates that might be in ISO (YYYY-MM-DD) or D-M-YYYY formats."""
    if pd.isna(x):
        return pd.NaT
    if isinstance(x, pd.Timestamp):
        return x
    s = str(x).strip()
    if not s:
        return pd.NaT

    # Try ISO parsing first (fast, avoids warning)
    try:
        return pd.to_datetime(s, format="%Y-%m-%d %H:%M:%S", errors="raise")
    except Exception:
        try:
            # Try again with dayfirst=True (handles 14-7-2025, 14/07/2025, etc.)
            return pd.to_datetime(s, dayfirst=True, errors="raise")
        except Exception:
            # Last resort: dateutil parser (very flexible but slower)
            try:
                return parser.parse(s, dayfirst=True)
            except Exception:
                return pd.NaT


def fill_missing_from_dbtr(merged_df, dbtr_df, target_col, dbtr_col, id_col='PRECOG_ID', dbtr_id_col='PRED_DBTR', sentinel_values = ['999']):
    """
    Fill missing values in merged_df[target_col] using values from dbtr_df[dbtr_col] matched on IDs.

    Parameters
    ----------
    merged_df : pd.DataFrame
        The main DataFrame with missing values.
    dbtr_df : pd.DataFrame
        The reference DataFrame to look up missing values.
    target_col : str
        Column in merged_df to fill.
    dbtr_col : str
        Column in dbtr_df to take values from.
    id_col : str
        Column in merged_df that identifies subjects (default 'PRECOG_ID').
    dbtr_id_col : str
        Column in dbtr_df that identifies subjects (default 'PRED_DBTR').

    Returns
    -------
    merged_df : pd.DataFrame
        DataFrame with missing values filled where possible.
    """
    
    merged_df[target_col] = merged_df[target_col].replace(sentinel_values, np.nan)
    print(f"{target_col} initially (nans): {merged_df[target_col].isna().sum()}")

    # Find missing rows
    missing = merged_df[merged_df[target_col].isna()]

    # Look them up in dbtr
    found = dbtr_df[dbtr_df[dbtr_id_col].isin(missing[id_col])]
    found = found[[dbtr_id_col, dbtr_col]]
    print(f"Found in dbtr (nans in lookup column): {found[dbtr_col].isna().sum()}")

    # Merge with main df
    merged_df = merged_df.merge(found, left_on=id_col, right_on=dbtr_id_col, how='left')

    # Fill missing values
    merged_df[target_col] = merged_df[target_col].fillna(merged_df[dbtr_col])

    # Drop extra columns
    merged_df.drop(columns=[dbtr_id_col, dbtr_col], inplace=True)

    print(f"{target_col} after filling from dbtr (nans): {merged_df[target_col].isna().sum()}")
    
    return merged_df


def check_date_column(df, col):
    print(f"NaNs: {df[col].isna().sum()}")
    parsed = pd.to_datetime(df[col], errors="coerce", dayfirst=True)
    bad = df.loc[parsed.isna() & df[col].notna(), col]
    if bad.empty:
        print("All non-null values are valid dates")
    else:
        print("Invalid date entries:")
        print(bad)

####  Master validation runner ####
def validate_merged_df(df: pd.DataFrame, key_col: str = "subject_id", numeric_cols: list[str] = None) -> None:
    print("=== Validation Report ===")
    test_no_duplicate_subjects(df, key_col)
    test_no_duplicate_column_names(df)
    test_conflicting_suffix_columns(df)
    test_missing_key_values(df, key_col)
    test_null_summary(df)
    if numeric_cols:
        test_numeric_ranges(df, numeric_cols)
    print("=== End of Report ===")


def add_excluded_subjects(excluded_dict, reason, subject_series):
    """
    Add excluded subjects to the tracker dictionary.

    Parameters:
    - excluded_dict: the dictionary tracking exclusions
    - reason: str, reason for exclusion
    - subject_series: pandas Series or list of subject IDs to exclude
    """
    subject_ids = list(subject_series)
    excluded_dict[reason] = {
        'count': len(subject_ids),
        'subjs_id': subject_ids
    }
    
def test_no_duplicate_subjects(df: pd.DataFrame, key_col: str = "subject_id") -> None:
    """Checks for duplicate subject IDs."""
    dups = df[df.duplicated(subset=key_col, keep=False)]
    if not dups.empty:
        print(f"[FAIL] Found {len(dups)} duplicate rows based on '{key_col}'.")
    else:
        print("[PASS] No duplicate subject IDs found.")

def test_no_duplicate_column_names(df: pd.DataFrame) -> None:
    """Checks if there are duplicate column names."""
    duplicates = df.columns[df.columns.duplicated()]
    if len(duplicates) > 0:
        print(f"[FAIL] Duplicate column names found: {list(duplicates)}")
    else:
        print("[PASS] No duplicate column names.")

def test_conflicting_suffix_columns(df: pd.DataFrame) -> None:
    """Checks for columns created by merge conflicts (like *_x and *_y)."""
    base_cols = [c[:-2] for c in df.columns if c.endswith(('_x', '_y'))]
    conflicts = sorted(set([c for c in base_cols if (c + '_x' in df.columns and c + '_y' in df.columns)]))
    if conflicts:
        print(f"[FAIL] Conflicting column pairs found (x/y): {conflicts}")
    else:
        print("[PASS] No conflicting x/y columns found.")

def test_missing_key_values(df: pd.DataFrame, key_col: str = "subject_id") -> None:
    """Checks for missing IDs in the key column."""
    missing = df[key_col].isna().sum()
    if missing > 0:
        print(f"[FAIL] Found {missing} missing values in '{key_col}'.")
    else:
        print(f"[PASS] No missing values in '{key_col}'.")

def test_null_summary(df: pd.DataFrame) -> None:
    """Prints a summary of null counts per column."""
    null_counts = df.isna().sum()
    print("[INFO] Null value counts per column:")
    print(null_counts[null_counts > 0])

def test_numeric_ranges(df: pd.DataFrame, numeric_cols: list[str]) -> None:
    """Checks basic numeric ranges for given columns."""
    desc = df[numeric_cols].describe()
    print("[INFO] Numeric column summary:")
    print(desc.T)
