#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 10:57:31 2025

Data loading script: In this script we defne all the functions that are needed for loading the data 
and creating the dataframe of the different analyses. 

@author: mlzimmermann
"""

__author__ = "Mona Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "22-10-2025"   ### Date it was created
__status__ = "Finished" ### Production = still being developed. Else: Concluded/Finished.



####################
# Review History   #
####################

# Reviewed by Sebastien Dam, 22-04-2026 ### 

####################
# Libraries        #
####################

# Standard imports  ### (Put here built-in libraries - https://docs.python.org/3/library/)


# Third party imports ### (Put here third-party libraries e.g. pandas, numpy)
import numpy as np 
import pandas as pd
import scipy.io
import seaborn as sns
import glob
import os
import matplotlib.pyplot as plt
import nibabel as nib
from scipy.stats import spearmanr
import re
import pyreadstat


#########################################
#####        ACTIVITY PREPROC       #####
#########################################
#pipeline to standradrize the activity metrics on the HCs and add them to the main dataframe

def reshape_activity_data(df, subject_id):
    """
    Function to reshape the activity dataframes into a sub/roi x activity dataframe. 

    Parameters
    ----------
    df : pd.DataFrame,
        pandas dataframe containing the activity data for one subject 
        (output of the alternative_compute_broadbandpower.py)
        
    subject_id : str,
        string of the subject_id of the current subject. Will be added as a column 
        to the reshaped dataframe.

    Returns
    -------
    df_transposed : pd.DataFrame, 
        reshaped dataframe with subject_id, roi and activity metric columns. Each row 
        represents an roi.

    """
    #safety: first copy dataframe
    df_local = df.copy()
    
    #transpose the dataframe 
    df_transposed = df_local.T
    
    #add rois and subject_ID
    df_transposed[ "roi"] = np.arange(1,211)
    df_transposed["sub"] = subject_id

    return(df_transposed)


def loop_reshaping_data(subjects="all"):
    """
    Function to loop over all subject-specific activity dfs and concatenating 
    them into one big one containing all of the subjects activity data in long format.

    Parameters
    ----------
    subjects : str,
        defines which subjects will be analysed. Depending on value the correct 
        input directory is chosen where data is stored. Default is all. Options are 
        all and HCs.

    Returns
    -------
    df_all_subs : pd.DataFrame, 
        dataframe containing all subjects data in long format.

    """
    base_path = "/path/to/bbp/files/"
    
    if subjects == "all":
        
        input_dir = f"{base_path}/06_all_patients/"
    
    elif subjects == "HCs":
        
        input_dir = "/pat/to/HC/files/"
    
    else: print("Specify subjects properly! Options: all or HCs!")
        
    #get a list of all files in the directory that contain the bbp, offset, slope data
    all_activity_files = glob.glob(os.path.join(input_dir, "*bbp_offset_slope.csv"))
    print(all_activity_files)
    
    all_subjects = []
    
    for file in all_activity_files: 
        
        #extract filename to make distinction between patients and HCs to better extract the subject ID
        filename = os.path.basename(file)
        
        if "sub" in filename:
            #for patients
            sub_id = re.search(r"(sub-\d{4})", filename)[0]
            print(sub_id)
            
        else:
            #for HCs
            sub_id = re.search(r"([A-z]{4}_\d{3,4})", filename)[0]
            print(sub_id)
        
        #load in activity dataframe
        df = pd.read_csv(file, delimiter=",", header=0, index_col=0)
        print(df.shape)
        print(df.head())
        
        #bring dataframe in right format
        reshaped_df = reshape_activity_data(df, sub_id)
        
        #store reshaped df in list for concatenation later
        all_subjects.append(reshaped_df)
    
    #concatenate all subjects dfs into one
    df_all_subs = pd.concat(all_subjects, ignore_index = True)
    print(df_all_subs.head())
    
    return(df_all_subs)



def standardize_activity(df_all_patients, df_all_HCs):
    """
    Function to regionally standardize the activity metrics based on the HCs.

    Parameters
    ----------
    df_all_patients : pd.DataFrame,
        dataframes of the patients containing the activity data.
    df_all_HCs : pd.DatFrame,
        dataframe of the HC containing the activity data on which we want to standardize.

    Returns
    -------
    df_patients_standardized, df_HCs_standardized : 
        standardized dataframes with activity values

    """
    #safety first: copy all input data 
    df_all_patients_local = df_all_patients.copy()
    df_all_HCs_local = df_all_HCs.copy()
    
    ### Step 1: Get the HC stats (regional mean and std) ### 
    df_HCs_stats = (df_all_HCs_local.groupby("roi").agg({
       "Slope" : ["mean", "std"], 
       "Offset" : ["mean", "std"],
       "Broadband Power" : ["mean", "std"]})).reset_index()
    
    
    df_HCs_stats.columns = ['roi', 'HC_slope_mean', 'HC_slope_std', 
                         'HC_offset_mean', 'HC_offset_std', 
                         'HC_bbp_mean', 'HC_bbp_std']
    
    assert (df_HCs_stats['HC_bbp_std'] != 0).all(), "Standard deviation for one or more ROIs is zero for bbp in the HC group!"
    assert (df_HCs_stats['HC_offset_std'] != 0).all(), "Standard deviation for one or more ROIs is zero for offset in the HC group!"
    assert (df_HCs_stats['HC_slope_std'] != 0).all(), "Standard deviation for one or more ROIs is zero for slope in the HC group!"
    
    ### Step 2: Standardize patient data ###
    
    #merge patient df with HCs stats 
    df_patients_standardized = pd.merge(df_all_patients_local, df_HCs_stats, on="roi", how="left", validate = "many_to_one") #how=left to make sure the patient dataframe stays preserved
    
    #standardize the three activity metrics
    df_patients_standardized['Offset_z'] = (
        (df_patients_standardized['Offset'] - df_patients_standardized['HC_offset_mean']) /
        df_patients_standardized['HC_offset_std']
    )
    
    df_patients_standardized['BBP_z'] = (
        (df_patients_standardized['Broadband Power'] - df_patients_standardized['HC_bbp_mean']) /
        df_patients_standardized['HC_bbp_std']
    )
    
    df_patients_standardized['Slope_z'] = (
        (df_patients_standardized['Slope'] - df_patients_standardized['HC_slope_mean']) /
        df_patients_standardized['HC_slope_std']
    )
    
    ### Step 3: Standardize HC data ### 
    #merge HC df with HCs stats
    df_HCs_standardized = pd.merge(df_all_HCs_local, df_HCs_stats, on="roi", how="left")
    
    #standardize the three activity metrics
    df_HCs_standardized['Offset_z'] = (
        (df_HCs_standardized['Offset'] - df_HCs_standardized['HC_offset_mean']) /
        df_HCs_standardized['HC_offset_std']
    )
    
    df_HCs_standardized['BBP_z'] = (
        (df_HCs_standardized['Broadband Power'] - df_HCs_standardized['HC_bbp_mean']) /
        df_HCs_standardized['HC_bbp_std']
    )
    
    df_HCs_standardized['Slope_z'] = (
        (df_HCs_standardized['Slope'] - df_HCs_standardized['HC_slope_mean']) /
        df_HCs_standardized['HC_slope_std']
    )
    
    return(df_patients_standardized, df_HCs_standardized)



def create_activity_overlaps_dataframe(df_patients_standardized):
    """
    Creates the dataframe that contains the standardized activity and tumor rim connectivity
    data.

    Parameters
    ----------
    
    df_patients_standardized : pd.DataFrame,
        dataframe that contains the standardized brain activity data.

    Returns
    -------
    df_activity_overlaps: pd.DataFrame,
        dataframe containing the standardized activity and 
        tumor overlaps per subject in long format. Each row one subject and roi.
    

    """
    #safety first: copy patients dataframe
    df_patients_standardized_local = df_patients_standardized.copy()
    
    #add tumor overlaps dataframe to dataframe so that know which regions are tumor per patients and should be excluded
    df_overlaps = pd.read_csv("/path/to/overlaps/file", index_col=False) #tumor defined as 20%
    
    #merge the activity, connectivity and overlaps dataframes together
    df_activity_overlaps = pd.merge(df_patients_standardized_local, df_overlaps, on = ['sub', 'roi'], how = 'inner')
    
    
    return(df_activity_overlaps)



###################################################################
###             NORMATIVE L-TDI & TDM PREPROC                   ###
###################################################################

def load_TDI(matlab_file, mat_name):
    """
    Function to load a matlab object and extract the needed variable.

    Parameters
    ----------
    matlab_file : str,
        path to matlab file 
    mat_name : str,
        name of object in matlab file that want to extract (e.g."pli_{freq}_full_raw")

    Returns
    -------
    extracted_var : np.arr, 
        array containing each subjects L-TDI value

    """
    
    matlab_data = scipy.io.loadmat(matlab_file)
    extracted_mat = matlab_data[mat_name]
    
    return(extracted_mat)

def create_L_TDI_df():
    """
    Function to extract the L-TDI values for each subject and 
    store all subjects data in a dataframe with sub_id and L-TDI column.

    Parameters
    ----------
    None

    Returns
    -------
    df_final : pd.DataFrame, 
        dataframe containing the L-TDI values for all patients

    """
    path = "/path/to/files/"
    filepaths = glob.glob(os.path.join(path, "sub-*", "*_L_TDI_MNI152NLin2009bAsym_res-05mm.mat"))
    print("Found files:", filepaths)
    
    li_dfs = []
    for i, file in enumerate(filepaths): 
        
        sub_id = os.path.basename(os.path.dirname(file))
        print(f"Processing sub number {i}, {sub_id}")
        
        LTDI = load_TDI(file, "L_TDI")[0]
        print(LTDI)
        
        df_LTDI = pd.DataFrame(LTDI,columns = ["L-TDI"])
        
        df_LTDI["sub"]  = sub_id
        
        #temp store dataframe for later concatenation
        li_dfs.append(df_LTDI)
        
    df_final = pd.concat(li_dfs, axis = 0, ignore_index = True)
    
    return(df_final)


def create_tumor_activity_LTDI_df(df_activity_overlaps, df_LTDI, avg_metric, tumor_def, split=1):
    """
    Creates a dataframe that contains the bbp_z meaned over the tumoral region
    and the LTDI value per subject. Also returns membership of subject to low vs. high 
    tumoral activity group based on median split of the mean tumoral activity.

    Parameters
    ----------
    df_activity_overlaps : pd.DataFrame,
        contains the standardized brain activity values and the tumor overlaps in long format.
    df_LTDI : pd.DataFrame,
        contains for every subject the L-TDI value.
    avg_metric : str, 
        determines whether the mean or the median tumoral activity is used to determine the 
        median split.
    tumor_def : int, 
        determines how much % overlap with the tumor mask is needed to define a region as tumoral.
    split : str or float, 
        determines based on what the tumor activity groups should be defined. median for median split. 
        Any float for split based on number. Default is 1.

    Returns
    -------
    df_tumor_activity_LTDI: pd.DataFrame,
        contains the standardized bbp (bbp_z) values averaged over the tumoral area 
        and the LTDI value per subject. Also contains a column indicating whether 
        subject belongs to low or high tumoral activity group based on median split 
        of the tumoral activity.

    """
    #safety first: copy original dataframes
    df_activity_overlaps_local = df_activity_overlaps.copy()
    df_LTDI_local = df_LTDI.copy()
    
    #extract only tumoral regions (i.e. regions with a n% overlap with the tumor)
    if tumor_def == "any":
        
        df_tumor_activity = df_activity_overlaps_local[df_activity_overlaps_local["perc_filt"] > 0].copy()
        print("averaged over regions with any overlap with the tumor")
        
    else:
        df_tumor_activity = df_activity_overlaps_local[df_activity_overlaps_local["perc_filt"] >= tumor_def].copy()
        print(f"averaged over regions with {tumor_def} overlap with the tumor")
    
    #get number of rois that are seen as tumoral
    counts_tumor_rois = df_tumor_activity.groupby('sub').size().reset_index(name="count_tumor_rois")
    
    #average tumor activity
    df_tumor_activity_avg = df_tumor_activity.groupby("sub").agg(mean_tumoral_bbp_z = ("BBP_z", "mean"),
                                                                   median_tumoral_bbp_z = ("BBP_z", "median")).reset_index()
    
    #add counts of how many rois are tumoral 
    df_tumor_activity_avg_merge = pd.merge(df_tumor_activity_avg, counts_tumor_rois, on="sub", how="left")
    
    #merge L-TDI dataframe with tumor activity dataframe 
    df_tumor_activity_LTDI = pd.merge(df_tumor_activity_avg_merge, df_LTDI_local,  on = "sub", how = "left")
    
    #safety: copy dataframe
    df_tumor_activity_LTDI_local =  df_tumor_activity_LTDI.copy()
    
    if split == "median":
        
        #determine median based on mean bbp_z to create a median split
        split_used = np.median(df_tumor_activity_LTDI_local[f"{avg_metric}_tumoral_bbp_z"])
        print(f"{avg_metric} tumoral activity for median split is: {split_used}")
        
        df_tumor_activity_LTDI_local["split_tumor_activity"] = np.where(
            df_tumor_activity_LTDI_local["mean_tumoral_bbp_z"] > split_used, "High", "Low")
        
    else: 
        df_tumor_activity_LTDI_local["split_tumor_activity"] = np.where(
            df_tumor_activity_LTDI_local["mean_tumoral_bbp_z"] > split, "High", "Low")
        
        print(f"{avg_metric} tumoral activity for split is: {split}")
        
          
    return(df_tumor_activity_LTDI_local)


def create_final_df_LTDI_analysis(avg_metric, tumor_def, split=1):
    """
    Wrapper function that creates the final dataframe for the activity LTDI analysis 
    using the functions above.

    Parameters
    ----------
    avg_metric : str,
        determines whether we will use the mean or the median tumoral activity to determine the median 
        split for the creation of high vs. low tumoral activity groups.
    tumor_def : int, 
        determines how much % overlap with the tumor mask is needed to define a region as tumoral.
    split : str or float, 
        determines based on what the tumor activity groups should be defined. median for median split. 
        Any float for split based on number. Default is 1.

    Returns
    -------
    df_tumor_activity_LTDI : pd.DataFrame, 
         the final dataframe for the LTDI analysis. Contains the average tumoral activity (mean and median)
         and the LTDI per patient as well as the grouping variable that describes whether a patient has a low or 
         high tumoral activity (based on the mean or median).
    """
    
    #Import raw activity data    
    df_patients = loop_reshaping_data(subjects="all") 
    df_HCs = loop_reshaping_data(subjects="HCs") #using mumo data
    
    #Regionally standardize activity data (bbp)
    df_patients_standardized, df_HCs_standardized = standardize_activity(df_patients, df_HCs)
    
    #Include the tumor overlaps in the dataframe
    df_activity_overlaps = create_activity_overlaps_dataframe(df_patients_standardized)
    
    #Create the LTDI dataframe 
    df_LTDI = create_L_TDI_df()
    
    #Put together the LTDI and activity dataframes, average the activity over the tumoral area and indicate which subjects have high or low 
    #tumoral activity (based on median split)
    df_tumor_activity_LTDI = create_tumor_activity_LTDI_df(df_activity_overlaps, df_LTDI, avg_metric, tumor_def, split)
    
    return(df_tumor_activity_LTDI)


##################################################
###     PATIENT DTI CONNECTIVITY PREPROC       ###
##################################################


def extract_rim_and_reorder(mat):
  

    """
    Function to extract cortical regions of tumor rim column (column 225) and reorder to original BNA ordering (left, right, left, right etc.)

    Parameters
    ----------
    mat : np.array,
        connectivity matrix. For the patients DTI, these are ordered in the modified way (all left, all right)

    Returns
    -------
    tumor_rim_cort : np.array, 
        1D array containing the cortical tumor rim connectivities (modified order, all left, all right)
        
    reordered_tumor_rim : np.array, 
        1D array containing the reordered cortical rim connectivities (original BNA order, left, right, left, right)

    """
    tumor_rim = mat[:, 224] #extract only rim connectivities (last column)
    tumor_rim_cort = tumor_rim[np.r_[0:105, 112:217]] #extract only cortical regions
    
    #reorder based on left, right, left right to be able to add to activity dataframe (that is also ordered left right left right)
    n = len(tumor_rim_cort) // 2
    left = tumor_rim_cort[:n]
    right = tumor_rim_cort[n:]
    reordered_tumor_rim = [val for pair in zip(left, right) for val in pair]

    return(tumor_rim_cort, reordered_tumor_rim)


def loop_extract_tumor_rim():
    """
    Function to create dataframe that contains 1.) all connectivities per subject and 2.) the tumor rim connectivities
    for every subjects to the 210 cortical regions of the BNA. 

    
    Returns
    -------
    df_whole_brain_conn : pd.DataFrame, 
        contains all connectivities per subject of the upper triangle of the DTI matrix (conns x sub)
        
    df_tumor_rim : pd. Dataframe, 
        dataframe containing the tumor rim connectivities to the 210 cortical rois of the BNA for all subjects 
        in long format.

    """
    
    
    subs = ["sub-xxxx", "sub-xxxx", "sub-xxxx","sub-xxxx", "..."]
    
    #create an empty list to later make a dataframe containing all subjects data
    all_subjects_whole_brain = [] #to store dataframes with all connectivities per patient (cols: subs)
    all_subjects_rim = [] #to store the dataframes with all patients tumor rims
    
    
    for sub in subs: 
        
        print(sub)
        
        #load in DTI matrix
        input_file = f"/path/to/connectivity/file/{sub}_ses-T1_acq-024_run-1_atlas-BNA_tumor_desc-streams_connmatrix.csv"
        matrix = np.loadtxt(input_file, delimiter=",")
        
        # --- All connectivities --- # 
        
        #obtain only cortical regions of the BNA 
        #filter out subcortical regions and tumor rim (i.e. keep only regions 1-105 and 112-217)
        rows_to_keep = np.r_[0:105, 112:217]
        cols_to_keep = np.r_[0:105, 112:217]
        
        mat_cort = matrix[rows_to_keep][:, cols_to_keep]
        print(mat_cort.shape)
        
        #extract upper triangle (without diagonal)
        mat_cort_triu = mat_cort[np.triu_indices(len(mat_cort), k=1)]
        
        df_all_conns = pd.DataFrame(mat_cort_triu, columns = [f"{sub}"])
        
        #store all the dataframe in a list for later overall dataframe construction (all subs)
        all_subjects_whole_brain.append(df_all_conns)
        
        
        # --- Tumor rim connectivities --- # 
        
        #extract (cortical) tumor rim and reorder (to left, right, left, right etc.)
        tumor_rim_cort, tumor_rim_cort_reordered = extract_rim_and_reorder(matrix)
    
        #save sub, roi and connectivity in dataframe
        df_sub_tumor_rim = pd.DataFrame({
            "sub" : [sub] *210,
            "roi" : np.arange(1,211),
            "tumor_conn" : tumor_rim_cort_reordered,#reordered to original BNA ordering (left, right, left, right)
            "tumor_conn_modified_ordering" : tumor_rim_cort #modified ordering (all left, all right)
            }
            )
        
        all_subjects_rim.append(df_sub_tumor_rim)
        
    # put together whole-brain connectivities in one dataframe
    df_whole_brain_conn = pd.concat(all_subjects_whole_brain, axis = 1)
    
    #put together tumor rim connectivities in one dataframe
    df_tumor_rim = pd.concat(all_subjects_rim, ignore_index = True, axis = 0)
    
    return(df_whole_brain_conn, df_tumor_rim)


############################################################################
###    CREATE FULL DATAFRAME  ACTIVITY & TUMOR CONNECTIVITY/OVERLAPS     ###
############################################################################

def create_activity_conn_df(df_tumor_rim, df_patients_standardized):
    """
    Creates the dataframe that contains the activity and tumor rim connectivity
    information.

    Parameters
    ----------
    df_tumor_rim : pd.DataFrame,
        dataframe that contains the tumor rim connectivities per subject.
    df_patients_standardized : pd.DataFrame,
        dataframe that contains the standardized brain activity data.

    Returns
    -------
    df_activity_conn_overlaps: pd.DataFrame,
        dataframe containing the standardized activity, tumor rim connectivity and 
        tumor overlaps per subject in long format. Each row one subject and roi.
    

    """
    base_overlaps = "/path/to/dataframes/"
    
    #put together dataframe for activity and tumor rim connectivity
    df_activity_conn = pd.merge(df_tumor_rim, df_patients_standardized, on = ['sub', 'roi'], how = 'inner')
   
   
    path_overlaps = f"{base_overlaps}20240628_tumor_overlaps_20_perc.csv"
        
        
    #add tumor overlaps dataframe to dataframe so that know which regions are tumor per patients and should be excluded
    df_overlaps = pd.read_csv(path_overlaps, index_col=False) #tumor defined as 20%
    
    #merge the activity, connectivity and overlaps dataframes together
    df_activity_conn_overlaps = pd.merge(df_activity_conn, df_overlaps, on = ['sub', 'roi'], how = 'inner')
    
    return(df_activity_conn_overlaps)



def create_tumor_activity_df(df_activity_conn_overlaps, avg_metric, split=1, tumor_def=20):
    """
    Creates dataframe of mean and median over tumoral regions (n% overlap with tumor)
    per subject and grouping as low vs. high tumoral activity based on median split 
    (that is created using the mean of median tumoral activity).

    Parameters
    ----------
    df_activity_conn_overlaps : pd.DataFrame,
        Contains the activity and overlaps information per subject. Each row is a subject/roi.
    avg_metric : str,
        determines whether the mean or the median tumoral activity is used to determine the 
        median split.
    split : float, 
        activity value used to define low and high tumor activity groups. Default is 1.
    tumor_def: float, optional
       defines the percentage overlap a region should have with the tumor mask to be considered a tumor region.
       Default is 20 (%).
             

    Returns
    -------
   df_sub_tumor_activity: pd.DataFrame,
       contains the meaned and median tumoral activity and information whether the subject falls into the 
       high or low tumoral activity group.

    """
    #safety first: copy input dataframe
    df_activity_conn_overlaps_local = df_activity_conn_overlaps.copy()
    
    if tumor_def == "any":
        
        df_tumor_activity = df_activity_conn_overlaps_local[df_activity_conn_overlaps_local["perc_filt"] > 0].copy()
        print("averaged over regions with any overlap with the tumor")
        
    else:
        df_tumor_activity = df_activity_conn_overlaps_local[df_activity_conn_overlaps_local["perc_filt"] >= tumor_def].copy()
        print(f"averaged over regions with {tumor_def} overlap with the tumor")

    #get subject mean over tumor regions
    df_sub_tumor_activity = df_tumor_activity.groupby("sub").agg(
        mean_tumoral_BBP_z = ("BBP_z", "mean"),
        median_tumoral_BBP_z = ("BBP_z", "median")).reset_index()
    
    
    #assign group membership according to tumoral activity split
    df_sub_tumor_activity["split_tumor_activity"] = np.where(
        df_sub_tumor_activity[f"{avg_metric}_tumoral_BBP_z"] > split, "High", "Low")
    
    print(f"{avg_metric} tumoral activity for split is: {split}")
    
    return(df_sub_tumor_activity)


def create_avg_activity_df(df_activity_conn_overlaps, df_sub_tumor_activity, threshold=1000, tumor_def=20):
    """
    Creates a dataframe that is averaged over the connected and not connected regions separately
    per subjects and also contains the information whether a subject belongs to the high or 
    low tumoral activity group.

    Parameters
    ----------
    df_activity_conn_overlaps : pd.DataFrame,
        Contains the activity and overlaps information per subject. Each row is a subject/roi.
        
    df_sub_tumor_activity : pd.DataFrame,
        contains the meaned and median tumoral activity and information whether the subject falls into the 
        high or low tumoral activity group.
        
    threshold : int, optional
        threshold detemining whether a region is considered to be connected to the tumor. 
        The default is 1000.
        
    tumor_def : float,optional
        defines the percentage overlap a region should have with the tumor mask to be considered a tumor region.
        Default is 20 (%).

    Returns
    -------
    df_subject_avg : pd.DataFrame, 
        contains the mean and median BBP_z over the connected and non-connected connected regions per subject and indicates whether a subject
        belongs to the low or high activity group.

    """
    

    #merge grouping information based on tumoral information back to original dataframe
    df_activity_conn_groups = pd.merge(df_activity_conn_overlaps, df_sub_tumor_activity[['sub', 'split_tumor_activity']], on = 'sub', how = 'left').reset_index(drop = True)
    print(df_sub_tumor_activity['sub'].is_unique)
    print(df_activity_conn_groups.head())

    #code whether a region is connected or not to the tumor (this still also includes tumor regions themselves)
    df_activity_conn_groups['tumor_conn_binary'] = np.where(
        df_activity_conn_groups['tumor_conn'] >= threshold,  
        'Connected',                                         
        'Not_connected'                                           
    )
    
    #make copy to be sure 
    df_activity_conn_groups_local = df_activity_conn_groups.copy()
    
    #remove all tumoral regions (regions with n% overlap with the tumor)
    if tumor_def == "any":
        df_activity_conn_groups_no_tumor = df_activity_conn_groups_local[df_activity_conn_groups_local['perc_filt'] == 0].copy()
        
        print("regions with any overlap to the tumor were removed")
    else:
        df_activity_conn_groups_no_tumor = df_activity_conn_groups_local[df_activity_conn_groups_local['perc_filt'] < tumor_def].copy()
        print(f"regions with {tumor_def} overlap to the tumor were removed")

    #create a subject average (mean or median) over the low and high tumor connected regions
    df_subject_avg = df_activity_conn_groups_no_tumor.groupby(
        ["sub", "split_tumor_activity", "tumor_conn_binary"]).agg(
        mean_BBP_z = ("BBP_z", "mean"),
        median_BBP_z = ("BBP_z", "median")).reset_index()
  
            
    return(df_activity_conn_groups_no_tumor, df_subject_avg)    
    


def create_final_df_tumorconn_analysis(split_at=1,
                                       tumor_def=20):
    """
    Wrapper function to create the final dataframe for the tumor connected region analysis.
    This dataframe contains the average activity (mean and median) over the weakly and strongly
    connected regions (to the tumor) and indicates whether a patient belongs to the low or high
    tumoral activity group.

    Parameters
    ----------
    split_at : float, 
        value at which high and low tumor activity groups should be split. Default is 1.
    tumor_def : str or int,
        determines how we define tumoral regions. If the keyword "any" is given,then regions with 
        any overlap with the tumor are seen as tumor regions. If a number is given (e.g. 20) then 
        regions with that % overlap or higher are considered tumoral regions.
    

    Returns
    -------
    df_tumor_connected_rois_avg : pd.DataFrame, 
        contains the average activity over the high and low tumor connected regions (excluding the tumor
        regions themselves) per patient and information on whether a patient belongs to the high or low 
        tumor activity group.

    """
    
    #Import raw activity data    
    df_patients = loop_reshaping_data()
    df_HCs = loop_reshaping_data("HCs") #using mumo data
    
    #Regionally standardize activity data (bbp)
    df_patients_standardized, df_HCs_standardized = standardize_activity(df_patients, df_HCs)
    
    #get the tumor rim data
    _, df_tumor_rim = loop_extract_tumor_rim()
    
    #create df of tumor rim, roi activity and tumor overlaps
    df_activity_conn_overlaps = create_activity_conn_df(df_tumor_rim, df_patients_standardized)
    
    #create the tumor activity averaged df that is used to determine the median split for grouping into high vs. low tumoral activity
    df_tumor_activity_split = create_tumor_activity_df(df_activity_conn_overlaps, "mean", split = split_at, tumor_def=tumor_def)
    
    #create final dataframe where activity is averaged across rois that are strongly or weakly
    #connected to the tumor rim, also indicates whether subjects belong to high or low tumoral activity 
    #group (based on split of regions with [tumor_def]% overlap with the tumor)
    df_activity_conn_no_tum, df_tumor_connected_rois_avg = create_avg_activity_df(df_activity_conn_overlaps, df_tumor_activity_split, 
                                                         threshold = 1000, tumor_def=tumor_def)
    
    return(df_tumor_connected_rois_avg, df_tumor_activity_split, df_activity_conn_no_tum)



################################################
###      TDM & TUMOR OCCURRENCE              ###
################################################


def load_tdm_tumor_occurrence():
    """
    Loads the tumor occurrence map and the tdm map and flatten them to bring them
    in right format for voxelwise correlation and plotting.

    Returns
    -------
    df : pd.DataFrame, 
        contains the tumor ccurrence and tdm voxelwise values as columns.

    """
    #load TDM
    tdm = nib.load("/path/to/tdm/nifti.nii.gz").get_fdata()
    print(tdm.shape)
    
    #load brainmaps
    tumor_occurrence_map_boston_mgh = nib.load("path/to/tumor/occurrence/map1").get_fdata()
    print(tumor_occurrence_map_boston_mgh.shape)
    tumor_occurrence_map_tcga = nib.load("/path/to/tumor/occurrence/map2").get_fdata()
    print(tumor_occurrence_map_tcga.shape)
    
    #flatten the maps to arrays for scatterplot and correlation
    flat_tdm = tdm.flatten()
    flat_occurrence_map_boston_mgh = tumor_occurrence_map_boston_mgh.flatten().astype(int)
    flat_occurrence_map_tcga = tumor_occurrence_map_tcga.flatten().astype(int)
    
    #put in dataframe
    data_boston_mgh = {"tumor_occurrence" : flat_occurrence_map_boston_mgh, "tdm" : flat_tdm}
    data_tcga = {"tumor_occurrence" : flat_occurrence_map_tcga, "tdm" : flat_tdm}
    
    df_boston_mgh = pd.DataFrame(data_boston_mgh)
    print(df_boston_mgh["tumor_occurrence"].value_counts())
    
    df_tcga = pd.DataFrame(data_tcga)
    print(df_tcga["tumor_occurrence"].value_counts())
    
    return(df_boston_mgh, df_tcga)

##################################################
###     PREP CLINICAL ANALYSIS                 ###
##################################################

def prep_clinical_analysis(df_LTDI_analysis, df_PATNET, df_tumor_conn_analysis):
    """
    Function that prepares the dataframes (LTDI, tumor activity and activity in tumor connected rois) for the 
    survival/clinical analysis. Loads in the dataframes for the LTDI and tumor connectivity analyses and 
    the patients information on progression and the covariates.
    Filters for only relevant columns and prepares the variables to be in the right format 
    for the coph analysis (dummy).

    Parameters
    ----------
    df_LTDI_analysis : pd.DataFrame,
        dataframe of the LTDI and tumoral activity analysis.
    df_PATNET : pd.DataFrame, 
        dataframe of the PATNET analysis.
    df_tumor_conn_analysis : pd.DataFrame,
        dataframe of the tumor connectivity analysis.

    Returns
    -------
    df_dict : dictionary, 
        dictionary of dataframes that will be used in the cox proportional hazards analyses.

    """
    # --- Clinical Info --- # 
    #load in SPSS files 
    info = pyreadstat.read_sav('/path/to/spss/info')[0] 
    info["sub"] = info["Case_ID"].astype(int).astype(str).str.zfill(4).apply(lambda x: f"sub-{x}")
    
    info_final = info.copy()
    
    #adapt covariates
    #turn sex into a dummy variable (reference group: males)
    info_final["sex_binary"] = info_final["sex"].replace({1: 0, 2: 1})
    assert not info_final["sex"].isna().any(), "NaN values found in 'sex' column!"
    
    #turn kps into a dummy variable (reference group: low KPS (80 or lower))
    info_final["kps_binary"] = np.where(info_final["kps_total"]<= 80, 0, 1)
    info_final["kps_groups"] = np.where(info_final["kps_total"]<= 80, 'Low', 'High')
    assert not info_final["kps_total"].isna().any(), "NaN values found in 'kps_total' column!"
    
    # --- LTDI --- #
    #copy dataframe just to be sure 
    df_LTDI_analysis_local = df_LTDI_analysis.copy()
    
    #prep and merge the dataframes 
    #L-TDI dataframe
    df_LTDI_merged = pd.merge(df_LTDI_analysis_local, info_final, how='left', on='sub').reset_index()
    
    median_LTDI = np.median(df_LTDI_merged['L-TDI'])
    print(f'Median L-TDI is: {median_LTDI}')
    df_LTDI_merged['split_LTDI'] = np.where(df_LTDI_merged['L-TDI']> median_LTDI,'High', 'Low')

    # --- PATNET --- #
    # copy dataframe just to be sure 
    df_PATNET_local = df_PATNET.copy()
    
    median_PATNET = np.median(df_PATNET_local['Connected'])
    print(f'Median PATNET is: {median_PATNET}')
    
    df_PATNET_local['median_split_PATNET'] = np.where(df_PATNET_local['Connected']> median_PATNET,'Strong', 'Weak')
    
    df_PATNET_merged = pd.merge(df_PATNET_local, info_final, how='left', on='sub').reset_index()

    # --- Tumor connected regions --- #
    #Tumor-conn dataframe
    #copy dataframe just to be sure
    df_tumor_conn_analysis_local = df_tumor_conn_analysis.copy()
    
    #filter for only high tumor connected regions
    df_tumor_conn_filt = df_tumor_conn_analysis_local[df_tumor_conn_analysis_local["tumor_conn_binary"] == "Connected"].copy()
    
    df_tumor_conn_merged = pd.merge(df_tumor_conn_filt, info_final, how='left', on='sub').reset_index()
    
    df_tumor_conn_merged['split_high_tumor_conn_activity'] = np.where(df_tumor_conn_merged['mean_BBP_z'] > 1,'High', 'Low')
    
    # --- FINAL PREP --- #
    #extract only relevant columns to prepare anaysis
    relevant_cols_info = ['sub', 'age_at_diagnosis', 'sex_binary', 'progressie', 'duration_progression_OK', 'status_death',
                          'duration_death_OK', 'kps_binary', 'kps_total', 'kps_groups','epilepsy_dich','tumor_volume_mL']
                          
   
    df_LTDI_final = df_LTDI_merged[relevant_cols_info + ['L-TDI', 'split_LTDI', 'mean_tumoral_bbp_z', 'split_tumor_activity']]
    df_tumor_conn_final = df_tumor_conn_merged[relevant_cols_info + ['mean_BBP_z', 'split_high_tumor_conn_activity']]
    df_PATNET_final = df_PATNET_merged[relevant_cols_info + ['Connected', 'median_split_PATNET']]
    
    df_dict = {'L-TDI_clinical' : [df_LTDI_final, 'L-TDI'],
               'High_tumor_conn_clinical' : [df_tumor_conn_final, 'mean_BBP_z'],
               'Tumor_activity_clinical' : [df_LTDI_final, 'mean_tumoral_bbp_z'], 
               'PATNET_clinical' : [df_PATNET_final, 'Connected']}
    
    return(df_dict)


##################################################################
###     TUMOR ACTIVITY AND ACTIVITY IN HEALTH                 ###
##################################################################

def prep_tumor_HCs_activity(tumor_def=20, split=1):
    """
    Function to prep a dataframe for relating tumoral activity with activity in 
    those regions in HCs.

    Parameters
    ----------
    tumor_def : str or int,
           defines the percentage overlap a region should have with the tumor mask to be considered a tumor region.
           Default is 20 (%). If str any is given it averages over all regions with any overlap.
    split : float, 
        value at which high and low tumor activity groups should be split. Default is 1.

    Returns
    -------
    sub_avg : subject averaged dataframe over tumoral regions. Also outputs the average BBP
            for those regions in HCs. 

    """
    
    #Import raw activity data    
    df_patients = loop_reshaping_data("all") 
    df_HCs = loop_reshaping_data("HCs") #using mumo data

    #Regionally standardize activity data (bbp)
    df_patients_standardized, df_HCs_standardized = standardize_activity(df_patients, df_HCs)

    #Include the tumor overlaps in the dataframe
    df_activity_overlaps = create_activity_overlaps_dataframe(df_patients_standardized)

    #Filter for only tumoral regions
    #extract only tumoral regions (i.e. regions with a tumor_def% overlap with the tumor)
    if tumor_def == "any":
        
        df_tumor_activity = df_activity_overlaps[df_activity_overlaps["perc_filt"] > 0].copy()
        print("averaged over regions with any overlap with the tumor")
        
    else:
        df_tumor_activity = df_activity_overlaps[df_activity_overlaps["perc_filt"] >= tumor_def].copy()
        print(f"averaged over regions with {tumor_def} overlap with the tumor")
    
    #get the average bbp_z and HC BBP over those tumoral regions per patient
    sub_avg = df_tumor_activity.groupby("sub").agg(mean_HC_bbp = ("HC_bbp_mean", "mean"),
                                                   mean_tumoral_bbp_z = ("BBP_z", "mean")).reset_index()
    
    
    #split tumor activity into high and low 
    #assign group membership according to tumoral activity split
    sub_avg["split_tumor_activity"] = np.where(
        sub_avg["mean_tumoral_bbp_z"] > split, "High", "Low")
    
    print(f"mean tumoral activity for split is: {split}")

    return(sub_avg)

################################################
###                 PATNET                   ###
################################################

def PATNET(df, threshold=1000, tumor_def=20):
    """
    Function to determine the number of strong and weak connections to the tumor.
    Excludes tumor regions.

    Parameters
    ----------
    df : pd.DataFrame,
        dataframe containing information on the tumor overlaps and the tumor rim connectivity.
    threshold : int, optional
        Threshold used to indicate whether a connection is considered strong. Anything higher or
        above that number is considered strong.
        The default is 1000.
    tumor_def : int, optional
        Determines what we consider a tumoral regiopn, i.e. what % overlap the tumor mask 
        has to have with a region. The default is 20.

    Returns
    -------
    counts_connectivity : pd.DataFrame,
        Shows per subject the amount of regions that are strongly or weakly connected to the tumor.

    """
    
    df_local = df.copy()
    
    
    df_local['tumor_conn_binary'] = np.where(
        df_local['tumor_conn'] >= threshold,  
        'Connected',                                         
        'Not_connected'                                           
    )
    #remove all tumoral regions (regions with n% overlap with the tumor)
    if tumor_def == "any":
        df_no_tumor = df_local[df_local['perc_filt'] == 0].copy()
        
        print("regions with any overlap to the tumor were removed")
    else:
        df_no_tumor = df_local[df_local['perc_filt'] < tumor_def].copy()
        print(f"regions with {tumor_def} overlap to the tumor were removed")

    
    #Get the number of strongly and weakly connected regions
    categories = ['Connected', 'Not_connected']
    
    counts_connectivity = (
        df_no_tumor
        .groupby('sub')['tumor_conn_binary']
        .value_counts()
        .unstack(fill_value=0)  # fill missing categories with 0
        .reindex(columns=categories, fill_value=0)  # ensure column order
        .reset_index())

    
    return(counts_connectivity)


def PATNET_activity(thr, split=1, tum_def=20):
    """
    Function to create the dataframe needed for the tumor integration/tumor activity
    analysis. Dataframe will include information on number of strong and weak connections of
    tumor and tumor activity for every patient.
    

    Parameters
    ----------
    thr : int,
        Threshold to consider a connection strong.
    split : float, 
        value at which high and low tumor activity groups should be split. Default is 1.
    tum_def : int, OPTIONAL
        Determines what we consider a tumoral regiopn, i.e. what % overlap the tumor mask 
        has to have with a region. The default is 20.

    Returns
    -------
    df_tumor_integration_activity : pd.DataFrame, 
        Dataframe will include information on number of strong and weak connections of
        tumor and tumor activity for every patient.

    """
    
    #Import raw activity data    
    df_patients = loop_reshaping_data("all") 
    df_HCs = loop_reshaping_data("HCs") #using mumo data
    
    #Regionally standardize activity data (bbp)
    df_patients_standardized, df_HCs_standardized = standardize_activity(df_patients, df_HCs)
    
    #get the tumor rim data
    _, df_tumor_rim = loop_extract_tumor_rim()
    
    #create df of tumor rim, roi activity and tumor overlaps
    df_activity_conn_overlaps = create_activity_conn_df(df_tumor_rim, df_patients_standardized)

    #filter for only tumoral rois
    if tum_def == "any":
        
        df_tumor_activity = df_activity_conn_overlaps[df_activity_conn_overlaps["perc_filt"] > 0].copy()
        print("averaged over regions with any overlap with the tumor")
        
    else:
        df_tumor_activity = df_activity_conn_overlaps[df_activity_conn_overlaps["perc_filt"] >= tum_def].copy()
        print(f"averaged over regions with {tum_def} overlap with the tumor")
    
    #get number of rois that are seen as tumoral
    counts_tumor_rois = df_tumor_activity.groupby('sub').size().reset_index(name="count_tumor_rois")
   
    
    #get subject mean over tumor regions to obtain tumoral activity per subject
    df_sub_tumor_activity = df_tumor_activity.groupby("sub").agg(
        mean_tumoral_BBP_z = ("BBP_z", "mean"),
        median_tumoral_BBP_z = ("BBP_z", "median")).reset_index()
    
    #indicate which subjects show high and low tumoral activity 
    if split == "median":
        
        #determine median based on mean bbp_z to create a median split
        split_used = np.median(df_sub_tumor_activity["mean_tumoral_BBP_z"])
        print("mean tumoral activity for median split is: {split_used}")
        
        df_sub_tumor_activity["split_tumor_activity"] = np.where(
            df_sub_tumor_activity["mean_tumoral_BBP_z"] > split_used, "High", "Low")
        
    else: 
        df_sub_tumor_activity["split_tumor_activity"] = np.where(
            df_sub_tumor_activity["mean_tumoral_BBP_z"] > split, "High", "Low")
        
        print(f"mean tumoral activity for split is: {split}")

    #prep tumor integration dataframe (i.e. indicates how many strong/weak connections the tumor has)
    df_tumor_integration = PATNET(df_activity_conn_overlaps, thr, tum_def)
    
    #merge tumor activity and tumor integration dataframes
    df_tumor_integration_activity = pd.merge(df_sub_tumor_activity, df_tumor_integration, on = "sub", how= "left")
    
    #add column telling how many rois were counted as tumoral 
    df_tumor_integration_activity_final = pd.merge(df_tumor_integration_activity, counts_tumor_rois, on="sub", how="left")
   
    
    return(df_tumor_integration_activity_final)



def PATNET_tum_conn_analysis(df_tumor_integration, df_tumor_conn_analysis):
    """
    Prepare dataframe for analysis where look at differences in activity in tumor connected regions
    for patients with an strongly integrated tumor and weakly integrated tumor.

    Parameters
    ----------
    df_tumor_integration : pd.DataFrame,
        includes information on the number of connected regions per patient. 
    df_tumor_conn_analysis : pd.DataFrame, 
        dataframe that contains averaged activity over connected and non-connected regions per patient.

    Returns
    -------
    df_merged : pd.DataFrame,
        contains information on whether a patient has a strongly or weakly integrated tumor based on the 
        nr of connected regions.

    """
    
    #copy input dataframes just to be sure 
    df_tum_int_local = df_tumor_integration.copy()
    df_tum_conn_local = df_tumor_conn_analysis.copy()
    
    #merge the counts of the tumor integration (high and low) and the dataframe for the tumor conn analysis
    df_merged = pd.merge(df_tum_int_local, df_tum_conn_local, on="sub", how = "right")
    
    #make categorization for strongly integrated tumors and weakly integrated tumors
    #based on median nr. of strong connections and on 50% of all connections (i.e. 105 strongly connected regions)
    median = np.median(df_merged["Connected"])
    print(f"The median for the nr. of strong tumor connections is {median}")
    
    #categorize based on median number of strongly connected regions
    df_merged["median_split_PATNET"] = np.where(df_merged["Connected"]> median, "Strong", "Weak")
    
    #categorize based on median percentage of strongly connected regions
    df_merged["perc_strong"] = (df_merged["Connected"]/(210-df_merged["count_tumor_rois"])) * 100
    median_perc = np.median(df_merged["perc_strong"])
    print(f"The median for the % of strong tumor connections is {median_perc}")
    
    df_merged["median_split_PATNET_perc"] = np.where(df_merged["perc_strong"]> median_perc, "Strong", "Weak")

    return(df_merged)  

def PATNET_TDI_analysis(df_PATNET, df_LTDI_analysis):
    """
    Function to put together the PATNET and LTDI dataframes to be able to correlate PATNET with 
    LTDI.

    Parameters
    ----------
    df_PATNET : pd.DataFrame, 
        dataframe with information on PATNET.
    df_LTDI_analysis : pd.DataFrame, 
        dataframe with information on LTDI.

    Returns
    -------
    df_filt : pd.DataFrame, 
        dataframe containing both PATNET and LTDI.

    """
    
    #copy input dfs to be sure
    df_PATNET_local = df_PATNET.copy()
    df_LTDI_analysis_local = df_LTDI_analysis.copy()
    
    #merge the two dataframes on subject
    df_merged = pd.merge(df_PATNET_local, df_LTDI_analysis_local[["sub", "L-TDI"]], on="sub", how="left")
    
    #only extract elevant columns
    df_filt = df_merged[["sub", "L-TDI", "Connected"]].copy()
    
    return(df_filt)

#############################################################
########## EUCLIDEAN DISTANCES 
#############################################################

def compute_tumor_distances(
    tumor_overlaps: pd.DataFrame,
    bna_path: str = "/path/to/BNA/coords.txt",
    perc_filt_threshold: float = 20.0,
) -> pd.DataFrame:
    """
    Compute minimum and median distance from each non-tumor ROI to tumor ROIs,
    using MNI BNA coordinates for comparability across subjects.

    Parameters
    ----------
    tumor_overlaps : pd.DataFrame
        Columns: sub, roi (1-indexed integers 1-210), perc_filt.
    bna_path : str
        Path to BNA MNI coordinate file. Rows are ordered ROI 1-210.
    perc_filt_threshold : float
        Minimum overlap percentage to classify a ROI as tumoral. Default = 20.

    Returns
    -------
    pd.DataFrame
        One row per subject per ROI (all 210), with columns:
        sub, roi, is_tumor, median_tumor_dist, min_tumor_dist.
        Tumor ROIs have NaN for distance columns.
        Subjects with no tumor ROIs are excluded.
    """

    # -- Load BNA MNI coordinates, explicitly label ROI numbers ----------------
    bna_coords    = np.loadtxt(bna_path)[:210, :]
    bna_df        = pd.DataFrame(bna_coords, columns=["x_mni", "y_mni", "z_mni"])
    bna_df["roi"] = np.arange(1, 211)

    # -- Build distance matrix, indexed by roi number --------------------------
    diff         = bna_coords[:, np.newaxis, :] - bna_coords[np.newaxis, :, :]
    roi_labels   = np.arange(1, 211)
    bna_distance = pd.DataFrame(
        np.sqrt((diff ** 2).sum(axis=-1)),
        index   = roi_labels,                               # row    = roi number
        columns = roi_labels,                               # column = roi number
    )

    # -- Identify tumor ROIs ---------------------------------------------------
    tumor_regions = tumor_overlaps[tumor_overlaps["perc_filt"] >= perc_filt_threshold].copy()

    # -- Loop over subjects ----------------------------------------------------
    results  = []
    subjects = tumor_overlaps["sub"].unique()

    for sub in subjects:

        sub_tumor = tumor_regions[tumor_regions["sub"] == sub]["roi"].values

        if len(sub_tumor) == 0:
            print(f"[compute_tumor_distances] Skipping {sub}: no tumor ROIs found")
            continue

        # Slice by roi number directly  no index arithmetic
        tumor_dists = bna_distance.loc[sub_tumor, :]       # (n_tumor, 210)

        median_dist = tumor_dists.median(axis=0)            # Series, index = roi
        min_dist    = tumor_dists.min(axis=0)               # Series, index = roi

        # Mask tumor ROIs
        median_dist.loc[sub_tumor] = np.nan
        min_dist.loc[sub_tumor]    = np.nan

        for roi in roi_labels:
            results.append({
                "sub":               sub,
                "roi":               roi,
                "is_tumor":          roi in sub_tumor,
                "median_tumor_dist": median_dist.loc[roi],
                "min_tumor_dist":    min_dist.loc[roi],
            })

    return pd.DataFrame(results)

def prepare_eucl_distance_df(df_activity_conn_no_tumor):
    """
    Function to prepare the dataframe for the euclidean distances analysis. 
    In the end includes the information on whether a roi is connected or not and the minimum distance to the 
    centroids of the rois of the tumoral regions.

    Parameters
    ----------
    df_activity_conn_no_tumor : pd.DataFrame,
        dataframe with activity data and information if a roi is connected or not. 
        Excludes tumoral rois.

    Returns
    -------
    df_merged : pd.DataFrame, 
        dataframe with all rois and their minimum distance to the centroids of the tumor rois.
    
    df_merged_only_conn : pd.DataFrame, 
        dataframe with only tumor connected rois and their minimum distance to the centroids of the tumor rois.
    

    """
    
    #copy input df just to be sure 
    df_input_local = df_activity_conn_no_tumor.copy()
    
    #load in tumor overlaps
    tumor_overlaps = pd.read_csv("/path/to/tumor/overlaps.csv", index_col=False)

    #compute the median and minimum distance of every rois centroid to the tumor centroids 
    tumor_distances = compute_tumor_distances(tumor_overlaps)
    
    #merge activity and connectivity information (excluding tumor rois) with distances 
    df_merged = pd.merge(df_input_local, tumor_distances, on=["sub", "roi"], how = "left")
    
    #only extract connected rois (i.e. higher than connection strength 1000)
    df_merged_only_conn = df_merged[df_merged['tumor_conn_binary']=='Connected'].copy()
    
    return(df_merged, df_merged_only_conn)
