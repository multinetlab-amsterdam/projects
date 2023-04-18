#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create master dataframe

Scripts to create a dataframe with all offset and network metrics for different densities and frequencies, with tumor overlaps and raw values in one #

"""
__author__ = "Mona Lilo Margarethe Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "2022/04/04"
__status__ = "Finished"

####################
# Review History   #
####################

# Reviewed by Eduarda Centeno 20230202

# %%
####################
# Libraries        #
####################

# Standard imports  ###
import glob
import re


# Third party imports ###
import numpy as np # version 1.23.5
import pandas as pd  # version 1.1.5
from scipy.stats import zscore # version 1.9.3

#%%
########################
# FUNCTION STANDARDIZATION OF REGIONAL CLUSTERING AND CENTRALITY
########################


def standardize_metrics(
    input_dir, input_HC_dir, output_dir, freq_li, HC_subs, HC=False
):
    """
    Function used to standardize network metrics based on regional mean and std of HCs,
    to aquire a regional measure of deviation from HCs. These standardized
    metrics are stored in the original dataframe.


    Parameters
    ----------
    input_dir : str,
        directory where df with network metrics is stored per frequency band (patients)
    input_HC_dir : str,
        directory where df with network metrics is stored per frequency band (HCs)
    output_dir : str,
        directory where dfs including the standardized CC and EC metrics will be stored.
    freq_li: list, 
        list of frequencies that want to investigate
    HC_subs: pd.DatFrame, 
        contains HCs subjects (str) that want to standardize patient metrics on.Will extract Case_ID of that dataframe to get HCs.
    HC, bool:
        True if want to standardize HCs on themselves, default = False

    Returns
    -------
    None.

    """

    # loop through frequencies and densities
    for freq in freq_list:
        print(freq)

        for d in [0.2, 0.3]: #loop through densities 20%, 30% used for this analysis
            print(d)

            # get network metric file for particular frequency and density
            HC_file = glob.glob(input_HC_dir + freq + "/*" + str(d) + "*")[0]
            print(HC_file)
            file = glob.glob(input_dir + freq + "/*" + str(d) + "*")[0]
            print(file)

            # load in HC CC_EC dataframe for freq and density
            HC_df = pd.read_csv(HC_file, index_col=0)
            HC_df["sub"] = HC_df["sub"].str.lstrip("mumo_0") #remove subject identifiers to be able to match with included HC (are differently coded in the different df)
            HC_df["sub"] = HC_df["sub"].str.lstrip("Case_0")
            
            
            # only include matched HCs (on sex and age)
            HC_df = HC_df[HC_df["sub"].isin(HC_subs["Case_ID"])]
            
           
            # determine if want to standardize HCs (on themselves) or patients & load in dataframe
            if HC == True:
                subs_df = HC_df
            else:
                subs_df = pd.read_csv(file, index_col=0)
                
            # --- Step 1.) Get regional mean and std of HCs EC and CC --- #
            ### Clustering Coefficient ###
            HC_mean_CC = HC_df.groupby("rois")["CC"].mean()
            HC_std_CC = HC_df.groupby("rois")["CC"].std()
            
            
            # store means and std in dataframe
            df_HC_CC = HC_mean_CC.to_frame()
            df_HC_CC["std"] = HC_std_CC
            df_HC_CC["rois"] = HC_std_CC.index
            df_HC_CC.rename(
                columns={"CC": f"CC_mean_{freq}_{str(d)}_HC", "std": f"CC_std_{freq}_{str(d)}_HC"}, inplace=True
            )
            
            ### Eigenvector centrality ###
            HC_mean_EC = HC_df.groupby("rois")["EC"].mean()
            HC_std_EC = HC_df.groupby("rois")["EC"].std()

            # store means and std in dataframe
            df_HC_EC = HC_mean_EC.to_frame()
            df_HC_EC["std"] = HC_std_EC
            df_HC_EC["rois"] = HC_std_EC.index
            df_HC_EC.rename(
                columns={"EC": f"EC_mean_{freq}_{str(d)}_HC", "std": f"EC_std_{freq}_{str(d)}_HC"}, inplace=True
            )

            # --- Step 2.) add means and std repeatedly to rows (i.e. rows) of subs dataframe to allow for easy standardization later--- #

            ### Clustering Coefficient ###
            subs_df[f"CC_mean_{freq}_{str(d)}_HC"] = subs_df["rois"].apply(
                lambda x: df_HC_CC[f"CC_mean_{freq}_{str(d)}_HC"].loc[x]
            )
            subs_df[f"CC_std_{freq}_{str(d)}_HC"] = subs_df["rois"].apply(
                lambda x: df_HC_CC[f"CC_std_{freq}_{str(d)}_HC"].loc[x]
            )

            ### Eigenvector centrality ###
            subs_df[f"EC_mean_{freq}_{str(d)}_HC"] = subs_df["rois"].apply(
                lambda x: df_HC_EC[f"EC_mean_{freq}_{str(d)}_HC"].loc[x]
            )
            subs_df[f"EC_std_{freq}_{str(d)}_HC"] = subs_df["rois"].apply(
                lambda x: df_HC_EC[f"EC_std_{freq}_{str(d)}_HC"].loc[x]
            )


            # --- Step 3.) Standardize CC and EC and add as new columns to dataframe --- #

            ### Clustering Coefficient ###
            subs_df[f"CC_{freq}_{str(d)}_z"] = (subs_df["CC"] - subs_df[f"CC_mean_{freq}_{str(d)}_HC"]) / subs_df[
                f"CC_std_{freq}_{str(d)}_HC"
            ]

            ### Eigenvector centrality ###
            subs_df[f"EC_{freq}_{str(d)}_z"] = (subs_df["EC"] - subs_df[f"EC_mean_{freq}_{str(d)}_HC"]) / subs_df[
                f"EC_std_{freq}_{str(d)}_HC"
            ]

            # --- Step 4.) Save df for current frequency and density ---#
            subs_df.to_csv(f"{output_dir}/df_std_EC_CC_{freq}_{str(d)}.csv")


# %%
########################
# STANDARDIZE NETWORK METRICS
########################

input_dir = "/path/to/network_metric/folder/patients/"
input_HC_dir = "/path/to/network_metric/folder/HCs/"

output_dir = "/path/to/folder/where/output/should/be/stored/"
output_dir_HC = "/path/to/folder/where/output/should/be/stored/"

#read in list of included HCs 
HC_subs_df = pd.read_csv("/path/to/csv/with/HC/IDs/", index_col=0)


#%%
freq_list = ["delta", "theta", "lower_alpha"]

standardize_metrics(
    input_dir, input_HC_dir, output_dir, freq_list, HC_subs_df
)  # standardize patients metrics

standardize_metrics(
    input_HC_dir, input_HC_dir, output_dir_HC, freq_list, HC_subs_df, HC=True
)  # standardize HCs metrics on themselves


# %%
########################
# FUNCTION STANDARDIZATION ACTIVITY METRIC
########################


def standardize_activity_metrics(
    input_dir, input_HC_dir, output_dir, HC_subs, HC=False
):
    """
    Function used to standardize activity metrics based on regional mean and std of HCs,
    to aquire a regional measure of deviation from HCs. These standardized
    metrics are stored in new dataframes.


    Parameters
    ----------
    input_dir : str,
        directory where dfs with activity is stored
    input_HC_dir : str,
        directory where df with activity is stored for the HCs
    output_dir : str,
        directory where dfs including the standardized BB, offset and slope will be stored.
    HC_subs: pd.DataFrame,
        contains one column with Ids of included HCs
    HC: bool, default = False,
        set to True when standardizing data of HC.


    Returns
    -------
    None.

    """

    # get list of all files containing offset,broadband power and slope
    HC_files = sorted(glob.glob(input_HC_dir + "/*_bbp*"))
    sub_files = sorted(glob.glob(input_dir + "/*_bbp*"))

    li_HC = []
    
    # iterate through all files (subjects)
    for file in HC_files:
       
        # load in file
        HC_df = pd.read_csv(file, index_col=0).T.iloc[0:210]

        # search subject_id and add as column to dataframe
        HC_df["sub"] = re.search("([A-z]{4}_\d{3,4})", file)[0]

        # make a column for the rois
        HC_df["roi"] = np.arange(1, 211)

        # append each dataframe to list
        li_HC.append(HC_df)
        
    # make dataframe for HC data
    HC_df = pd.concat(li_HC)
    HC_df["sub"] = HC_df["sub"].str.lstrip("mumo_0") #remove subject identifiers to be able to match with included HC (are differently coded in the different df)
    HC_df["sub"] = HC_df["sub"].str.lstrip("Case_0")
    print(HC_df["sub"].dtypes)

    
    # only include matched HCs (on sex and age)
    HC_df = HC_df[HC_df["sub"].isin(HC_subs["Case_ID"])]
    

    # determine if want to standardize HCs (on themselves) or patients & load in dataframe
    if HC == True:
        pat_df = HC_df
    else:
        li_pat = []
        
        for file in sub_files:
            # load in dataframe containing activity metrics for every subject
            pat_df = pd.read_csv(file, index_col=0).T.iloc[
                0:210
            ]  # dataframe is turned to have activity metrics as columns instead of rows

            # search subject_id and add as column to dataframe
            pat_df["sub"] = re.search("(sub-\d{4})", file)[0]

            # make a column for the rois
            pat_df["roi"] = np.arange(1, 211)

            # store every patients dataframe
            li_pat.append(pat_df)
            
        # construct big dataframe containing data from all patients
        pat_df = pd.concat(li_pat)
        
    # --- Step 1.) Get regional mean and std of HC power, offset and slope --- #
    li_mean = []
    for measure in ["Broadband Power", "Offset", "Slope"]:
        
        # obtain regional mean and std of particular activity metric
        HC_mean = HC_df.groupby("roi")[measure].mean()
        HC_std = HC_df.groupby("roi")[measure].std() 

        # store means and std in dataframe
        df_mean = HC_mean.to_frame()
        df_mean[measure + "_std_HC"] = HC_std
        df_mean["rois"] = HC_std.index
        df_mean.rename(columns={measure: measure + "_mean_HC"}, inplace=True)
        df_mean.reset_index(drop=True, inplace=True)

        # store dataframes of every subject in list
        li_mean.append(df_mean)
        
    # make a dataframe with all activity metrics for all subjects
    df_mean = pd.concat(li_mean, axis=1)
    

    # drop duplicate columns (here: rois)
    df_mean = df_mean.loc[:, ~df_mean.columns.duplicated()]
    df_mean.set_index("rois", inplace=True)

    # --- Step 2.) add means and std repeatedly to rows (i.e. rois) of subs dataframe to allow for easy standardization later--- #

    for measure in ["Broadband Power", "Offset", "Slope"]:
        pat_df[measure + "_mean_HC"] = pat_df["roi"].apply(
            lambda x: df_mean[measure + "_mean_HC"].loc[x]
        )
        pat_df[measure + "_std_HC"] = pat_df["roi"].apply(
            lambda x: df_mean[measure + "_std_HC"].loc[x]
        )
    # --- Step 3.) Standardize BB, offset and slope and add as new columns to dataframe --- #

    for name, measure in {
        "BB_welch": "Broadband Power",
        "offset": "Offset",
        "slope": "Slope",
    }.items():
        pat_df[name + "_z"] = (pat_df[measure] - pat_df[measure + "_mean_HC"]) / pat_df[
            measure + "_std_HC"
        ]
        

    # --- Step 4.) Save df ---#
    if HC == True:
        pat_df.to_csv(output_dir + "/df_activity_standardized_HC.csv")
    else:
        pat_df.to_csv(output_dir + "/df_activity_standardized_patients.csv")
        
    
    return(sub_files)



# %%
########################
# STANDARDIZE ACTIVITY METRICS
########################

input_dir_bbos = "/path/to/activity/folder/patients/""
input_HC_dir_bbos = "/path/to/activity/folder/patients/"

#read in list of included HCs
HC_subs = pd.read_csv("/path/to/HCs/csv/list", index_col=0)


output_dir = "/path/to/folder/where/output/should/be/stored/"
output_dir_HC = "/path/to/folder/where/output/should/be/stored/"


standardize_activity_metrics(
    input_dir_bbos, input_HC_dir_bbos, output_dir, HC_subs, HC=False
)  # standardize patients metrics

standardize_activity_metrics(
    input_HC_dir_bbos, input_HC_dir_bbos, output_dir_HC, HC_subs, HC=True
)  # standardize HCs metrics on themselves



#%%


def construct_network_offset_df(input_dir_net, input_file_offset, output_dir, HC = False):
    """
    Function to conactenate all raw and standardized network and activity metrics in one dataframe.
    Contains the raw and standardized measures and also the mean and std of the HCs on which the data was standardized. 
    Does NOT contain tumor overlaps yet.

    Parameters
    ----------
    input_dir_net : str,
        path to directory containing dfs with standardized network metrics for 
        the different freqs and densities.
    input_file_offset : str,
        path to file containing the standardized activity metrics
    output_dir: str,
        path to directory where df with standardized and raw offset and networkmetrics should be stored
    HC : bool, optional
        boolean to set whether construct df for patients of HCs. Default = False.

    Returns
    -------
    None.

    """
    #concatenate network metrics
    li = []
    for file in glob.glob(input_dir_net): #input dir contains the df with raw, std , mean and std (HC) EC and CC for a the different frequency bands and densities
    
        if HC == False:
            freq_den = file[172:-4] #extract the frequency and density that were looking at based on the specific name of the file (string)
         
        else:
            freq_den = file[167:-4]

        #load in df containing network measures (incl. std measures) for every freq band and density 
        df_net = pd.read_csv(file, index_col = 0)
        df_net.sort_values(["sub", "rois"], inplace = True) 
        df_net.rename(columns = {"EC": f"EC_{freq_den}", "CC": f"CC_{freq_den}"}, inplace = True) #get frequency and density out of name to add to the EC/CC raw values
        
        li.append(df_net)
    
    #concatenate all networkmetrics
    df_all_network = pd.concat(li, axis = 1)
    df_all_network = df_all_network.loc[:, ~df_all_network.columns.duplicated()]
    df_all_network.reset_index(inplace = True, drop = True)
    
    #load in offset dataframe
    df_offset = pd.read_csv(input_file_offset, index_col = 0)
    df_offset.sort_values(["sub", "roi"], inplace = True)
    df_offset.reset_index(inplace = True, drop = True)
    
    #concatenate network df with offset df
    df_final = pd.concat([df_offset, df_all_network], axis = 1)
    df_final = df_final.loc[:, ~df_final.columns.duplicated()]
    df_final.reset_index(inplace = True, drop = True)
    
    if HC == False: 
        df_final.to_csv(f"{output_dir}df_std_network_offset_patients.csv")
    else: 
        df_final.to_csv(f"{output_dir}df_std_network_offset_HCs.csv")
        
    return(df_final)


#%%
########################
# Construct dataframes (without tumor overlaps)
########################

########################
# PATIENTS
########################

input_networks_pat = "path/to/folder/where/standardized/network/metrics/stored/*"
input_file_offset_pat = "path/to/csv_file/where/standardized/activity/stored.csv"
output_dir= "/path/to/folder/where/output/should/be/stored/"

df_final_pat = construct_network_offset_df(input_networks_pat, input_file_offset_pat, output_dir, False)

#%%
########################
# HCs
########################

input_networks_HC = "path/to/folder/where/standardized/network/metrics/stored/*"
input_file_offset_HC = "path/to/csv_file/where/standardized/activity/stored.csv"
output_dir = "/path/to/folder/where/output/should/be/stored/"

df_final_HC = construct_network_offset_df(input_networks_HC, input_file_offset_HC, output_dir, True)



# %%
########################
# TUMOR OVERLAPS TO DEFINE CONTRALATERAL HOMOLOGUE
########################


# load in the dataframe containing percentage overlap of tumor with the BNA regions
df_overlaps = pd.read_csv("path/to/csv_file/where/tumor_BNA_overlaps/stored.csv", index_col = 0)

#filter the overlaps based on the subjects that include in current analysis 
df_overlaps = df_overlaps[df_overlaps["sub"].isin(df_final_pat["sub"])]
df_overlaps.reset_index(inplace = True, drop = True)

#%%

# --- Step 1: Define the contralateral homologues --- #

# Get the homologue contralateral non-peritumoral rois for even (right) and uneven (left) peritumoral rois:
#   - if the peritumoral region is even (i.e. right hemisphere): get the roi from row before (i.e. uneven, left hemisphere)
#   - if the peritumroal region is uneven (i.e. left hemisphere): get the next roi (i.e. even, right hemisphere)

df_overlaps["partner_roi_right"] = (
    df_overlaps["roi"]
    .shift(-1)
    .where((df_overlaps["peritumor"] == 1) & (df_overlaps["roi"] % 2 != 0))
)
df_overlaps["partner_roi_left"] = (
    df_overlaps["roi"]
    .shift(1)
    .where((df_overlaps["peritumor"] == 1) & (df_overlaps["roi"] % 2 == 0))
)

# --- Step 2: Indicate which area the ROI belongs to (contralateral homologue?) --- #

# Group the df by subjects and get the list of (uneven (left) and even (right)) non-peritumoral rois per subject for both. Then match those to the index of the ROIs and set a 1,
# to indicate that that ROI is part of the non-peritumoral, contralateral area.First this is stored in lists per subject (in the groupby object).
# Then it is exploded  to put it back into the dataframe.

# Left hemisphere
grouped_left = pd.DataFrame(
    df_overlaps.groupby("sub").apply(
        lambda x: np.where(x["roi"].isin(x["partner_roi_left"]), 1, 0)
    )
)
grouped_left_exploded = grouped_left[0].explode().reset_index(drop=False)

df_overlaps.reset_index(drop=True, inplace=True)
df_overlaps["non_peritumoral_left"] = grouped_left_exploded[0]

# Right hemisphere
grouped_right = pd.DataFrame(
    df_overlaps.groupby("sub").apply(
        lambda x: np.where(x["roi"].isin(x["partner_roi_right"]), 1, 0)
    )
)
grouped_right_exploded = grouped_right[0].explode().reset_index(drop=False)

df_overlaps["non_peritumoral_right"] = grouped_right_exploded[0]

# Merge the left and right columns to one non-peritumoral area column
df_overlaps["non_peritumoral_hom"] = (
    df_overlaps["non_peritumoral_left"] + df_overlaps["non_peritumoral_right"]
)
#%%
# Save df_overlaps for future reference
df_overlaps.to_csv("path/to/csv_file/where/new_overlaps?should/be/stored.csv")

# %%
# --- Step 3: Add all overlaps to master_df --- #

# load in master dataframe of patients
#df_master_pat = pd.read_pickle("/path/to/master/df.pkl")

# sort df_overlaps based on sub and roi
df_overlaps = df_overlaps.sort_values(["sub", "roi"])
df_overlaps.reset_index(drop=True, inplace=True) #double chekc that indices are the same
df_final_pat.reset_index(drop=True, inplace=True) #double check that indices are the same 

# concatenate df_overlaps with master df of patients
df_master_pat_final = pd.concat([df_final_pat, df_overlaps], axis=1)
df_master_pat_final = df_master_pat_final.loc[:, ~df_master_pat_final.columns.duplicated()]


#%%                                       
# save master df of patients
df_master_pat_final.to_csv("path/to/csv_file/where/master_df/should/be/stored.csv")


#%%
#construct specific dataframes that will use in analysis:
# - df_peritumoral
# - df_homologue
# - df_non_tum 

#df_peritumor will only include subjects/regions that: 
#   - have tumor with 12% overlap with regions
#   - have unilateral tumor (exclude anyone that has tumor left and right)

#df_homologue will only include subjects/regions that: 
#   - have tumor with 12% overlap with regions
#   - have unilateral tumor (exclude anyone that has tumor left and right)
  

subs_no_overlap = df_master_pat_final.groupby('sub').sum() #get the sum of col peritumor for every subject (if the sum of peritumor is 0 then no 12% overlap)
subs_no_overlap.reset_index(inplace = True) 
li_subs_no_overlap = subs_no_overlap.loc[subs_no_overlap['peritumor'] == 0, 'sub'] # N = 5 subjects that do not have a 12% region of tumor with any region

#%%
#remove all subjects that do not have 12% overlap of tumor with any region from master dataframe of patients
df_master_only_overlaps = df_master_pat_final.loc[~df_master_pat_final['sub'].isin(li_subs_no_overlap)]
df_master_only_overlaps.reset_index(inplace = True, drop = True)

#%%
#remove subjects that have tumor in both hemispheres 
#subjects that have both a value at partner roi left and partner roi right

subs_loc_groupby = df_master_only_overlaps.groupby('sub').sum() #get again summary of sum of all patients. If a patient has non-zero number for partner roi left and right then means that they have a bilateral tumor.
subs_loc_groupby.reset_index(inplace = True)

li_subs_bilateral = subs_loc_groupby.loc[(subs_loc_groupby['partner_roi_left'] != 0) & (subs_loc_groupby['partner_roi_right'] != 0), 'sub'] #N = 12 subjects have bilateral tumors

#%%
#remove bilateral tumors
df_master_unilateral = df_master_only_overlaps.loc[~df_master_only_overlaps['sub'].isin(li_subs_bilateral)] #N = 67 
df_master_unilateral.reset_index(inplace = True, drop = True)

#%%
#Make peritumoral and homologue dataframes --> only contains peritumoral regions

#Peritumoral 
df_peritumor = df_master_unilateral[df_master_unilateral['peritumor'] == 1]
df_peritumor.reset_index(inplace = True, drop = True)

df_homologue = df_master_unilateral[df_master_unilateral['non_peritumoral_hom'] == 1]
df_homologue.reset_index(inplace = True, drop = True)
#%%

#save dataframes
df_peritumor.to_csv("path/to/csv_file/where/peritumoral_df/should/be/stored.csv")
df_homologue.to_csv("path/to/csv_file/where/homologue_df/should/be/stored.csv")

#%%

#df_non_tum (rest of the brain) will include subjects/regions that: 
# - do not have any tumor (0% overlap) 
# - also in subjects that do not have a 12% overlap! 

df_non_tum = df_master_pat_final[df_master_pat_final["perc"] ==0]
df_non_tum.reset_index(inplace = True, drop = True)
df_non_tum.to_csv("path/to/csv_file/where/non_tum_df/should/be/stored.csv")


#%% 

########################
# Make subtype dataframes#
########################

#IDH-mutant, 1p19q codeleted
IDH_codel_peri = df_peritumor[df_peritumor["sub"].str.lstrip('sub-0').astype(int).isin(IDH_codeleted["Case_ID"])] 
IDH_codel_hom = df_homologue[df_homologue["sub"].str.lstrip('sub-0').astype(int).isin(IDH_codeleted["Case_ID"])] 
IDH_codel_non_tum = df_non_tum[df_non_tum["sub"].str.lstrip('sub-0').astype(int).isin(IDH_codeleted["Case_ID"])] 


#%%

#IDH-mutant, 1p19q non-codeleted
IDH_noncodel_peri = df_peritumor[df_peritumor["sub"].str.lstrip('sub-0').astype(int).isin(IDH_noncodeleted["Case_ID"])] 
IDH_noncodel_hom = df_homologue[df_homologue["sub"].str.lstrip('sub-0').astype(int).isin(IDH_noncodeleted["Case_ID"])] 
IDH_noncodel_non_tum = df_non_tum[df_non_tum["sub"].str.lstrip('sub-0').astype(int).isin(IDH_noncodeleted["Case_ID"])] 


#%%

#IDH-wt
IDH_wt_peri = df_peritumor[df_peritumor["sub"].str.lstrip('sub-0').astype(int).isin(IDH_wt["Case_ID"])] 
IDH_wt_hom = df_homologue[df_homologue["sub"].str.lstrip('sub-0').astype(int).isin(IDH_wt["Case_ID"])] 
IDH_wt_non_tum = df_non_tum[df_non_tum["sub"].str.lstrip('sub-0').astype(int).isin(IDH_wt["Case_ID"])] 


#%%    
     
          
########################
# Z-scored dataframes for beta analysis 
########################


def z_score_df(df_original, output_dir, name):
    """
    Function to z-score numeric columns on themselves in given dataframe. Saves
    new dataframe in given output directory.

    Parameters
    ----------
    df_original : pd.DataFrame,
        dataframe containing theraw values that should be z-scored
    output_dir : str,
        output path where new dataframe should be stored
    name : str,
        name that dataframe should have

    Returns
    -------
    None.

    """
    
    numeric_columns = df_original.select_dtypes("number").columns #extract only numeric columns
    numeric_columns = numeric_columns.drop(["roi", "rois"])#remove rois column from numeric colums index to be able to retain them later
    
    df_original[numeric_columns] = df_original[numeric_columns].apply(zscore, ddof = 1)
    
    df_original.to_csv(f"{output_dir}20231003_df_{name}_z.csv")
    



#%%
output_dir = "path/to/folder/where/z_scored_df/should/be/stored.csv""

########################
#  Main dataframes
########################

z_score_df(df_peritumoral, output_dir, "peri")
z_score_df(df_homologue, output_dir, "hom")
z_score_df(df_non_tum, output_dir, "non_tum")

########################
#  Subtype dataframes
########################

z_score_df(IDH_codel_non_tum, output_dir, "IDH_codel_non_tum")
z_score_df(IDH_noncodel_non_tum, output_dir, "IDH_noncodel_non_tum")
z_score_df(IDH_wt_non_tum, output_dir, "IDH_wt_non_tum")

#%%
########################
#  HC dataframe
########################

z_score_df(df_HC, output_dir,"HC" )



