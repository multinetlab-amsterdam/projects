#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


In this script I prepare the dataframes for the main statistical analysis.


"""

__author__ = "Mona Lilo Margarethe Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "2023/06/16", "first personal review MZ, 18/09/24"
__status__ = "Fina"

####################
# Review History   #
####################

# Reviewed by 

# %%
####################
# Libraries        #
####################

# Standard imports  ###
import numpy as np
import pandas as pd  # version 1.1.5
import matplotlib.pyplot as plt
import seaborn as sns
import pyreadstat

# Third party imports ###


#%%

#Load in all subjects that included in the sudy (N = 37), including bilateral tumors 

li_subs = pd.read_csv("/path/to/list/subs.csv")

#%%

#### MAKE OVERLAP DFs #### 

#script copied from 2022_activity_network project masks.py and adapted 


def calc_perc_overlap(all_subs, path, name_overlap_file,name_volume_file, area):
    """
    Function to calculate the percentage overlap of the tumor and the peritumoral area
    with the 210 gray matter BNA regions for every patient. Importantly it takes the lateralization into account, 
    thereby controling for mistakes the registration made. It also returns whether a regions fall into
    the tumoral area (if the tumor overlaps with 12% of the region it is defined as being part of the tumor area).

    Parameters
    ----------
    all_subs : pd.DataFrame,
        dataframe containing a list of all the subs that should be included in the study
    path : str,
        path to the directory where the overlaps and volume files are stored 
    name_overlap_file : str,
        name of the overlap file
    name_volume_file : str,
        ame of the volume file
    area : str,
        determine if you want to look at the tumoral area or peritumoral area (dialated mask), important for name of columns.

    Returns
    -------
    pd.DataFrame, 
        contains the overlap per region, also corrected for lateralization of the tumor. 

    """
    
    df_li = []
    for sub in all_subs["Case_ID"]:
        print(sub)
        print(type(all_subs))
  
        # get overlaps and volumes (first cols: number of voxels, second cols: cubic volume) of tumor from overlap file
        file_overlap = f"{path}{sub}/{name_overlap_file}"
        df_overlap = pd.read_csv(file_overlap, sep=" ", header=None)
        df_overlap.columns = ["voxels_overlap", "volumes_overlap", "Nan"]
        df_overlap = df_overlap.iloc[0:210]
   
        # get volumes of all regions
        file_volumes = f"{path}{sub}/{name_volume_file}"
        df_volumes = pd.read_csv(file_volumes, sep=" ", header=None)
        df_volumes.columns = ["voxels_all", "volumes_all", "Nan"]
        df_volumes = df_volumes.iloc[0:210]
    
        # concatenate the two dataframes to be able to compute percentage overlap of tumor with particular region
        df_full = pd.concat(
            [df_overlap["voxels_overlap"], df_volumes["voxels_all"]], axis=1
        )
        df_full["roi"] = np.arange(1, 211)
        df_full["sub"] = sub
    
    
        # calculate the percentage overlap of the tumor voxels with the brain region
        df_full["perc"] = (
            df_overlap["voxels_overlap"].div(df_volumes["voxels_all"].values) * 100
        )
        
        print(df_full.shape)
  
    
        #get the lateralization of the tumor 
        lateralization = all_subs.loc[all_subs["Case_ID"] == sub,"lateralization"].values[0]
        
        df_full["lateralization"] = lateralization
        
        ### LATERALIZATION OF TUMOR #### (important for step when creating peritumoral df) 
    
        #filter the overlap based on the lateralization of the tumor, to filter out wrong overlaps
        df_full["perc_filt"] = df_full["perc"]
        
        #define condition that checks if roi is even (right) or uneven (left)
        is_even = df_full["roi"] % 2 == 0
            
        #make all percentage overlap 0 for right regions (even number roi) if the lateralization is left 
        if (df_full["lateralization"] == "left").all():
            
            df_full.loc[is_even, "perc_filt"] =  0
            
        #make all percentage overlap 0 for left regions (odd number roi) if the lateralization if right 
        elif (df_full["lateralization"] == "right").all():
            
            df_full.loc[~is_even, "perc_filt"] = 0
        
        ### Area definition ### 
        
        if area == "peritumor":
            # define peritumoral regions as regions that have 12% or more overlap with the tumor
            df_full[f"{area}"] = df_full["perc_filt"].apply(lambda x: 1 if x >= 12 else 0)
        else:
            #for all other regions take all overlaps
            df_full[f"{area}"] = "NA" #because dont define regions like this for cavity
        
        # store subjects overlap information in dataframe in list
        df_li.append(df_full[["voxels_overlap","voxels_all", "roi", "perc", "sub", "lateralization", "perc_filt", f"{area}"]])
        
    # make dataframe containing all subjects overlap information
    df_all = pd.concat(df_li)
    
    # change name of columns, st know whether tumoral area or tumoral area plus peritumoral area (tumor mask dialated 2 times)
    df_all.columns = [f"voxels_overlap_{area}",f"all_voxels_region_{area}", "roi", f"perc_{area}", "sub", "lateralization", "perc_filt",  f"{area}"]
    
    return(df_all)

#%%
all_subs = li_subs
path = "/path/to/tumormasks/"
name_overlap_file = "bna_gm_TumorOverlapVolume.txt"
name_volume_file = "bna_gm_RoiVolumes.txt"
area = "tumor" 

df_tumor_overlaps = calc_perc_overlap(all_subs, path, name_overlap_file,name_volume_file, area)


#%%
all_subs = li_subs
path = "/path/to/tumormasks/"
name_overlap_file = "bna_gm_TumorOverlapVolume_dialated_2.txt"
name_volume_file = "bna_gm_RoiVolumes_dialated_2.txt"
area = "peritumor" 

df_peri_overlaps = calc_perc_overlap(all_subs, path, name_overlap_file,name_volume_file, area)
#%%
#resection cavity
all_subs = li_subs
path = "/path/to/tumormasks/"
name_overlap_file = "bna_gm_surrounding_cavity_OverlapVolume.txt"
name_volume_file = "bna_gm_RoiVolumes.txt"
area = "cavity" 

df_cavity_overlaps = calc_perc_overlap(all_subs, path, name_overlap_file,name_volume_file, area)

#%%
#enhancing tumor
all_subs = li_subs
all_subs_enhancing_tumor = all_subs.loc[all_subs["Case_ID"].isin(["sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxxx"])]
path = "/path/to/tumormasks/"
name_overlap_file = "bna_gm_surrounding_enhancingtumor_OverlapVolume.txt"
name_volume_file = "bna_gm_RoiVolumes.txt"
area = "enhancing_tumor" 

df_enhancing_tumor_overlaps = calc_perc_overlap(all_subs_enhancing_tumor, path, name_overlap_file,name_volume_file, area)


#%% 

### HISTOGRAM ### 
#Plot histogramm to check how much overlap the tumor hsa with any region --> to see how much we should use to define a region as tumoral (for the peritumoral analysis)

sns.histplot(data = df_tumor_overlaps, x = "perc_filt", bins = 40)
plt.ylim(0,100)
plt.xlim(0, 50)
plt.title("Perc overlap tumor")

#%%
sns.histplot(data = df_peri_overlaps, x = "perc_filt", bins = 40)
plt.ylim(0,100)
plt.xlim(0, 50)
plt.title("Perc overlap peri-tumor")

#%%

### EXTRACT IPSILATERAL FPN (EXCLUDING TUMORAL/PERITUMORAL AREA) ###

#Define FPN regions 
left_FPN = [17, 19, 21, 29, 31, 99, 137, 147, 177]
right_FPN = [16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 100, 138, 142]

#%%
#load in standardized dataframes from baseline and FU
df_patients_baseline_std = pd.read_csv("/path/to/df_activity_standardized_patients_baseline.csv")
df_patients_baseline_std.reset_index(drop = True, inplace = True)


df_patients_FU_std = pd.read_csv("/path/to/df_activity_standardized_patients_FU.csv")
df_patients_FU_std.reset_index(drop = True, inplace = True)

#%%

def make_lateralized_df(df_activity, df_overlaps, area):
    """
    Function to concatenate the activity and overlaps dataframes and return only the rows 
    belongig to the ipsilateral FPN.Also removes the regions that fall into the tumoral areas
    to purely obtain ipsilateral FPN regions.
    

    Parameters
    ----------
    df_activity : pd.DataFrame,
        contains all activity data 
    df_overlaps : pd.DataFrame,
        contains the tumor overlaps and lateralization data
    area: string,
        do you want to make the contralateral or ipsilateral FPN dataframe?

    Returns
    -------
    None.

    """

    ### Dataframe preparation ###
    #just to be sure, reset_index of loaded dataframes and sort by Case_ID and roi
    df_activity.reset_index(drop = True, inplace = True)
    df_activity.sort_values(by = ["sub", "roi"], inplace = True)
    
    df_overlaps.reset_index(drop = True, inplace = True)
    df_overlaps.sort_values(by = ["sub", "roi"], inplace = True)
    
    #filter out patients with bilateral tumors from overlaps dataframe
    df_overlaps_filt = df_overlaps.drop(df_overlaps[df_overlaps["lateralization"] == "bilateral"].index)
    
    #merge the overlaps dataframe and standardized dataframes (containing the activity data), dropping rows that are not in both dataframes
    df_concat = pd.merge(df_activity, df_overlaps_filt, on= ["sub", "roi"], how = "inner")
    
    print(df_activity.shape)
    print(df_overlaps_filt.shape)
    print(df_concat.shape)
    print(df_concat.head())
 
    
    ### Construct ipsilateral or contralateral dataframes ###
    # This part was coded with the help of Chatgpt 
    #define the masks for the lateralization
    left_mask = df_concat["lateralization"] == "left"
    
    right_mask = df_concat["lateralization"] == "right"
    
    
    
    if area == "ipsilateral":
        
        #apply filters to obtain desired rows 
        left_df = df_concat[left_mask].loc[df_concat["roi"].isin(left_FPN)]
        left_df.reset_index(drop = True, inplace = True)
        
        right_df = df_concat[right_mask].loc[df_concat["roi"].isin(right_FPN)]
        right_df.reset_index(drop = True, inplace = True)
        
        
    elif area == "contralateral":
        
        #apply filters to obtain desired rows 
        left_df = df_concat[left_mask].loc[df_concat["roi"].isin(right_FPN)]#if subject has left tumor, retain only right areas
        left_df.reset_index(drop = True, inplace = True)
        
        right_df = df_concat[right_mask].loc[df_concat["roi"].isin(left_FPN)]#if subject has right tumor, retain only left areas
        right_df.reset_index(drop = True, inplace = True)
        
    #put dataframes back together to obtain dataframe with only ipsilateral/contralateral rois
    df = pd.concat([left_df, right_df])
    df.sort_values(["sub", "roi"], inplace = True)
    df.reset_index(drop = True, inplace = True)
    
    #remove rows from FPN that have any overlap with tumoral/peritumoral area to obtain pure FPN dataframes
    df_pure = df.drop(df[df["perc_filt"] != 0].index)
    df_pure.reset_index(drop = True, inplace = True)

    return(df_concat, df, df_pure, df_overlaps_filt)


#%%

#Make ipsilateral and contralateral dataframe incl. removing overlaps of peritumoral area 
#BASELINE 
df_concat_ipsi, df_ipsi, df_ipsi_pure, df_overlaps_filt = make_lateralized_df(df_patients_baseline_std, df_peri_overlaps, "ipsilateral")
df_concat_contra, df_contra, df_contra_pure, df_overlaps_filt_contra = make_lateralized_df(df_patients_baseline_std, df_peri_overlaps, "contralateral")


#%%
#1-YEAR FU
df_concat_FU, df_ipsi_FU, df_ipsi_pure_FU, _ = make_lateralized_df(df_patients_FU_std, df_peri_overlaps, "ipsilateral")
df_concat_FU, df_contra_FU, df_contra_pure_FU, _ = make_lateralized_df(df_patients_FU_std, df_peri_overlaps, "contralateral")



#%%
#Average over the ipsilateral/contralateral regions to get one average FPN value per patient for the contralateral and the ipsilateral hemisphere
#BASELINE
df_ipsi_averaged = df_ipsi_pure.groupby("sub", as_index=False).mean(numeric_only = True)
df_contra_averaged = df_contra_pure.groupby("sub", as_index=False).mean(numeric_only = True)


#%%
#FU
df_ipsi_averaged_FU = df_ipsi_pure_FU.groupby("sub").mean(numeric_only = True)
df_contra_averaged_FU = df_contra_pure_FU.groupby("sub").mean(numeric_only = True) 



#%%
### Construct dataframe for the TUMORAL/PERITUMORAL/CAVITY/ENHANCING TUMOR AREA ###

#also includes bilateral tumors 

#concatenate baseline std activity dataframe (make sure that includes also sub-xxx that bilateral) and overlaps dataframe 

def make_area_df(df_activity, df_overlaps, area):
    """
    Function to make dataframe only including the peritumoral, cavity or enhancing tumor areas.

    Parameters
    ----------
    df_activity : pd.DataFrame,
        dataframe with activity data
    df_overlaps : pd.DataFrame
        dataframe with overlaps data
    area: str
        determines which area we look at

    Returns
    -------
    pd.DatFrame, 
        contains the data filtered for only peritumoral data

    """
    #preprocess dataframes
    df_activity.reset_index(drop = True, inplace = True)
    df_activity.sort_values(by = ["sub", "roi"], inplace = True)
    
    df_overlaps.reset_index(drop = True, inplace = True)
    df_overlaps.sort_values(by = ["sub", "roi"], inplace = True)
    
    #merge the overlaps dataframe and standardized dataframes (containing the activity data), dropping rows that are not in both dataframes (sub-xxx excluded for peritumor bc no tumor overlap)
    df_concat = pd.merge(df_activity, df_overlaps, on = ["sub", "roi"], how = "inner")
    
    
    if area == "peritumor":
        #select only data of the peritumoral area
        df_area_filt = df_concat[df_concat["peritumor"] == 1]
    else:
        #select all regions that have some overlap with the area
        df_area_filt = df_concat[df_concat["perc_filt"]>0]
    
    return(df_area_filt)


#%%
### PERITUMORAL DATAFRAME ###
#BASELINE
df_peri = make_area_df(df_patients_baseline_std, df_peri_overlaps,"peritumor" )
#FU 
df_peri_FU = make_area_df(df_patients_FU_std, df_peri_overlaps, "peritumor")

#%%
#average for every subject over the rois
#BASELINE
df_peri_averaged = df_peri.groupby("sub", as_index=False).mean(numeric_only = True)

#FU
df_peri_FU_averaged = df_peri_FU.groupby("sub", as_index=False).mean(numeric_only = True)


#%% 
### Create delta dataframes, for longitudinal change analysis ### 

df_peri_combined = pd.merge(df_peri_averaged, df_peri_FU_averaged, on='sub', suffixes=('_T1', '_T2'))

#%%
df_ipsi_combined = pd.merge(df_ipsi_averaged, df_ipsi_averaged_FU, on='sub', suffixes=('_T1', '_T2'))

#%%%
df_contra_combined = pd.merge(df_contra_averaged, df_contra_averaged_FU, on='sub', suffixes=('_T1', '_T2'))

#%%
### CAVITY DATAFRAMES ###

#Baseline
df_cavity = make_area_df(df_patients_baseline_std, df_cavity_overlaps, 'cavity')

#%%
# remove sub-xxx from resection cavity overlaps as will only include them in enhancing tumor (has an enhancing tumor rim)
df_cavity = df_cavity[df_cavity["sub"]!= "sub-xxx"]
#%%
#FU
df_cavity_FU = make_area_df(df_patients_FU_std, df_cavity_overlaps, 'cavity')

# remove sub-xxx from resection cavity overlaps as will only include them in enhancing tumor (has an enhancing tumor rim)
df_cavity_FU = df_cavity_FU[df_cavity_FU["sub"]!= "sub-xxx"]

#%%
#Average for every subject over the rois to get one value per subject
df_cavity_avg = df_cavity.groupby("sub", as_index=False).mean(numeric_only = True)

df_cavity_FU_avg = df_cavity_FU.groupby("sub", as_index=False).mean(numeric_only = True)

#%%
#create delta df 
df_cavity_combined = pd.merge(df_cavity_avg, df_cavity_FU_avg, on="sub", suffixes = ('_T1', '_T2'))

#%%
### ENHANCING TUMOR DATAFRAMES ###

#Baseline
df_enhancing_tumor = make_area_df(df_patients_baseline_std, df_enhancing_tumor_overlaps, "enhancing")

#FU
df_enhancing_tumor_FU = make_area_df(df_patients_FU_std, df_enhancing_tumor_overlaps, "enhancing")

#%%
#Average for every subject over the rois to get one value per subject
df_enhancing_tumor_avg = df_enhancing_tumor.groupby("sub", as_index=False).mean()

df_enhancing_tumor_FU_avg = df_enhancing_tumor_FU.groupby("sub", as_index=False).mean()

#%%
#create delta df 
df_enhancing_tumor_combined = pd.merge(df_enhancing_tumor_avg, df_enhancing_tumor_FU_avg, on="sub", suffixes = ('_T1', '_T2'))


#%%
###########################################
#### PREPARE DFS FOR POST-HOC ANALYSES ####
###########################################

##############################
### PERI- RESECTION CAVITY ###
##############################

#Focus on cavity and enhancing tumor 
df_cavity_avg["MM"] = "T1_Baseline"
df_cavity_FU_avg["MM"] = "T2_FU"


#%%
### Prepare clinical covariates (progression, epilepsy and molecular info) 

### Progression ###

#Subjects that have progression before FU MEG or within 4 months after (N = 10)
df_cavity_cov_prog = pd.concat([df_cavity_avg, df_cavity_FU_avg], axis =0)
df_cavity_cov_prog.reset_index(inplace= True)
 #%%
li_prog = ["sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx",  "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx"]
df_cavity_cov_prog["progression"] = df_cavity_cov_prog["sub"].apply(lambda x: 'progression' if x in li_prog else 'no_progression')

#%%
### Epilepsy ###

df_cavity_cov_epi= pd.concat([df_cavity_avg, df_cavity_FU_avg], axis =0)
df_cavity_cov_epi.reset_index(inplace= True)
li_no_epi = ["sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx"]
df_cavity_cov_epi["epilepsy_aed"] = df_cavity_cov_epi["sub"].apply(lambda x: 'no_epilepsy' if x in li_no_epi else 'epilepsy_aed')

#For epilepsy analysis remove subs that no AEDs or switched from having epilepsy to no epilepsy 
df_cavity_cov_epi_final = df_cavity_cov_epi[~df_cavity_cov_epi["sub"].isin(["sub-xxx", "sub-xxx"])]

#%%
### Molecular status ###
df_cavity_cov_mol= pd.concat([df_cavity_avg, df_cavity_FU_avg], axis =0)
df_cavity_cov_mol.reset_index(inplace= True)

spss_file = "/path/to/file/with/cohort/info.sav"
df,meta = pyreadstat.read_sav(spss_file, apply_value_formats = True)
df_cavity_cov_mol["Case_ID"] = df_cavity_cov_mol["sub"].str.extract("(\d+)").astype(float)

df_cavity_cov_mol_final = pd.merge(df[["Case_ID", "IDH_1p19q"]],df_cavity_cov_mol, on = "Case_ID", how = "inner" )

#%%

###  Post hoc covariates: RTH and tumor grade ###

df_cavity_cov_rth_grade= pd.concat([df_cavity_avg, df_cavity_FU_avg], axis =0)
df_cavity_cov_rth_grade.reset_index(inplace= True)

spss_file = "/path/to/file/with/cohort/info.sav"

df,meta = pyreadstat.read_sav(spss_file, apply_value_formats = True)
df_cavity_cov_rth_grade["Case_ID"] = df_cavity_cov_rth_grade["sub"].str.extract("(\d+)").astype(float)

#add grade and dummy variable whether RTH at all
df_cavity_cov_rth_grade_final = pd.merge(df[["Case_ID", "Syntax_dummy_RTH_all", "graad"]],df_cavity_cov_rth_grade, on = "Case_ID", how = "inner" )




#%%

##############################
###     ENHANCING TUMOR    ###
##############################

#prepare dataframe for anova --> long format 
df_enhancing_tumor_avg["MM"] = "T1_Baseline"
df_enhancing_tumor_FU_avg["MM"] = "T2_FU"

#concat the two dataframes 
df_enhancing_tumor_combined_long = pd.concat([df_enhancing_tumor_avg, df_enhancing_tumor_FU_avg], axis =0)
df_enhancing_tumor_combined_long.reset_index(inplace= True)


#%%
#Prepare clinical covariates (progression, epilepsy and molecular info) 

#progression
#Subjects that have progression befor FU MEG or within 4 months after (N = 10)
li_prog = ["sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx"]
df_enhancing_tumor_combined_long["progression"] = df_enhancing_tumor_combined_long["sub"].apply(lambda x: 'progression' if x in li_prog else 'no_progression')

#epilepsy
li_no_epi = ["sub-xxx", "sub-xxx", "sub-xxx",  "sub-xxx"]
df_enhancing_tumor_combined_long["epilepsy_aed"] = df_enhancing_tumor_combined_long["sub"].apply(lambda x: 'no_epilepsy' if x in li_no_epi else 'epilepsy_aed')

#For epilepsy analysis remove subs that no AEDs or switched from having epilepsy to no epilepsy (forgot to do that previously) --> save in different dataframe)
df_enhanc_tumor_filtered = df_enhancing_tumor_combined_long[~df_enhancing_tumor_combined_long["sub"].isin(["sub-xxx", "sub-xxx"])]

#%%
#concat the two dataframes 
df_enhancing_tumor_combined_long = pd.concat([df_enhancing_tumor_avg, df_enhancing_tumor_FU_avg], axis =0)
df_enhancing_tumor_combined_long.reset_index(inplace= True)

#molecular status
spss_file = "/path/to/file/with/cohort/info.sav"
df,meta = pyreadstat.read_sav(spss_file, apply_value_formats = True)
df_enhancing_tumor_combined_long["Case_ID"] = df_enhancing_tumor_combined_long["sub"].str.extract("(\d+)").astype(float)

df_enhancing_tumor_combined_long_mol = pd.merge(df[["Case_ID", "IDH_1p19q"]],df_enhancing_tumor_combined_long, on = "Case_ID", how = "inner" )


#%%%

##############################
###     PERIUMORAL AREA    ###
##############################
### Dataframe has to be in long format with columns indicating the timpoint and progression ### 
df_peri_averaged["MM"] = "T1_Baseline"
df_peri_FU_averaged["MM"] = "T2_FU"


#%%
# Progression
#subjects that have progression before FU MEG or within 4 months after (N = 11) #
#concat the two dataframes to put into long format
df_peri_cov_prog = pd.concat([df_peri_averaged, df_peri_FU_averaged], axis =0)
df_peri_cov_prog.reset_index(inplace= True)

li_prog = ["sub-xxx", "sub-xxx", "sub-xxx","sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx"]
df_peri_cov_prog["progression"] = df_peri_cov_prog["sub"].apply(lambda x: 'progression' if x in li_prog else 'no_progression')

#%%%
# Epilepsy
#prepare dataframe where it indicates which subjects have epilepsy and take AEDs at baseline and the one-year FU. 
df_peri_cov_epi = pd.concat([df_peri_averaged, df_peri_FU_averaged], axis =0)
df_peri_cov_epi.reset_index(inplace= True)

li_no_epi = ["sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx",  "sub-xxx"]
df_peri_cov_epi["epilepsy_aed"] = df_peri_cov_epi["sub"].apply(lambda x: 'no_epilepsy' if x in li_no_epi else 'epilepsy_aed')

#remove subjects that no epilepsy at FU (sub-xxx, sub-xxx)
df_peri_cov_epi_final = df_peri_cov_epi[~df_peri_cov_epi["sub"].isin(["sub-xxx", "sub-xxx"])]

#%%
# Molecular statis
#add molecular status to dataframe that includes all subjects (before filtering the aeds)
df_peri_cov_mol = pd.concat([df_peri_averaged, df_peri_FU_averaged], axis =0)
df_peri_cov_mol.reset_index(inplace= True)

spss_file = "/path/to/file/with/cohort/info.sav"
df,meta = pyreadstat.read_sav(spss_file, apply_value_formats = True)
df_peri_cov_mol["Case_ID"] = df_peri_cov_mol["sub"].str.extract("(\d+)").astype(float)

df_peri_cov_mol_final = pd.merge(df[["Case_ID", "IDH_1p19q"]],df_peri_cov_mol, on = "Case_ID", how = "inner" )

#%%
###  Post hoc covariates: RTH and tumor grade ###

df_peri_cov_rth_grade= pd.concat([df_peri_averaged, df_peri_FU_averaged], axis =0)
df_peri_cov_rth_grade.reset_index(inplace= True)

spss_file = "/path/to/file/with/cohort/info.sav"
df,meta = pyreadstat.read_sav(spss_file, apply_value_formats = True)
df_peri_cov_rth_grade["Case_ID"] = df_peri_cov_rth_grade["sub"].str.extract("(\d+)").astype(float)

#add grade and dummy variable whether RTH at all
df_peri_cov_rth_grade_final = pd.merge(df[["Case_ID", "Syntax_dummy_RTH_all", "graad"]],df_peri_cov_rth_grade, on = "Case_ID", how = "inner" )


#%% 
####### PREPARATION COX PFS AND SURVIVAL ANALYSIS #########
#Prepare a dataframe for the different areas including
#delta absolute change between T1 and T2 
#percentage change between T1 and T2 

#load in dataframes, same as df_peri_combined and df_cavity_combined 
df_peri_delta = pd.read_csv("/path/to/20230628_dataframe_peritumoral_averaged_baseline_FU_combined.csv")
df_cavity_delta = pd.read_csv("/path/to/20240214_dataframe_cavity_averaged_baseline_FU_combined.csv")


#%%
def prepare_df_COX(df, area):
    
    #BB_welch_z
    #calculate absolute difference between T1 and T2
    df["BB_welch_z_delta"] = df["BB_welch_z_T2"] - df["BB_welch_z_T1"]
    
    #calculate percentage change between T1 and T2
    df["BB_welch_z_perc_change"] = (df["BB_welch_z_delta"]/df["BB_welch_z_T1"])*100

    #offset_z
    #calculate absolute difference between T1 and T2
    df["offset_z_delta"] = df["offset_z_T2"] - df["offset_z_T1"]
 
    #calculate percentage change between T1 and T2
    df["offset_z_perc_change"] = (df["offset_z_delta"]/df["offset_z_T1"])*100
    
    return(df)

#%%
df_peri_change = prepare_df_COX(df_peri_delta, "peritumor")
df_cavity_change = prepare_df_COX(df_cavity_delta, "cavity")

#%%
### further prepare dataframe for COX hazards analysis 

#%%
#exclude patients with progression between T1 and T2 or censored date on T2 moment
li_prog = ["sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx", "sub-xxx","sub-xxx"]

df_cavity_change_excl_prog = df_cavity_change[~df_cavity_change["sub"].isin(li_prog)]

#%%
#add covariates (age at diagnosis, sex, IDH-1p19q status, KPS)
#load in spss 
spss_file = "/path/to/file/with/cohort/info.sav"

df,meta = pyreadstat.read_sav(spss_file, apply_value_formats = True)

#%%
#prepare df_cavity_change for merging
df_cavity_change_excl_prog["Case_ID"] = df_cavity_change_excl_prog["sub"].str.extract("(\d+)").astype(float)

#add covariates
df_cavity_change_met_covariates = pd.merge(df[["Case_ID", "IDH_1p19q", "kps_total", "age_at_diagnosis", "sex"]],df_cavity_change_excl_prog, on = "Case_ID", how = "inner")

#%%
#add PFS from T2 (plus 4 months) and event variable (did progression occur or not?)
df_cavity_change_met_covariates_and_PFS = pd.merge(df[["Case_ID", "progr", "PFS_weeks_MEG_FU"]], df_cavity_change_met_covariates,  on = "Case_ID", how = "inner")


