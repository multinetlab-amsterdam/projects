#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 20:58:00 2025

@author: mlzimmermann
"""

"""
Script for the post-hoc analysis to investigate whether the tumor location matters for FPN activity. 
We will investigate two covariates in that sense: The presence of the tumor in teh FPN and whether the tumor is frontal
or not. 

"""

__author__ = "Mona Lilo Margarethe Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "2024/04/04"
__status__ = "Final"

####################
# Review History   #
####################

# Reviewed by 

# %%
####################
# Libraries        #
####################

# Standard imports  ###


# Third party imports ###
import numpy as np
import pandas as pd  # version 1.1.5
import matplotlib.pyplot as plt
import seaborn as sns
import pyreadstat
from scipy.stats import ttest_ind, levene, mannwhitneyu
import statsmodels.api as sm
import matplotlib.pyplot as plt


##########################
#    DATAFRAME PREP      #
##########################

#What do I need? 
#- standardized activity dataframe with all regions (remove FPN regions per subject)
#- indicate per subject whether they have a tumor in the FPN or not/have a frontal tumor or not
#- compare the FPN averaged activity between the groups (using t-test e.g.)

#%%
df_activity = pd.read_csv("/path/to/df_activity_standardized_patients_baseline.csv")
df_activity.reset_index(drop = True, inplace = True)

df_overlaps = pd.read_csv("/path/to/20250404_post_hoc_peritumor_overlaps.csv")

#%%
#merge the two dataframes 
df_activity_overlap = pd.merge(df_activity, df_overlaps, on = ["sub", "roi"], how = "inner")

#%%
#Define the FPN regions
left_FPN = [17, 19, 21, 29, 31, 99, 137, 147, 177]
right_FPN = [16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 100, 138, 142]

FPN_all = left_FPN + right_FPN

#%%
#- define whether a subject has their peritumoral area in the FPN (get a list of patients)
#- remove FPN regions that have any overlap with tumor) 
#- average FPN regions per subject
#- add column saying whether tumor in FPN or not and column whether tumor frontal or not

#%%
# Step 1.) 
# Filter df for the rois that are in the FPN and peritumoral
filtered = df_activity_overlap[(df_activity_overlap['roi'].isin(FPN_all)) & (df_activity_overlap['perc_filt'] != 0)]

# Find which subjects have their tumor in the FPN (N = 18)
subjects_tumor_in_FPN = filtered['sub'].unique()


#%%
# Step 2.) 
# only extract FPN regions and remove FPN regions that have any overlap with the tumor 
df_FPN = df_activity_overlap[df_activity_overlap['roi'].isin(FPN_all)]
df_FPN_tumor_removed = df_FPN[df_FPN['perc_filt'] == 0]

#%%
#Step 3.) 
# average the FPN regions per subject to get one value 
df_FPN_tumor_removed_avg = df_FPN_tumor_removed.groupby('sub', as_index = False).agg(
    BBP_z_mean = ('BB_welch_z', 'mean'),
    Offset_z_mean = ('offset_z', 'mean')

    )

#%%
#Step 4.) 
# add column that indicates if a subject had their tumor in the FPN or not
df_FPN_tumor_removed_avg['tumor_FPN'] = df_FPN_tumor_removed_avg['sub'].apply(lambda sub: 'tumor_in_FPN' if sub in subjects_tumor_in_FPN else 'tumor_not_in_FPN')

#%%
#Step 5.) 
# add column indicating whether tumor is frontal  or not
#load in SPSS file with frontal tumor or not dummy variable and add to dataframe

spss_file = "/path/to/info/file.sav"
df_info,meta = pyreadstat.read_sav(spss_file, apply_value_formats = True)

df_FPN_tumor_removed_avg["Case_ID"] = df_FPN_tumor_removed_avg["sub"].str.extract("(\d+)").astype(float)

#%%
df_FPN_final = pd.merge(df_FPN_tumor_removed_avg, df_info[["Case_ID","Dummy_frontal_or_not_20250627_final"]], on = "Case_ID", how = "inner")

#%%
#descriptive results 
descriptives_in_FPN = df_FPN_final.groupby("tumor_FPN").describe()
descriptives_frontal = df_FPN_final.groupby("Dummy_frontal_or_not_20250627_final").describe()
#%%
##########################
#      STATISTICS        #
##########################

# Step 1.) Prepare dataframes for the 4 groups (2 activity, 2 covariates)
### Tumor in FPN or not
#BBP
BB_z_tumor_in_FPN = df_FPN_final[df_FPN_final["tumor_FPN"] == "tumor_in_FPN"]["BBP_z_mean"]
BB_z_tumor_not_in_FPN = df_FPN_final[df_FPN_final["tumor_FPN"] == "tumor_not_in_FPN"]["BBP_z_mean"]

#Offset
offset_z_tumor_in_FPN = df_FPN_final[df_FPN_final["tumor_FPN"] == "tumor_in_FPN"]["Offset_z_mean"]
offset_z_tumor_not_in_FPN = df_FPN_final[df_FPN_final["tumor_FPN"] == "tumor_not_in_FPN"]["Offset_z_mean"]

#%%
### Tumor Frontal or not
#BBP
BB_z_tumor_frontal= df_FPN_final[df_FPN_final["Dummy_frontal_or_not_20250627_final"] == 1]["BBP_z_mean"]
BB_z_tumor_not_frontal = df_FPN_final[df_FPN_final["Dummy_frontal_or_not_20250627_final"] == 0]["BBP_z_mean"]

#Offset
offset_z_tumor_frontal = df_FPN_final[df_FPN_final["Dummy_frontal_or_not_20250627_final"] == 1]["Offset_z_mean"]
offset_z_tumor_not_frontal = df_FPN_final[df_FPN_final["Dummy_frontal_or_not_20250627_final"] == 0]["Offset_z_mean"]

#%%
# Step 2.) 
U_BBP_FPN, p_BBP_FPN = mannwhitneyu(BB_z_tumor_in_FPN, BB_z_tumor_not_in_FPN, alternative="two-sided")
t_offset_FPN, p_offset_FPN = ttest_ind(offset_z_tumor_in_FPN, offset_z_tumor_not_in_FPN)

#%%
U_BBP_frontal, p_BBP_frontal = mannwhitneyu(BB_z_tumor_frontal, BB_z_tumor_not_frontal, alternative = "two-sided")
t_offset_frontal, p_offset_frontal = ttest_ind(offset_z_tumor_frontal, offset_z_tumor_not_frontal)

#%%
#Assumptions for independent t-test
#Normality

from scipy.stats import shapiro

stat1, p1 = shapiro(BB_z_tumor_in_FPN) # --> significant so data maybe not normal --> use Mann-Whitney U test
stat2, p2 = shapiro(BB_z_tumor_not_in_FPN)

stats3, p3 = shapiro(offset_z_tumor_in_FPN)
stats4, p4 = shapiro(offset_z_tumor_not_in_FPN)
#%%
stats5, p5 = shapiro(BB_z_tumor_frontal) # --> significant so data maybe not normal --> use Mann-Whitney U test
stats6, p6 = shapiro(BB_z_tumor_not_frontal)

stats7, p7 = shapiro(offset_z_tumor_frontal)
stats8, p8 = shapiro(offset_z_tumor_not_frontal)

#%%
#Normality check for group 1 and group 2 --> as shapiro test might indicate that are not normal
sm.qqplot(BB_z_tumor_in_FPN, line='s')
sm.qqplot(BB_z_tumor_not_in_FPN, line='s')

#%%
# Homogeneity of Variances
levene_stat_FPN, levene_p_FPN = levene(BB_z_tumor_in_FPN, BB_z_tumor_not_in_FPN)
levene_stat_FPN_off, levene_p_FPN_off = levene(offset_z_tumor_in_FPN, offset_z_tumor_not_in_FPN)

#%%
levene_stat_frontal, levene_p_frontal = levene(BB_z_tumor_frontal, BB_z_tumor_not_frontal)
levene_stat_frontal_off, levene_p_frontal_off = levene(offset_z_tumor_frontal, offset_z_tumor_not_frontal)

#--> variances all good!