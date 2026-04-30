#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 15:32:26 2025

Main script to run data preparation, analyses and visualizations. 

@author: mlzimmermann
"""

__author__ = "Mona Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "22-10-2025"   ### Date it was created
__status__ = "Finished" ### Production = still being developed. Else: Concluded/Finished.

#Created with the help of chatgpt & claude

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
import sys
import statsmodels.formula.api as smf

sys.path.append(os.path.join(os.path.dirname(__file__), 'github'))
# Own imports 
from src.data_loading import *
from src.statistics import * 

#%%
#################################################################################
###             NORMATIVE L-TDI & TUMORAL ACTIVITY ANALYSIS                   ###
#################################################################################
#In this analysis we relate the normative L-TDI to tumoral activity.

# --- load the data --- #
df_LTDI_analysis = create_final_df_LTDI_analysis(avg_metric = 'mean', tumor_def = 20, split=1)

#%%
# --- do statistics --- #
# --- independent t-test --- #
#Test difference in L-TDI between low and high tumoral activity group
out = '/path/to/output_file/ind_ttest_results.xlsx'
results_L_TDI = independent_ttest(df = df_LTDI_analysis, groups_to_test = 'split_tumor_activity', outcome_var= 'L-TDI', output_file = out)

#%%
# --- Spearmans rho --- #
#Test correlation L-TDI and tumoral activity
out_path = '/path/to/output_file/spearman_r_results.xlsx'
results_spearmanr = spearman(df = df_LTDI_analysis, x='L-TDI', y='mean_tumoral_bbp_z', output_file=out_path)

#%%
#################################################################################
###             ACTIVITY TUMOR CONNECTED REGIONS ANALYSIS                     ###
#################################################################################
#In this analysis we investigate the activity in strongly connected regions depending on the level of tumoral activity.

# --- load the data --- # 
#load in the data using three different options of tumor definition to then use the dataframe in R
df_tumorconn_analysis, df_tumor_activity_split,df_activity_conn_no_tum = create_final_df_tumorconn_analysis(split_at = 1, tumor_def = 20)

#%%
################################################
###                 PATNET                    ###
################################################
#In this analysis we investigate PATNET (i.e. the integration of the tumor in the patients own structural connectome).
#We investigate its relation to tumoral activity, KPS and the activity we see in tumor connected regions. 

# --- load data --- #
df_PATNET = PATNET_activity(thr=1000, tum_def=20) 
df_PATNET_tum_conn = PATNET_tum_conn_analysis(df_PATNET, df_tumorconn_analysis)

#%%
# --- do statistics --- #
out_tum_int = '/path/to/output_file/tumor_integration/PATNET_tum_activity_spearman_rho.xlsx'
spearman(df_PATNET,x='mean_tumoral_BBP_z', y='Connected', output_file=out_tum_int)

#%%
#Test difference in PATNET between low and high tumoral activity group
out = '/path/to/output_file/ttest_results_PATNET_tum_activity.xlsx'
independent_ttest(df = df_PATNET, groups_to_test = 'split_tumor_activity', outcome_var= 'Connected', output_file = out)

#%%
################################################
###                 TDI & PATNET             ###
################################################
#In this analysis we investigate the relationship between PATNET and L-TDI.

# --- load data --- #
df_PATNET_TDI = PATNET_TDI_analysis(df_PATNET, df_LTDI_analysis)

# --- do statistics --- #
# --- Spearman rho --- #
out_corr = '/path/to/output_file/spearman_rho_LTDI_PATNET.xlsx'
spearman(df_PATNET_TDI,x='L-TDI', y='Connected', output_file=out_corr)
#%%
########################################################
###          TDM & TUMOR OCCURRENCE  ANALYSIS        ###
########################################################
#In this analysis we relate the tract density map of the normative tractorgram to tumor occurrence, 
#to see whether tumors more often occur in regions with stronger structural embedding.

# --- load the data --- #
df_boston, df_tcga = load_tdm_tumor_occurrence()

#%%
# --- do statistics --- #
# --- Spearman rho --- #

#TCGA dataset
out_path_tdm_tcga = '/path/to/output_file/tdm_tumor_occurrence/spearman_r_results_tcga.xlsx'
results_spearmanr_tdm_tcga = spearman(df=df_tcga, x='tumor_occurrence', y='tdm', output_file=out_path_tdm_tcga)

#%%
##################################################
###             CLINICAL ANALYSES              ###
##################################################
#In this analysis we investigate the clinical application of L-TDI and activity in tumor connected regions.

# --- load the data --- #
dict_clinical_dfs = prep_clinical_analysis(df_LTDI_analysis, df_PATNET, df_tumorconn_analysis)

#%%
# --- do statistics --- #
# --- COX PH ANALYSIS --- #
output_path_cox_progr = '/path/to/output_file/clinical_analysis/cox_hazards_model_no_covs_with_assumptions_PFS.xlsx'
write_COX_to_excel(dict_clinical_dfs, 'duration_progression_OK', 'progressie', output_path_cox_progr, covs_incl=False)

output_path_cox_death = '/path/to/output_file/clinical_analysis/20260112_cox_hazards_model_no_covs_with_assumptions_OS.xlsx'
write_COX_to_excel(dict_clinical_dfs, 'duration_death_OK', 'status_death', output_path_cox_death, covs_incl=False)

#%%
# --- KPS & L-TDI ANALYSIS --- #
out_kps = '/path/to/output_file/clinical_analysis/kps_LTDI_ttest.xlsx'
independent_ttest(dict_clinical_dfs['L-TDI_clinical'][0], groups_to_test = 'kps_groups', outcome_var='L-TDI',output_file = out_kps)

#%%
# --- KPS & PATNET ANALYSIS --- #
out_kps_PATNET = '/path/to/output_file/clinical_analysis/kps_PATNET_ttest.xlsx'
independent_ttest(dict_clinical_dfs['PATNET_clinical'][0], groups_to_test = 'kps_groups', outcome_var='Connected',output_file = out_kps_PATNET)

#%%
# --- KPS & mean bbpz in high tumor connected regions --- #
out_kps_tumor_conn = '/path/to/output_file/clinical_analysis/kps_high_tumor_conn_activity_ttest.xlsx'
independent_ttest(dict_clinical_dfs['High_tumor_conn_clinical'][0], groups_to_test = 'kps_groups', outcome_var='mean_BBP_z',output_file = out_kps_tumor_conn)

#%%
# --- Sex, Epilepsy, tum_vol & tumoral activity --- #
out_sex = '/path/to/output_file/clinical_analysis/indttest_tum_activity_sex.xlsx'
independent_ttest(dict_clinical_dfs['Tumor_activity_clinical'][0], groups_to_test = 'sex_binary', outcome_var='mean_tumoral_bbp_z',output_file = out_sex, activity_kps_split=False)

out_epilepsy = '/path/to/output_file/clinical_analysis/indttest_tum_activity_epilepsy.xlsx'
independent_ttest(dict_clinical_dfs['Tumor_activity_clinical'][0], groups_to_test = 'epilepsy_dich', outcome_var='mean_tumoral_bbp_z',output_file = out_epilepsy, activity_kps_split = False)

out_tum_vol = '/path/to/output_file/clinical_analysis/spearman_tum_activity_tum_vol.xlsx'
spearman(df=dict_clinical_dfs['Tumor_activity_clinical'][0], x='mean_tumoral_bbp_z', y='tumor_volume_mL', output_file = out_tum_vol)

out_age = '/path/to/output_file/clinical_analysis/spearman_tum_activity_age.xlsx'
spearman(df=dict_clinical_dfs['Tumor_activity_clinical'][0], x='mean_tumoral_bbp_z', y='age_at_diagnosis', output_file = out_age)
#%%
# --- Sex, Epilepsy, tum_vol & PATNET --- #
out_sex = '/path/to/output_file/clinical_analysis/indttest_PATNET_sex.xlsx'
independent_ttest(dict_clinical_dfs['PATNET_clinical'][0], groups_to_test = 'sex_binary', outcome_var='Connected',output_file = out_sex, activity_kps_split=False)

out_epilepsy = '/path/to/output_file/clinical_analysis/indttest_PATNET_epilepsy.xlsx'
independent_ttest(dict_clinical_dfs['PATNET_clinical'][0], groups_to_test = 'epilepsy_dich', outcome_var='Connected',output_file = out_epilepsy, activity_kps_split = False)

out_tum_vol = '/path/to/output_file/clinical_analysis/spearman_PATNET_tum_vol.xlsx'
spearman(df=dict_clinical_dfs['PATNET_clinical'][0], x='Connected', y='tumor_volume_mL', output_file = out_tum_vol)

out_age = '/path/to/output_file/clinical_analysis/spearman_PATNET_age.xlsx'
spearman(df=dict_clinical_dfs['PATNET_clinical'][0], x='Connected', y='age_at_diagnosis', output_file = out_age)
#%%
# --- Sex, Epilepsy, tum_vol & tumm_conn activity in connected regions --- #
out_sex = '/path/to/output_file/clinical_analysis/indttest_high_tum_conn_sex.xlsx'
independent_ttest(dict_clinical_dfs['High_tumor_conn_clinical'][0], groups_to_test = 'sex_binary', outcome_var='mean_BBP_z',output_file = out_sex, activity_kps_split=False)

out_epilepsy = '/path/to/output_file/clinical_analysis/indttest_high_tum_conn_epilepsy.xlsx'
independent_ttest(dict_clinical_dfs['High_tumor_conn_clinical'][0], groups_to_test = 'epilepsy_dich', outcome_var='mean_BBP_z',output_file = out_epilepsy, activity_kps_split = False)

out_tum_vol = '/path/to/output_file/clinical_analysis/spearman_high_tum_conn_tum_vol.xlsx'
spearman(df=dict_clinical_dfs['High_tumor_conn_clinical'][0], x='mean_BBP_z', y='tumor_volume_mL', output_file = out_tum_vol)

out_age = '/path/to/output_file/clinical_analysis/spearman_high_tum_conn_age.xlsx'
spearman(df=dict_clinical_dfs['High_tumor_conn_clinical'][0], x='mean_BBP_z', y='age_at_diagnosis', output_file = out_age)
#%% 
# --- Sex, Epilepsy, tum_vol & L-TDI --- #
out_sex = '/path/to/output_file/clinical_analysis/indttest_LTDI_sex.xlsx'
independent_ttest(dict_clinical_dfs['L-TDI_clinical'][0], groups_to_test = 'sex_binary', outcome_var='L-TDI',output_file = out_sex, activity_kps_split=False)

out_epilepsy = '/path/to/output_file/clinical_analysis/indttest_LTDI_epilepsy.xlsx'
independent_ttest(dict_clinical_dfs['L-TDI_clinical'][0], groups_to_test = 'epilepsy_dich', outcome_var='L-TDI',output_file = out_epilepsy, activity_kps_split = False)

out_tum_vol = '/path/to/output_file/clinical_analysis/spearman_LTDI_tum_vol.xlsx'
spearman(df=dict_clinical_dfs['L-TDI_clinical'][0], x='L-TDI', y='tumor_volume_mL', output_file = out_tum_vol)

out_age = '/path/to/output_file/clinical_analysis/spearman_LTDI_age.xlsx'
spearman(df=dict_clinical_dfs['L-TDI_clinical'][0], x='L-TDI', y='age_at_diagnosis', output_file = out_age)



#%%

##################################################################
###     TUMOR ACTIVITY AND ACTIVITY IN HEALTH                 ###
##################################################################
#In this analysis we relate tumoral activity to activity seen in HCs at the tumor location.
# --- load the data --- #
df_tum_HC_activity = prep_tumor_HCs_activity(tumor_def=20, split=1)

#%%
# --- do statistics --- #
out_pat_HC = '/path/to/output_file/clinical_analysis/tum_activity_HC_activity_spearman_rho.xlsx'
spearman(df_tum_HC_activity,x='mean_HC_bbp', y='mean_tumoral_bbp_z', output_file=out_pat_HC)

#%%
#correlation only in high activity group
out_pat_HC_high = '/path/to/output_file/clinical_analysis/tum_activity_HC_activity_high_spearman_rho.xlsx'
df_tum_HC_activity_high = df_tum_HC_activity[df_tum_HC_activity['split_tumor_activity']=='High'].copy()
spearman(df_tum_HC_activity_high,x='mean_HC_bbp', y='mean_tumoral_bbp_z', output_file=out_pat_HC_high)


#%%
################################################
###        EUCLIDEAN DISTANCES               ###
################################################
#In this analysis we see whether the euclidean distance to the tumor matters for activity levels of tumor connected and unconnected regions

# --- load data --- #
df_eucl_distances_all, df_eucl_distances_only_conn = prepare_eucl_distance_df(df_activity_conn_no_tum)


#%%
#--- do statistics --- #

model = smf.mixedlm(
    formula="BBP_z ~ min_tumor_dist",
    data=df_eucl_distances_only_conn,
    groups=df_eucl_distances_only_conn["sub"]
)

result = model.fit()
print(result.summary())


#%%
# --- extract descriptive statistics  ---

# Per subject: median min_tumor_dist for tumor connected and unconnected regions
per_subject = (
    df_eucl_distances_all.groupby(['sub', 'tumor_conn_binary'])['min_tumor_dist']
    .median()
    .unstack('tumor_conn_binary')
)

# Overall: median min_tumor_dist for tumor connected and unconnected regions
overall = (
    df_eucl_distances_all.groupby('tumor_conn_binary')['min_tumor_dist']
    .median()
)

print("Per subject:")
print(per_subject)

print("\nOverall:")
print(overall)


overall_range = df_eucl_distances_all.groupby('tumor_conn_binary')['min_tumor_dist'].agg(['min', 'max'])


