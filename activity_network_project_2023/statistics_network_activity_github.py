#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Relation regional activity and network metrics

Script for main analysis: regional activity in relation to network characteristics


"""
__author__ = "Mona Lilo Margarethe Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "2022/03/29"
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
from datetime import datetime

# Third party imports ###
import numpy as np # version 1.23.5
import pandas as pd  # version 1.1.5
import matplotlib.pyplot as plt # version 3.6.2
import statsmodels.formula.api as smf # version 0.13.2
import statsmodels.stats.multitest # version 0.13.2
from scipy.stats import pearsonr # version 1.9.3
from scipy.stats import ttest_1samp, shapiro, wilcoxon, ttest_ind, mannwhitneyu # version 1.9.3

#%%

####################
# DATAFRAME IMPORT #
####################

####################
# Patients #
####################

#All patients (N = 84), all regions 
df_all = pd.read_csv("path/to/csv_file.csv")

#N = 67 patients (unilateral tumors, 12% overlap), tumoral regions
df_peritumoral = pd.read_csv("path/to/csv_file.csv")

#N = 67 patients (unilateral tumors, 12% overlap), non-tumoral homologue regions
df_homologue = pd.read_csv("path/to/csv_file.csv")

#N = 84 patients, non-umoral regions (0% overlap)
df_non_tum = pd.read_csv("path/to/csv_file.csv")

#%%

####################
# HCs #
####################

#HCs (N = 61) --> matched on sex and age on the patients 
df_HC = pd.read_csv("path/to/csv_file.csv")

#%%

####################
# Patients subgroups #
####################

IDH_codel_peri = pd.read_csv("path/to/csv_file.csv")
IDH_codel_hom = pd.read_csv("path/to/csv_file.csv")
IDH_codel_non_tum = pd.read_csv("path/to/csv_file.csv")

IDH_noncodel_peri = pd.read_csv("path/to/csv_file.csv")
IDH_noncodel_hom = pd.read_csv("path/to/csv_file.csv")
IDH_noncodel_non_tum = pd.read_csv("path/to/csv_file.csv")

IDH_wt_peri = pd.read_csv("path/to/csv_file.csv")
IDH_wt_hom = pd.read_csv("path/to/csv_file.csv")
IDH_wt_non_tum = pd.read_csv("path/to/csv_file.csv")


# %%
####################
# Standradized dataframes #
####################

#These dataframes were standardized on themselves to be able to run the LMM analysis to get std betas (as a means of effect size)

####################
# Main dataframes #
####################

df_peri_z = pd.read_csv("path/to/csv_file.csv")
df_hom_z = pd.read_csv("path/to/csv_file.csv")
df_non_tum_z = pd.read_csv("path/to/csv_file.csv")
df_HC_z = pd.read_csv("path/to/csv_file.csv")

####################
# Patients subgroups #
####################

df_IDH_codel_non_tum_z = pd.read_csv("path/to/csv_file.csv")
df_IDH_noncodel_non_tum_z = pd.read_csv("path/to/csv_file.csv")
df_IDH_wt_non_tum_z = pd.read_csv("path/to/csv_file.csv")



#%%

#########################
# --- LINEAR MIXED MODELS TO TEST RELATIONSHIPS --- #
#########################

def mixed_model(df, metric_activity, output_dir, area, group):
    """
    Function to statistically investigate the relationship between functional
    network metrics (EC and CC) and neural activity with a linear mixed model.
    This model fits a random intercept for subjects.
    The output for all frequencies and denisties is saved in a csv file.

    Parameters
    ----------
    df : pd.DataFrame,
        dataframe containing the data on network and activity.
    metric_activity : str,
        determines the activity metric that should be investigated (BB, BB_welch, offset, slope)
    output_dir:str, 
        path to directory where the results should be stored
    area : str,
        determines the are that is being investigated (whole brain, peritumoral, non-peritumoral homologue)
    group : str,
        determines the group that is being investigated. (patients, HCs)

    Returns
    -------
    None.

    """
    #now = datetime.now()
    li = []
    for freq in ["lower_alpha", "delta", "theta"]:
        for den in ["02", "03"]:
            print(f"{freq}_{den}")
            # fit a linear mixed model with subjects as random intercept
            model = smf.mixedlm(
                metric_activity + "_z ~ EC_" + freq + "_" + den + "_z + CC_"+freq + "_" + den + "_z",
                df,
                groups=df["sub"]
            )
            results = model.fit()

            # get the summary table of the results
            results_summary = results.summary()

            # extract the model specs, parameters (coefficients) and p-values
            df_model_specs = results_summary.tables[0]
            df_results = results_summary.tables[1]
            df_results.reset_index(inplace=True)

            df_pvalues = pd.DataFrame(results.pvalues)
            df_pvalues.reset_index(inplace=True)
            df_pvalues.columns = ["tmp", "p"]

            # save everything in a dataframe
            df_all = pd.concat([df_model_specs, df_results, df_pvalues], axis=1)
            df_all["Frequency"] = freq
            df_all["Density"] = str(int(den[1]) * 10) + "%"

            li.append(df_all)
    # construct a dataframe with all model results from all frequencies and densities
    df_final = pd.concat(li, axis=0)

    # save df
    df_final.to_csv(f"{output_dir}mixed_model_{area}_{group}_{metric_activity}.csv")
    


# %%
########################
# PATIENTS ANALYSIS
########################

output_dir = "path/to/coutput_folder"

# Peritumoral
mixed_model(df_peritumoral, "offset",output_dir, "peritumoral", "patients")

# Non-peritumoral homologue
mixed_model(df_homologue, "offset", output_dir,"homologue", "patients")

# Rest of the brain (all non-tumoral areas)
mixed_model(df_non_tum, "offset", output_dir, "non_tum", "patients")

# %%

########################
# HCs ANALYSIS
########################

mixed_model(df_HC, "offset", output_dir, "whole_brain", "HC")


# %%

########################
# SUBTYPE ANALYSIS
########################

# IDH wildtype
mixed_model(IDH_wt_non_tum, "offset", output_dir, "rest", "IDH_wt")
mixed_model(IDH_codel_non_tum, "offset", output_dir, "rest", "IDH_codeleted")
mixed_model(IDH_codel_non_tum, "offset", output_dir, "rest", "IDH_noncodeleted")


# %%

########################
# STANDARDIZED ANALYSIS (TO OBTAIN BETA COEFFICIENTS)
########################

output_dir = "path/to/output_folder/"

### Effect size 1: Standardization of coefficients (betas)###
mixed_model(df_non_tum_z, "offset",output_dir, "non_tum", "patients_z")
mixed_model(df_HC_z, "offset", output_dir, "whole_brain", "HC_z")
mixed_model(df_IDH_wt_non_tum_z, "offset", output_dir, "rest", "IDH_wt_z")
mixed_model(df_IDH_codel_non_tum_z, "offset", output_dir, "rest", "IDH_codeleted_z")
mixed_model(df_IDH_noncodel_non_tum_z, "offset", output_dir, "rest", "IDH_noncodeleted_z")

# %%


def mixed_model_interactions(df, metric_activity, output_dir, area):
    """
    Function to statistically investigate the relationship between functional
    network metrics (EC and CC) and neural activityand differences between the groups
    with a linear mixed model.
    This model fits a random intercept subjects.
    The output for all frequencies and denisties is saved in a csv file.

    Parameters
    ----------
    df : pd.DataFrame,
        dataframe containing the data on network and activity.
    metric_activity : str,
        determines the activity metric that should be investigated (BB, BB_welch, offset, slope)
    output_dir:str, 
        path to directory where the results should be stored
    area : str,
        determines the are that is being investigated (whole brain, peritumoral, non-peritumoral homologue)


    Returns
    -------
    None.

    """
    now = datetime.now()
    li = []
    for freq in ["delta", "theta", "lower_alpha"]:
        for den in ["02", "03"]:
            # mixed model with random intercept for subjects; interaction terms are added to investigate differences between patients and HCs
            model = smf.mixedlm(
                f"{metric_activity}_z ~ EC_{freq}_{den}_z + CC_{freq}_{den}_z + C(group) + C(group)*EC_{freq}_{den}_z + C(group)*CC_{freq}_{den}_z",
                df,
                groups="sub",
            )
            results = model.fit()

            # get the summary table of the results
            results_summary = results.summary()

            # extract the model specs, parameters (coefficients) and p-values
            df_model_specs = results_summary.tables[0]
            df_results = results_summary.tables[1]
            df_results.reset_index(inplace=True)

            df_pvalues = pd.DataFrame(results.pvalues)
            df_pvalues.reset_index(inplace=True)
            df_pvalues.columns = ["tmp", "p"]

            # save everything in a dataframe
            df_all = pd.concat([df_model_specs, df_results, df_pvalues], axis=1)
            df_all["Frequency"] = freq
            df_all["Density"] = str(int(den[1]) * 10) + "%"

            li.append(df_all)
    # construct a dataframe with all model results from all frequencies and densities
    df_final = pd.concat(li, axis=0)

    # save df
    df_final.to_csv(
        f"{output_dir}{area}_interaction_{metric_activity}.csv")


# %%

# Make dataframes containing both patients and HC data in preparation for
# the interaction analysis
df_HC["group"] = "HCs"
df_HC["sub"] = "Case_" + df_HC["sub"].astype(str)

# --- Peritumoral dataframes --- #
df_peritumoral["group"] = "patients"
df_peritumoral_all = pd.concat([df_peritumoral, df_HC], axis=0)
df_peritumoral_all.reset_index(drop=True, inplace=True)

# --- Non-Peritumoral Homologue dataframes --- #
df_homologue["group"] = "patients"
df_homologue_all = pd.concat([df_homologue, df_HC], axis=0)
df_homologue_all.reset_index(drop=True, inplace=True)

# --- All Non-peritumoral areas --- #
df_non_tum["group"] = "patients"
df_non_tum_all = pd.concat([df_non_tum, df_HC], axis=0)
df_non_tum_all.reset_index(drop=True, inplace=True)


# %%
########################
# INTERACTION ANALYSIS
########################
output_dir = "path/to/output_folder/"

# --- Peritumoral areas --- #
mixed_model_interactions(df_peritumoral_all, "offset", output_dir, "peritumor")

#%%
# --- Homologue areas --- #
mixed_model_interactions(df_homologue_all, "offset", output_dir, "homologue")

#%%
# --- Nontumoral areas --- #
mixed_model_interactions(df_non_tum_all, "offset", output_dir,"non-tumoral")

# %%
########################
# FDR CORRECTION
########################

#For all FDR corrections in this manuscript, we additionally used this website to check the corrections:https://tools.carbocation.com/FDR

import glob
import statsmodels.stats.multitest # not used?

# path only to non-interaction models --> this is for non-interaction models
files = glob.glob("path/to/folder/with/results/*.csv")#comment/uncomment depending on whether you want to correct the norma model or the interaction model

# path to interaction analyses
#files = glob.glob("path/to/folder/with/interaction/results/*.csv")

for file in files:
    ##print(name)
    df = pd.read_csv(file)

    df_nan = df.loc[df["p"].isna()]
    df_not_nan = df.loc[~df["p"].isna()]

    df_IVs = df_not_nan.loc[df_not_nan['index'].str.startswith(('CC', 'EC'))] #uncomment when looking at normal analysis
    # df_IVs = df_not_nan.loc[
    #     df_not_nan["index"].str.startswith(("CC", "EC", "C(group)[T.patients]:")) #uncomment when looking at intercation analysis
    # ]

    # correct only with p-values from IVs (EC, CC)
    df_IVs["corrected_only_IVs"] = statsmodels.stats.multitest.fdrcorrection(
        df_IVs["p"]
    )[1]

    df_not_nan = df_not_nan.loc[~df_not_nan["index"].str.startswith(("CC", "EC"))]
    df_sorted = pd.concat([df_nan, df_not_nan, df_IVs]).sort_index()

    # correct with all pvalues from model (including intercept)
    df_nan = df_sorted.loc[df_sorted["p"].isna()]
    df_not_nan = df_sorted.loc[~df_sorted["p"].isna()]

    df_not_nan["corrected"] = statsmodels.stats.multitest.fdrcorrection(
        df_not_nan["p"]
    )[1]

    df_final = pd.concat([df_nan, df_not_nan]).sort_index()

    df_final.to_csv(file)
# %%
########################
# WITHIN-SUBJECT CORRELATION ANALYSIS: PEARSONS CORRELATIONS
########################
# In this alternative analysis we correlate EC/CC and activity metrics per subject to understand the within-subject relationships between these measures.
# This is an alternative to the Linear mixed models used above, to check for replicability. We use the Pearson correlation to calculate the correlation
# per subject.We calculate the Pearson correlation, Fisher z-score the correlations (per frequency band and density) and do a t-test between groups.

# %%


def make_corrs(df, group, metric, metric_activity, output_dir):
    """
    Description?


    Parameters
    ----------
    df : pd.DataFrame,
        Dataframe comtaining the data that should be correlated.
    group : str,
        Describes the group that is being investigated
    metric : str,
        network metric (e.g. EC, CC) that should be correlated to activity (e.g. offset)
    metric_activity : str,
        activity metric (e.g. offset) that should be correlated to network metrics (e.g. EC, CC)
    output_dir : str, 
        path to directory where correlations should be stored 

    Returns
    -------
    df_corrs : pd.DataFrame,
        dataframe containing the pearsons correlations

    df_corrs_z :pd.DataFrame,
        dataframe containing the Fisher z transformed correlations

    """

    # group the dataframe by subject
    df_grouped = df.groupby("sub")
    li = []
    for freq in ["delta", "theta", "lower_alpha"]:
        for d in ["02", "03"]:
            # per subject, correlate the activity metric with the network metric and extract the correlation coefficient
            corrs = pd.DataFrame(
                df_grouped.apply(
                    lambda x: pearsonr(
                        x[f"{metric}_{freq}_{d}_z"], x[f"{metric_activity}_z"]
                    )[0]))
                
            
            corrs.columns = [f"corrs_{freq}_{d}"]

            li.append(corrs)
            
    # put all correlations into one dataframe (all frequencies and densities) and save it
    df_corrs = pd.concat(li, axis=1)
    df_corrs.to_csv(
        f"{output_dir}20230803_df_{metric}_{metric_activity}_within_sub_corrs_{group}.csv"
    )

    # Fisher z-transform the correlations and save the dataframe
    df_corrs_z = pd.DataFrame(df_corrs.apply(lambda x: np.arctanh(x)))
    df_corrs_z.to_csv(
        f"{output_dir}20230803_df_{metric}_{metric_activity}_within_sub_corrs_z_{group}.csv"
    )

    return (df_corrs, df_corrs_z)


# %%
########################
# WHOLE GROUP PATIENTS ANALYSIS: Non-tumoral areas
########################
output_dir = "path/to/output_folder/"

df_corrs_patient_rest_CC, df_corrs_z_patient_rest_CC = make_corrs(
    df_non_tum, "patients_non_tum", "CC", "offset", output_dir
)
df_corrs_patient_rest_EC, df_corrs_z_patient_rest_EC = make_corrs(
    df_non_tum, "patients_non_tum", "EC", "offset", output_dir
)


########################
# HCs analysis
########################
df_corrs_HCs_CC, df_corrs_z_HCs_CC = make_corrs(df_HC, "HCs", "CC", "offset", output_dir)
df_corrs_HCs_EC, df_corrs_z_HCs_EC = make_corrs(df_HC, "HCs", "EC", "offset", output_dir)

# %%
########################
# SUBTYPE ANALYSIS
########################

# --- IDH-wildtype --- #
df_corrs_wt_rest_CC, df_corrs_z_wt_rest_CC = make_corrs(
    IDH_wt_non_tum, "IDH_wt", "CC", "offset", output_dir
)
df_corrs_wt_rest_EC, df_corrs_z_wt_rest_EC = make_corrs(
    IDH_wt_non_tum, "IDH_wt", "EC", "offset", output_dir
)

# --- IDH-mutant, 1p19q-codeleted --- #
df_corrs_codel_rest_CC, df_corrs_z_codel_rest_CC = make_corrs(
    IDH_codel_non_tum, "IDH_mut_codeleted", "CC", "offset", output_dir
)
df_corrs_codel_rest_EC, df_corrs_z_codel_rest_EC = make_corrs(
    IDH_codel_non_tum, "IDH_mut_codeleted", "EC", "offset", output_dir
)

# --- IDH-mutant, 1p19q-noncodeleted --- #
df_corrs_noncodel_rest_CC, df_corrs_z_noncodel_rest_CC = make_corrs(
    IDH_noncodel_non_tum, "IDH_mut_noncodeleted", "CC", "offset", output_dir
)
df_corrs_noncodel_rest_EC, df_corrs_z_noncodel_rest_EC = make_corrs(
    IDH_noncodel_non_tum, "IDH_mut_noncodeleted", "EC", "offset", output_dir
)


# %%
########################
# GROUP LEVEL CORRELATION PREPARATION
########################

# To obtain group-level correlation coefficients we used this tool: https://www.psychometrica.de/correlation.html (Section 8)
# Here, we calculate the number of regions involved in the correlation per patients (for the weighting)

# --- Non-tumoral areas --- #
nr_rois_non_tum = pd.DataFrame(df_non_tum["sub"].value_counts())
nr_rois_non_tum.reset_index(inplace=True)
nr_rois_non_tum.columns = ["sub", "nr_rois"]
nr_rois_non_tum.sort_values("sub", inplace=True)

#%%
# --- IDH-wildtype --- #
nr_rois_IDH_wt = pd.DataFrame(IDH_wt_non_tum["sub"].value_counts())
nr_rois_IDH_wt.reset_index(inplace=True)
nr_rois_IDH_wt.columns = ["sub", "nr_rois"]
nr_rois_IDH_wt.sort_values("sub", inplace=True)

# --- IDH-mutant, 1p19q-codeleted --- #
nr_rois_IDH_mut_codel = pd.DataFrame(IDH_codel_non_tum["sub"].value_counts())
nr_rois_IDH_mut_codel.reset_index(inplace=True)
nr_rois_IDH_mut_codel.columns = ["sub", "nr_rois"]
nr_rois_IDH_mut_codel.sort_values("sub", inplace=True)

# --- IDH-mutant, 1p19q-noncodeleted --- #
nr_rois_IDH_mut_noncodel = pd.DataFrame(IDH_noncodel_non_tum["sub"].value_counts())
nr_rois_IDH_mut_noncodel.reset_index(inplace=True)
nr_rois_IDH_mut_noncodel.columns = ["sub", "nr_rois"]
nr_rois_IDH_mut_noncodel.sort_values("sub", inplace=True)

# %%
########################
# SIGNIFICANCE OF CORRELATIONS
########################

# --- Calculate siginificance of Fisher transformed correlations against 0 (1 sample t-test)--- #

########################
# WHOLE GROUP PATIENTS ANALYSIS: Non-tumoral areas
########################

CC_t_non_tum, CC_p_non_tum = wilcoxon(df_corrs_z_patient_rest_CC)
df_CC_corr_against_0_non_tum = pd.concat(
    [pd.DataFrame(CC_t_non_tum), pd.DataFrame(CC_p_non_tum)], axis=1
)
df_CC_corr_against_0_non_tum.columns = ["t", "p"]

EC_t_non_tum, EC_p_non_tum = wilcoxon(df_corrs_z_patient_rest_EC)
df_EC_corr_against_0_non_tum = pd.concat(
    [pd.DataFrame(EC_t_non_tum), pd.DataFrame(EC_p_non_tum)], axis=1
)
df_EC_corr_against_0_non_tum.columns = ["t", "p"]

#%%
test = wilcoxon(df_corrs_z_patient_rest_CC["corrs_delta_02"])
# %%

########################
# HCs analysis
########################

CC_t_HCs, CC_p_HCs = wilcoxon(df_corrs_z_HCs_CC)
df_CC_corr_against_0_HCs = pd.concat(
    [pd.DataFrame(CC_t_HCs), pd.DataFrame(CC_p_HCs)], axis=1
)
df_CC_corr_against_0_HCs.columns = ["t", "p"]

EC_t_HCs, EC_p_HCs = wilcoxon(df_corrs_z_HCs_EC)
df_EC_corr_against_0_HCs = pd.concat(
    [pd.DataFrame(EC_t_HCs), pd.DataFrame(EC_p_HCs)], axis=1
)
df_EC_corr_against_0_HCs.columns = ["t", "p"]

# %%
########################
# SUBTYPE ANALYSIS
########################

# --- IDH-wildtype --- #
CC_t_wt_rest, CC_p_wt_rest = wilcoxon(df_corrs_z_wt_rest_CC)
df_CC_corr_against_0_wt = pd.concat(
    [pd.DataFrame(CC_t_wt_rest), pd.DataFrame(CC_p_wt_rest)], axis=1
)
df_CC_corr_against_0_wt.columns = ["t", "p"]

EC_t_wt_rest, EC_p_wt_rest = wilcoxon(df_corrs_z_wt_rest_EC)
df_EC_corr_against_0_wt = pd.concat(
    [pd.DataFrame(EC_t_wt_rest), pd.DataFrame(EC_p_wt_rest)], axis=1
)
df_EC_corr_against_0_wt.columns = ["t", "p"]


# --- IDH-mutant, 1p19q-codeleted --- #
CC_t_codel_rest, CC_p_codel_rest = wilcoxon(df_corrs_z_codel_rest_CC)
df_CC_corr_against_0_codel = pd.concat(
    [pd.DataFrame(CC_t_codel_rest), pd.DataFrame(CC_p_codel_rest)], axis=1
)
df_CC_corr_against_0_codel.columns = ["t", "p"]

EC_t_codel_rest, EC_p_codel_rest = wilcoxon(df_corrs_z_codel_rest_EC)
df_EC_corr_against_0_codel = pd.concat(
    [pd.DataFrame(EC_t_codel_rest), pd.DataFrame(EC_p_codel_rest)], axis=1
)
df_EC_corr_against_0_codel.columns = ["t", "p"]

# --- IDH-mutant, 1p19q-noncodeleted --- #
CC_t_noncodel_rest, CC_p_noncodel_rest = wilcoxon(df_corrs_z_noncodel_rest_CC)
df_CC_corr_against_0_noncodel = pd.concat(
    [pd.DataFrame(CC_t_noncodel_rest), pd.DataFrame(CC_p_noncodel_rest)], axis=1
)
df_CC_corr_against_0_noncodel.columns = ["t", "p"]

EC_t_noncodel_rest, EC_p_noncodel_rest = wilcoxon(df_corrs_z_noncodel_rest_EC)
df_EC_corr_against_0_noncodel = pd.concat(
    [pd.DataFrame(EC_t_noncodel_rest), pd.DataFrame(EC_p_noncodel_rest)], axis=1
)
df_EC_corr_against_0_noncodel.columns = ["t", "p"]

# %%
########################
# FDR CORRECTION
########################



output_dir = "path/to/output_folder/"

########################
# WHOLE GROUP PATIENTS ANALYSIS: Non-tumoral areas
########################
df_CC_corr_against_0_non_tum[
    "FDR_corrected_p"
] = statsmodels.stats.multitest.fdrcorrection(df_CC_corr_against_0_non_tum["p"])[1]
df_EC_corr_against_0_non_tum[
    "FDR_corrected_p"
] = statsmodels.stats.multitest.fdrcorrection(df_EC_corr_against_0_non_tum["p"])[1]
df_CC_corr_against_0_non_tum.to_csv(f"{output_dir}df_CC_patients_against_0.csv")
df_EC_corr_against_0_non_tum.to_csv(f"{output_dir}df_EC_patients_against_0.csv")


########################
# HCs ANALYSIS
########################
df_CC_corr_against_0_HCs["FDR_corrected_p"] = statsmodels.stats.multitest.fdrcorrection(
    df_CC_corr_against_0_HCs["p"]
)[1]
df_EC_corr_against_0_HCs["FDR_corrected_p"] = statsmodels.stats.multitest.fdrcorrection(
    df_EC_corr_against_0_HCs["p"]
)[1]
df_CC_corr_against_0_HCs.to_csv(f"{output_dir}df_CC_HCs_against_0.csv")
df_EC_corr_against_0_HCs.to_csv(f"{output_dir}df_EC_HCs_against_0.csv")

# %%

########################
# SUBTYPE ANALYSIS
########################

# --- IDH wildtype --- #
df_CC_corr_against_0_wt["FDR_corrected_p"] = statsmodels.stats.multitest.fdrcorrection(
    df_CC_corr_against_0_wt["p"]
)[1]
df_EC_corr_against_0_wt["FDR_corrected_p"] = statsmodels.stats.multitest.fdrcorrection(
    df_EC_corr_against_0_wt["p"]
)[1]
df_CC_corr_against_0_wt.to_csv(f"{output_dir}df_CC_wt_against_0.csv")
df_EC_corr_against_0_wt.to_csv(f"{output_dir}df_EC_wt_against_0.csv")

# --- IDH mutant, 1p19q-codeleted --- #
df_CC_corr_against_0_codel[
    "FDR_corrected_p"
] = statsmodels.stats.multitest.fdrcorrection(df_CC_corr_against_0_codel["p"])[1]
df_EC_corr_against_0_codel[
    "FDR_corrected_p"
] = statsmodels.stats.multitest.fdrcorrection(df_EC_corr_against_0_codel["p"])[1]
df_CC_corr_against_0_codel.to_csv(f"{output_dir}df_CC_codel_against_0.csv")
df_EC_corr_against_0_codel.to_csv(f"{output_dir}df_EC_codel_against_0.csv")

# --- IDH mutant, 1p19q-noncodeleted --- #
df_CC_corr_against_0_noncodel[
    "FDR_corrected_p"
] = statsmodels.stats.multitest.fdrcorrection(df_CC_corr_against_0_noncodel["p"])[1]
df_EC_corr_against_0_noncodel[
    "FDR_corrected_p"
] = statsmodels.stats.multitest.fdrcorrection(df_EC_corr_against_0_noncodel["p"])[1]
df_CC_corr_against_0_noncodel.to_csv(
    f"{output_dir}df_CC_noncodel_against_0.csv"
)
df_EC_corr_against_0_noncodel.to_csv(
   f"{output_dir}df_EC_noncodel_against_0.csv"
)


# %%
########################
# DIFFERENCE BETWEEN FISHERS Z TRANSFORMED CORRELATIONS OF PATIENTS AND HCS
########################

t_CC_comp_HCs_patients_rest, p_CC_comp_HCs_patients_rest = mannwhitneyu(
    df_corrs_z_patient_rest_CC, df_corrs_z_HCs_CC
)
df_CC_corr_comp = pd.concat(
    [
        pd.DataFrame(t_CC_comp_HCs_patients_rest),
        pd.DataFrame(p_CC_comp_HCs_patients_rest),
    ],
    axis=1,
)
df_CC_corr_comp.columns = ["t", "p"]

t_EC_comp_HCs_patients_rest, p_EC_comp_HCs_patients_rest = mannwhitneyu(
    df_corrs_z_patient_rest_EC, df_corrs_z_HCs_EC
)
df_EC_corr_comp = pd.concat(
    [
        pd.DataFrame(t_EC_comp_HCs_patients_rest),
        pd.DataFrame(p_EC_comp_HCs_patients_rest),
    ],
    axis=1,
)
df_EC_corr_comp.columns = ["t", "p"]

# %%
########################
# FDR CORRECTION
########################
df_CC_corr_comp["FDR_corrected_p"] = statsmodels.stats.multitest.fdrcorrection(
    df_CC_corr_comp["p"]
)[1]
df_EC_corr_comp["FDR_corrected_p"] = statsmodels.stats.multitest.fdrcorrection(
    df_EC_corr_comp["p"]
)[1]
df_CC_corr_comp.to_csv("path/to/file/where/fdr/correction/should/be/stored.csv")
df_EC_corr_comp.to_csv("path/to/file/where/fdr/correction/should/be/stored.csv")

# %%
########################
# RELATION BETWEEN PERITUMORAL ACTIVITY AND PEARSONS CORRELATIONS IN NON-TUMORAL REGIONS
########################

# import correlation dataframes
df_CC_corrs = pd.read_csv(
    "path/to/file/containing/correlations.csv"
)
df_EC_corrs = pd.read_csv(
    "path/to/file/containing/correlations.csv"
)


# %%


def make_corrs_offset_corr(df_corrs, df_offset, metric, output_dir, plot=True):
    """ Description?


    Parameters
    ----------
    df_corrs : pd.DataFrame,
        dataframe containing the correlations between networkmetrics (e.g. CC or EC) and offset
    df_offset : pd.DataFrame,
        dataframe containing the activity data (e.g. offset)
    metric : str,
        described which network metric was correlated with offset, for saving of dataframe
    plot : bool, optional
        Plotting of relationship, The default is True.

    Returns
    -------
    None.

    """
    #COMMENT OR UNCOMMENT (when want to investigate entire peritumoral area)
    # average offset_z over all regions in peritumoral area
    # df_offset_avg = pd.DataFrame(
    #     df_offset.groupby("sub")["offset_z"].mean()
    # )  

    #COMMENT OR UNCOMMENT (when want to investigate 3 highest values)
    # average only over 3 highest offset_z values
    df_offset_3_largest = pd.DataFrame(
        df_offset.groupby("sub")["offset_z"].nlargest(3)
    )  # comment/uncomment depending on which are you want to average the offset over
    
    df_offset_avg = pd.DataFrame(
        df_offset_3_largest.groupby("sub")["offset_z"].mean()
    )  # comment/uncomment depending on which are you want to average the offset over

    corrs_dict = {}
    for freq in ["lower_alpha"]:
        for d in ["02", "03"]:
            # correlate peritumoral offset with the correlation between network metric and offset
            corr, p = pearsonr(df_offset_avg["offset_z"], df_corrs[f"corrs_{freq}_{d}"])
            corrs_dict[f"corrs_{freq}_{d}"] = corr, p  # store correlation in dictionary

            if plot == True:
                plt.scatter(df_corrs[f"corrs_{freq}_{d}"], df_offset_avg["offset_z"])
                plt.xlabel(f"corrs_{freq}_{d}_{metric}")
                plt.ylabel("offset_z")
                plt.show()
                
    # store all correlations in one dataframe
    df_corrs = pd.DataFrame.from_dict(corrs_dict).T
    df_corrs.columns = ["corr", "p"]
    df_corrs.to_csv(f"/{output_dir}20231003_df_{metric}_offset_corr.csv")

    return df_corrs


# %%
#####################################
# IMPORT AND PREPARE DATAFRAMES CORRs
#####################################

#PATIENTS
#import correlation dataframes (patients non_tumoral)

#[pd.read_csv...]


#%% 
#As the above dataframes are based on the non-tumoral dataframe that includes all subjects 
#(also subjects with bilateral tumor and no 12% overlaps) we need to filter them to 
#the subjects that are in the peritumoral dataframe (unilateral, only subjects that have at least one region with 12% overlaps)

#PATIENTS
df_CC_corrs = df_CC_corrs[df_CC_corrs['sub'].isin(df_peritumoral['sub'])]
df_EC_corrs = df_EC_corrs[df_EC_corrs['sub'].isin(df_peritumoral['sub'])]

#PATIENTS SUBGROUP 
#IDH-codel peritumoral N = 11
df_CC_corrs_IDH_codel = df_CC_corrs_IDH_codel[df_CC_corrs_IDH_codel['sub'].isin(IDH_codel_peri['sub'])]
df_EC_corrs_IDH_codel = df_EC_corrs_IDH_codel[df_EC_corrs_IDH_codel['sub'].isin(IDH_codel_peri['sub'])]

#IDH-noncodel peritumoral N = 23
df_CC_corrs_IDH_noncodel = df_CC_corrs_IDH_noncodel[df_CC_corrs_IDH_noncodel['sub'].isin(IDH_noncodel_peri['sub'])]
df_EC_corrs_IDH_noncodel = df_EC_corrs_IDH_noncodel[df_EC_corrs_IDH_noncodel['sub'].isin(IDH_noncodel_peri['sub'])]

#IDH-wt peritumoral N = 24
df_CC_corrs_IDH_wt = df_CC_corrs_IDH_wt[df_CC_corrs_IDH_wt['sub'].isin(IDH_wt_peri['sub'])]
df_EC_corrs_IDH_wt = df_EC_corrs_IDH_wt[df_EC_corrs_IDH_wt['sub'].isin(IDH_wt_peri['sub'])]


#%%

########################
# WHOLE GROUP PATIENTS ANALYSIS: Non-tumoral areas
########################
output_dir = 'path/to/output/directory'

df_EC_corrs_offset = make_corrs_offset_corr(df_EC_corrs, df_peritumoral, "EC", output_dir, True) 
df_CC_corrs_offset = make_corrs_offset_corr(df_CC_corrs, df_peritumoral, "CC", output_dir, True)


# %%
########################
# SUBTYPE ANALYSIS
########################

# --- IDH-wildtype --- #
df_EC_corrs_offset_wt = make_corrs_offset_corr(
    df_EC_corrs_IDH_wt, IDH_wt_peri, "EC_wt", output_dir
)
df_CC_corrs_offset_wt = make_corrs_offset_corr(
    df_CC_corrs_IDH_wt, IDH_wt_peri, "CC_wt", output_dir
)

# --- IDH-mutant, 1p19q-codeleted --- #
df_EC_corrs_offset_mut_codel = make_corrs_offset_corr(
    df_EC_corrs_IDH_codel, IDH_codel_peri, "EC_mut_codel", output_dir
)
df_CC_corrs_offset_mut_codel = make_corrs_offset_corr(
    df_CC_corrs_IDH_codel, IDH_codel_peri, "CC_mut_codel", output_dir
)

# --- IDH mutant, 1p19q-noncoldeleted --- #
df_EC_corrs_offset_mut_noncodel = make_corrs_offset_corr(
    df_EC_corrs_IDH_noncodel, IDH_noncodel_peri, "EC_mut_noncodel", output_dir
)
df_CC_corrs_offset_mut_noncodel = make_corrs_offset_corr(
    df_CC_corrs_IDH_noncodel, IDH_noncodel_peri, "CC_mut_noncodel", output_dir
)

# %%

########################
# THREE HIGHEST ACTIVITY PERITUMORAL AREAS
########################

########################
# WHOLE GROUP PATIENTS ANALYSIS: Non-tumoral areas
########################
df_EC_corrs_offset = make_corrs_offset_corr(df_EC_corrs, df_peritumoral, "EC_3", output_dir)
df_CC_corrs_offset = make_corrs_offset_corr(df_CC_corrs, df_peritumoral, "CC_3", output_dir)


########################
# SUBTYPE ANALYSIS
########################

# --- IDH-wildtype --- #
df_EC_corrs_offset_wt = make_corrs_offset_corr(
    df_EC_corrs_IDH_wt, IDH_wt_peri, "EC_wt_3", output_dir
)
df_CC_corrs_offset_wt = make_corrs_offset_corr(
    df_CC_corrs_IDH_wt, IDH_wt_peri, "CC_wt_3", output_dir
)

# --- IDH-mutant, 1p19q-codeleted --- #
df_EC_corrs_offset_mut_codel = make_corrs_offset_corr(
    df_EC_corrs_IDH_codel, IDH_codel_peri, "EC_codel_3", output_dir
)
df_CC_corrs_offset_mut_codel = make_corrs_offset_corr(
    df_CC_corrs_IDH_codel, IDH_codel_peri, "CC_codel_3", output_dir
)

# --- IDH mutant, 1p19q-noncoldeleted --- #
df_EC_corrs_offset_mut_noncodel = make_corrs_offset_corr(
    df_EC_corrs_IDH_noncodel, IDH_noncodel_peri, "EC_noncodel_3", output_dir
)
df_CC_corrs_offset_mut_noncodel = make_corrs_offset_corr(
    df_CC_corrs_IDH_noncodel, IDH_noncodel_peri, "CC_noncodel_3", output_dir
)



#%%

########################
# PREPARATION SPINTEST
########################
def prepare_spintest(df, output_dir, name): 
    """
    Function to calculate the regional mean of the raw offset and network metrics.
    These are stored in dataframes per frequency and density in goal directory. 
    These Dataframes will be used in the Spin-test analysis.

    Parameters
    ----------
    df : pd.DataFrame,
        Dataframe storing the raw offset and network data
    output_dir : str,
        path to directory where dataframes with regional means should be stored
    name : str, 
        used to describe which group was analysed (make output file name clear)

    Returns
    -------
    None.

    """
    
    #calculate the regional average
    grouped_by_roi = df.groupby("roi").mean()  
    
    for freq in ["delta", "theta", "lower_alpha"]:
        
        for d in ["02", "03"]:
            
            # extract the needed columns per frequency and density and store in new dataframe
            # first col is offset, second col is network metric
            df_CC = grouped_by_roi[["Offset", f"CC_{freq}_{d}"]]
            df_CC.to_csv(
                f"{output_dir}df_CC_{freq}_{d}_spintest_{name}.csv",
                header=False,
                index=False,
            ) 
            
            # extract the needed columns per frequency and density and store in new dataframe
            # first col is offset, second col is network metric
            
            df_EC = grouped_by_roi[["Offset", f"EC_{freq}_{d}"]]
            df_EC.to_csv(
                f"{output_dir}df_EC_{freq}_{d}_spintest_{name}.csv",
                header=False,
                index=False,
            )  
            
#%%
output_dir = '/path/to/output/dir'

prepare_spintest(df_non_tum, output_dir, "patients_non_tum")
prepare_spintest(df_HC, output_dir, "HCs")