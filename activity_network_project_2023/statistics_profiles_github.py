#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Profiles of activity and network metrics

Script for main analysis: statistical analysis concerning the profiles of (disturbances in) network metrics and activity. 



"""

__author__ = "Mona Lilo Margarethe Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "2022/02/15"
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
from datetime import datetime
import re

# Third party imports ###
import pandas as pd  # version 1.1.5
import matplotlib.pyplot as plt # version 3.6.2
import statsmodels.api as sm # version 0.13.2
import statsmodels.formula.api as smf # version 0.13.2
import scipy.stats as stats # version 1.9.3
from scipy.stats import wilcoxon # version 1.9.3
from scipy.stats import mannwhitneyu # version 1.9.3
import seaborn as sns # version 0.11.0
#%%
####################
# DATAFRAME IMPORT #
####################

#PATIENTS

#import correlation dataframes (patients non_tumoral)

#[pd.read_csv...]
#df_all, df_peritumoral, df_homologue, df_non_tum, df_HC, Subtype dataframes


# %%
###################
# --- Statistical analysis of averaged areas in patients vs. HC averaged whole brain --- #
###################


def calc_area_vs_HC(df_pat, df_HC, metric, area, output_dir, activity=False):
    """
    Function to compare the peritumoral/non-peritumoral homologue and rest of
    the brain with the whole brain HC network/activity metrics using a Mann Whitney U test.
    These metrics are averaged over the areas: peritumoral area, non-peritumoral
    contralateral homologue, all non-peritumoral areas and the whole brain (in HC).
    Saves the test statistics and p-values subsequently in a dataframe.

    Parameters
    ----------
    df_pat : pd.DataFrame,
        df contining relevant data for patients
    df_HC : pd.DataFrame,
        df contining relevant data for HCs
    metric : str,
        network metric to be calculates (CC, clustering coefficient or EC, eigenvector centrality)
    area : str,
        the area that should be compared to whole brain HC metrics
    output_dir : str,
        path to directory where statistics should be saved
    activity : boo,
        boolean to set whether activity (e.g. BB) data is being analysed. Default is False.

    Returns
    -------
    None.

    """
    now = datetime.now()

    # Network comparisons
    if activity == False:
        dict_stat = {}

        for freq in ["theta", "delta", "lower_alpha"]:
            print(freq)
            for d in ["02", "03"]:
                print(d)
                # average the network metrics per subject
                avg_pat = df_pat.groupby("sub")[f"{metric}_{freq}_{d}_z"].mean()
                avg_HC = df_HC.groupby("sub")[f"{metric}_{freq}_{d}_z"].mean()

                # Mann-Whitney U test
                stat, p = mannwhitneyu(avg_pat, avg_HC)
                dict_stat[f"{freq}_{d}"] = [stat, p]  # store results
    
    # Activity comparisons
    else:
        dict_stat = {}

        # average the activity metric per subject
        avg_pat = df_pat.groupby("sub")[f"{metric}_z"].mean()
        avg_HC = df_HC.groupby("sub")[f"{metric}_z"].mean()

        # Mann-Whitney U test
        stat, p = mannwhitneyu(avg_pat, avg_HC)
        dict_stat[area] = [stat, p]  # store results
    
    # store all results in dataframe
    df_stats = pd.DataFrame.from_dict(dict_stat, orient="index")
    # return(df_stats)

    # save dataframe
    df_stats.to_csv(
        output_dir
        + now.strftime("%Y%m%d")
        + "_pat_vs_HC_"
        + metric
        + "_"
        + area
        + ".csv"
    )


# %%
########################
# WHOLE GROUP ANALYSIS
########################

######
# Comparison activity (offset_z)
######

output_dir = "/path/to/output_folder/"

calc_area_vs_HC(
    df_peritumoral, df_HC, "offset", "peritumor", output_dir, activity=True
)  # peritumoral area vs. whole brain HCs

calc_area_vs_HC(
    df_homologue, df_HC, "offset", "homologue", output_dir, activity=True
)  # contralateral homologue vs. whole brain HCs

calc_area_vs_HC(
    df_non_tum, df_HC, "offset", "non_tum", output_dir, activity=True
)  # rest of the brain (all non-tumoral areas) vs. wholebrain HCs

# %%
######
# Comparison network metrics (EC_z, CC_z)
######
calc_area_vs_HC(
    df_peritumoral, df_HC, "CC", "peritumor", output_dir, activity=False
)  # peritumoral area vs. whole brain HCs
calc_area_vs_HC(
    df_peritumoral, df_HC, "EC", "peritumor", output_dir, activity=False
)  # peritumoral area vs. whole brain HCs

calc_area_vs_HC(
    df_homologue, df_HC, "CC", "homologue", output_dir, activity=False
)  # contralateral homologue vs. whole brain HCs
calc_area_vs_HC(
    df_homologue, df_HC, "EC", "homologue", output_dir, activity=False
)  # contralateral homologue vs. whole brain HCs

calc_area_vs_HC(
    df_non_tum, df_HC, "CC", "non_tum", output_dir, activity=False
)  # rest of the brain (all non-tumoral areas) vs. wholebrain HCs
calc_area_vs_HC(
    df_non_tum, df_HC, "EC", "non_tum", output_dir, activity=False
)  # rest of the brain (all non-tumoral areas) vs. wholebrain HCs

# %%
########################
# SUBTYPE ANALYSIS
########################

output_dir = '/path/to/output/folder'


######
# Comparison activity (offset_z)
######

# IDH wildtype
calc_area_vs_HC(
    IDH_wt_peri,
    df_HC,
    "offset",
    "peritumor_IDH_wt",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_wt_hom,
    df_HC,
    "offset",
    "homologue_IDH_wt",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_wt_non_tum, df_HC, "offset", "rest_IDH_wt", output_dir, activity=True
)

# IDH mutant codeleted
calc_area_vs_HC(
    IDH_codel_peri,
    df_HC,
    "offset",
    "peritumor_IDH_mut_codel",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_codel_hom,
    df_HC,
    "offset",
    "homologue_IDH_mut_codel",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_codel_non_tum,
    df_HC,
    "offset",
    "rest_IDH_mut_codel",
    output_dir,
    activity=True,
)

# IDH mutant non-codeleted
calc_area_vs_HC(
    IDH_noncodel_peri,
    df_HC,
    "offset",
    "peritumor_IDH_mut_noncodel",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_noncodel_hom,
    df_HC,
    "offset",
    "homologue_IDH_mut_noncodel",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_noncodel_non_tum,
    df_HC,
    "offset",
    "rest_IDH_mut_noncodel",
    output_dir,
    activity=True,
)

#%%
######
# Comparison network metrics (EC_z, CC_z)
######

#IDH wildtype


calc_area_vs_HC(
    IDH_wt_non_tum, df_HC, "EC", "rest_IDH_wt", output_dir, activity=False
)


calc_area_vs_HC(
    IDH_wt_non_tum, df_HC, "CC", "rest_IDH_wt", output_dir, activity=False
)

# IDH mutant codeleted
calc_area_vs_HC(
    IDH_codel_non_tum,
    df_HC,
    "EC",
    "rest_IDH_mut_codel",
    output_dir,
    activity=False,
)


calc_area_vs_HC(
    IDH_codel_non_tum,
    df_HC,
    "CC",
    "rest_IDH_mut_codel",
    output_dir,
    activity=False,
)


# IDH mutant non-codeleted

calc_area_vs_HC(
    IDH_noncodel_non_tum,
    df_HC,
    "EC",
    "rest_IDH_mut_noncodel",
    output_dir,
    activity=False,
)

calc_area_vs_HC(
    IDH_noncodel_non_tum,
    df_HC,
    "CC",
    "rest_IDH_mut_noncodel",
    output_dir,
    activity=False,
)


# %%
### --- Wilcoxon signed rank test to test differences in local clustering and eigenvector centrality and activity metrics between peritumoral and contralateral homologue --- ###
# on averaged data --> paired data of patientsn (wilcoxon signed rank test); as in Clustering paper


def calc_peri_vs_hom(
    df_pat_peri, df_pat_hom, metric, output_dir, name, activity=False
):  
    """
    Function to statistically test the difference between a network or activity
    metric in the peritumoral and peritumoral homologue area (averaged over these areas).
    Test statistics are being saved in a csv file.

    Parameters
    ----------
    df_pat_peri : pd.DataFrame,
        df containing data of the peritumoral area of patients
    df_pat_hom : pd.DataFrame,
        df containing data of the peritumoral contralateral homologue area of patients
    metric : str,
        determines which metric should be investigated (CC,EC,offset)
    output_dir : str,
        path to directory in which statistics should be stored.
    name: str,
        name that is added to the saved csv
    activity : bool, optional
        set to determine if activity or network metric is being investigated. The default is False (i.e. network metrics will be investigated).

    Returns
    -------
    None.

    """

    now = datetime.now()
    dict_stat = {}
    if activity == False:
        for freq in ["theta", "delta", "lower_alpha"]:
            print(freq)
            for d in ["02", "03"]:
                print(d)

                # average the network metric per subject
                avg_peritumoral = df_pat_peri.groupby("sub")[
                    f"{metric}_{freq}_{d}_z"
                ].mean()
                avg_homologue = df_pat_hom.groupby("sub")[
                    f"{metric}_{freq}_{d}_z"
                ].mean()

                # Wilcoxon-signed rank test
                stat, p = wilcoxon(avg_peritumoral, avg_homologue)

                dict_stat[freq + "_" + d] = [stat, p]  # store results
    else:
        # average the activity metric per subject
        avg_peritumoral = df_pat_peri.groupby("sub")[
            metric + "_z"
        ].mean()  # changed for new analysis with subtypes
        avg_homologue = df_pat_hom.groupby("sub")[
            metric + "_z"
        ].mean()  # changed for new analysis with subtypes

        # Wilcoxon-signed rank test
        stat, p = wilcoxon(avg_peritumoral, avg_homologue)

        dict_stat[metric] = [stat, p]  # store results
    # store all results in dataframe
    df_stats = pd.DataFrame.from_dict(dict_stat, orient="index")

    # save dataframe
    df_stats.to_csv(
        output_dir
        + now.strftime("%Y%m%d")
        + "_peri_vs_hom_"
        + name
        + metric
        + ".csv"
    )


# %%
output_dir = "/path/to/output_folder/"

# %%

########################
# WHOLE GROUP ANALYSIS
########################

######
# Comparison activity (offset_z) and network metrics between peritumoral area and homologue
######

calc_peri_vs_hom(df_peritumoral, df_homologue, "CC", output_dir, "_", activity=False)
calc_peri_vs_hom(df_peritumoral, df_homologue, "EC", output_dir, "_", activity=False)
calc_peri_vs_hom(df_peritumoral, df_homologue, "offset", output_dir, "_", activity=True)


# %%
########################
# SUBTYPE ANALYSIS
########################

######
# Comparison activity (offset_z) between peritumoral area and homologue for the different tumor subtypes
######

output_dir = 'path/to/output/folder/'

# IDH wt
calc_peri_vs_hom(
    IDH_wt_peri, IDH_wt_hom, "offset", output_dir, "IDH_wt_", activity=True
)

# IDH mutant codeleted
calc_peri_vs_hom(
    IDH_codel_peri,
    IDH_codel_hom,
    "offset",
    output_dir,
    "IDH_mut_codeleted_",
    activity=True,
)

# IDH mutant non-codeleted
calc_peri_vs_hom(
    IDH_noncodel_peri,
    IDH_noncodel_hom,
    "offset",
    output_dir,
    "IDH_mut_noncodeleted_",
    activity=True,
)


#%%

########################
# MEANS AND STDs 
########################

def get_means_std(df_orig, output_dir, name):
    """
    Function to get the descriptive stats for every column in the dataframes. 
    These are stored in a pd.DataFrame and saved in the desired directory.

    Parameters
    ----------
    df_orig : pd.DataFrame,
        dataframe for which descriptives should be calculated
    output_dir : str,
        path to directory where descriptives pd.DataFrames should be stored
    name : str,
        name for the file to identify which df was analysed

    Returns
    -------
    None.

    """
    
    descriptives = df_orig.describe()
    
    descriptives.to_csv(f"{output_dir}{name}_descriptives.csv")
    

#%%

output_dir = "/path/to/output_folder/"

########################
# Main dataframes
########################

get_means_std(df_peritumoral, output_dir, "df_peri")
get_means_std(df_homologue, output_dir, "df_hom")
get_means_std(df_non_tum, output_dir, "df_non_tum")
get_means_std(df_HC, output_dir, "df_HC")

########################
# Patients subtypes
########################

#IDH-mut codeleted
get_means_std(IDH_codel_peri, output_dir, "df_IDH_codel_peri")
get_means_std(IDH_codel_hom, output_dir, "df_IDH_codel_hom")
get_means_std(IDH_codel_non_tum, output_dir, "df_IDH_codel_non_tum")

#IDH-mut noncodel
get_means_std(IDH_noncodel_peri, output_dir, "df_IDH_noncodel_peri")
get_means_std(IDH_noncodel_hom, output_dir, "df_IDH_noncodel_hom")
get_means_std(IDH_noncodel_non_tum, output_dir, "df_IDH_noncodel_non_tum")

#IDH-wt
get_means_std(IDH_wt_peri, output_dir, "df_IDH_wt_peri")
get_means_std(IDH_wt_hom, output_dir, "df_IDH_wt_hom")
get_means_std(IDH_wt_non_tum, output_dir, "df_IDH_wt_non_tum")







