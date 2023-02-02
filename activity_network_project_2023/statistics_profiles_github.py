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

# Third party imports ###
import pandas as pd  # version 1.1.5
import numpy as np  ## not used?
import matplotlib.pyplot as plt
import re
import statsmodels.api as sm
import statsmodels.formula.api as smf
import scipy.stats as stats
from scipy.stats import wilcoxon
from scipy.stats import mannwhitneyu
import seaborn as sns

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

        for freq in ["theta", "delta", "alpha"]:
            print(freq)
            for d in ["02", "03"]:
                print(d)
                # average the network metrics per subject
                avg_pat = df_pat.groupby("sub")[metric + "_z_" + freq + "_" + d].mean()
                avg_HC = df_HC.groupby("sub")[metric + "_z_" + freq + "_" + d].mean()

                # Mann-Whitney U test
                stat, p = mannwhitneyu(avg_pat, avg_HC)
                dict_stat[freq + "_" + d] = [stat, p]  # store results
    # Activity comparisons
    else:
        dict_stat = {}

        # average the activity metric per subject
        avg_pat = df_pat.groupby("sub")[metric + "_z"].mean()
        avg_HC = df_HC.groupby("sub")[metric + "_z"].mean()

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
        + "_desired_file_name_"
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

output_dir = "/path/to/output/folder/"

calc_area_vs_HC(
    df_peritumoral, df_master_HC, "offset", "peritumor", output_dir, activity=True
)  # peritumoral area vs. whole brain HCs

calc_area_vs_HC(
    df_homologue, df_master_HC, "offset", "homologue", output_dir, activity=True
)  # contralateral homologue vs. whole brain HCs

calc_area_vs_HC(
    df_non_tum, df_master_HC, "offset", "non_tum", output_dir, activity=True
)  # rest of the brain (all non-tumoral areas) vs. wholebrain HCs

# %%
######
# Comparison network metrics (EC_z, CC_z)
######
calc_area_vs_HC(
    df_peritumoral, df_master_HC, "CC", "peritumor", output_dir, activity=False
)  # peritumoral area vs. whole brain HCs
calc_area_vs_HC(
    df_peritumoral, df_master_HC, "EC", "peritumor", output_dir, activity=False
)  # peritumoral area vs. whole brain HCs

calc_area_vs_HC(
    df_homologue, df_master_HC, "CC", "homologue", output_dir, activity=False
)  # contralateral homologue vs. whole brain HCs
calc_area_vs_HC(
    df_homologue, df_master_HC, "EC", "homologue", output_dir, activity=False
)  # contralateral homologue vs. whole brain HCs

calc_area_vs_HC(
    df_non_tum, df_master_HC, "CC", "non_tum", output_dir, activity=False
)  # rest of the brain (all non-tumoral areas) vs. wholebrain HCs
calc_area_vs_HC(
    df_non_tum, df_master_HC, "EC", "non_tum", output_dir, activity=False
)  # rest of the brain (all non-tumoral areas) vs. wholebrain HCs

# %%
########################
# SUBTYPE ANALYSIS
########################

# output_dir = '/path/to/output/folder'

######
# Comparison activity (offset_z)
######

# IDH wildtype
calc_area_vs_HC(
    IDH_wt_peritumoral_all,
    df_master_HC,
    "offset",
    "peritumor_IDH_wt",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_wt_homologue,
    df_master_HC,
    "offset",
    "homologue_IDH_wt",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_wt_rest, df_master_HC, "offset", "rest_IDH_wt", output_dir, activity=True
)

# IDH mutant codeleted
calc_area_vs_HC(
    IDH_mut_codeleted_peritumoral_all,
    df_master_HC,
    "offset",
    "peritumor_IDH_mut_codel",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_mut_codeleted_homologue,
    df_master_HC,
    "offset",
    "homologue_IDH_mut_codel",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_mut_codeleted_rest,
    df_master_HC,
    "offset",
    "rest_IDH_mut_codel",
    output_dir,
    activity=True,
)

# IDH mutant non-codeleted
calc_area_vs_HC(
    IDH_mut_noncodeleted_peritumoral_all,
    df_master_HC,
    "offset",
    "peritumor_IDH_mut_noncodel",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_mut_noncodeleted_homologue,
    df_master_HC,
    "offset",
    "homologue_mut_noncodel",
    output_dir,
    activity=True,
)
calc_area_vs_HC(
    IDH_mut_noncodeleted_rest,
    df_master_HC,
    "offset",
    "rest_IDH_mut_noncodel",
    output_dir,
    activity=True,
)


######
# Comparison network metrics (EC_z, CC_z)
######

# IDH wildtype
calc_area_vs_HC(
    IDH_wt_peritumoral_all,
    df_master_HC,
    "EC",
    "peritumor_IDH_wt",
    output_dir,
    activity=False,
)
calc_area_vs_HC(
    IDH_wt_homologue, df_master_HC, "EC", "homologue_IDH_wt", output_dir, activity=False
)
calc_area_vs_HC(
    IDH_wt_rest, df_master_HC, "EC", "rest_IDH_wt", output_dir, activity=False
)

calc_area_vs_HC(
    IDH_wt_peritumoral_all,
    df_master_HC,
    "CC",
    "peritumor_IDH_wt",
    output_dir,
    activity=False,
)
calc_area_vs_HC(
    IDH_wt_homologue, df_master_HC, "CC", "homologue_IDH_wt", output_dir, activity=False
)
calc_area_vs_HC(
    IDH_wt_rest, df_master_HC, "CC", "rest_IDH_wt", output_dir, activity=False
)


# IDH mutant codeleted
calc_area_vs_HC(
    IDH_mut_codeleted_peritumoral_all,
    df_master_HC,
    "EC",
    "peritumor_IDH_mut_codel",
    output_dir,
    activity=False,
)
calc_area_vs_HC(
    IDH_mut_codeleted_homologue,
    df_master_HC,
    "EC",
    "homologue_IDH_mut_codel",
    output_dir,
    activity=False,
)
calc_area_vs_HC(
    IDH_mut_codeleted_rest,
    df_master_HC,
    "EC",
    "rest_IDH_mut_codel",
    output_dir,
    activity=False,
)

calc_area_vs_HC(
    IDH_mut_codeleted_peritumoral_all,
    df_master_HC,
    "CC",
    "peritumor_IDH_mut_codel",
    output_dir,
    activity=False,
)
calc_area_vs_HC(
    IDH_mut_codeleted_homologue,
    df_master_HC,
    "CC",
    "homologue_IDH_mut_codel",
    output_dir,
    activity=False,
)
calc_area_vs_HC(
    IDH_mut_codeleted_rest,
    df_master_HC,
    "CC",
    "rest_IDH_mut_codel",
    output_dir,
    activity=False,
)


# IDH mutant non-codeleted
calc_area_vs_HC(
    IDH_mut_noncodeleted_peritumoral_all,
    df_master_HC,
    "EC",
    "peritumor_IDH_mut_noncodel",
    output_dir,
    activity=False,
)
calc_area_vs_HC(
    IDH_mut_noncodeleted_homologue,
    df_master_HC,
    "EC",
    "homologue_mut_noncodel",
    output_dir,
    activity=False,
)
calc_area_vs_HC(
    IDH_mut_noncodeleted_rest,
    df_master_HC,
    "EC",
    "rest_IDH_mut_noncodel",
    output_dir,
    activity=False,
)

calc_area_vs_HC(
    IDH_mut_noncodeleted_peritumoral_all,
    df_master_HC,
    "CC",
    "peritumor_IDH_mut_noncodel",
    output_dir,
    activity=False,
)
calc_area_vs_HC(
    IDH_mut_noncodeleted_homologue,
    df_master_HC,
    "CC",
    "homologue_mut_noncodel",
    output_dir,
    activity=False,
)
calc_area_vs_HC(
    IDH_mut_noncodeleted_rest,
    df_master_HC,
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
):  # df_pat, metric, output_dir, activity = False): #changed for analysis with subtypes
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
        for freq in ["theta", "delta", "alpha"]:
            print(freq)
            for d in ["02", "03"]:
                print(d)

                # average the network metric per subject
                avg_peritumoral = df_pat_peri.groupby("sub")[
                    f"{metric}_z_{freq}_{d}"
                ].mean()
                avg_homologue = df_pat_hom.groupby("sub")[
                    f"{metric}_z_{freq}_{d}"
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
        + "_desired_file_name_"
        + name
        + metric
        + ".csv"
    )


# %%
output_dir = "/path/to/output/folder/"

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

# output_dir = 'path/to/output/folder/'

# IDH wt
calc_peri_vs_hom(
    IDH_wt_peritumoral, IDH_wt_homologue, "offset", output_dir, "IDH_wt_", activity=True
)

# IDH mutant codeleted
calc_peri_vs_hom(
    IDH_mut_codeleted_peritumoral,
    IDH_mut_codeleted_homologue,
    "offset",
    output_dir,
    "IDH_mut_codeleted_",
    activity=True,
)

# IDH mutant non-codeleted
calc_peri_vs_hom(
    IDH_mut_noncodeleted_peritumoral,
    IDH_mut_noncodeleted_homologue,
    "offset",
    output_dir,
    "IDH_mut_noncodeleted_",
    activity=True,
)


# %%
########################
# FDR CORRECTION
########################

# Peritumor vs. Homologue
import glob
import statsmodels.stats.multitest


# files = glob.glob('/path/to/folder/results_to_correct/*')
files = glob.glob("/path/to/folder/results_to_correct/*.csv")

for file in files:
    df = pd.read_csv(file, header=None)
    df.columns = ["metric", "statistic", "p-value"]

    # FDR corrected based on the column p-value
    df["FDR_corrected"] = statsmodels.stats.multitest.fdrcorrection(df["p-value"])[1]

    # FDR correction with only 2 densities
    list_2_den = [
        "delta_02",
        "delta_03",
        "theta_02",
        "theta_03",
        "alpha_02",
        "alpha_03",
    ]
    df_2_den = df.loc[df["metric"].isin(list_2_den)]
    df_2_den["FDR_corrected_only_2_den"] = statsmodels.stats.multitest.fdrcorrection(
        df_2_den["p-value"]
    )[1]

    # put both corrections together in one dataframe
    concatenated = pd.concat([df, df_2_den["FDR_corrected_only_2_den"]], axis=1)
    concatenated.to_csv(file)
