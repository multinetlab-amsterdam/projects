#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 10:10:05 2025

Functions used for the statistical analysis.

@author: mlzimmermann
"""

__author__ = "Mona Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "28-10-2025"   ### Date it was created
__status__ = "Production" ### Production = still being developed. Else: Concluded/Finished.

#Created with the help of chatgpt

####################
# Review History   #
####################

# Reviewed by Name Date ### 

####################
# Libraries        #
####################

# Standard imports  ### (Put here built-in libraries - https://docs.python.org/3/library/)


# Third party imports ### (Put here third-party libraries e.g. pandas, numpy)
import numpy as np 
import pandas as pd
import scipy.io
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
import nibabel as nib
from scipy.stats import spearmanr
import re
import pingouin as pg
from lifelines import CoxPHFitter
import io
import statsmodels.api as sm
import statsmodels.stats.api as sms
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy import stats
import contextlib


#%%

###################################################################
###             GENERAL STATISTICS  FUNCTIONS                   ###
###################################################################

def independent_ttest(df, groups_to_test, outcome_var, output_file, activity_kps_split=True):
    """
    Independent t-test. Tests the difference in the outcome_var between patients with high vs. low 
    in the group to test. If Normality is violated defaults to the Mann-Whitney U test. In case of Unequal variances
    switches to Welchs test.

    Parameters
    ----------
    df : pd.DataFrame, 
        dataframe conatining the median split grouping (low vs. high)
        and the outcome values per patient.
    output_file : str, optional
        path to location where the output file should be stored.
    activity_kps_split : bool, 
        determines whether we do analysis for activity or KPS differences as there groups to test are coded as Low or High.
        For other (clinical) analyses variables are often coded as 0 or 1, then we need to give False. Default = True.
        

    Returns
    -------
    results_df : pd.DataFrame,
        dataframe with the results. Also shows results of normality and equal variance test.

    """
    #safety measure: copy the dataframe
    df_local = df.copy()
    
    # Split the data based on tumor activity groups
    if activity_kps_split == True: 
        low = df_local[df_local[f"{groups_to_test}"] == "Low"][f"{outcome_var}"].copy()
        high = df_local[df_local[f"{groups_to_test}"] == "High"][f"{outcome_var}"].copy()
    else: 
        low = df_local[df_local[f"{groups_to_test}"] == 0][f"{outcome_var}"].copy()
        high = df_local[df_local[f"{groups_to_test}"] == 1][f"{outcome_var}"].copy()
        
    # Calculate descriptive statistics for both groups
    low_mean = low.mean()
    high_mean = high.mean()
    low_median = low.median()
    high_median = high.median()
    low_std = low.std()
    high_std = high.std()
    # median_split  = df['mean_bbp_z'].median()

    # Collect the results in a dictionary
    results = {}
    results['Low Sample Size'] = len(low)
    results['High Sample Size'] = len(high)
    results['Low Mean'] = low_mean
    results['High Mean'] = high_mean
    results['Low Median'] = low_median
    results['High Median'] = high_median
    results['Low Std Dev'] = low_std
    results['High Std Dev'] = high_std
    #results['Median whole group (median split at)'] = median_split
    
    # 1. Check for normality using the Shapiro-Wilk test
    low_normality = stats.shapiro(low)
    high_normality = stats.shapiro(high)
    results['Low Normality p-value'] = low_normality.pvalue
    results['High Normality p-value'] = high_normality.pvalue
    
    # 2. Check for equal variances using Levene's test
    levene_test = stats.levene(low, high)
    results['Levene’s Test p-value'] = levene_test.pvalue
    
    # Determine if normality is violated
    normality_violated = low_normality.pvalue < 0.05 or high_normality.pvalue < 0.05
    
    # 3. Conduct the appropriate test based on normality
    if normality_violated:
        # Use Mann-Whitney U test (non-parametric)
        u_stat, p_value = stats.mannwhitneyu(low, high, alternative='two-sided', method='exact')
        results['Test Type'] = 'Mann-Whitney U Test'
        results['U-statistic'] = u_stat
        results['U-test p-value'] = p_value
        results['Significant Difference'] = 'Yes' if p_value < 0.05 else 'No'
        
    else:
        # Conduct the two-sample t-test (parametric)
        
        # Determine if equal variance can be assumed (p > 0.05 for Levene's test)
        equal_variance_assumed = levene_test.pvalue > 0.05
        
        if equal_variance_assumed:
            test_type = 'Two-Sample T-Test (Equal Variance)'
            assumption_violated = 'None'
        else:
            test_type = 'Welch’s T-Test (Unequal Variance)'
            assumption_violated = 'Homogeneity of Variance'
            
        t_stat, p_value = stats.ttest_ind(low, high, equal_var=equal_variance_assumed, alternative='two-sided')
        
        results['Test Used'] = test_type
        results['Test Statistic'] = t_stat
        results['P-Value'] = p_value
        results['Assumption Violated'] = assumption_violated

    # Save the results to a CSV file (append if file already exists)
    results_df = pd.DataFrame([results])
    
    
    # --- Write all to Excel ---
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        results_df.to_excel(writer, sheet_name='independent t-test', index=False)

    print(f"Results saved to {output_file}")
    
    return results_df
    

def spearman(df,x,y, output_file):
    """
    Spearmanns rho. Tests the correlation between two vectors and saves the results in 
    a excel file. Used in normative LTDI analysis, TDM/tumor occurrence analysis,
    

    Parameters
    ----------
    df : pd.DataFrame, 
        dataframe conatining the tumoral activity and L-TDI per patient.
    x : str,
        name of column for the x variable to correlate.
    y : str,
        name of column for the y variable to correlate.
    output_file : str,
        path to file where results should be stored.Includes path and filename.
        

    Returns
    -------
    results_df : pd.DataFrame,
        contains the results of the spearman rho correlation.

    """
    
    #safety: copy input dataframe
    df_local = df.copy()
    
    x = df_local[x]
    print(x)
    y = df_local[y]
    
    #Collect the results in a dictionary
    results = {}
    
    #calculate the spearmanr
    r_corr, p_value = spearmanr(x, y)
    print(p_value)
    print(f"Exact p-value: {p_value:.15f}")
    
    #store results in dictionary
    results["Test_Type"] = "Spearmanns r"
    results["rho"] = r_corr
    results["p_value"] = str(p_value)
    results['Significant Difference'] = 'Yes' if p_value < 0.05 else 'No'
    
    
    
    # Save the results to a CSV file (append if file already exists)
    results_df = pd.DataFrame([results])

    print(f"Results saved to {output_file}")
    
    # --- Write all to Excel   ---
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        results_df.to_excel(writer, sheet_name='Spearmans rho', index=False)
    
    return results_df
    


   
############################################################################
###                        SURVIVAL ANALYSIS                             ###
############################################################################
def COXPH(df, metric, duration, event, covariates = False):
    """
    Function to model the cox proportional hazard for known covariates and LTDI or high tumor connected regions. 
    Covariates are included if covariates is True,(sex, age at diagnosis and KPS), which are dummy coded or prepared for analysis.
    Also included is the whole-brain excitability (mean_slope_z). For the analysis we recode the variable with 
    more positive values indicating higher excitability in this script for better interpretability of results.

    Parameters
    ----------
    df : pd.DataFrame,
        dataframe with relevant data. In function, we only retain the relevant columns.
    duration : str, 
        name of the duration column (e.g. duration progression_OK)
    event : str,
        name of the event column (e.g. progressie)
    covariates : bool, 
        determine if want to include covariates sex, age and kps in analysis or not. 
        Default is False.

    Returns
    -------
    cph : cox proportional hazards model, 
        summary of cox-proportional hazards analysis.

    """
    #safety: copy the original dataframe
    df_local = df.copy()
    
    df_local.reset_index(drop= True, inplace = True)
    print(df_local.head)
    
    # --- prepare covariates --- #
    if covariates == True: 
        print("Include covariates age, sex and kps in analysis!")
        
        covariates = ["sex_binary", "age_at_diagnosis", "kps_binary"]
        
        
        # --- prepare df --- #
        df_cox = df_local[[duration, event, f"{metric}"] + covariates].copy()
        
        
        # --- fit model --- #
        cph = CoxPHFitter()
        cph.fit(df_cox, duration_col=duration, event_col=event)
        cph.print_summary()
        
    elif covariates == False:
        print("No covariates included in analysis!")
         #--- prepare df --- #
        df_cox = df_local[[duration, event,f"{metric}"]].copy()
        
        
        # --- fit model --- #
        cph = CoxPHFitter()
        cph.fit(df_cox, duration_col=duration, event_col=event)
        cph.print_summary()
       
    else: 
        print("!!! Should covariates be included? Not properly defined !!!")
        
   
    # --- Assumption tests --- #
    print("\n--- Proportional hazards assumption tests ---")
    buffer = io.StringIO()
    try:
        with contextlib.redirect_stdout(buffer):
            cph.check_assumptions(df_cox, p_value_threshold=0.05, show_plots=False)
        assumption_text = buffer.getvalue()
    except Exception as e:
        assumption_text = f"Assumption test failed:\n{e}"


    # Convert text to simple DataFrame for Excel storage
    assumption_df = pd.DataFrame({"assumption_output": assumption_text.split("\n")})
    
    
    return(cph, assumption_df)

def safe_sheet_name(name):
    # Excel allows max 31 chars and forbids : \ / ? * [ ]
    cleaned = "".join(c for c in name if c not in r'[]:*?/\\')
    return cleaned[:31]



def write_COX_to_excel(dfs, duration, event, output_path, covs_incl=False):
    """
    Function that takes a dictionary of different dataframes that want to analyse
    using the cox proportional hazards model. Returns an excel sheet with different 
    sheets for the results per dataframe.

    Parameters
    ----------
    dfs : dict,
        dictionary containing the identifiers and dataframes.
    covs_incl : bool, 
        determine if want to include covariates sex, age and kps in analysis or not. 
        Default is False.
    output_path : str,
        path to file that should be saved.

    Returns
    -------
    None.

    """
    
    results = {}
    

    # Run Cox model for each dataset
    for cohort_name, (df, metric) in dfs.items():
        print(f"Running Cox model for {cohort_name}...")
        cph, assumption_df = COXPH(df, metric, duration, event, covs_incl)
        results[cohort_name] = {
            "summary": cph.summary,
            "assumptions": assumption_df
        }
        
        
    
    # Write all results to one Excel file with multiple tabs
    output_path = f"{output_path}"
    
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for cohort_name, content in results.items():

            model_sheet = safe_sheet_name(f"{cohort_name}_model")
            assum_sheet = safe_sheet_name(f"{cohort_name}_assumptions")

            content["summary"].to_excel(writer, sheet_name=model_sheet)
            content["assumptions"].to_excel(writer, sheet_name=assum_sheet, index=False)

            
    print(f"Results saved to {output_path}")    
    



  