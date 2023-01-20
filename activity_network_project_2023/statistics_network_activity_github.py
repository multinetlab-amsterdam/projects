#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:58:55 2022

Script for main analysis: regional network characteristics in relation to activity

@author: m.zimmermann
"""

#---IMPORTS---# 
import numpy as np 
import pandas as pd
from datetime import datetime

import matplotlib.pyplot as plt
import seaborn as sns

import statsmodels.api as sm
import statsmodels.formula.api as smf 
import scipy.stats as stats
from statsmodels.stats.diagnostic import het_white 
import statsmodels.stats.multitest

from scipy.stats import zscore
from scipy.stats import pearsonr


from scipy.stats import ttest_1samp
from scipy.stats import ttest_ind



#%%

#########################
# --- LINEAR MIXED MODELS TO TEST RELATIONSHIPS --- # 
#########################

#%%


def mixed_model(df, metric_activity, area, group):
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
    area : str,
        determines the are that is being investigated (whole brain, peritumoral, non-peritumoral homologue)
    group : str,
        determines the group that is being investigated. (patients, HCs)

    Returns
    -------
    None.

    """
    now = datetime.now()
    li = []
    for freq in ['alpha','delta','theta']: 
        for den in ['02','03']:
            
            #fit a linear mixed model with subjects as random intercept 
            model = smf.mixedlm(metric_activity + '_z ~ EC_z_'+ freq +'_' + den + '+ CC_z_' + freq + '_' + den, df, groups = df['sub'])
            results = model.fit()
            
            #get the summary table of the results
            results_summary = results.summary()
    
            #extract the model specs, parameters (coefficients) and p-values
            df_model_specs = results_summary.tables[0]
            df_results = results_summary.tables[1]
            df_results.reset_index(inplace = True)
            
            df_pvalues = pd.DataFrame(results.pvalues)
            df_pvalues.reset_index(inplace = True)
            df_pvalues.columns = ['tmp','p']
            
            #save everything in a dataframe
            df_all = pd.concat([df_model_specs, df_results, df_pvalues], axis = 1)
            df_all['Frequency'] = freq
            df_all['Density'] = str(int(den[1])*10) + '%'
            
            li.append(df_all)
            
    #construct a dataframe with all model results from all frequencies and densities
    df_final = pd.concat(li, axis = 0)
    
    
    #save df 
    df_final.to_csv('/path/to/goal/directory/' 
                    + now.strftime('%Y%m%d') + '_desired_folder_name/' + now.strftime('%Y%m%d') + '_' + area + '_' + group + '_' + metric_activity +'.csv')


#%%
########################
#PATIENTS ANALYSIS
########################

#Peritumoral
mixed_model(df_peritumoral, 'offset', 'peritumoral', 'patients')

#Non-peritumoral homologue
mixed_model(df_homologue, 'offset', 'homologue', 'patients')

#Rest of the brain (all non-tumoral areas)
mixed_model(df_non_tum, 'offset', 'non_tum', 'patients')

#%%

########################
#HCs ANALYSIS
########################

mixed_model(df_master_HC, 'offset', 'whole_brain', 'HC')


#%%

########################
#SUBTYPE ANALYSIS
########################

#IDH wildtype 
mixed_model(IDH_wt_rest, 'offset', 'rest', 'IDH_wt')
mixed_model(IDH_mut_codeleted_rest, 'offset', 'rest', 'IDH_codeleted')
mixed_model(IDH_mut_noncodeleted_rest, 'offset', 'rest', 'IDH_noncodeleted')


#%% 

########################
#STANDARDIZED ANALYSIS (TO OBTAIN BETA COEFFICIENTS)
########################
### Effect size 1: Standardization of coefficients (betas)### 
mixed_model(df_non_tum_z, 'offset', 'non_tum', 'patients')
mixed_model(df_master_HC_z, 'offset', 'whole_brain', 'HC')
mixed_model(df_IDH_wt_z, 'offset', 'rest', 'IDH_wt')
mixed_model(df_IDH_mut_codel_z, 'offset', 'rest', 'IDH_codeleted')
mixed_model(df_IDH_mut_noncodel_z, 'offset', 'rest', 'IDH_noncodeleted')

#%%

def mixed_model_interactions(df, metric_activity, area):
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
    area : str,
        determines the are that is being investigated (whole brain, peritumoral, non-peritumoral homologue)


    Returns
    -------
    None.

    """
    now = datetime.now()
    li = []
    for freq in ['delta', 'theta', 'alpha']: 
        for den in ['02', '03']:
            
            #mixed model with random intercept for subjects; interaction terms are added to investigate differences between patients and HCs
            model = smf.mixedlm(f'{metric_activity}_z ~ EC_z_{freq}_{den} + CC_z_{freq}_{den} + C(group) + C(group)*EC_z_{freq}_{den} + C(group)*CC_z_{freq}_{den}', df, groups = 'sub')
            results = model.fit()
            
            #get the summary table of the results
            results_summary = results.summary()
    
            #extract the model specs, parameters (coefficients) and p-values
            df_model_specs = results_summary.tables[0]
            df_results = results_summary.tables[1]
            df_results.reset_index(inplace = True)
            
            df_pvalues = pd.DataFrame(results.pvalues)
            df_pvalues.reset_index(inplace = True)
            df_pvalues.columns = ['tmp','p']
            
            #save everything in a dataframe
            df_all = pd.concat([df_model_specs, df_results, df_pvalues], axis = 1)
            df_all['Frequency'] = freq
            df_all['Density'] = str(int(den[1])*10) + '%'
            
            li.append(df_all)
            
    
    #construct a dataframe with all model results from all frequencies and densities
    df_final = pd.concat(li, axis = 0)
    
    
    #save df 
    df_final.to_csv('/path/to/goal/directory/' 
                    + now.strftime('%Y%m%d') + '_desired_folder_name/' + now.strftime('%Y%m%d') + '_' + area + '_interaction_' + metric_activity +'.csv')
    

#%% 

#Make dataframes containing both patients and HC data in preparation for
#the interaction analysis

# --- Peritumoral dataframes --- #
df_peritumoral['group'] = 'patients'
df_peritumoral_all = pd.concat([df_peritumoral, df_master_HC], axis = 0)
df_peritumoral_all.reset_index(drop = True, inplace = True)

# --- Non-Peritumoral Homologue dataframes --- #
df_homologue['group'] = 'patients'
df_homologue_all = pd.concat([df_homologue, df_master_HC], axis = 0)
df_homologue_all.reset_index(drop = True, inplace = True)

# --- All Non-peritumoral areas --- # 
df_non_tum['group'] = 'patients'
df_master_HC['group'] = 'HCs'
df_non_tum_all = pd.concat([df_non_tum, df_master_HC], axis = 0)
df_non_tum_all.reset_index(drop = True, inplace = True)


#%%
########################
#INTERACTION ANALYSIS
########################

# --- Peritumoral areas --- # 
mixed_model_interactions(df_peritumoral_all, 'offset', 'peritumor')

# --- Homologue areas --- # 
mixed_model_interactions(df_homologue_all, 'offset', 'homologue')

# --- Nontumoral areas --- # 
mixed_model_interactions(df_non_tum_all, 'offset', 'non-tumoral')

#%%
########################
#FDR CORRECTION
########################

import glob
import statsmodels.stats.multitest

#path only to non-interaction models --> this is for non-interaction models 
#files = glob.glob('/path/to/files/with/model/results/*.csv')#comment/uncomment depending on whether you want to correct the norma model of the 

#path to interaction analyses
files = glob.glob('/path/to/files/with/interaction/results/*.csv')

for file in files:
    ##print(name)
    df = pd.read_csv(file)

    df_nan = df.loc[df['p'].isna()]
    df_not_nan = df.loc[~df['p'].isna()]
    
    #df_IVs = df_not_nan.loc[df_not_nan['index'].str.startswith(('CC', 'EC'))]
    df_IVs = df_not_nan.loc[df_not_nan['index'].str.startswith(('CC', 'EC', 'C(group)[T.patients]:'))]
    
    #correct only with p-values from IVs (EC, CC)
    df_IVs['corrected_only_IVs'] = statsmodels.stats.multitest.fdrcorrection(df_IVs['p'])[1]

    df_not_nan = df_not_nan.loc[~df_not_nan['index'].str.startswith(('CC', 'EC'))]
    df_sorted = pd.concat([df_nan, df_not_nan, df_IVs]).sort_index()
    
    #correct with all pvalues from model (including intercept)
    df_nan = df_sorted.loc[df_sorted['p'].isna()]
    df_not_nan = df_sorted.loc[~df_sorted['p'].isna()]
    
    df_not_nan['corrected'] = statsmodels.stats.multitest.fdrcorrection(df_not_nan['p'])[1]
    
    df_final = pd.concat([df_nan, df_not_nan]).sort_index()
    
    df_final.to_csv(file)
    
    
#%%
########################
#WITHIN-SUBJECT CORRELATION ANALYSIS: PEARSONS CORRELATIONS
########################
# In this alternative analysis we correlate EC/CC and activity metrics per subject to understand the within-subject relationships between these measures. 
# This is an alternative to the Linear mixed models used above, to check for replicability. We use the Pearson correlation to calculate the correlation 
# per subject.We calculate the Pearson correlation, Fisher z-score the correlations (per frequency band and density) and do a t-test between groups. 

#%%

def make_corrs(df, group, metric, metric_activity):
    """
    

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

    Returns
    -------
    df_corrs : pd.DataFrame, 
        dataframe containing the pearsons correlations
        
    df_corrs_z :pd.DataFrame, 
        dataframe containing the Fisher z transformed correlations

    """
    
    #group the dataframe by subject
    df_grouped = df.groupby('sub')
    li = []
    for freq in ['delta', 'theta', 'alpha']: 
        for d in ['02', '03']: 
          
            #per subject, correlate the activity metric with the network metric and extract the correlation coefficient
            corrs = pd.DataFrame(df_grouped.apply(lambda x: pearsonr(x[f'{metric}_z_{freq}_{d}'], x[f'{metric_activity}_z'])[0]))
            corrs.columns = [f'corrs_{freq}_{d}']
            
            
            li.append(corrs)
    #put all correlations into one dataframe (all frequencies and densities) and save it
    df_corrs = pd.concat(li, axis = 1)
    df_corrs.to_csv(f'_desired_folder_name/df_{metric}_{metric_activity}_within_sub_corrs_{group}.csv')

    #Fisher z-transform the correlations and save the dataframe
    df_corrs_z = pd.DataFrame(df_corrs.apply(lambda x: np.arctanh(x)))
    df_corrs_z.to_csv(f'_desired_folder_name/df_{metric}_{metric_activity}_within_sub_corrs_z_{group}.csv')
    
    return(df_corrs,df_corrs_z)

#%%
########################
#WHOLE GROUP PATIENTS ANALYSIS: Non-tumoral areas
########################
df_corrs_patient_rest_CC, df_corrs_z_patient_rest_CC = make_corrs(df_non_tum,'patients_non_tum','CC','offset')
df_corrs_patient_rest_EC, df_corrs_z_patient_rest_EC = make_corrs(df_non_tum,'patients_non_tum', 'EC','offset')


########################
#HCs analysis
########################
df_corrs_HCs_CC, df_corrs_z_HCs_CC = make_corrs(df_master_HC,'HCs','CC','offset')
df_corrs_HCs_EC, df_corrs_z_HCs_EC = make_corrs(df_master_HC,'HCs','EC','offset')


#%%
########################
#SUBTYPE ANALYSIS
######################## 

# --- IDH-wildtype --- # 
df_corrs_wt_rest_CC, df_corrs_z_wt_rest_CC = make_corrs(IDH_wt_rest,'IDH_wt','CC','offset')
df_corrs_wt_rest_EC, df_corrs_z_wt_rest_EC = make_corrs(IDH_wt_rest,'IDH_wt', 'EC','offset')

# --- IDH-mutant, 1p19q-codeleted --- # 
df_corrs_codel_rest_CC, df_corrs_z_codel_rest_CC = make_corrs(IDH_mut_codeleted_rest,'IDH_mut_codeleted','CC','offset')
df_corrs_codel_rest_EC, df_corrs_z_codel_rest_EC = make_corrs(IDH_mut_codeleted_rest,'IDH_mut_codeleted', 'EC','offset')

# --- IDH-mutant, 1p19q-noncodeleted --- # 
df_corrs_noncodel_rest_CC, df_corrs_z_noncodel_rest_CC = make_corrs(IDH_mut_noncodeleted_rest,'IDH_mut_noncodeleted','CC','offset')
df_corrs_noncodel_rest_EC, df_corrs_z_noncodel_rest_EC = make_corrs(IDH_mut_noncodeleted_rest,'IDH_mut_noncodeleted', 'EC','offset')


#%%
########################
#GROUP LEVEL CORRELATION PREPARATION
########################

#To obtain group-level correlation coefficients we used this tool: https://www.psychometrica.de/correlation.html (Sectie 8)
#Here, we calculate the number of regions involved in the correlation per patients (for the weighting)

# --- Non-tumoral areas --- #
nr_rois_non_tum = pd.DataFrame(df_non_tum['sub'].value_counts())
nr_rois_non_tum.reset_index(inplace = True)
nr_rois_non_tum.columns = ['sub', 'nr_rois']
nr_rois_non_tum.sort_values('sub', inplace = True)

# --- IDH-wildtype --- # 
nr_rois_IDH_wt = pd.DataFrame(IDH_wt_rest['sub'].value_counts())
nr_rois_IDH_wt.reset_index(inplace = True)
nr_rois_IDH_wt.columns = ['sub', 'nr_rois']
nr_rois_IDH_wt.sort_values('sub', inplace = True)

# --- IDH-mutant, 1p19q-codeleted --- # 
nr_rois_IDH_mut_codel = pd.DataFrame(IDH_mut_codeleted_rest['sub'].value_counts())
nr_rois_IDH_mut_codel.reset_index(inplace = True)
nr_rois_IDH_mut_codel.columns = ['sub', 'nr_rois']
nr_rois_IDH_mut_codel.sort_values('sub', inplace = True)

# --- IDH-mutant, 1p19q-noncodeleted --- # 
nr_rois_IDH_mut_noncodel = pd.DataFrame(IDH_mut_noncodeleted_rest['sub'].value_counts())
nr_rois_IDH_mut_noncodel.reset_index(inplace = True)
nr_rois_IDH_mut_noncodel.columns = ['sub', 'nr_rois']
nr_rois_IDH_mut_noncodel.sort_values('sub', inplace = True)

#%%
########################
#SIGNIFICANCE OF CORRELATIONS
########################

# --- Calculate siginificance of Fisher transformed correlations against 0 (1 sample t-test)--- # 

########################
#WHOLE GROUP PATIENTS ANALYSIS: Non-tumoral areas
########################

CC_t_non_tum, CC_p_non_tum = ttest_1samp(df_corrs_z_non_tum_CC, 0)
df_CC_corr_against_0_non_tum = pd.concat([pd.DataFrame(CC_t_non_tum), pd.DataFrame(CC_p_non_tum)], axis = 1)
df_CC_corr_against_0_non_tum.columns = ['t', 'p']

EC_t_non_tum, EC_p_non_tum = ttest_1samp(df_corrs_z_non_tum_EC, 0)
df_EC_corr_against_0_non_tum = pd.concat([pd.DataFrame(EC_t_non_tum), pd.DataFrame(EC_p_non_tum)], axis = 1)
df_EC_corr_against_0_non_tum.columns = ['t', 'p']

#%%

########################
#HCs analysis
########################

CC_t_HCs, CC_p_HCs = ttest_1samp(df_corrs_z_HCs_CC, 0)
df_CC_corr_against_0_HCs = pd.concat([pd.DataFrame(CC_t_HCs), pd.DataFrame(CC_p_HCs)], axis = 1)
df_CC_corr_against_0_HCs.columns = ['t', 'p']

EC_t_HCs, EC_p_HCs = ttest_1samp(df_corrs_z_HCs_EC, 0)
df_EC_corr_against_0_HCs = pd.concat([pd.DataFrame(EC_t_HCs), pd.DataFrame(EC_p_HCs)], axis = 1)
df_EC_corr_against_0_HCs.columns = ['t', 'p']

#%%
########################
#SUBTYPE ANALYSIS
########################

# --- IDH-wildtype --- # 
CC_t_wt_rest, CC_p_wt_rest = ttest_1samp(df_corrs_z_wt_rest_CC, 0)
df_CC_corr_against_0_wt = pd.concat([pd.DataFrame(CC_t_wt_rest), pd.DataFrame(CC_p_wt_rest)], axis = 1)
df_CC_corr_against_0_wt.columns = ['t', 'p']

EC_t_wt_rest, EC_p_wt_rest = ttest_1samp(df_corrs_z_wt_rest_EC, 0)
df_EC_corr_against_0_wt = pd.concat([pd.DataFrame(EC_t_wt_rest), pd.DataFrame(EC_p_wt_rest)], axis = 1)
df_EC_corr_against_0_wt.columns = ['t', 'p']


# --- IDH-mutant, 1p19q-codeleted --- #
CC_t_codel_rest, CC_p_codel_rest = ttest_1samp(df_corrs_z_codel_rest_CC, 0)
df_CC_corr_against_0_codel = pd.concat([pd.DataFrame(CC_t_codel_rest), pd.DataFrame(CC_p_codel_rest)], axis = 1)
df_CC_corr_against_0_codel.columns = ['t', 'p']

EC_t_codel_rest, EC_p_codel_rest = ttest_1samp(df_corrs_z_codel_rest_EC, 0)
df_EC_corr_against_0_codel = pd.concat([pd.DataFrame(EC_t_codel_rest), pd.DataFrame(EC_p_codel_rest)], axis = 1)
df_EC_corr_against_0_codel.columns = ['t', 'p']

# --- IDH-mutant, 1p19q-noncodeleted --- #
CC_t_noncodel_rest, CC_p_noncodel_rest = ttest_1samp(df_corrs_z_noncodel_rest_CC, 0)
df_CC_corr_against_0_noncodel = pd.concat([pd.DataFrame(CC_t_noncodel_rest), pd.DataFrame(CC_p_noncodel_rest)], axis = 1)
df_CC_corr_against_0_noncodel.columns = ['t', 'p']

EC_t_noncodel_rest, EC_p_noncodel_rest = ttest_1samp(df_corrs_z_noncodel_rest_EC, 0)
df_EC_corr_against_0_noncodel = pd.concat([pd.DataFrame(EC_t_noncodel_rest), pd.DataFrame(EC_p_noncodel_rest)], axis = 1)
df_EC_corr_against_0_noncodel.columns = ['t', 'p']

#%% 
########################
#FDR CORRECTION
########################


########################
#WHOLE GROUP PATIENTS ANALYSIS: Non-tumoral areas
########################
df_CC_corr_against_0_non_tum['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_CC_corr_against_0_non_tum['p'])[1]
df_EC_corr_against_0_non_tum['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_EC_corr_against_0_non_tum['p'])[1]
df_CC_corr_against_0_non_tum.to_csv('desired_folder_name/df_CC_patients_against_0.csv')
df_EC_corr_against_0_non_tum.to_csv('_desired_folder_name/df_EC_patients_against_0.csv')


########################
#HCs ANALYSIS
########################
df_CC_corr_against_0_HCs['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_CC_corr_against_0_HCs['p'])[1]
df_EC_corr_against_0_HCs['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_EC_corr_against_0_HCs['p'])[1]
df_CC_corr_against_0_HCs.to_csv('desired_folder_name/df_CC_HCs_against_0.csv')
df_EC_corr_against_0_HCs.to_csv('desired_folder_name/df_EC_HCs_against_0.csv')

#%%

########################
#SUBTYPE ANALYSIS
########################

# --- IDH wildtype --- # 
df_CC_corr_against_0_wt['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_CC_corr_against_0_wt['p'])[1]
df_EC_corr_against_0_wt['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_EC_corr_against_0_wt['p'])[1]
df_CC_corr_against_0_wt.to_csv('/desired_folder_name/df_CC_wt_against_0.csv')
df_EC_corr_against_0_wt.to_csv('/desired_folder_name/df_EC_wt_against_0.csv')

# --- IDH mutant, 1p19q-codeleted --- #
df_CC_corr_against_0_codel['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_CC_corr_against_0_codel['p'])[1]
df_EC_corr_against_0_codel['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_EC_corr_against_0_codel['p'])[1]
df_CC_corr_against_0_codel.to_csv('/desired_folder_name/df_CC_codel_against_0.csv')
df_EC_corr_against_0_codel.to_csv('/desired_folder_name/df_EC_codel_against_0.csv')

# --- IDH mutant, 1p19q-noncodeleted --- #
df_CC_corr_against_0_noncodel['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_CC_corr_against_0_noncodel['p'])[1]
df_EC_corr_against_0_noncodel['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_EC_corr_against_0_noncodel['p'])[1]
df_CC_corr_against_0_noncodel.to_csv('/desired_folder_name/df_CC_noncodel_against_0.csv')
df_EC_corr_against_0_noncodel.to_csv('/desired_folder_name/df_EC_noncodel_against_0.csv')


#%% 
########################
#DIFFERENCE BETWEEN FISHERS Z TRANSFORMED CORRELATIONS OF PATIENTS AND HCS
########################

t_CC_comp_HCs_patients_rest, p_CC_comp_HCs_patients_rest = ttest_ind(df_corrs_z_non_tum_CC, df_corrs_z_HCs_CC)
df_CC_corr_comp = pd.concat([pd.DataFrame(t_CC_comp_HCs_patients_rest), pd.DataFrame(p_CC_comp_HCs_patients_rest)], axis = 1)
df_CC_corr_comp.columns = ['t', 'p']

t_EC_comp_HCs_patients_rest, p_EC_comp_HCs_patients_rest = ttest_ind(df_corrs_z_non_tum_EC, df_corrs_z_HCs_EC)
df_EC_corr_comp = pd.concat([pd.DataFrame(t_EC_comp_HCs_patients_rest), pd.DataFrame(p_EC_comp_HCs_patients_rest)], axis = 1)
df_EC_corr_comp.columns = ['t', 'p']

#%%
########################
# FDR CORRECTION
########################
df_CC_corr_comp['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_CC_corr_comp['p'])[1]
df_EC_corr_comp['FDR_corrected_p'] = statsmodels.stats.multitest.fdrcorrection(df_EC_corr_comp['p'])[1]
df_CC_corr_comp.to_csv('/desired_folder_name/df_CC_comp.csv')
df_EC_corr_comp.to_csv('/desired_folder_name/df_EC_comp.csv')

#%%
########################
#RELATION BETWEEN PERITUMORAL ACTIVITY AND PEARSONS CORRELATIONS IN NON-TUMORAL REGIONS
########################

#import correlation dataframes
df_CC_corrs = pd.read_csv('/desired_folder_name/df_CC_offset_within_sub_corrs_patients_non_tum.csv')
df_EC_corrs = pd.read_csv('/desired_folder_name/df_EC_offset_within_sub_corrs_patients_non_tum.csv')


#%%

def make_corrs_offset_corr(df_corrs, df_offset,  metric, plot = True):
    """
    

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
    
    #average offset_z over all regions in peritumoral area 
    df_offset_avg = pd.DataFrame(df_offset.groupby('sub')['offset_z'].mean())#comment/uncomment depending on which are you want to average the offset over
   
    
    #average only over 3 highest offset_z values
    df_offset_3_largest = pd.DataFrame(df_offset.groupby('sub')['offset_z'].nlargest(3))#comment/uncomment depending on which are you want to average the offset over
    df_offset_avg = pd.DataFrame(df_offset_3_largest.groupby('sub')['offset_z'].mean())#comment/uncomment depending on which are you want to average the offset over
    
    
    corrs_dict = {}
    for freq in ['alpha']: 
        for d in ['02', '03']: 
          
            #correlate peritumoral offset with the correlation between network metric and offset
            corr, p = pearsonr(df_offset_avg['offset_z'], df_corrs[f'corrs_{freq}_{d}'])
            corrs_dict[f'corrs_{freq}_{d}'] = corr, p #store correlation in dictionary
            
            if plot == True:
            
                plt.scatter(df_corrs[f'corrs_{freq}_{d}'], df_offset_avg['offset_z'])
                plt.xlabel(f'corrs_{freq}_{d}_{metric}')
                plt.ylabel('offset_z')
                plt.show()
    
    #store all correlations in one dataframe
    df_corrs = pd.DataFrame.from_dict(corrs_dict).T
    df_corrs.columns = ['corr', 'p']
    df_corrs.to_csv(f'/desired_folder_name/df_{metric}_offset_corr.csv')


        
    return(df_corrs)

#%%
########################
#ALL PERITUMORAL AREAS
########################

########################
#WHOLE GROUP PATIENTS ANALYSIS: Non-tumoral areas
########################

df_EC_corrs_offset = make_corrs_offset_corr(df_EC_corrs, df_peritumoral_all, 'EC')
df_CC_corrs_offset = make_corrs_offset_corr(df_CC_corrs, df_peritumoral_all, 'CC')


#%%
########################
#SUBTYPE ANALYSIS
########################

# --- IDH-wildtype --- #
df_EC_corrs_offset_wt = make_corrs_offset_corr(df_corrs_wt_rest_EC, IDH_wt_peritumoral_all, 'EC_wt')
df_CC_corrs_offset_wt = make_corrs_offset_corr(df_corrs_wt_rest_CC, IDH_wt_peritumoral_all, 'CC_wt')

# --- IDH-mutant, 1p19q-codeleted --- #
df_EC_corrs_offset_mut_codel = make_corrs_offset_corr(df_corrs_codel_rest_EC, IDH_mut_codeleted_peritumoral_all, 'EC_mut_codel')
df_CC_corrs_offset_mut_codel = make_corrs_offset_corr(df_corrs_codel_rest_CC, IDH_mut_codeleted_peritumoral_all, 'CC_mut_codel')

# --- IDH mutant, 1p19q-noncoldeleted --- #
df_EC_corrs_offset_mut_noncodel = make_corrs_offset_corr(df_corrs_noncodel_rest_EC, IDH_mut_noncodeleted_peritumoral_all, 'EC_mut_noncodel')
df_CC_corrs_offset_mut_noncodel = make_corrs_offset_corr(df_corrs_noncodel_rest_CC, IDH_mut_noncodeleted_peritumoral_all, 'CC_mut_noncodel')

#%%

########################
#THREE HIGHEST ACTIVITY PERITUMORAL AREAS
########################

########################
#WHOLE GROUP PATIENTS ANALYSIS: Non-tumoral areas
########################
df_EC_corrs_offset = make_corrs_offset_corr(df_EC_corrs, df_peritumoral_all, 'EC_3')
df_CC_corrs_offset = make_corrs_offset_corr(df_CC_corrs, df_peritumoral_all, 'CC_3')


########################
#SUBTYPE ANALYSIS
########################

# --- IDH-wildtype --- # 
df_EC_corrs_offset_wt = make_corrs_offset_corr(df_corrs_wt_rest_EC, IDH_wt_peritumoral_all, 'EC_wt_3')
df_CC_corrs_offset_wt = make_corrs_offset_corr(df_corrs_wt_rest_CC, IDH_wt_peritumoral_all, 'CC_wt_3')

# --- IDH-mutant, 1p19q-codeleted --- #
df_EC_corrs_offset_mut_codel = make_corrs_offset_corr(df_corrs_codel_rest_EC, IDH_mut_codeleted_peritumoral_all, 'EC_mut_codel_3')
df_CC_corrs_offset_mut_codel = make_corrs_offset_corr(df_corrs_codel_rest_CC, IDH_mut_codeleted_peritumoral_all, 'CC_mut_codel_3')

# --- IDH mutant, 1p19q-noncoldeleted --- #
df_EC_corrs_offset_mut_noncodel = make_corrs_offset_corr(df_corrs_noncodel_rest_EC, IDH_mut_noncodeleted_peritumoral_all, 'EC_mut_noncodel_3')
df_CC_corrs_offset_mut_noncodel = make_corrs_offset_corr(df_corrs_noncodel_rest_CC, IDH_mut_noncodeleted_peritumoral_all, 'CC_mut_noncodel_3')


#%%

########################
#SPINTEST PREPARATIONS
########################

### Load in dataframes containing raw (non-standardized) data ###
df_patients_raw = pd.read_pickle('/desired_folder_name/df_master_pat_no_distance_no_edema_z_and_raw.pkl')
df_patients_raw = df_patients_raw.loc[:, ~df_patients_raw.columns.duplicated()] #remove duplicate columns
df_patients_raw.reset_index(inplace = True)


########################
#PATIENTS RAW DATAFRAMES
########################

#df only containing subjects with tumor mask (so tumors with overlap to brain regions >12%)
df_patients_raw = df_patients_raw.loc[:, ~df_patients_raw.columns.duplicated()]
subs_no_overlap = df_patients_raw.groupby('sub').sum()
subs_no_overlap.reset_index(inplace = True)
li_subs_no_overlap = subs_no_overlap.loc[subs_no_overlap['peritumor'] == 0, 'sub']
df_master_only_overlaps = df_patients_raw.loc[~df_patients_raw['sub'].isin(li_subs_no_overlap)]

#subjects that have tumor that reaches into both hemispheres (N = 11) 
subs_hom = np.unique(df_master_only_overlaps.loc[(df_master_only_overlaps['non_peritumoral_hom'] == 1)&(df_master_only_overlaps['peritumor'] == 1), 'sub'])

#df containining subjects that do not have tumor that reaches to both hemispheres
df_master_uni = df_master_only_overlaps.loc[~df_master_only_overlaps['sub'].isin(subs_hom)]


#df containing only rest of the brain (non-tumoral data) filtered on subjects that are included in the study
patients_included = sorted(df_non_tum['sub'].unique())
df_non_tum_raw = df_master_only_overlaps.loc[df_master_only_overlaps['perc'] == 0]
df_non_tum_raw = df_non_tum_raw.loc[df_non_tum_raw['sub'].isin(patients_included)] #only include patients that are included in the main analysis

#%%
########################
#HCs RAW DATAFRAMES
########################
df_HC_raw = pd.read_pickle('/desired_folder_name/df_HC_master_z_and_raw.pkl')
df_HC_raw = df_HC_raw.loc[:, ~df_HC_raw.columns.duplicated()]

#%%
########################
#AVERGAE PER REGION 
########################

#obtain one average per brain region over all subjects
grouped_by_roi = df_non_tum_raw.groupby('roi').mean()#comment/uncomment depending on which group you want to average 
#grouped_by_roi = df_HC_raw.groupby('roi').mean()

group = 'patients_rest_non_tum'
for freq in ['delta', 'theta', 'lower_alpha']: 
    for d in ['2', '3']: 
        #extract the network metric columnsper frequency and density and store in dataframe
        df_CC = grouped_by_roi[['Offset', f'CC_{freq}_{d}']]
        df_CC.to_csv(f'/desired_folder_name/df_CC_{freq}_{d}_{group}.csv', header = False, index = False)#change name of file depending on which group you averaged for
        
        #extract the network metric columnsper frequency and density and store in dataframe
        df_EC = grouped_by_roi[['Offset', f'EC_{freq}_{d}']]
        df_EC.to_csv(f'/desired_folder_name/df_EC_{freq}_{d}_{group}.csv', header = False, index = False)#change name of file depending on which group you averaged for
             
