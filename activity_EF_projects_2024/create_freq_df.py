#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 08:10:16 2024

In this script I prepare the dataframes for the post-hoc analysis to compare 
baseline and follow-up relative power in the different frequency bands.  

@author: mlzimmermann

"""

__author__ = "Mona Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl" 
__date__ = "23/06/24"  
__status__ = "Final" 


####################
# Review History   #
####################

# Reviewed by 



## Script for preparation of relative frequency bands analysis ##

####################
# LIBRARIES       #
####################

# Standard imports ###
import pandas as pd
import numpy as np 


#%%

def standardize(row, mean, std):
    return (row-mean)/std

#%%

def make_rel_power_df(df_pat, df_HC, df_overlaps, sub_ids, freq, MM):
    
    ### Standardize pat rel power based on HCs ###
    
    #calculate the mean and std of HCs
    mean_HC = df_HC.mean(axis = 1)
    std_HC = df_HC.std(axis = 1)
    
    #standardize the patient df based on the regional mean and std of HCs
    standardized_df = df_pat.apply(lambda row: standardize(row, mean_HC[row.name], std_HC[row.name]), axis=1)
    
    #add subject ids to dataframe
    standardized_df.columns = sub_ids
    
    #reshape dataframe to bring into same format as the masks dataframe
    reshaped_df = pd.melt(standardized_df, var_name='sub', value_name=f'{freq}_z')
    
    #add rois repeatedly to dataframe
    rois = np.tile(np.arange(1, 211), int(len(reshaped_df) / 210))
    
    # Add the sequence as a new column 'rois'
    reshaped_df['rois'] = rois
    
    
    
    ### Get peritumoral mean power for every pat ###
    
    # extract the roi, peritumor and sub column from the df_peri_overlaps column
    df_overlaps_filt = df_overlaps[['roi', 'sub', 'peritumor']]
    
    
    #concatenate the dataframes
    #reset indices just to be sure
    reshaped_df.reset_index(inplace = True)
    df_overlaps_filt.reset_index(inplace = True)
    
    #concatenate dataframes 
    df_overlaps_final = pd.concat([reshaped_df, df_overlaps_filt], axis = 1)
    df_overlaps_final_dropped = df_overlaps_final.loc[:, ~df_overlaps_final.columns.duplicated()]
    
    # filter rows where 'peritumor' column is 1
    peritumor_df = df_overlaps_final_dropped[df_overlaps_final_dropped['peritumor'] == 1]
    print(peritumor_df)
   
    
    # Calculate the mean of 'value' column for each 'sub'
    df_peri_mean = peritumor_df.groupby('sub')[f'{freq}_z'].mean()
        
    
    return(df_peri_mean)

#%%

sub_ids = ['sub-xxx', 'sub-xxx', ... ]
MM = 'T4'

#DEFINE INPUTS
df_overlaps = pd.read_csv('/path/to/peri_overlaps.csv')
df_overlaps = df_overlaps.loc[df_overlaps['sub'].isin(sub_ids)] #only include the subjects that want to include for the timepoint
print(df_overlaps.shape) 
freq_bands = ['alpha1', 'alpha2', 'beta', 'delta', 'gamma', 'theta']


#%%
all_dfs = []
for freq in freq_bands:
    print(freq)
    
    #load in dataframes
    df_pat = pd.read_csv(f'/path/to/20240515_avg_rel_power_{freq}_{MM}.csv', header= None)
    
    #special case T4 -->dataframe is stored in different orientation when run or only one subject (1x210 instead of 210x1)
    #df_pat = df_pat.T
    
    
    df_HC = pd.read_csv(f'/path/to/20240516_avg_rel_power_{freq}_HC.csv', header= None)
    print(df_pat.shape)

    #standardize and calculate peritumoral value 
    df_peri_mean = make_rel_power_df(df_pat, df_HC, df_overlaps, sub_ids, freq, MM)
    print(df_peri_mean.shape)

    #store in list
    all_dfs.append(df_peri_mean)
    
### Merge all dataframes ### 

# Initial merge with the first DataFrame
df_all = all_dfs[0]

# Iterate over the remaining DataFrames and merge on the common column 'sub'
for df in all_dfs[1:]:
    df_all = pd.merge(df_all, df, on='sub', how='inner')


#%%

#Load in all standardized rel power peritumoral dataframes for baseline and the FUs

### Baseline
df_peri_baseline = pd.read_csv('/path/to/20240517_std_rel_power_baseline.csv')

### FUs
df_peri_T2 = pd.read_csv('/path/to/20240517_std_rel_power_T2.csv')
df_peri_T3 = pd.read_csv('/path/to/20240517_std_rel_power_T3.csv')
df_peri_T4 = pd.read_csv('/path/to/20240517_std_rel_power_T4.csv')

#concatenate to get one FU dataframe
df_peri_FU = pd.concat([df_peri_T2, df_peri_T3, df_peri_T4], axis = 0)

#%%
### Prepare df for Plotting ### 

df_peri_baseline['MM'] = 'Baseline'
df_peri_FU['MM'] = 'FU'

df_both_timepoints = pd.concat([df_peri_baseline, df_peri_FU], axis = 0)
    
