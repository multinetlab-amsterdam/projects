#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create master dataframe

Scripts to create a dataframe with all offset and network metrics for different densities and frequencies, with tumor overlaps and raw values in one #

"""
__author__ = "Mona Lilo Margarethe Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "2022/02/01"   
__status__ = "Finished"

####################
# Review History   #
####################

# Reviewed by Eduarda Centeno 20230202

#%%
####################
# Libraries        #
####################

# Standard imports  ###
import glob
from functools import reduce

# Third party imports ###
import numpy as np 
import pandas as pd #version 1.1.5
import re


#%%
########################
#FUNCTION STANDARDIZATION ACTIVITY METRIC
########################


def standardize_activity_metrics(input_dir, input_HC_dir, output_dir, HC_subs, HC = False): 
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
    
    
        
    #get list of all files containing offset,broadband power and slope
    HC_files = glob.glob(input_HC_dir + '/*_bbp*' )
    sub_files = glob.glob(input_dir + '/*_bbp*' )
    
    li_HC = []
    #iterate through all files (subjects)
    for file in HC_files: 
        
        #load in file
        HC_df = pd.read_csv(file, index_col = 0).T.iloc[0:210]
        
        #search subject_id and add as column to dataframe
        HC_df['sub'] = re.search('([A-z]{4}_\d{3,4})', file)[0] 
        
        #make a column for the rois
        HC_df['roi'] = np.arange(1,211)
         
        #append each dataframe to list
        li_HC.append(HC_df)
        
    
    #make dataframe 
    HC_df = pd.concat(li_HC) 
    
    #only include matched HCs (on sex and age)
    HC_df = HC_df[HC_df['sub'].isin(HC_subs['sub'])]
    
    #determine if want to standardize HCs (on themselves) or patients & load in dataframe
    if HC == True: 
        pat_df = HC_df 
        
    else: 
        li_pat = []
        for file in sub_files: 
            
            #load in dataframe containing activity metrics for every subject
            pat_df = pd.read_csv(file,index_col = 0).T.iloc[0:210] #dataframe is turned to have activity metrics as columns instead of rows
            
            #search subject_id and add as column to dataframe
            pat_df['sub'] = re.search('(sub-\d{4})', file)[0] 
            
            #make a column for the rois
            pat_df['roi'] = np.arange(1,211)
            
            #store every patients dataframe
            li_pat.append(pat_df)
         
        #construct big dataframe containing data from all patients    
        pat_df = pd.concat(li_pat)  
    
    
    # --- Step 1.) Get regional mean and std of HC power, offset and slope --- #
    li_mean = []
    for measure in ['Broadband Power', 'Offset','Slope']: 
        
        #obtain regional mean and std of particular activity metric
        HC_mean = HC_df.groupby('roi')[measure].mean()
        HC_std = HC_df.groupby('roi')[measure].std()  
        
        #store means and std in dataframe 
        df_mean = HC_mean.to_frame()
        df_mean[measure +'_std_HC'] = HC_std 
        df_mean['rois'] = HC_std.index
        df_mean.rename(columns = {measure:measure + '_mean_HC'}, inplace = True)
        df_mean.reset_index(drop = True,inplace = True)
        
        #store dataframes of every subject in list
        li_mean.append(df_mean)
    
    #make a dataframe with all activity metrics for all subjects
    df_mean = pd.concat(li_mean, axis = 1)
    
    #drop duplicate columns (here: rois)
    df_mean = df_mean.loc[:, ~df_mean.columns.duplicated()]
    df_mean.set_index('rois', inplace = True)
    
    
    # --- Step 2.) add means and std repeatedly to rows (i.e. rois) of subs dataframe to allow for easy standardization later--- #
    
    for measure in ['Broadband Power', 'Offset','Slope']:
        pat_df[measure + '_mean_HC'] = pat_df['roi'].apply(lambda x: df_mean[measure + '_mean_HC'].loc[x])
        pat_df[measure + '_std_HC'] = pat_df['roi'].apply(lambda x: df_mean[measure + '_std_HC'].loc[x])
    
    
    # --- Step 3.) Standardize BB, offset and slope and add as new columns to dataframe --- # 
    
    for name, measure in {'BB_welch':'Broadband Power', 'offset':'Offset','slope':'Slope'}.items(): 
        pat_df[name + '_z'] = (pat_df[measure]- pat_df[measure + '_mean_HC'])/pat_df[measure + '_std_HC']
        

    # --- Step 4.) Save df ---# 
    if HC == True: 
        pat_df.to_csv(output_dir + '/name_of_file_HC.csv')
        

    else: 
        pat_df.to_csv(output_dir + '/name_of_file.csv')
        

    
    
            
#%%         
########################
#STANDARDIZE ACTIVITY METRICS
########################

input_dir_bbos = '/folder/with/unstandardized/activity/measures/patients/'
input_HC_dir_bbos = '/folder/with/unstandardized/activity/measures/HCs/'

HC_subs_li = pd.read_csv('/path/to/HC/subject/list.csv', index_col = 0)
HC_subs_li.columns = ['sub']

output_dir = '/folder/where/standardized/dataframes/shouldbe/stored/'


standardize_activity_metrics(input_dir_bbos, input_HC_dir_bbos, output_dir, HC_subs_li, HC = False) #standardize patients metrics
standardize_activity_metrics(input_HC_dir_bbos, input_HC_dir_bbos, output_dir,HC_subs_li, HC = True)#standardize HCs metrics on themselves

#%% 
########################
#CONSTRUCT ONE MAIN DATAFRAME CONTAINING ALL METRICS
########################
 
#       - z_scored BB welch 
#       - z_scored offset 
#       - z_scored slope 
#       - z_scoredd CC & EC  
#       - tumor overlap --> peritumoral/non-peritumoral & regions contralateral homologue of peritumoral area
#%%


#%%
########################
#FUNCTION TO CONSTRUCT DATAFRAME CONTAINING ALL ACTIVITY AND NETWORK METRICS
########################

def concat_power_network(df_BB_welch, dict_net_df, HC = False):
    """
    Function to concatenate z_scored power data and network markers. 

    Parameters
    ----------
    df_BB_fft : pd.DataFrame, 
        df containing BB power caluclated with the Fast Fourier Transform
        
    df_BB_welch : pd.DataFrame,
        df containing BB power caluclated with the welch method
    dict_net_df : TYPE
        DESCRIPTION.
    HC : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """

    # --- Step 1.) Prepare activity dataframe --- #
    
    #extract columns needed from welch BB dataframe
    df_activity = df_BB_welch[['BB_welch_z',  'offset_z', 'slope_z', 'sub', 'roi']].sort_values(['sub', 'roi'])
    df_activity.reset_index(inplace = True, drop = True)
   
    #remove sub-0002 & sub-9028 (exclusion)
    df_activity = df_activity[~df_activity['sub'].isin(['sub-0002', 'sub-9028'])]
    df_activity.reset_index(inplace = True, drop = True )
    
    
    # --- Step 2.) Prepare network dataframe --- #
    
    li = []
    
    #loop through the network dataframes
    for key, df in dict_net_df.items(): 
        df = df[['CC_z', 'EC_z', 'sub', 'rois']].sort_values(['sub', 'rois'])
        df.reset_index(inplace = True, drop = True )
        
        #rename columns per frequency and density
        df.columns = ['CC_z' + key, 'EC_z' + key, 'sub', 'rois']
        
        if HC == False: 
            #remove sub-0002 & sub-9028 (exclusion) from the patient dataframes
            df = df[~df['sub'].isin(['sub-0002', 'sub-9028'])
            df.reset_index(inplace = True, drop = True )
            
            #store dataframe without excluded subjects
            li.append(df)
            
        else: 
            li.append(df)
    
    #construct one dataframe containing all network metrics
    df_network =  pd.concat(li, axis =1) 
    
    #add network metric dataframe to activity dataframe
    df_master = pd.concat([df_activity, df_network], axis = 1)
    
    #drop duplicate columns rois so that only have one column roi 
    df_master = df_master.loc[:, ~df_master.columns.duplicated()]
    df_master.drop(['rois'], axis = 1)
    
    return(df_master)

#%% 
########################
#LOAD IN NEEDED DATAFRAMES
########################

# --- BB, offset, slope --- #
df_BB_welch = pd.read_csv('/folder/to/activity_dataframe/patients.csv')
df_BB_welch_HC = pd.read_csv('/folder/to/activity_dataframe/HCs.csv')

# --- CC & EC --- #
#load in all dataframes needed for network metrics and put them in dictionary 

dict_dfs = {'_theta_01': df_theta_01, '_theta_02': df_theta_02, '_theta_03': df_theta_03, 
                '_delta_01': df_delta_01, '_delta_02': df_delta_02, '_delta_03': df_delta_03, 
                '_alpha_01': df_alpha_01, '_alpha_02': df_alpha_02, '_alpha_03': df_alpha_03}


dict_dfs_HC = {'_theta_01': df_theta_01_HC, '_theta_02': df_theta_02_HC, '_theta_03': df_theta_03_HC, 
                '_delta_01': df_delta_01_HC, '_delta_02': df_delta_02_HC, '_delta_03': df_delta_03_HC, 
                '_alpha_01': df_alpha_01_HC, '_alpha_02': df_alpha_02_HC, '_alpha_03': df_alpha_03_HC}

#%% 

########################
#CONSTRUCT MASTER DATAFRAME
########################

df_master_pat = concat_power_network(df_BB_welch, dict_dfs, HC = False)
df_master_HC = concat_power_network(df_BB_welch, dict_dfs_HC, HC = True)


#%%
########################
#TUMOR OVERLAPS TO DEFINE CONTRALATERAL HOMOLOGUE
########################


#load in the dataframe containing percentage overlap of tumor with the BNA regions
df_overlaps = pd.read_csv('/path/to/overlap/file.csv', index_col = 0)

# --- Step 1: Define the contralateral homologues --- #
    
#Get the homologue contralateral non-peritumoral rois for even (right) and uneven (left) peritumoral rois:
#   - if the peritumoral region is even (i.e. right hemisphere): get the roi from row before (i.e. uneven, left hemisphere)
#   - if the peritumroal region is uneven (i.e. left hemisphere): get the next roi (i.e. even, right hemisphere)

df_overlaps['partner_roi_right'] = df_overlaps['roi'].shift(-1).where((df_overlaps['peritumor'] == 1) & (df_overlaps['roi'] % 2 != 0))
df_overlaps['partner_roi_left'] = df_overlaps['roi'].shift(1).where((df_overlaps['peritumor'] == 1) & (df_overlaps['roi'] % 2 == 0))

# --- Step 2: Indicate which area the ROI belongs to (contralateral homologue?) --- #

#Group the df by subjects and get the list of (uneven (left) and even (right)) non-peritumoral rois per subject for both. Then match those to the index of the ROIs and set a 1, 
#to indicate that that ROI is part of the non-peritumoral, contralateral area.First this is stored in lists per subject (in the groupby object). 
#Then it is exploded  to put it back into the dataframe.

#Left hemisphere 
grouped_left = pd.DataFrame(df_overlaps.groupby('sub').apply(lambda x: np.where(x['roi'].isin(x['partner_roi_left']),1,0)))
grouped_left_exploded = grouped_left[0].explode().reset_index(drop = False)

df_overlaps.reset_index(drop = True, inplace = True)
df_overlaps['non_peritumoral_left'] = grouped_left_exploded[0]

#Right hemisphere 
grouped_right = pd.DataFrame(df_overlaps.groupby('sub').apply(lambda x: np.where(x['roi'].isin(x['partner_roi_right']),1,0)))
grouped_right_exploded = grouped_right[0].explode().reset_index(drop = False)

df_overlaps['non_peritumoral_right'] = grouped_right_exploded[0]

#Merge the left and right columns to one non-peritumoral area column
df_overlaps['non_peritumoral_hom'] = df_overlaps['non_peritumoral_left'] + df_overlaps['non_peritumoral_right']

#Save df_overlaps for future reference 
df_overlaps.to_csv('/path/to/store/new/overlaps/df.csv')

#%%
# --- Step 3: Add all overlaps to master_df --- #

#load in master dataframe of patients
df_master_pat = pd.read_pickle('/path/to/master/df.pkl')

#sort df_overlaps based on sub and roi 
df_overlaps = df_overlaps.sort_values(['sub', 'roi'])
df_overlaps.reset_index(drop = True, inplace = True)
df_master_pat.reset_index(drop = True, inplace = True)

#concatenate df_overlaps with master df of patients 
df_master_pat_final = pd.concat([df_master_pat, df_overlaps], axis = 1)

#save master df of patients 
df_master_pat_final.to_pickle('/path/to/store/new/master/df.csv')


#%%
########################
#PREPARATION SPIN-TEST: CREATE RAW VALUES DATAFRAME
########################

#Contains: 
# - BB_welch, offset
# - CC and EC 


#%%
########################
#PATIENTS
########################

#load in master dataframe
df_master = pd.read_csv('/path/to/master/df.csv')

#get list of all included subjects
df_pat_subs = np.unique(df_master['sub'])

#define path to folder where all data was stored (zscored and raw, in sperate files for the different frequencies and densities)
input_dir_networks = '/path/to/raw/data/networks/'
input_dir_offset = '/path/to/raw/data/activity/'
#%%

# --- Step 1.) Network metrics --- #

for freq in ['delta', 'theta', 'lower_alpha']:
    
    #access frequency specific folder and get all files (for the different densities)
    df_path = f'{input_dir_networks}/{freq}/*.csv'
    df_list = glob.glob(df_path)
    
    
    for file in df_list: 
        print(file)
        
        #load in frequency and density specific dataframe
        df = pd.read_csv(file)
        
        #search which density we currently look at
        d = re.search('(_\d.\d)', file)[1][-1]
       
        
        #filter for patients that want to include 
        df = df.loc[df['sub'].isin(df_pat_subs)]
        df.sort_values(['sub', 'rois'], inplace = True)
        df.reset_index(drop = True,inplace = True)
        
        #make a new column in master_df containing raw network metrics
        df_master[f'CC_{freq}_{d}'] = df['CC']
        df_master[f'EC_{freq}_{d}'] = df['EC']
        
        
#%%

# --- Step 2.) Activity metrics --- #

#get activity files for every subject
sub_files = glob.glob(input_dir_offset + '/*_bbp*' )

li_pat = []
for file in sub_files: 
    
    #read in the activity file as dataframe
    pat_df = pd.read_csv(file, index_col = 0).T.iloc[0:210]
    
    #determine which subject the file belongs to and add to dataframe
    pat_df['sub'] = file[87:95] 
    
    #add rois to dataframe
    pat_df['roi'] = np.arange(1,211)
    
    #store dataframe
    li_pat.append(pat_df)
    
#construct a dataframe containing all raw activity values 
patients_df = pd.concat(li_pat) 

#only include subjects that want to include 
patients_df = patients_df[patients_df['sub'].isin(df_pat_subs)]

patients_df.sort_values(['sub', 'roi'], inplace = True)
patients_df.reset_index(drop = True,inplace = True)

#%%
# --- Step 3.) Make one dataframe of raw activity and network metrics  --- #      
        
df_master_raw = pd.concat([df_master, patients_df], axis = 1)
df_master_raw.to_pickle('/path/to/store/raw/and/zscore/df.csv')       
        
#%%

########################
#HCs
########################

#load in master dataframe for HCs
df_HC = pd.read_csv('/path/to/master/df.csv')

#get list of all included HCs subjects
df_HC_subs = np.unique(df_HC['sub'])

#define path to folder where all data was stored (zscored and raw, in sperate files for the different frequencies and densities)
input_dir_networks = '/path/to/raw/data/networks/HCs'
input_dir_offset = '/path/to/raw/data/activity/HCs'

#%%

# --- Step 1.) Network metrics --- #

for freq in ['delta', 'theta', 'lower_alpha']:
    
    #access frequency specific folder and get all files (for the different densities)
    df_path = f'{input_dir_networks}/{freq}/*.csv'
    df_list = glob.glob(df_path)
    
    for file in df_list: 
        
        #load in frequency and density specific dataframe
        df = pd.read_csv(file)
        d = re.search('(_\d.\d)', file)[1][-1]
        
        
        #filter for patients that want to include 
        df = df.loc[df['sub'].isin(df_HC_subs)]
        df.sort_values(['sub', 'rois'], inplace = True)
        df.reset_index(drop = True,inplace = True)
        
        #make a new column in master_df containing raw network metrics
        df_HC[f'CC_{freq}_{d}'] = df['CC']
        df_HC[f'EC_{freq}_{d}'] = df['EC']
        
    
        df_HC.sort_values(['sub', 'rois'], inplace = True)
        df_HC.reset_index(drop = True,inplace = True)
 
#%%

# --- Step 2.) Activity metrics --- #

#get activity files for every subject
HC_files = glob.glob(input_dir_offset + '/*_bbp*' )

li_HC = []
for file in HC_files: 
    
    #read in the activity file as dataframe
    HC_df = pd.read_csv(file, index_col = 0).T.iloc[0:210]
    
    #determine which subject the file belongs to and add to dataframe
    HC_df['sub'] = re.search('([A-z]{4}_\d{3,4})', file)[0] 
    
    #add rois to dataframe
    HC_df['roi'] = np.arange(1,211)
     
    #store dataframe
    li_HC.append(HC_df)
    

#construct a dataframe containing all raw activity values 
HC_df = pd.concat(li_HC)  

#only include subjects that want to include
HC_df = HC_df[HC_df['sub'].isin(df_HC_subs)]

HC_df.sort_values(['sub', 'roi'], inplace = True)
HC_df.reset_index(drop = True,inplace = True)

#%%

# --- Step 3.) Make one dataframe of raw activity and network metrics  --- #

df_HC_raw = pd.concat([df_HC, HC_df], axis = 1)
df_HC.to_pickle('/path/to/store/raw/and/zscore/df.csv')       
        
    
