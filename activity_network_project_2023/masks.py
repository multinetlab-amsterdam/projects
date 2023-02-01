#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Define the peritumoral area


Scripts to calculate the percentage overlap of the tumor with brain regions of the BNA to determine the peritumoral area.


"""

__author__ = "Mona Lilo Margarethe Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "2021/07/05"   
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
import os 
from glob import glob
from shutil import copy 

# Third party imports ###
import numpy as np 
import pandas as pd #version 1.1.5
import matplotlib.pyplot as plt


#%%

########################
#COPY OVERLAP AND VOLUMES FILES TO DIRECTORY 
########################

# --- Step 1.) Preparations --- #

#get a list of all subjects included in study 
df_main_patients = pd.read_pickle('/path/to/main/df.csv')
all_subs = np.unique(df_main_patients['sub']) 

#specify directory where all overlap and volumes files should be stored
goal_dir_overlap = '/desired/output/folder/overlaps/'
goal_dir_volumes = '/desired/output/folder/volumes/'

#%%

# --- STEP 2.) Copy overlap and volumes txt files to goal directory --- #

subs_no_masks = []
for sub in all_subs:
    print(sub)
    
    #get overlap file per subject
    path_overlap = glob('/path/to/folder/Scans/' + sub + '/T1/anat/' + sub + '_sienax/bna_gm_TumorOverlapVolume.txt')   
    print('Found files overlap:')
    print(len(path_overlap))
    
    #get ovolumes file per subject
    path_volumes = glob('/path/to/folder/Scans/' + sub + '/T1/anat/' + sub + '_sienax/bna_gm_RoiVolumes.txt') 
    print('Found files volume:')
    print(path_volumes)
    
    if len(path_overlap) > 0:
        print('Got here')
        
        #create directories for every subject at goal location
        goal_dir_o = os.path.join(goal_dir_overlap, sub)
        os.makedirs(os.path.join(goal_dir_overlap, sub), exist_ok = True)
        
        goal_dir_v = os.path.join(goal_dir_volumes, sub)
        os.makedirs(os.path.join(goal_dir_volumes, sub), exist_ok = True)
        
        #copy overlap and volume files to goal directory
        copy(path_overlap[0], goal_dir_o)
        copy(path_volumes[0], goal_dir_v)  
    
        
   
    else:
        
        #all subjects that do not have overlap/volume files will be added to the list as check
        subs_no_masks.append(sub)
    

#%%
# --- STEP 3.) Determine which regions fall into the peritumoral area --- #

#get overlaps per subject: which areas does the tumor cover? 
df_li = []
for sub in all_subs: 
    print(sub)
    
    
    
    #get overlaps and volumes (first cols: number of voxels, second cols: cubic volume) of tumor 
    file_overlap = glob('/path/to/overlaps/' + sub + '/*.txt')[0]
    df_overlap = pd.read_csv(file_overlap, sep = ' ', header = None)
    df_overlap.columns = ['voxels_overlap', 'volumes_overlap','Nan']
    df_overlap = df_overlap.iloc[0:210]
    
    #get volumes of all regions
    file_volumes = glob('/path/to/volumes/of_regions/' + sub + '/*.txt')[0]
    df_volumes = pd.read_csv(file_volumes, sep = ' ', header = None)
    df_volumes.columns = ['voxels_all', 'volumes_all','Nan']
    df_volumes = df_volumes.iloc[0:210]
    
    #concatenate the two dataframes to be able to compute percentage overlapof tumor with particular region
    df_full = pd.concat([df_overlap['voxels_overlap'], df_volumes['voxels_all']], axis = 1)
    
    
    #which areas show 12% overlap with tumor?
    df_full = df_full.iloc[0:210]#filter on the 210 cortical BNA regions
    
    #calculate the percentage overlap of the tumor voxels with the brain region
    df_full['perc'] = df_overlap['voxels_overlap'].div(df_volumes['voxels_all'].values)*100 
    
    #define peritumoral regions ase regions that have 12% or more overlap with the tumor
    df_full['peritumor'] = df_full['perc'].apply(lambda x: 1 if x >= 12 else 0)
    df_full['roi']= np.arange(1,211)
    df_full['sub'] = sub
    
    
    #store subjects overlap information in dataframe in list 
    df_li.append(df_full[['roi','perc','peritumor', 'sub']])

#make dataframe containing all subjects overlap information
df_all = pd.concat(df_li) 



#%%

#save dataframe 
df_all.to_csv('/desired/path/to/folder/df_overlap.csv')

#%% 
########################
#DATA EXPLORATION TO SET THRESHOLD FOR PERITUMORAL AREA (I.E. 12%)
########################

#plot histogram of overlaps to determine threshold for peritumoral area 

n, bins, patches = plt.hist(df_all['perc'], bins = 100)
plt.ylim(ymin = 0, ymax = 100)

