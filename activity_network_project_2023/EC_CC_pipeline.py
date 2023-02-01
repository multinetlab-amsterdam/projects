#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Network metrics pipeline.

   Script to threshold the networks, calculate network metriccs and standadrize those based on regional data from HCs.

"""

__author__ = "Mona Lilo Margarethe Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "2022/02/11"   
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
import os
import warnings
warnings.filterwarnings('ignore')
import pickle
from datetime import datetime

# Third party imports ###
import numpy as np
import pandas as pd #version 1.1.5
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import re
import networkx as nx #version 2.3




#%%
########################
#FUNCTIONS NETWORK METRIC CALCULATION
########################

def densthr(d, i,DIAGNOSTIC=False):

    """Creating a binarized graph with a specific density

    Parameters
    ---------
    d: float
        density value

    i: numpy matrix
        connectivity matrix

    Returns
    -------
    finaldensity: float
        final density value

    G1: networkx graph
        graph with the specified density
        
    Author
    -------
    Fernando Nobrega Santos

    """

    np.fill_diagonal(i,0)
    # Will flatten it and rank corr values (stronger comes first).
    temp = sorted(i.ravel(), reverse=True)
    size = len(i)
   
    cutoff = np.ceil(d * (size * (size-1)))
    
    tre = temp[int(cutoff)]
    G0 = nx.from_numpy_matrix(i)
    G0.remove_edges_from(list(nx.selfloop_edges(G0)))
    G1 = nx.from_numpy_matrix(i)
    for u,v,a in G0.edges(data=True):
        if (a.get('weight')) <= tre:
            G1.remove_edge(u, v)
    finaldensity = nx.density(G1)
    if DIAGNOSTIC==True: # for sanity check :D
        print(finaldensity)

    return G1




def calculate_network_metrics(PLI_dir, output_dir, freq_dict, HC = False): 
    """
    Function to calculate the eigenvector centrality and clustering coefficient 
    of the thresholded (density 20,30%) networks based on PLI adjacency matrices. 
    After creating the thresholded networks, the eigenvector centrality and 
    clustering coefficients are calculated on the binarized (unweighted network).
    These are stored per frequency band and density in a pandas DataFrame.

    Parameters
    ----------
    PLI_dir : str,
        directory where the PLIs are stored
    output_dir : str,
        directory where output dfs with EC and CC should be stored.
    freq_dict: dict or list, 
        dictionary of list containing the frequencies that one is intereseted in. 
        Example: freq_dict = {'delta': [0.5, 4],'theta':[4, 8],'lower_alpha':[8, 10]}
    HC : bool, optional
        boolean to set whether the data of HCs is being analysed. The default is False.

    Returns
    -------
    None.

    """
    # --- Step 1.) create thresholded graphs from PLIs --- #
    for key, value in freq_dict.items():

        #access right directory were PLI is stored    
        pli_directory = PLI_dir + key + '/*'
        print(pli_directory)

        #iterate through densities 
        for d in [0.1,0.2,0.3]:
            
            #initialize list to store all measures for one density 
            dfs_density = []
            
            #loop through PLI files  
            for pli in glob.glob(pli_directory):
                print(pli)
                
                #determine which subject we are working with
                if HC == True: 
                    sub = re.search('([A-z]{4}_\d{3,4})', pli)[0]
                
                else:
                    sub = re.search('(sub-\d{4})', pli)[0]
                
                #load in the PLI matrix
                adj_mat = pd.read_csv(pli).to_numpy()
                
                
                #threshold the matrix based on the given density
                G = densthr(d, adj_mat)
                
                # --- Step 2: Compute network metrics --- # 
                print('Computing network metrics for sub')
                print(sub)
                
                ### Eigenvector centrality ###
                eigen = nx.eigenvector_centrality(G)
                df_EC = pd.DataFrame(eigen.items(), columns = ['rois', 'EC'])
                df_EC['sub'] = sub
                df_EC['rois'] = np.arange(1, 211)
    
                ### Clustering coefficient ###
                clustering = nx.clustering(G)
                df_CC = pd.DataFrame(clustering.items(), columns = ['rois', 'CC'])
                df_CC['sub'] = sub
                df_CC['rois'] = np.arange(1, 211)
                
                df = pd.concat([df_EC, df_CC], axis = 1)
                
                
                #save subs df in density list
                dfs_density.append(df)
                
            # --- Step 3.) make df with all subs EC and CC for current frequency and density --- #
            df_all_density = pd.concat(dfs_density, axis = 0)
            
            #drop duplicate columns (rois, sub) 
            df_all_density = df_all_density.loc[:, ~df_all_density.columns.duplicated()]
            
            
            # --- Step 4.) Save df ---# 
            df_all_density.to_csv(output_dir + key + '/desired_file_name_' + str(d) + '.csv')


#%%
########################
#CALCULATE NETWORK METRICS
########################

freq_dict = {'delta': [0.5, 4],'theta':[4, 8],'lower_alpha':[8, 10]}
PLI_dir = '/folder/where/PLIs/are/stored'
output_dir = '/desired/output/folder'

calculate_network_metrics(PLI_dir, output_dir, freq_dict, HC=True)
calculate_network_metrics(PLI_dir, output_dir, freq_dict, HC=False)



#%%

########################
# FUNCTION STANDARDIZATION OF REGIONAL CLUSTERING AND CENTRALITY
########################

def standardize_metrics(input_dir, input_HC_dir, output_dir, freq_dic, HC_subs_list, HC = False): 
    """
    Function used to standardize network metrics based on regional mean and std of HCs, 
    to aquire a regional measure of deviation from HCs. These standardized 
    metrics are stored in the original the dataframe.
    

    Parameters
    ----------
    input_dir : str,
        directory where df with network metrics is stored per frequency band (patients)
    input_HC_dir : str,
        directory where df with network metrics is stored per frequency band (HCs)
    output_dir : str,
        directory where dfs including the standardized CC and EC metrics will be stored.

    Returns
    -------
    None.

    """
    
    #loop through frequencies and densities
    for key, freq in freq_dict.items(): 
        print(freq)
        
        for d in [0.1,0.2,0.3]:
            print(d)
            
            #get network metric file for particular frequency and density
            HC_file = glob.glob(input_HC_dir + key + '/*' + str(d) + '*')[0]
            file = glob.glob(input_dir + key + '/*' + str(d) + '*')[0]
        
            #load in HC CC_EC dataframe for freq and density
            HC_df = pd.read_csv(HC_file, index_col=0)
            
            #only include matched HCs (on sex and age)
            HC_df = HC_df[HC_df['sub'].isin(HC_subs_li['sub'])]
            
            #determine if want to standardize HCs (on themselves) or patients & load in dataframe
            if HC == True:
                subs_df = HC_df
                
            else:
                subs_df = pd.read_csv(file, index_col=0)
            
            # --- Step 1.) Get regional mean and std of HCs EC and CC --- #
            ### Clustering Coefficient ###
            HC_mean_CC = HC_df.groupby('rois')['CC'].mean()
            HC_std_CC = subs_df.groupby('rois')['CC'].std()  
            
            #store means and std in dataframe 
            df_HC_CC = HC_mean_CC.to_frame()
            df_HC_CC['std'] = HC_std_CC
            df_HC_CC['rois'] = HC_std_CC.index
            df_HC_CC.rename(columns = {'CC':'CC_mean_HC', 'std':'CC_std_HC'}, inplace = True)
            
            ### Eigenvector centrality ###
            HC_mean_EC = HC_df.groupby('rois')['EC'].mean()
            HC_std_EC = subs_df.groupby('rois')['EC'].std()   
           
            #store means and std in dataframe 
            df_HC_EC = HC_mean_EC.to_frame()
            df_HC_EC['std'] = HC_std_EC
            df_HC_EC['rois'] = HC_std_EC.index
            df_HC_EC.rename(columns = {'EC':'EC_mean_HC', 'std':'EC_std_HC'}, inplace = True)
            
            
            
            # --- Step 2.) add means and std repeatedly to rows (i.e. rows) of subs dataframe to allow for easy standardization later--- #
            
            ### Clustering Coefficient ###
            subs_df['CC_mean_HC'] = subs_df['rois'].apply(lambda x: df_HC_CC['CC_mean_HC'].loc[x])
            subs_df['CC_std_HC'] = subs_df['rois'].apply(lambda x: df_HC_CC['CC_std_HC'].loc[x])
    
            ### Eigenvector centrality ###
            subs_df['EC_mean_HC'] = subs_df['rois'].apply(lambda x: df_HC_EC['EC_mean_HC'].loc[x])
            subs_df['EC_std_HC'] = subs_df['rois'].apply(lambda x: df_HC_EC['EC_std_HC'].loc[x])
      
                
            # --- Step 3.) Standardize CC and EC and add as new columns to dataframe --- # 
            
            ### Clustering Coefficient ###
            subs_df['CC_z'] = (subs_df['CC']- subs_df['CC_mean_HC'])/subs_df['CC_std_HC']
            
            ### Eigenvector centrality ###
            subs_df['EC_z'] = (subs_df['EC']- subs_df['EC_mean_HC'])/subs_df['EC_std_HC']
    
            # --- Step 4.) Save df for current frequency and density ---# 
            subs_df.to_csv(output_dir + key + '/desired_file_name_' + str(d) + '.csv')
            
#%%    
########################
#STANDARDIZE NETWORK METRICS
########################
    
input_dir = '/folder/where/raw/values/are/stored'
input_HC_dir = '/folder/where/raw/values/HCs/are/stored'

output_dir = '/folder/where/standardized/dataframes/should be stored'
output_dir_HC = '/folder/where/standardized/dataframes/HCs/should be stored'

HC_subs_li = pd.read_csv('/path/to/HCs/subject/list.csv', index_col = 0)
HC_subs_li.columns = ['sub']

freq_dict = {'delta': [0.5, 4],'theta':[4, 8],'alpha':[8, 10]}

standardize_metrics(input_dir, input_HC_dir, output_dir, freq_dict, HC_subs_li) #standardize patients metrics
standardize_metrics(input_HC_dir, input_HC_dir, output_dir_HC, freq_dict, HC_subs_li, HC = True) #standardize HCs metrics on themselves