#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Network metrics pipeline.

   Script to threshold the networks, calculate network metriccs and standadrize those based on regional data from HCs.

"""

__author__ = "Mona Lilo Margarethe Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "2022/04/04"
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
import re

# Third party imports ###
import numpy as np # version 1.23.5
import pandas as pd  # version 1.1.5
import networkx as nx  # version 2.3


# %%
########################
# FUNCTIONS NETWORK METRIC CALCULATION
########################


def densthr(d, i, DIAGNOSTIC=False):
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

    np.fill_diagonal(i, 0)
    # Will flatten it and rank corr values (stronger comes first).
    temp = sorted(i.ravel(), reverse=True)
    size = len(i)

    cutoff = np.ceil(d * (size * (size - 1)))

    tre = temp[int(cutoff)]
    G0 = nx.from_numpy_matrix(i)
    G0.remove_edges_from(list(nx.selfloop_edges(G0)))
    G1 = nx.from_numpy_matrix(i)
    for u, v, a in G0.edges(data=True):
        if (a.get("weight")) <= tre:
            G1.remove_edge(u, v)
    finaldensity = nx.density(G1)
    if DIAGNOSTIC == True:  # for sanity check :D
        print(finaldensity)
    return G1


def calculate_network_metrics(PLI_dir, output_dir, freq_dict, dens, HC=False):
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
    dens: list, 
        list of densities that we want to look at. For example [0.2, 0.3] to calculate for 20% and 30% density.
    HC : bool, optional
        boolean to set whether the data of HCs is being analysed. The default is False.

    Returns
    -------
    None.

    """
    # --- Step 1.) create thresholded graphs from PLIs --- #
    for key, value in freq_dict.items():
        # access right directory were PLI is stored
        pli_directory = f"{PLI_dir}{key}/*"
        print(pli_directory)

        # iterate through densities
        for d in dens: 
            # initialize list to store all measures for one density
            dfs_density = []

            # loop through PLI files
            for i, pli in enumerate(sorted(glob.glob(pli_directory))):
                

                # determine which subject we are working with
                if HC == True:
                    sub = re.search("([A-z]{4}_\d{3,4})", pli)[0]
                else:
                    sub = re.search("(sub-\d{4})", pli)[0]
                
                # load in the PLI matrix
                adj_mat = pd.read_csv(pli).to_numpy()
                

                # threshold the matrix based on the given density
                G = densthr(d, adj_mat)

                # --- Step 2: Compute network metrics --- #
                print("Computing network metrics for sub %d" %i )
                print(sub)

                ### Eigenvector centrality ###
                eigen = nx.eigenvector_centrality(G)
                df_EC = pd.DataFrame(eigen.items(), columns=["rois", "EC"])
                df_EC["sub"] = sub
                df_EC["rois"] = np.arange(1, 211) #number of rois

                ### Clustering coefficient ###
                clustering = nx.clustering(G)
                df_CC = pd.DataFrame(clustering.items(), columns=["rois", "CC"])
                df_CC["sub"] = sub
                df_CC["rois"] = np.arange(1, 211) #number of rois

                df = pd.concat([df_EC, df_CC], axis=1)

                # save subs df in density list
                dfs_density.append(df)
                
            # --- Step 3.) make df with all subs EC and CC for current frequency and density --- #
            df_all_density = pd.concat(dfs_density, axis=0)

            # drop duplicate columns (rois, sub)
            df_all_density = df_all_density.loc[:, ~df_all_density.columns.duplicated()]

            #--- Step 4.) Save df ---#
            df_all_density.to_csv(f"{output_dir}/{key}/EC_CC_{key}_den_{d}.csv")
              

# %%
########################
# CALCULATE NETWORK METRICS
########################

#specify the frequency that you want to look at 
freq_dict = {"delta": [0.5, 4], "theta": [4, 8], "lower_alpha": [8, 10]}

PLI_dir_pat = "/path/to/PLIs/folder/patients/"
PLI_dir_HC = "/path/to/PLIs/folder/HCs"


output_dir_pat = "/path/to/folder/where/output/should/be/stored/"
output_dir_HC = "/path/to/folder/where/output/should/be/stored/"


#run function to calculate the CC and EC for patients and HCs
calculate_network_metrics(PLI_dir_pat, output_dir_pat, freq_dict, [0.2, 0.3], HC=False)
calculate_network_metrics(PLI_dir_HC, output_dir_HC, freq_dict, [0.2, 0.3], HC=True)

