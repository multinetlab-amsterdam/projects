# MuMoBrain
MuMo is a multimodal project in which MRI (resting-state, diffusion, spectroscopy), MEG, EEG, and NPA data was collected from healthy human subjects to explore the multimodal correlates of cognitive functioning through the framework of multilayer networks. 

## mumo_script_clean_v1.m
This script is used to convert pre-processed neuroimaging data from different modalities (fMRI, DWI, MEG) into connectivity matrices, and to perform monolayer network analyses on and construct supra-adjacency matrices for multilayer analyses from these matrices. It consists of four parts:

1.	Calculating connectivity
2.	Randomizing data
3.	Supra-adjacency matrices
4.	Network measures

A more detailed breakdown of these four parts is provided below. The following functions, to be found in the /functions folder, are used in the script:

* pli_matteo.m
* fft_filt_BNA.m
* remove_empty_regions.m
* kruskal_algorithm.m

It further uses several functions from the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/)


### Part 1: Calculating connectivity
In Part 1, functional connectivity for MEG and fMRI timeseries is determined and a DWI structural connectivity matrix is constructed. Empty regions (i.e. missing signal) are then removed from these matrices, data is normalized, and minimum spanning trees are constructed.

### Part 2: Randomizing data
In Part 2, the normalized matrices from Part 1 are shuffled while maintaining degree-, weight-, and strength- distributions to create randomized data. Minimum spanning trees are then constructed from these randomized data.

### Part 3: Supra-adjacency matrices
In Part 3, four supra-adjacency matrices are constructed for later (external) multilayer analysis: 

* *supra_weighted*, with the normalized weighted connectivity matrices;
* *rand_supra_weighted*, with the randomized weighted connectivity matrices;
* *supra_mst*, with the minimum spanning tree matrices; and
* *rand_supra_mst*, with the minimum spanning tree matrices constructed from the randomized data.

## Part 4: Network measures
In Part 4, nodal eigenvector centrality (EC) is calculated for each of the monolayers, and subsequently EC of the fronto-parietal network (FPN) is extracted and averaged. A dataset containing multilayer nodal EC is then loaded into the workspace, for which again EC of the FPN is extracted and averaged. Finally, the multilayer EC of the FPN and the multiple monolayer EC of the FPN are written to a table and exported to a .csv file for easy importation of the network measures in SPSS.

