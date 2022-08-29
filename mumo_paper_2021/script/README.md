# mumo_script_clean_v1.m
This script is used to convert pre-processed neuroimaging data from different modalities (dMRI, fMRI, MEG) into connectivity matrices, and to perform monolayer network analyses on and construct supra-adjacency matrices for multilayer analyses from these matrices. It consists of four parts:

1.	Calculating connectivity
2.	Randomizing data
3.	Constructing supra-adjacency matrices
4.	Computing network measures

A more detailed breakdown of these four parts is provided below. The following functions, to be found in the /projects/mumo/functions folder, are used in the script:

* pli_matteo.m
* fft_filt_BNA.m
* remove_empty_regions.m
* kruskal_algorithm.m

It further uses several functions from the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/), version 2019_03_03.


### Part 1: Calculating connectivity
In Part 1, functional connectivity for MEG and fMRI timeseries is determined and a dMRI structural connectivity matrix is constructed. Empty regions (i.e. missing signal) are then removed from these matrices, matrices are normalized, and minimum spanning trees are constructed.

### Part 2: Randomizing data
In Part 2, the normalized matrices from Part 1 are shuffled while maintaining degree-, weight-, and strength- distributions to create randomized data. From the BCT website: _"Values of many network measures are greatly influenced by basic network characteristics, such as the number of nodes and links, and the degree distribution. Consequently, the significance of network statistics should often be established by comparison with statistics calculated on null network models."_ Additionally, minimum spanning trees are constructed from these randomized data.

### Part 3: Constructing supra-adjacency matrices
In Part 3, four supra-adjacency matrices are constructed for later (external) multilayer analysis: 

* *supra_weighted*, with the normalized weighted connectivity matrices;
* *rand_supra_weighted*, with the randomized weighted connectivity matrices;
* *supra_mst*, with the minimum spanning tree matrices; and
* *rand_supra_mst*, with the minimum spanning tree matrices constructed from the randomized data.

### Part 4: Computing network measures
In Part 4, nodal eigenvector centrality (EC) is calculated for each of the monolayers, and subsequently EC of the fronto-parietal network (FPN) is extracted and averaged. Here, regions belonging to the FPN are based on the classical seven-network parcellation by [Yeo et al. (2011)](https://doi.org/10.1152/jn.00338.2011). A dataset containing multilayer nodal EC is then loaded into the workspace, for which again EC of the FPN is extracted and averaged. Finally, the multilayer EC of the FPN and the multiple monolayer EC of the FPN are written to a table and exported to a .csv file for easy importation of the network measures in SPSS.
