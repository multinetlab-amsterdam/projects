This is a project that has resulted in the paper "Cellular Substrates of Functional Network Integration and Memory in Temporal Lobe Epilepsy", published in Cerebral Cortex in 2021 (open access): https://academic.oup.com/cercor/article/32/11/2424/6375261

The repo contains:
- multiscale_integration_memory_tle_final.m --> the main code to do unilayer network analysis, statistical analysis and create plots
- full_output_final.mat --> the full output generated of the main code when ran on the files within the data/ folder
- data/micro_clinical_data.mat --> contains all clinical data (identical to table 1 in the paper) and micro-scale variables
- data/aal_fmri_adjmats.mat --> contains the AAL fMRI adjacency matrices
- data/aal_meg_adjmats_t.mat --> contains the AAL MEG adjacency matrices
- data/additional_analyses.mat --> contains the data for all additional analyses (additional frequency bands, additional atlas)

In addition to these files, we used:
- The Brain Connectivity Toolbox version 2017-15-01 https://sites.google.com/site/bctnet/Home/updates/version2017-15-01majorupdate 
- Our Python Multilayer code https://github.com/multinetlab-amsterdam/data_analysis/tree/Multilayer/Multilayer 
- A slightly adjusted implementation of the Minimum Spanning Tree algorithm based on https://nl.mathworks.com/matlabcentral/fileexchange/13457-kruskal-algorithm

