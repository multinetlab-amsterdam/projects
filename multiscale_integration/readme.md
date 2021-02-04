This is a project that has resulted in the paper "Neuronal morphology and physiology, functional brain networks and memory in temporal lobe epilepsy", which can be found on bioRxiv. 

biorxiv.org/content/10.1101/2021.01.31.428369v1 

The repo contains:
- multiscale_integration_memory_tle.m --> the main code to do unilayer network analysis, statistical analysis and create plots
- full_output.mat --> the full output generated of the main code when ran on the files within the data/ folder
- data/micro_clinical_data.mat --> contains all clinical data (identical to table 1 in the paper) and micro-scale variables
- data/aal_fmri_adjmats.mat --> contains the fMRI adjacency matrices
- data/aal_meg_adjmats_t.mat --> contains the MEG adjacency matrices

In addition to these files, we used:
- The Brain Connectivity Toolbox version 2017-15-01 https://sites.google.com/site/bctnet/Home/updates/version2017-15-01majorupdate 
- Our Python Multilayer code https://github.com/multinetlab-amsterdam/data_analysis/tree/Multilayer/Multilayer 

