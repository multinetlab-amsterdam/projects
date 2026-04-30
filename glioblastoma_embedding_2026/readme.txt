This repository includes all the code used in the analysis pipeline of the project "Structural network embedding governs peritumor and distant pathological brain activity in glioblastoma" by Zimmermann et al., 2026. 

The script main.py includes the loading of the data and most statistical analyses. 

The functions used to create all metrics, to load the data (i.e. to create the dataframes used in the analysis in main.py) and to run statistical tests can be found in the folder scr. Within this folder you can find several scripts:

-  dTOR_compute_fiber_weights.m, tumor_tract_density_indices.m,whole_brain_tdm.m --> scripts to extract the normative embedding metrics (TDM, LTDI) 

- data_loading.py --> script with functions to create the dataframes used in the analyses (putting together brain activity, tumor embedding and 		      clinical data). Broadband power was extracted using our publicly available script at  https://github.com/multinetlab-		              amsterdam/data_analysis/tree/master/fooof

- statistics.py --> script with functions for the statistical analyses

- mixed_anova.r --> script for the mixed ANOVA analysis (type 3).

