#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Spin-test binomial confidence interval


Script to calculate how many of the spin-test correlations (from the permutations) are above the final correlation, 
to calculate binomial p-value on website: https://statpages.info/confint.html


"""

__author__ = "Mona Lilo Margarethe Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "2022/12/02"   
__status__ = "Finished" 

####################
# Review History   #
####################

# Reviewed by Eduarda Centeno 20230202

#%%
####################
# Libraries        #
####################

# Third party imports  ###
import pandas as pd #version 1.1.5
import numpy as np
import scipy.io as sio



#%%
#Get number of positives spin-test 
results = '/path/to/spin-test/result/results_df_EC_lower_alpha_3_patients_rest_non_tum.mat'
test = sio.loadmat(results)

#get results tab out of dictionary
results_data = test['results']

#get third element (i.e. the permutation results, correlation of every iteration)
results = results_data[:,3]
permutation_results = results[0]

#test how many of the permutation results are above correlation value 
positives = sum(abs(permutation_results[:,0]) > 0.27588) #change correlation value depending on which correlation you are looking at 

