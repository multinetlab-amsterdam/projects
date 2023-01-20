#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 15:02:47 2022

Script to calculate how many of the spin-test correlations (from the permutations) are above correlation, 
to calculate binomial p-value on website: https://statpages.info/confint.html

@author: m.zimmermann
"""
import pandas as pd
import numpy as np
import scipy.io as sio




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

