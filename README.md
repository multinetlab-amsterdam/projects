#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  6 11:48:36 2026

@author: ekoderman
"""
This is the official code repository for the paper: Predicting long-term postoperative functional status in contrast-enhancing glioma
Preregistration: https://osf.io/f29xm/overview
DOI: TBD

PROJECT OVERVIEW
This project predicts the long-term postoperative Karnofsky performance score (KPS) using only preoperative data. It uses a population of patients with contrast-enhancing glioma undergoing resection for the first time. The model predicts a three level outcome at 12 months postoperative: 1) death (KPS = 0) 2) functional dependence (KPS = 10 - 60) 3) functional independence (KPS = 70 - 100). The input features are based on the following 3 types: clinical, radiomics, and tumor volumetrics features. To obtain them, the following tools and databases were used:
1) clinical: Dutch Brain Tumor Registry in combination with available in-house research based databases. Additionally, for missing main outcome, a trained clinical researcher assigned the postoperative KPS based on the information available in the patient health records.
2) radiomics: two toolboxes were used in this process. First, the PICTURE toolbox was used to obtain the tumor masks. Second, the GSI-RADS toolbox was used to obtain the automatic reports of these tumor masks which contain quantifiable tumor related metrics.
3) tumor volumetrics: the tumor mask as obtained by the PICTURE toolbox produced 3 tumor components: necrotic core, enhancing component, and T2 hyperintensity. Then, a bespoke script was used to calculate the volumes of each of these components in the MNI space.

REPOSITORY STRUCTURE
data/raw: raw files containing clinical information. The MRI scans are stored separetely on the server.
data/processed: processed dataframes after the three groups of input features have been merged and harmonized.
envs: the two yaml files (one for preprocessing and one for modeling) for creating python virtual environments with all the listed dependendencies.
src_modeling: scripts used for the modeling steps (model training, testing, and evaluation)
src_preprocessing: scripts used for the preprocessing steps of the 3 groups of input features
models: stored information regarding models trained using the src_modeling scripts (model results,  hyperparameter tuning specifications, ...)
subj_ids: stored information regarding subject ids (excluded subjects, subjects overlap between clinical and brain imaging data availability, and the train_test_ID split of the 70% development set and the 30% test set).

DATA DESCRIPTION
Important folders have additional README.md explanations (data/raw, data/processed, subj_ids). 
Privacy concerns: This project involves sensitive data. Raw data (such as date of birth or MRI scans) are not included in this repository. All analyses were performed on anonymized or pseudonymized datasets stored in secure, access-controlled environments. Access to the original data requires appropriate approvals. No direct identifiers, access credentials, or linkage files are shared here.

DEPENDENCIES
This project uses two YAML configuration files:
preprocessing_env.yaml: dependencies required for running scripts located in /src_preprocessing
modeling_env.yaml: dependencies required for running scripts located in /src_modeling
Both files must be present before running the code.

To create a virtual environment that provides reproducible results, run the following:
conda env create -f preprocessing_env.yaml
conda activate preprocessing_env

HOW TO RUN
1. Setup virtual environments
3. The four files should be present in the data/raw --> see data/raw/README.md
2. /src_preprocessing: follow script order --> interim files stored in /data/processed and /subjs_ids
3. /src_modeling: follow script order
4. final model results: models/ordinal_reg_xgb/

CITATIONS/LICENSES
CC BY 4.0: Anyone can use, modify, redistribute, but must give credit.
Citation: TO BE ADDED ONCE KNOWN

CONTACT
e.koderman@amsterdamumc.nl