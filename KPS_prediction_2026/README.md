## KPS prediction in contrast-enhancing glioma
This is the official code repository for the paper: **Predicting long-term postoperative functional status in contrast-enhancing glioma**

**Preregistration:** [View on OSF](https://osf.io/f29xm/overview)

**DOI:** TBD

<img width="640" height="588" alt="image" src="https://github.com/user-attachments/assets/3ade13f4-1184-42b3-86b3-543ba9c34937" />

## PROJECT OVERVIEW
This project predicts the long-term postoperative Karnofsky performance score (KPS) using only preoperative data. It uses a population of patients with contrast-enhancing glioma undergoing resection for the first time. The model predicts a three level outcome at 12 months postoperative: 1) death (KPS = 0) 2) functional dependence (KPS = 10 - 60) 3) functional independence (KPS = 70 - 100). The input features are based on the following 3 types: clinical, radiomics, and tumor volumetrics features. To obtain them, the following tools and databases were used:
1) clinical: Dutch Brain Tumor Registry in combination with available in-house research based databases. Additionally, for missing main outcome, a trained clinical researcher assigned the postoperative KPS based on the information available in the patient health records.
2) radiomics: two toolboxes were used in this process. First, the PICTURE toolbox was used to obtain the tumor masks. Second, the GSI-RADS toolbox was used to obtain the automatic reports of these tumor masks which contain quantifiable tumor related metrics.
3) tumor volumetrics: the tumor mask as obtained by the PICTURE toolbox produced 3 tumor components: necrotic core, enhancing component, and T2 hyperintensity. Then, a bespoke script was used to calculate the volumes of each of these components in the MNI space.

## REPOSITORY STRUCTURE
```
├── data/
│   ├── raw/              # Raw files containing clinical information. MRI scans stored separately on server
│   └── processed/        # Processed dataframes after merging and harmonizing input features
├── envs/                 # YAML files for creating Python virtual environments with dependencies
│   ├── preprocessing_env.yml
│   └── modeling_env.yml
├── src_preprocessing/    # Scripts for preprocessing the 3 groups of input features
├── src_modeling/         # Scripts for model training, testing, and evaluation
├── models/               # Stored model results and hyperparameter tuning specifications
└── subj_ids/             # Subject IDs, excluded subjects, train/test split information
```

**Note:** For patient privacy reasons, the folders `subj_ids/`, `data/raw/`, and `data/processed/` are excluded from this repository.

## DEPENDENCIES
This project uses two YAML configuration files:
preprocessing_env.yaml: dependencies required for running scripts located in /src_preprocessing
modeling_env.yaml: dependencies required for running scripts located in /src_modeling
Both files must be present before running the code.

To create a virtual environment that provides reproducible results, run the following:
conda env create -f preprocessing_env.yaml
conda activate preprocessing_env

## HOW TO RUN
1. Setup virtual environments
2. /src_preprocessing: follow script order --> interim files stored in /data/processed and /subjs_ids
3. /src_modeling: follow script order
4. final model results: models/ordinal_reg_xgb/ or models/mcp

## CITATIONS/LICENSES
CC BY 4.0: Anyone can use, modify, redistribute, but must give credit.
Citation: TO BE ADDED ONCE KNOWN

## CONTACT
e.koderman@amsterdamumc.nl
