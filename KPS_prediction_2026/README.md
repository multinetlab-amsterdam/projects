## KPS prediction in contrast-enhancing glioma
This is the official code repository for the paper: **Predicting long-term postoperative functional status in contrast-enhancing glioma**.

**Preregistration:** [View on OSF](https://osf.io/f29xm/overview)

**Open-source app:** [App](https://glioma-kps-prediction.streamlit.app/)

**DOI:** TBD


## Project overview
This project develops a model to predict the long-term postoperative Karnofsky performance score (KPS) using only preoperative data. It uses a population of patients with contrast-enhancing glioma undergoing resection for the first time. The model predicts a three level outcome at 12 months postoperatively:

      1) mortality (KPS = 0)
      2) functional dependence (KPS = 10 - 60)
      3) functional independence (KPS = 70 - 100)

The input features are the following: clinical, radiomics, and tumor volumetrics features. The following toolboxes were used in the project:

2) **Radiomics:** First, the [PICTURE toolbox](https://gitlab.com/picture-production/picture-qni-robust-glioma-segmentation) was used to obtain the tumor masks. Second, the [GSI-RADS toolbox](https://github.com/SINTEFMedtek/GSI-RADS) was used to obtain the automatic reports of these tumor masks which contain quantifiable tumor related metrics.

3) **Tumor volumetrics:** the tumor mask as obtained by the PICTURE toolbox produced 3 tumor components: necrotic core, enhancing component, and T2 hyperintensity. Then, a bespoke script was used to calculate the volumes of each of these components in the MNI space (`06_tumor_volume_components.py`).

## Repository structure
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

**Note:** For patient privacy reasons, the folders `subj_ids/`, `models/`, `data/raw/`, `data/processed/` are excluded from this repository.

## Dependencies
This project uses two YAML configuration files: 

`preprocessing_env.yaml`: dependencies required for running scripts in `/src_preprocessing`

`modeling_env.yaml`: dependencies required for running scripts in `/src_modeling`

To create a virtual environment that provides reproducible results, run the following:

`conda env create -f preprocessing_env.yaml`

`conda activate preprocessing_env`

## How to run
1. Setup virtual environments
2. /src_preprocessing: follow script order --> interim files stored in /data/processed and /subjs_ids
3. /src_modeling: follow script order
4. final model results: stored in models/ordinal_reg_xgb/ or models/mcp

## Citations / licenses
CC BY 4.0: Anyone can use, modify, redistribute, but must give credit.

Citation: TBD

## Contact
e.koderman@amsterdamumc.nl
