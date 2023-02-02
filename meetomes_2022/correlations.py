#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 13:10:21 2023


This file includes all analyses performed as part of the "Meetomes" project. 
After identifying six interactional roles and counting how often each of the team members produced role-specific actions, 
we wanted to correlate the occurrence of such actions to network measures as well as meeting evaluations. For that purpose, we use this script.
"""

__author__ = "Marloes Bet"
__contact__ = "m.bet@amsterdamumc.nl"  # or l.douw@amsterdamumc.nl
__date__ = "2022"  ### Date it was created
__status__ = "Finished"

####################
# Review History   #
####################

####################
# Libraries        #
####################

# Standard imports
import os
import warnings

# Third party imports
import numpy as np # version 1.23.5
import pandas as pd # version 1.5.2
import statsmodels.api as sm # version 0.13.2
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from scipy import stats # version 1.9.3
import scipy.stats
from scipy.stats import pearsonr
import seaborn as sns # version 0.12.1
import matplotlib.pyplot as plt # version 3.6.2
from matplotlib.pyplot import figure

warnings.filterwarnings("ignore")

# %%

# 0. Set working directory
os.chdir("path/to/data/")

# %%

# 1. Load data (header=0 takes the first row as column names)
data = pd.read_excel("finaldata_def.xls", header=0)
data["local_eff"] = data["local_eff"].astype(float)

# Remove survey item on mood, was only for validation purposes
data = data.drop(columns='mood_norm')

# %%
# 2. Transform data
# Calculate total number of role-specific utterances
data["total_actions"] = (
    data["chair"]
    + data["skeptic"]
    + data["expert"]
    + data["clarifier"]
    + data["connector"]
    + data["practical"]
)

# Change counts of role-specific actions to proportions relative to the total count
roles = ["chair", "skeptic", "expert", "clarifier", "connector", "practical"]
for col in roles:
    data[col] = (
        data[col] / data["total_actions"]
    )  # calculate relative proportions of role-specific actions
    data[col] = data[col].fillna(0)  # omit NA

# Rename integers to letters from the alphabet; otherwise, analysis will treat variable "name" as integer
alphabet = [
    "AA",
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "U",
    "V",
    "W",
    "X",
    "Y",
    "Z",
    "ZZ",
]

counter = 0

for j in data["name"]:
    for i in range(28):
        print(i)
        if data.iloc[counter, 1] == i:  # refers to the name column
            data.iloc[counter, 1] = alphabet[i]
    counter += 1

# Save outcomes in following dataframe
cols = [
    "chair_stdB",
    "chair_p",
    "skeptic_stdB",
    "skeptic_p",
    "expert_stdB",
    "expert_p",
    "clarifier_stdB",
    "clarifier_p",
    "connector_stdB",
    "connector_p",
    "practical_stdB",
    "practical_p",
]
rows = [
    "local_eff",
    "nod_str",
    "partcoef",
    "wmdegree",
    "betwcent",
    "effect_norm",
    "worth_norm",
    "overall_norm",
    "motivation_norm",
    "education_norm",
    "translation_norm",
    "connection_norm",
    "understand_norm",
]
output = pd.DataFrame(columns=cols, index=rows)

# %%
# 3. Correlate proportions of role-specific actions to meeting evaluations

# Variables to correlate to each other
q_items = [
    "effect_norm",
    "worth_norm",
    "overall_norm",
    "motivation_norm",
    "education_norm",
    "translation_norm",
    "connection_norm",
    "understand_norm",
]

columns_for_analysis = ["name"] + roles + q_items
data_evals = data[columns_for_analysis].dropna()

# Reset indices
data_evals.reset_index(inplace=True, drop=True)

# Create dummies to indicate that the digits in the col 'name' are categories and not integers
dummies = pd.get_dummies(data_evals["name"])

# Merge the dummy df with the data df
data_evals = pd.concat([data_evals, dummies], axis=1)

# Rename dummy columns for clarity
for i in data_evals.iloc[:, 16:]:
    data_evals.rename(columns={str(i): "dummy_" + str(i)}, inplace=True)

# Loop over action proportions and correlate them to evaluations
for row_item in q_items:  # cells containing evaluations
    for col_item in roles:  # cells containing proportions
        string = "{} ~ {}".format(row_item, col_item)
        md = smf.mixedlm(string, data_evals, groups=data_evals["name"])
        mdf = md.fit()
        output.loc[row_item, col_item + "_p"] = mdf.pvalues[1]  # save p-value
        output.loc[row_item, col_item + "_stdB"] = mdf.params[1]  # save std Beta

        # Instantaneously report significant correlations
        if mdf.pvalues[1] <= 0.05:
            print(
                "Significant association between " + col_item + " and " + row_item + ":"
            )
            print("Coefficient: " + str(mdf.params[1]))
            print("P-value: " + str(mdf.pvalues[1]))
            print("Bonferroni-corrected P: " + str(mdf.pvalues[1].astype(float) * (len(q_items)*len(roles))))

# %%

# 4. Correlate proportions of role-specific actions to network measures
# Subset only the data that you need for these analysis and drop NA values
ntwrkmsrs = ["local_eff", "nod_str", "partcoef", "wmdegree", "betwcent"]
columns_for_analysis = ["name"] + roles + ntwrkmsrs
data_ntwrk = data[columns_for_analysis].dropna()

# Reset indices
data_ntwrk.reset_index(inplace=True, drop=True)

# Create dummies to indicate that the digits in the col 'name' are categories and not integers
dummies = pd.get_dummies(data_ntwrk["name"])

# Merge the dummy df with the data df
data_ntwrk = pd.concat([data_ntwrk, dummies], axis=1)

# Rename dummy columns for clarity
for i in data_ntwrk.iloc[:, 13:]:
    data_ntwrk.rename(columns={str(i): "dummy_" + str(i)}, inplace=True)

# Loop over action proportions and correlate them to network measures
for row_item in ntwrkmsrs:  # cells containing network measures
    for col_item in roles:  # cells containing action proportions
        string = "{} ~ {}".format(row_item, col_item)
        md = smf.mixedlm(
            string, data_ntwrk, groups=data_ntwrk["name"]
        )  # create linear mixed model
        mdf = md.fit()  # model fit parameters
        output.loc[row_item, col_item + "_p"] = mdf.pvalues[1]  # save p-value
        output.loc[row_item, col_item + "_stdB"] = mdf.params[
            1
        ]  # save std beta coefficient

        # Instantaneously report significant correlations
        if mdf.pvalues[1] <= 0.05:
            print(
                "Significant association between " + col_item + " and " + row_item + ":"
            )
            print("Std. B: " + str(mdf.params[1]))
            print("P-value: " + str(mdf.pvalues[1]))
            print("Bonferroni-corrected P: " + str(mdf.pvalues[1].astype(float) * (len(ntwrkmsrs)*len(roles))))

# Save output file
output.to_excel("output.xlsx")

# %%

# 5. Create a heatmap of the nodal network measure correlations
# Put the correlation coefficients that you want to plot in a matrix
r_values = [
    "chair_stdB",
    "skeptic_stdB",
    "expert_stdB",
    "clarifier_stdB",
    "connector_stdB",
    "practical_stdB",
]
mat = output.loc["local_eff":"betwcent", r_values].astype(float)

# Generate plot
netwmeasures = [
    "local efficiency",
    "nodal strength",
    "participation coefficient",
    "within-module degree",
    "betweenness centrality",
]
actiontypes = ["chair", "skeptical", "expert", "clarifying", "connecting", "practical"]
plot = sns.heatmap(
    mat, xticklabels=actiontypes, yticklabels=netwmeasures, center=0, cmap="magma"
)
plot.set(xlabel="Interactional role", ylabel="Network measure")
plot.set_title("Correlations between interactional roles and network measures")

# make sure that ticks are shown diagonally above the figure
plt.setp(plot.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.show()

# %%

# 6. Calculate correlations with global measures
globaldata = pd.read_excel("globaldata.xlsx")

# Make correlation matrix
globmat = globaldata.corr()


# Calculate p-values
def calculate_pvalues(df):
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how="outer")
    for r in df.columns:
        for c in df.columns:
            tmp = df[df[r].notnull() & df[c].notnull()]
            pvalues[r][c] = round(pearsonr(tmp[r], tmp[c])[1], 4)
    return pvalues


pvals = calculate_pvalues(globmat)
pvals_Bonfcor = calculate_pvalues(globmat) * (len(roles)*2)

# Reshape output
globmat = globmat.iloc[
    6:8, 0:6
]  # remove autocorrelations and duplicate correlation pairs
pvals = pvals.iloc[6:8, 0:6]  # do the same for the corresponding pvalues

globmat.loc["modularity_r"] = globmat.loc[
    "modularity"
]  # copy modularity coefs to a third row
globmat.rename(index={"modularity": "GE_norm_p"}, inplace=True)  # rename row
globmat.loc["modularity_p"] = pvals.loc[
    "modularity"
]  # add modularity p-values to global matrix
globmat.loc["GE_norm_p"] = pvals.loc[
    "GE_norm"
]  # transfer p-values to row underneath GE coefs

# Save output files
globmat.to_excel("globalmeasures.xlsx")
