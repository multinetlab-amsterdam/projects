#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
""" Example script on how raincloud plots were created for the modelling paper. 
    This script is used to make a raincloud plot for the data calculated with 
    the AEC over epochs but can be applied to each of the FC metrics that we 
    included. We made these rainclouds for each of the FC metrics and pasted 
    them in powerpoint to create figure 2. 

"""

__author__ = "Shanna Kulik & Lucas Breedt"
__contact__ = "l.douw@amsterdamumc.nl"
__date__ = "2020/12/7"  ### Date it was created
__status__ = "Finished"


####################
# Review History   #
####################

# Reviewed and updated by Eduarda Centeno 20220818


####################
# Libraries        #
####################

# Third party imports
import numpy as np  # version 1.19.1
import matplotlib.pyplot as plt  # version 3.3.0
import pandas as pd  # version 1.1.0
import seaborn as sns  # version 0.11.1
import matplotlib as mpl


#%% Load data

#### See parameter 'decimal' - this ensures the commas in your excel file are 
#### recognized as decimal points
df = pd.read_csv(
    "/path/to/data.csv",
    sep=";",
    decimal=",",
)

#%% Core figure

# Define pre-settings
w = 20
h = 20
title_size = 20
xlab_size = 15
ylab_size = 20

# Merge dataframe from wide to long for sns.pointplot
df_long = pd.melt(df, value_vars=["AEC epochs", "av SC AEC epochs"])
df_long.tail()


# Set the amount of jitter and create a dataframe containing the jittered x-axis values
jitter_2 = 0.05  # horizontal dispersion of dots
np.random.seed(3)
df_jitter_2 = pd.DataFrame(
    np.random.normal(loc=0, scale=jitter_2, size=df.values.shape),
    columns=df.columns,
)

# Update the dataframe with adding a number based on the length on the columns. 
# Otherwise all datapoints would be at the same x-axis location.
df_jitter_2 += np.arange(len(df.columns))

#### This is where we plot the raincloud plots!
# Create empty figure and plot the individual datapoints
mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
    color=["mediumblue", "teal"]
)  # #BDD7EE
fig, ax = plt.subplots(figsize=(15, 9))

for col in df:
    ax.plot(
        df_jitter_2[col],
        df[col],
        "o",
        alpha=0.8,
        zorder=2,
        ms=10,
        mew=1.5,
    )

for value in df_long:
    sns.violinplot(
        x="variable",
        y="value",
        data=df_long,
        hue="variable",
        split="True",
        inner="quartile",
        cut=1,
        dodge=True,
    )
    sns.boxplot(
        x="variable",
        y="value",
        data=df_long,
        hue="variable",
        dodge=True,
        width=0.2,
        fliersize=2,
    )

    # Additonal settings
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels((["AEC", "Average SC"]), size=xlab_size)
    ax.set_xlim(-1, len(df.columns))
    ax.set(ylim=(-0.05, 0.4))
    ax.set_ylabel(
        "Maximum correlation (simulated vs empirical FC)",
        size=ylab_size,
    )
    ax.set_xlabel("")
    ax.set_title('Average AEC per subject of full timeseries and SC', size = title_size)
    sns.despine()
    ax.legend_.remove()
    plt.setp(ax.collections, alpha=0.1)
