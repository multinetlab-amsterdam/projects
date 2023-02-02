"""PLI construction

Scripts to construct functional networks using the PLI for the different frequency bands and densities.

"""

__author__ = "Mona Lilo Margarethe Zimmermann"
__contact__ = "m.l.m.zimmermann@amsterdamumc.nl"
__date__ = "January 2022"
__status__ = "Finished"

####################
# Review History   #
####################

# Reviewed by Eduarda Centeno 20230202

# %%
####################
# Libraries        #
####################

# Standard imports  ###
import os
import time # not used?
from datetime import date


# Third party imports ###
import numpy as np
import pandas as pd  # version 1.1.5
from numpy import fft
from scipy.signal import hilbert
from scipy import fft # not used?
import scipy # not used? 
import statsmodels.api as sm # not used?
from tqdm import tqdm
import matplotlib.pyplot as plt # not used?
import seaborn as sns # not used?
import plotly.graph_objects as go # not used?
from plotly.offline import plot as plot_off # not used?
from scipy import signal
from scipy.fft import fft, irfft2, ifft2, irfftn, irfft, ifft # not used?
from scipy.fftpack import hilbert # not used?
import matplotlib as mpl # not used?


# Internal imports ###
from alternative_Compute_broadband import find_paths
from fft_filt import fft_filt
from make_csv import make_csv  # not used?


# %%
########################
# FUNCTIONS
########################


### Base PLI function ###
def pli_matteo(filtered_timeseries):
    """Calculate the phase lag index (ref) between all pairs of regions

    Parameters
    ----------
    filtered_timeseries: ndarray
        Rows are timepoints, columns are rois/electrodes

    Returns
    -------
    m: ndarray
        Array with PLI between channel i and j averaged over timebins

    """

    N = nch = filtered_timeseries.shape[1]  ###why double asignement?
    m = np.zeros([N, N])

    complex_a = signal.hilbert(filtered_timeseries, axis=0)

    for i in np.arange(nch):
        for j in np.arange(nch):
            if i < j:
                m[i, j] = abs(
                    np.mean(np.sign(((complex_a[:, i] / complex_a[:, j]).imag)))
                )
    m = m + np.transpose(m)
    return m


##### Construct connectivity matrices in loop #####
def loop_PLI(subjects, freq_value, fs, nr_rois):
    """Function to calculate the phase lag index (PLI) as a measure of functional
    connectivity in a given frequency band for every subject.

    Parameters
    ----------
    subjects : str,
        path to the csv including all subjects.
    freq_value : list,
        contains the lower and upper frequency bound (e.g. [0.5, 48]).
    fs : int,
        sample frequency.
    nr_rois : int,
        number of rois to include (e.g. BNA atlas 246).

    Returns
    -------
    all_pli_matrices: np.array,
        n x 1 array containing the PLI matrices for every subject (n)

    filtered_timeseries: np.array,
        n x m array containing the filtered timeseries for ROI (m) for every subject (n)

    """
    all_subs_pli_matrices = []

    # per subject get epoch files
    for index, row in pd.read_csv(subjects).iterrows():
        print(
            "\n\n//////////// Subject " + str(index) + " on subject_list ////////////"
        )
        files_list = find_paths(
            main_dir=row["Path"],
            subject=row["Case_ID"],
            extension=extension,
            timepoint=row["MM"],
            atlas=row["Atlas"],
            start=row["Start"],
            end=row["End"],
            selection=row["Selection"],
        )

        if len(files_list) == 0:
            print("No ASCIIs available for: " + row)
            continue
        else:
            sub_all_matrices = []

            for epoch_file, name in zip(range(len(files_list)), files_list):
                print(epoch_file)
                print(name)

                # load in single epoch´s timseries
                timeseries = pd.read_csv(
                    files_list[epoch_file], index_col=False, header=None, delimiter="\t"
                )

                # split the epoch in 4 pieces
                timeseries_split = np.array_split(timeseries, 4)

                for i in range(0, 3):  # why 3?
                    # extract the epoch piece
                    epoch_piece = timeseries_split[i]

                    # filter timeseries of epoch piece into frequency band
                    filt_piece = fft_filt(
                        epoch_piece,
                        freq_value[0],
                        freq_value[1],
                        fs=fs,
                        nr_rois=nr_rois,
                    )[0]

                    # construct PLI for epoch piece
                    pli_matrix = pli_matteo(filt_piece)

                    # store every PLI in list
                    sub_all_matrices.append(pli_matrix)
            # average the pli pieces per subject
            sub_avg_pli = np.mean(sub_all_matrices, axis=0)

            # store PLIs for all subjects
            all_subs_pli_matrices.append(sub_avg_pli)
        # turn into array
        all_pli_matrices = np.array(all_subs_pli_matrices)
    return all_pli_matrices


# %%
########################
# CONSTRUCT PLIS
########################

# dictionary of frequencies for which to construct the PLIs
freq_dict = {
    "BB": [0.5, 48],
    "delta": [0.5, 4],
    "theta": [4, 8],
    "lower_alpha": [8, 10],
}


subs_csv = "path/to/csv/file"


for freq, freq_value in tqdm(freq_dict.items()):
    print("This is frequency: ")
    print(freq)

    ##########################
    # Settings               #
    ###########################

    # set nice level to 10, especially FOOOF algorithm is heavy!
    os.nice(10)

    # 1. Create correctly your list of subjects you want to process
    # an example is given here: 'example_MEG_list.csv'
    subject_list = subs_csv ## why duplicate?

    # 2. Define the type of file extension your are looking for
    extension = ".asc"  # extension type

    # 3. Select which roi or rois you want to analyze
    # if you want to analyze 1 roi, specify its number (nr_rois = (10,))
    nr_rois = np.arange(
        210
    )  # (10,) if only run for roi 11 # note that python indexes at 0!
    # if you want to analyze multiple rois, create list with these rois
    # (for example nr_rois = np.arange(78) for all 78 cortical AAL rois)

    # 4. Set sample frequency (1250 Hz for Elekta data)
    Fs = 1250  # sample frequency

    # 5. Set frequency range you want to study
    freq_range = freq_value  # frequency range you want to analyze

    # 6. Give output directory
    dir_output_pli = "/path/to/store/output"

    # 7. Do you want to see the plots?
    # plot_choice = True

    # 7a. Do you want to save the output?
    save_output = True  # you can save output

    ###########################
    # Run analysis            #
    ###########################

    # output filtered timeseries?
    all_pli_matrices = loop_PLI(subject_list, freq_value, Fs, nr_rois)

    # save output
    if save_output == True:
        subjects = pd.read_csv(subject_list, delimiter=",", header=0)
        print("\n.....Saving PLIs ......")
        for index, row in subjects.iterrows():
            if all_pli_matrices.shape[0] > 1:
                # extract subject´s pli
                slice_pli = all_pli_matrices[index, :, :]

                # turn it into a dataframe for easier handling
                df_pli = pd.DataFrame(slice_pli)

                # set the roi numbers as column names
                df_pli.columns = [roi for roi in nr_rois]

                # save the pli
                df_pli.to_csv(
                    path_or_buf=dir_output_pli
                    + row["Case_ID"]
                    + "_"
                    + str(row["MM"])
                    + "_"
                    + str(row["Atlas"])
                    + "_"
                    + str(date.today().strftime("%Y%m%d"))
                    + "_"
                    + "pli"
                    + ".csv",
                    header=True,
                    index=False,
                )
            # if theres only one pli
            elif all_pli_matrices.shape[0] == 1:
                df_pli = pd.DataFrame(all_pli_matrices[0, :, :])

                df_pli.columns = [roi for roi in nr_rois]

                df_pli.to_csv(
                    path_or_buf=dir_output_pli
                    + row["Case_ID"]
                    + "_"
                    + str(row["MM"])
                    + "_"
                    + str(row["Atlas"])
                    + "_"
                    + str(date.today().strftime("%Y%m%d"))
                    + "_"
                    + "pli"
                    + ".csv",
                    header=True,
                    index=False,
                )
            else:
                print("Nothing to do here")
