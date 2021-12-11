#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 09:26:32 2018

@author: oudelaar, dowens
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from os import listdir, makedirs
from os.path import isfile, join



plt.rcParams['pdf.fonttype'] = 42

# This script is working to plot individual matrices outputted from
# 3_TileR_replicateAnalysis_iced_normMatrices_200327
# and
# merged matrices outputted from
# 6_TileR_subtractionAndSubsetR_DESeq_filt_200322

# this version will plot custom colours from Dark2 colour palette in RColourbrewer which is colourblind friendly

# D OWENS 1st APRIL 2020



###############################################################################
# CUSTOM CMAP FUNCTIONS

### from https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]



def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0, 1, len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp








###############################################################################

# Specify directories, file names, bin number

# set the root directory of matrices
my_dir = "C:/Users/Dominic/Desktop/Work/Paper/Bioinformatics/Tile-C/CTCF-KO/outMatrixFolder_indiReps_iced/"
#my_dir = "C:/Users/Dominic/Desktop/Work/Paper/Bioinformatics/Tile-C/CTCF-KO/outMatrixFolder_merged_iced/"


# get the list of directories (this is sample names
samples = [f for f in listdir(my_dir) if isfile(join(my_dir, f))]

# output folder (create if not exist)
my_dir_out = "C:/Users/Dominic/Desktop/Work/Paper/Bioinformatics/Tile-C/CTCF-KO/python/plots/merged_iced_varThresh_customCmap_indiReps/"
makedirs(my_dir_out, exist_ok=True)


## define which colours to use in hex format

Un_col = "#7570b3"
Flk1_col = "#d95f02"
CD41_col = "#1b9e77"


###############################################################################


# start the for sample loop
for sample in samples:

    # get the matrix file name
    mat_file = my_dir + sample

    # get the plot output name
    out_file = my_dir_out + sample + ".2kb"


    # bin number is different for whole and windowed matrices
    if "window" in sample:
        n_bins = 311 # standard window over Runx1 gene used in thesis
        #n_bins = 236 # EOMES smaller window
    else:
        n_bins = 1268 # 2kb res Runx1 region

    # can different thresholds for subtracts (increased contrast)
    if "subtract" in sample:
        thresh = 97
    else:
        thresh = 94

    # get the matrix
    matrix = np.zeros((n_bins, n_bins))

    with open(mat_file) as mat:
        for line in mat:
            x1, x2, count = line.split()
            bin1 = int(x1)
            bin2 = int(x2)
            matrix[bin1, bin2] = float(count)
            matrix[bin2, bin1] = float(count)

    mask = np.tri(matrix.shape[0], k=-1)
    matrix_half = np.ma.array(matrix, mask=mask)

    threshold = np.percentile(matrix_half, thresh)


    # set up some cmap basics that will be constant for all
    N = 256
    vals = np.ones((N, 4))
    vmin = 0.001

    # colour is different for normalised and subtracted matrices
    if "subtract" and "CD41_minus_Flk1" in sample:
        color = "CD41_minus_Flk1"
        vmin = -threshold
        hex1 = CD41_col
        hex2 = Flk1_col
        hex_list = [hex2, '#ffffff', hex1]
    elif "subtract" and "Flk1_minus_Un" in sample:
        color = "Flk1_minus_Un"
        vmin = -threshold
        hex1 = Flk1_col
        hex2 = Un_col
        hex_list = [hex2, '#ffffff', hex1]
    elif "Un" in sample:
        color = "Un"
        hex1 = Un_col
        hex_list = ['#ffffff', hex1, '#000000']
    elif "Flk1" in sample:
        color = "Flk1"
        hex1 = Flk1_col
        hex_list = ['#ffffff', hex1, '#000000']
    elif "CD41" in sample:
        color = "CD41"
        hex1 = CD41_col
        hex_list = ['#ffffff', hex1, '#000000']


    # make the custom cmap from hexlist
    newcmp = get_continuous_cmap(hex_list)



    ###############################################################################
    # PLOT

    plot = plt.imshow(matrix_half,
                      interpolation="nearest",
                      origin="upper", vmin=vmin,
                      vmax=threshold,
                      cmap=newcmp)

    # turn off the axis around the plot
    plt.axis('off')

    # turn on or off the colorbar
    #cbar = plt.colorbar()

    plt.savefig(out_file + "_" + color + "_" + str(threshold) + ".pdf",
                dpi=1000)

    ###############################################################################

    #report = my_dir_out + sample + "_2kb_report.txt"

    #r = open(report, "w")
    #r.write("Number of contacts:\t" + str(np.sum(matrix)))
    #r.close()
