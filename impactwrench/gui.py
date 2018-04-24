#!/usr/bin/env python3
import sys
print (sys.argv)
import os
import click
import os
import glob
from os import listdir
from plot  import *
from extract import *
from p_value_correction import *
from principal_component_analysis import *
from signal_processing import *
from perceptron import *
from logistic_regression import *
from gooey import Gooey
import argparse

@Gooey()
def main():
    parser = GooeyParser(description='LC-MS/MS Signal Processing')

    parser.add_argument(
    'required_field',
    metavar='Filepath',
    help='Enter a valid file path!')
    args = parser.parse_args()


    dicts ={}
    directory = args.Filepath
    keys = ["file_1", "file_2"]
    files_dir =  listdir(directory)
    filenames = []
    for i in keys:
        for names in files_dir:
            if names.endswith(".mzML"):
                filenames.append(names)
            dicts[i] = names

    run1 = pymzml.run.Reader(dicts['file_1'], MS1_Precision =5e-6, MSn_Precision = 20e-6)
    with click.progressbar(run1) as bar:
        for spectrum in bar:
            if spectrum["ms level"] == 1: 
                x =spectrum.mz
                y =spectrum.i
                z =run1['TIC'].peaks
            mass_to_charge =np.asarray(x)
            signal_intensity =np.asarray(y)

    visualization(x, y,z)
    plot_baseline_correction(x, y)
    smooth(y, 100)
    plot_resampling(z)
    plot_smoothing(x, y)
    plot_scatter(x, y)