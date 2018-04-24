import sys
print (sys.argv)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import pymzml
from os import listdir
import click
import ursgal


def remove_nan(x):
	return x.dropna(how = 'any')


def extract_runs_1(x):
	spectra_1 = []
	for spectrum in x:
		if spectrum["ms level"] == 1: 
				x =spectrum.mz
				y =spectrum.i
				z =x['TIC'].peaks
				spectra_1 = pd.DataFrame({'m/z':x, "Signal_Intensity":y})
				ID_1 = pd.Series([1] * (len(spectra_1)))
				spectra_1['Identity']  = ID_1.values
				spectra_1 = remove_nan(spectra_1)
		return spectra_1


def extract_runs_2(x):
	spectra_2 = []
	for spectrum in x:
		if spectrum["ms level"] == 1: 
				x =spectrum.mz
				y =spectrum.i
				#z =x['TIC'].peaks
				spectra_2 = pd.DataFrame({'m/z':x, "Signal_Intensity":y})
				ID_2 = pd.Series([2] * (len(spectra_2)))
				spectra_2['Identity']  = ID_2.values
				spectra_2 = remove_nan(spectra_2)
		return spectra_2
					
uc = ursgal.UController(
	params = {'database': 'fasta_file',
			'frag_mass_tolerance': 20,
			'frag_mass_tolerance_unit':'ppm',
			'modifications' : [
				'M,opt,any,Oxidation',        
				'C,fix,any,Carbamidomethyl',  
				'*,opt,Prot-N-term,Acetyl',],})

def psm_search_1(x):
	search_result_1 = uc.search(
		input_file       = 'x',
		engine           = 'omssa_2_1_9',
		force            = True,
		output_file_name = 'search_result_1.csv')
	psm_1 = pd.read_csv('search_result_1.csv')
	ID_1 = pd.Series([1] * (len(psm_1)))
	psm_1['Identity'] = ID_1.values
	return psm_1



def psm_search_2(x):
	search_result_2 = uc.search(
		input_file       = 'x',
		engine           = 'omssa_2_1_9',
		force            = True,
		output_file_name = 'search_result_2.csv')
	psm_2 = pd.read_csv('search_result_2.csv')
	ID_2 = pd.Series([2] * (len(psm_2)))
	psm_2['Identity'] = ID_2.values
	return psm_2

 

