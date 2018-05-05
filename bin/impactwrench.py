
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
from gooey import Gooey, GooeyParser
from gooey import *
import argparse
from folder_wide_processing import *

@Gooey(program_name = 'impactwrench', advanced = True)
def main():
	parser = GooeyParser(description='LC-MS/MS Signal Processing')

	parser.add_argument(
	'Filepath',
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
	
	click.echo('The next block of code performs LC-MS chemometric data evaluation?')

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
	

	click.echo('Are you satisfied with the defined signal processing parameters? [yn] ', nl=False)
	c = click.getchar()
	click.echo()
	if c == 'y':
		click.echo('The show must go on')
	elif c == 'n':
		click.echo('Abort!')
	else:
		click.echo('Invalid input :(')


	run1 = pymzml.run.Reader(dicts['file_1'], MS1_Precision =5e-6, MSn_Precision = 20e-6)
	x_1 = extract_runs_1(run1)

	run2 = pymzml.run.Reader(dicts['file_2'], MS1_Precision =5e-6, MSn_Precision = 20e-6)
	x_2 = extract_runs_2(run2)
	

	df = pd.concat([x_1, x_2])

	X = df.loc[:, ['m/z', 'Signal_Intensity']].values
	y = df.loc[:, ['Identity']].values
	y = np.ravel(y)

	

	X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=500,random_state=4)

	one_dimensional_decision_boundary(X_train, y_train)
	multiple_classifiers_decision_boundary(X_train, y_train)


# LC-MSMS Data extraction

   #fasta_file = click.prompt('Please give the name of the fasta file in directory ', type=str)

   #psm1 = psm_search_1(dicts['file_1'])
   #psm2 = psm_search_1(dicts['file_1'])

   #combined_psm = pd.concat([psm1, psm2])
   #X_1 = combined_psm.loc[:, ['OMSSA:evalue', 'Exp m/z', 'OMSSA:pvalue']].values
   #y_1 = combined_psm.loc[:, ['Identity']].values
   #y_1 = np.ravel(y_1)

   #p_value = combined_psm.loc[:, ['OMSSA:pvalue']]

#FDR Correction
   #corrected_pvalues = fdrcorrection(df9_pvalue, alpha=0.01, method='indep', is_sorted=True)

#splitting dataframe into train and test sets for LC-MS/MS Analysis

#X_train, X_test, y_train, y_test = train_test_split(X_1, y_1, train_size=,random_state=4)

#PCA Analysis
 #  plot_cumulative_sum(X_train)
  # best_pca_parameters(X_train)
   #singular_value_decomposition(X_train)
   #plot_pca(X_train)

#Classification and Regression
  # plot_decision_surface(X_train,y_train)
   #linear_svm(X_train, y_train)
   #non_linear_svm(X_train,y_train)
   #one_dimensional_decision_boundary(X_train, y_train)
   #classifier_comparison(X_train, y_train)
   #decision_slices(X_train, y_train)
   #decision_tree_stacking(X_train,y_train)
   #validation_error(X_train, y_train)
   #bagging_reggressor(X_train, y_train)


   #dict_models = batch_classify(X_train, Y_train, X_test, Y_test, no_classifiers = 8)
   #display_dict_models(dict_models)

# Rosenblatt Perceptron
  #perceptron(X_train, y_train)


# Bayesian Modelling and Inference


  #with pm.Model() as logistic_model:
   # pm.glm.GLM.from_formula('Identity ~ OMSSA:evalue + Exp m/z + OMSSA:pvalue', combined_psm, family=pm.glm.families.Binomial())
	#trace_logistic_model = pm.sample(2000, chains=1, tune=1000)

  #plot_traces(trace_logistic_model)
  #plot_distribution_of_factors(trace_logistic_model)

# Folder wide processing
	parser.add_argument(
	'Fasta',
	metavar='Fasta',
	help='Enter a valid file path to a fasta file!')
	args = parser.parse_args()
	decoy = generate_decoy(args.Fasta)

	parser = GooeyParser(description='LC-MS/MS Signal Processing')
	parser.add_argument('Folder',metavar='Fasta',help='Enter a valid file path to a fasta file!')
	args = parser.parse_args()
	folder = args.Folder

	parser = GooeyParser(description='LC-MS/MS Signal Processing')
	parser.add_argument('decoy',metavar='Fasta',help='Enter a valid file path to the decoy database')
	args = parser.parse_args()
	decoy_database = args.decoy

	parser = GooeyParser(description='LC-MS/MS Signal Processing')
	parser.add_argument('result_folder',metavar='Fasta',help='Enter a valid file path to the folder with PSM csvs')
	args = parser.parse_args()
	result_folder= args.result_folder

	combined_result_files = combine_results(result_folder)

	is_decoy = combined_result_files.loc[:, 'Is decoy']
	psm_count= combined_result_files.loc[:, 'PSM_counts']

	combined_result_files['FDR'] = calc_FDR(psm_count, is_decoy)

	combined_result_files['Category'] = np.sign(combined_result_files['FDR']-0.25)

	#using an RBF SVM to validate the protein probability scores
	X = feature_space.loc[:,['OMSSA:evalue', 'OMSSA:pvalue', 'PSM_counts', 'mass']].values
	y = feature_space.loc[:,['Category']].values
	y = y.astype(np.int64)
	y = np.ravel(y)
	combined_result_files['SVM Score'] = validate_results(X,y)



   # folder wide clustering, heatmaps and dendograms
   cluster_dendogram(Z,truncate_mode='lastp',p=12,leaf_rotation=90.,leaf_font_size=12.,show_contracted=True,annotate_above=10,  )
   plt.show()

   cluster_plots(X)
   cnf_matrix = confusion_matrix(y_test, y_pred)
   np.set_printoptions(precision=2)

   plot_confusion_matrix(cnf_matrix, classes=class_names,title='Confusion matrix, without normalization')
   plot_confusion_matrix(cnf_matrix, classes=class_names, normalize=True,title='Normalized confusion matrix')
   silhouette_analysis(X)
   affinity_propagation_clustering(X, 1)  # the integer value represents the desired number of centers


if __name__ == '__main__':
	main() 