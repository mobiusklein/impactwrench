import sys
print (sys.argv)
import os
import click
import os
import glob
import ursgal
import pymzml
from os import listdir
import csv
import pandas as pd
import numpy as np
import re
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.svm import SVC
from sklearn.cluster import KMeans
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
import itertools
import time
import warnings
from __future__ import print_function
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.cm as cm
from itertools import cycle

#generating target_decoy database

def generate_decoy(fasta_file_path):

	params = {
	    'enzyme': 'trypsin',
	    'decoy_generation_mode': 'reverse_protein',}

	uc = ursgal.UController(
	    params=params)

	new_target_decoy_db_name = uc.execute_misc_engine(
	    input_file=fasta_file_path,
	    engine='generate_target_decoy_1_0_0',
	    output_file_name='Etanercept_decoy.fasta',)
	return output_file_name


   #initialization of search engine and validation engines
    mzML_files = []
    for mzml in glob.glob(os.path.join('{0}'.format(folder), '*.mzML')):
    mzML_files.append(mzml)

    search_engines = ['omssa_2_1_9', 'myrimatch_2_2_140 ',]
    all_mods = ['C,fix,any,Carbamidomethyl','M,opt,any,Oxidation','N,opt, any, Deamidated', 'S,opt, any,Phospho','N, opt,anyHexNAc','M,opt,any,Oxidation','*,opt,Prot-N-term,Acetyl', ]

    params = {
        'database': decoy_database,
        'modifications': all_mods,
        'csv_filter_rules': [
            ['Is decoy', 'equals', 'false'],
            ['PEP', 'lte', 0.01],
        ]
    }

def search(x):
  '''This function uses the architucture from ursgal and open source search engines'''
	  uc = ursgal.UController(params = params,)

	  for search_engine in search_engines:
		unified_csvs = []
		for mzML_file in mzML_files:
			unified_search_engine_result_file = uc.search(
				input_file = mzML_file,
				engine = search_engine,)
	return	unified_csvs.append(unified_search_engine_result_file)

def combine_result_files(result_folder):
	csv_files = glob.glob(os.path.join(result_folder, '*csv'))

	combined = []

	''' Depending on the convention used to name files in the folder, the code below can be easily tailored to 
	cater to specific names'''

	Eta8 = []
	Eta9 = []
	for filename in csv_files:
    filename = str(filename)
    if filename.endswith('Eta8.csv'):
        Eta8.append(filename)
        Eta8 = pd.concat(pd.read_csv(f) for f in Eta8)
        Eta8['Identity'] = [1] * len(Eta8)
    elif filename.endswith('Eta9.csv'):
        Eta9.append(filename)
        Eta9= pd.concat(pd.read_csv(f) for f in Eta9)
        Eta9['Identity'] = [2] * len(Eta9) 
    return pd.concat([Eta8, Eta9])

def filter_results(x):
    '''This function adds two more columns to the dataframe for downstream processing'''
    x['PSM_counts'] = x.groupby('Sequence')['Sequence'].transform('count')
    x['Is decoy'] = x['Is decoy'].astype(int)
    x['mass']= ((df['Exp m/z'] * df['Charge'])-(df['Charge']-1)*1 )
    s = pd.Series(x['Sequence'])
	x['pep_len'] = s.str.len()
    columns = ['mass', 'pep_len', 'Charge', 'Retention Time (s)', 'Identity', 'Sequence', 'OMSSA:evalue', 'OMSSA:pvalue', 'Is decoy','PSM_counts','Accuracy(ppm)']
  	return pd.DataFrame(x, columns = columns)

def calc_FDR(PSM_count, false_positives):
    true_positives  = PSM_count - (2 * false_positives)
    if true_positives <= 0:  # prevent FDR above 1. Not sure if this is done?
        return 1.0
    FDR = false_positives / (false_positives + true_positives)
    return FDR

def validate_results(x, y):
 
	rbf_kernel_svm_clf = Pipeline((
	    ("scaler", StandardScaler()),
	    ("svm_clf", SVC(kernel = "rbf", gamma = 5, C = 0.001)) # gamma is the regularization hyperparameter))
	rbf_kernel_svm_clf.fit(X,y)))
	return rbf_kernel_svm_clf.score(X,y)  




def naive_bayes(prob_list):
    multiplied_probs = \
        functools.reduce(operator.mul, prob_list, 1)
    # multiplying the opposite of all probabilities: (1-a)*(1-b)*(1-c)
    multiplied_opposite_probs = \
        functools.reduce(operator.mul, (1-p for p in prob_list), 1)
    return multiplied_probs / (multiplied_probs + multiplied_opposite_probs)


 def cluster_dendogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata

def cluster_plots(X):
   ''' This function allows one to plots numerous clustering techniques'''


   #heirarchical clustering
	Z = linkage(X, 'ward')
	plt.figure(figsize=(25, 10))
	plt.title('Hierarchical Clustering Dendrogram')
	plt.xlabel('sample index')
	plt.ylabel('distance')dendrogram(Z,
	    leaf_rotation=90.,  # rotates the x axis labels
	    leaf_font_size=8.,)
	return plt.show()
    

    #3D clustering 
    np.random.seed(5)



	estimators = [('k_means_8', KMeans(n_clusters=8)),
	              ('k_means_3', KMeans(n_clusters=3)),
	              ('k_means_iris_bad_init', KMeans(n_clusters=3, n_init=1,
	                                               init='random'))]

	fignum = 1
	titles = ['8 clusters', '3 clusters', '3 clusters, bad initialization']
	for name, est in estimators:
	    fig = plt.figure(fignum, figsize=(4, 3))
	    ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
	    est.fit(X)
	    labels = est.labels_

	    ax.scatter(X[:, 1], X[:, 0], X[:, 2],
	               c=labels.astype(np.float), edgecolor='k')

    ax.w_xaxis.set_ticklabels([])
    ax.w_yaxis.set_ticklabels([])
    ax.w_zaxis.set_ticklabels([])
    ax.set_xlabel('mass')
    ax.set_ylabel('pvalue')
    ax.set_zlabel('evalue')
    ax.set_title(titles[fignum - 1])
    ax.dist = 12
    fignum = fignum + 1

	# Plot the ground truth
	fig = plt.figure(fignum, figsize=(4, 3))
	ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

	for name, label in [('Batch 1', 0),
	                    ('Batch 2', 1),
	                    ('Batch 3', 2)]:
	    ax.text3D(X[y == label, 1].mean(),
	              X[y == label, 0].mean(),
	              X[y == label, 2].mean() + 2, name,
	              horizontalalignment='center',
	              bbox=dict(alpha=.2, edgecolor='w', facecolor='w'))
	# Reorder the labels to have colors matching the cluster results
	y = np.choose(y, [1, 2, 0]).astype(np.float)
	ax.scatter(X[:, 1], X[:, 0], X[:, 2], c=y, edgecolor='k')

	ax.w_xaxis.set_ticklabels([])
	ax.w_yaxis.set_ticklabels([])
	ax.w_zaxis.set_ticklabels([])
	ax.set_xlabel('pvalue ')
	ax.set_ylabel('evalue')
	ax.set_zlabel('mass ')
	ax.set_title('Ground Truth')
	ax.dist = 12

	return fig.show()

def plot_confusion_matrix(cm, classes,normalize=False,title='Confusion matrix',cmap=plt.cm.Blues):
	class_names = ['Eta8', 'Eta9']
	if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    return plt.show()

def silhouette_analysis(X):
	range_n_clusters = [2, 3, 4, 5, 6]

for n_clusters in range_n_clusters:

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)
    ax1.set_xlim([-0.1, 1])
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])


    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)
    silhouette_avg = silhouette_score(X, cluster_labels)
    print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg)

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(n_clusters):

        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        y_lower = y_upper + 10  # 10 for the 0 samples

    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    colors = cm.spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                c=colors, edgecolor='k')
    centers = clusterer.cluster_centers_
    ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
                c="white", alpha=1, s=200, edgecolor='k')

    for i, c in enumerate(centers):
        ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
                    s=50, edgecolor='k')

    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("Feature space for the 1st feature")
    ax2.set_ylabel("Feature space for the 2nd feature")

    plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                  "with n_clusters = %d" % n_clusters),
                 fontsize=14, fontweight='bold')

    return plt.show()

def affinity_propagation_clustering(X, n):

	centers = [[n, n], [-n, -n], [n, -n]]
	labels= ['sample 1', 'sample 2']

	#computer clusters

	af = AffinityPropagation(preference=-50).fit(X)
	cluster_centers_indices = af.cluster_centers_indices_
	labels = af.labels_

	n_clusters_ = len(cluster_centers_indices)

	return ("Silhouette Coefficient: %0.3f"
	      % metrics.silhouette_score(X, labels, metric='sqeuclidean'))

	plt.close('all')
	plt.figure(1)
	plt.clf()

	colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
	for k, col in zip(range(n_clusters_), colors):
	    class_members = labels == k
	    cluster_center = X[cluster_centers_indices[k]]
	    plt.plot(X[class_members, 0], X[class_members, 1], col + '.')
	    plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
	             markeredgecolor='k', markersize=14)
	    for x in X[class_members]:
	        plt.plot([cluster_center[0], x[0]], [cluster_center[1], x[1]], col)

	plt.title('Estimated number of clusters: %d' % n_clusters_)
	return plt.show()










