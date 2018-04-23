import sys
print (sys.argv)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import pymzml
from scipy import signal, misc
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns
from mlxtend.plotting import plot_decision_regions
from sklearn import datasets
from sklearn.svm import SVR
import matplotlib.pyplot as plt
from mlxtend.plotting import plot_decision_regions
import matplotlib.gridspec as gridspec
import itertools
from scipy import signal, misc
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB 
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
def visualization(x, y,z):
	fig = plt.figure()
	plt.subplot(211)
	plt.plot(x, y)
	plt.title('Extracted Ion Chromatogram')
	plt.ylabel('Signal Intensity')
	plt.xlabel('Mass/Charge')

	plt.subplot(212)
	plt.plot(z)
	plt.title('Chromatogram')
	plt.ylabel('Signal Intensity')
	plt.xlabel('Retention Time')
	plt.show()

def signal_resampling(x):
	plt.figure(1)
	plt.subplot(211)
	plt.plot(signal.resample(run['TIC'].peaks, 1000))
	plt.title('Resampled TIC')
	plt.ylabel('Signal Intensity')
	plt.xlabel('Retention Time')

	plt.subplot(212)
	plt.plot(run['TIC'].peaks)
	plt.title('Original TIC')
	plt.ylabel('Signal Intensity')
	plt.xlabel('Retention Time')
	plt.show()

import numpy as np
from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve

def WhittakerSmooth(x,w,lambda_,differences=1):
	X=np.matrix(x)
	m=X.size
	i=np.arange(0,m)
	E=eye(m,format='csc')
	D=E[1:]-E[:-1] 
	W=diags(w,0,shape=(m,m))
	A=csc_matrix(W+(lambda_*D.T*D))
	B=csc_matrix(W*X.T)
	background=spsolve(A,B)
	return np.array(background)

def airPLS(x, lambda_=100, porder=1, itermax=15):
	m=x.shape[0]
	w=np.ones(m)
	for i in range(1,itermax+1):
		z=WhittakerSmooth(x,w,lambda_, porder)
		d=x-z
		dssn=np.abs(d[d<0].sum())
		if(dssn<0.001*(abs(x)).sum() or i==itermax):
			if(i==itermax): print ('WARNING max iteration reached!')
			break
		w[d>=0]=0 # d>0 means that this point is part of a peak, so its weight is set to 0 in order to ignore it
		w[d<0]=np.exp(i*np.abs(d[d<0])/dssn)
		w[0]=np.exp(i*(d[d<0]).max()/dssn) 
		w[-1]=w[0]
	return z   

def plot_baseline_correction(x, y):
	
	#preliminary processing
	_x = np.array(x)
	_y =np.array(y)
	_signal = _y
	
	#linear baseline
	_baseline_1= 5e-4*_x+ 0.2
	_baseline_2=0.2*np.sin(np.pi*_x/_x.max())
	_noise=np.random.normal(0, 1, 100)
	_signal_1 = _y 

	#correction of baselines
	_c1 =_signal_1-airPLS(_signal_1)
	
	fig, ax1 = plt.subplots()
	ax1.set_title('Adaptive Iterative Reweighted Partial Least Squares Baseline Correction')
	ax1.plot(_x,_signal_1, 'k')
	ax1.plot(_x,_c1, '-r')
	ax1.set_xlabel('mass/charge')
	ax1.set_ylabel('Intensity')
	plt.legend([_signal_1, _c1], ["Original", "Corrected"])
	plt.show()

def smooth(y, box_pts):
	box = np.ones(box_pts)/box_pts
	return  np.convolve(y, box_pts, mode='same')

def plot_smoothing(x, y):
	fig = plt.figure()
	plt.plot(x, y,'b')
	plt.plot(x, smooth(y,3), 'r-', lw=2)
	plt.plot(x ,smooth(y,3), 'g-', lw=2)
	plt.title('Smoothing by Convolution')
	plt.ylabel('Signal Intensity')
	plt.xlabel('Mass/Charge')
	plt.show()

def plot_resampling(z):
	f = signal.resample(z,1000)
	plt.figure(1)
	plt.subplot(211)
	plt.plot(f)
	plt.title('Resampled TIC')
	plt.ylabel('Signal Intensity')
	plt.xlabel('Retention Time')

	plt.subplot(212)
	plt.plot(run['TIC'].peaks)
	plt.title('Original TIC')
	plt.ylabel('Signal Intensity')
	plt.xlabel('Retention Time')
	plt.show()

#normalization
import matplotlib.colors as colors
def plot_normalization(x,y):
	_norm = colors.Normalize(0, 100, clip = False)
	for pix in np.nditer(y):
		_val = (norm(pix))
	
	fig = plt.figure()
	img=norm(y)
	plt.plot(x, img, 'b' )
	plt.show()

def plot_scatter(x, y):
	def add_titlebox(ax, text):
		ax.text(0.55, 0.8,text, horizontalalignment= 'center', transform = ax.transAxes, bbox= dict(facecolor= 'white', alpha = 0.6), fontsize= 12.5)
		return ax
	gridsize = (3, 2)
		
	fig = plt.figure(figsize=(12, 8))
	ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=2, rowspan=2)
	ax2 = plt.subplot2grid(gridsize, (2, 0))
	ax3 = plt.subplot2grid(gridsize, (2, 1))

	
	sctr =ax1.scatter(x , y , c = y, cmap = 'RdYlGn')
	plt.colorbar(sctr, ax=ax1, format = None)
	ax1.set_yscale('log')
	ax1.set_title('Scatter Plot M/z vs Intensity')
	ax2.hist(x, bins='auto')
	ax3.hist(y, bins='auto', log=True)

	add_titlebox(ax2, 'Histogram: mass to charge')
	add_titlebox(ax3, 'Histogram: signal intensity(log scl.)')	
	plt.show()



	# Loading some example data
def one_dimensional_decision_boundary(X_train, y_train):
		svr = SVR(C=0.5, kernel='linear')
		svr.fit(X_train, y_train)
		plot_decision_regions(X_train, y_train, clf=svr, legend=2)
			# Adding axes annotations
		plt.xlabel(' Intensity [units]')
		plt.title('SVR')
		plt.show()

from sklearn.svm import SVC

def multiple_classifiers_decision_boundary(x, y):
	# Initializing Classifiers
	clf1 = LogisticRegression(random_state=1)
	clf2 = RandomForestClassifier(random_state=1)
	clf3 = GaussianNB()
	clf4 = SVC()

	gs = gridspec.GridSpec(2, 2)

	fig = plt.figure(figsize=(10,8))

	labels = ['Logistic Regression', 'Random Forest', 'Naive Bayes', 'SVM']
	for clf, lab, grd in zip([clf1, clf2, clf3, clf4],
						 labels,
						 itertools.product([0, 1], repeat=2)):

		clf.fit(x, y)
		ax = plt.subplot(gs[grd[0], grd[1]])
		fig = plot_decision_regions(X=x, y=y, clf=clf, legend=2)
		plt.title(lab)

		plt.show()
		
def remove_nan(x):
    return x.dropna(how = 'any')

