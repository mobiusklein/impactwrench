# impactwrench

ImpactWrench is a Python CLI application for LC-MS and LC-MS/MS data analysis.
The pacakge allows one to conduct LC-MS Chemometric signal processing followed by application of embedded statistical methods and 
machine learning techniques for  classification, clustering and package recognition.


## LC MS Chemometric Signal Processing
In this module, pymzml is used as an mzML file parser and each mzML file can be processed using parameters that one desires and re-exported as an mzML file that is then subjected to peptide spectrum matching using the OMSSA search engine.

Chemometric Signal Processing should ideally proceed in the following fashion:
- Signal Resampling using architecture from scipy
- Signal Smoothing via Whittaker smoothing function/ Convolution
- Baseline Correction via the adaptive iterative reweighted partial least squares regression algorithm
-Signal Normalization
- Chromatographic Alignment via Dynamic Time Warping

## LC-MS Data Processing
This module is ideal for metabolomics as it allows one to exclusively carry out statistical analysis of LC-MS signals whilst disregarding LC-MS/MS signals. Ideally, one should carry out dimensionality reduction via included functions such as kernel_pca, rbf_pca before subjecting the resultant dataframes to machine learning methods.

## LC-MS/MS Data Processing
This module uses architecture from ursgal with parameters optimized for high resolution mass spectrometers (e.g Orbitrap Fusion Tribid Mass Spectrometer) and OMSSA search engine. In the case that one has additional open source or commercial search engines, this CLI has built in functionality to combine and correct for the FDR. For FDR calculation, Bonferroni, Benjamini-Hochberg and Benjamini-Yekutieli correction procedures are employed. Machine Learning methods available include  Support Vector Machines, Logistic Regression, Nearest Neighbours, Naive Bayes, Gradient Boosting Classifier and Random Forest. It should be known that different combinations of machine learning methods exist in teh script and can be easily changed to the user's liking.

Currently, with respect to neural networks, the package uses architecture from Theano and pymc3 to construct a simple Rosenblatt perceptron for classification and uses Bayes Inference to conduct Logistic Regression.d



## Basic setup

Install the requirements:
```
$ pip install -r requirements.txt
```

Run the application:
```
$ python -m impactwrench --help
```

To run the tests:
```
    $ pytest
```

To run the script:
```
$ python impactwrench.py
```
