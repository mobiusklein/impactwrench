import pandas as pd
import numpy as np
from padua import *
import csv
import matplotlib.pyplot as plt
import sys
from sklearn import neighbors
from biokit.viz import corrplot
from biokit.viz import hist2d
from biokit import viz
from biokit import ScatterHist
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm


def read_proteinGroups(p, header = 0):
	pG = pd.read_csv(f, delimiter='\t', header=header, index_col=index_col, **kwargs, low_memory = False)
    return pG 


def remove_potential_contaminants(df):
    return df.[df.Potential contaminant != "+"]

def remove_only_identified_by_site(df):
    return df.[df.Only identified by site !="+"]

def remove_reverse(df):
    return df.[df.Reverse !="+"]


def localization_probability_filter(df, threshold):
    return df[~(df['Localization Probability'] >= threshold)]  

def delta_score_threshold(df, threshold):
    return df[~(df['Delta Score'] >= threshold)] 


def localization_cut_off(df, threshold):
    return df[~(df['Localization Score'] >= threshold)] 

def log_2_transformation(x):
    return np.log2(x)

def guassian_missing_value_imputation(df, width = 0.3, downshift = -1.8, prefix = None):
	f = df.copy()

    imputed = df.isnull()  

    if prefix:
        mask = np.array([l.startswith(prefix) for l in df.columns.values])
        mycols = np.arange(0, df.shape[1])[mask]
    else:
        mycols = np.arange(0, df.shape[1])


    if type(width) is not list:
        width = [width] * len(mycols)

    elif len(mycols) != len(width):
        raise ValueError("Length of iterable 'width' does not match # of columns")

    if type(downshift) is not list:
        downshift = [downshift] * len(mycols)

    elif len(mycols) != len(downshift):
        raise ValueError("Length of iterable 'downshift' does not match # of columns")

    for i in mycols:
        data = df.iloc[:, i]
        mask = data.isnull().values
        mean = data.mean(axis=0)
        stddev = data.std(axis=0)

        m = mean + downshift[i]*stddev
        s = stddev*width[i]

        values = np.random.normal(loc=m, scale=s, size=df.shape[0])

        df.iloc[mask, i] = values[mask]

    return df, imputed

def anova(df):
    Intensity = df.loc[:,'Intensity']
    Peptides  = df.loc[:,'Peptides']

    model = ols('Intensity ~ C(Peptides)', df).fit()
    anovaResults = anova_lm(model)
    return anovaResults


def fit_knn(x, column):

    clf = neighbors.KNeighborsRegressor(n_neighbors=10)
    missing_idxes = np.where(pd.isnull(x[:, column]))[0]
    if len(missing_idxes) == 0:
        return None
    X_copy = np.delete(X, missing_idxes, 0)
    X_train = np.delete(X_copy, column, 1)
    col_mean = None
    if not is_categorical:
        col_mean = np.nanmean(X, 0)
    else:
        col_mean = np.nanmedian(X, 0)
    for col_id in range(0, len(col_mean) - 1):
        col_missing_idxes = np.where(np.isnan(X_train[:, col_id]))[0]
        if len(col_missing_idxes) == 0:
            continue
        else:
            X_train[col_missing_idxes, col_id] = col_mean[col_id]
    y_train = X_copy[:, column]
    clf.fit(X_train, y_train)
    return clf

def transform(x, column):
    clf = neighbors.KNeighborsRegressor(n_neighbors=10)
    missing_idxes = np.where(np.isnan(x[:, column]))[0]
    X_test = x[missing_idxes, :]
    X_test = np.delete(X_test, column, 1)
    col_mean = np.nanmedian(x, 0)
    for col_id in range(0, len(col_mean) - 1):
        col_missing_idxes = np.where(np.isnan(X_test[:, col_id]))[0]
        # if no missing values for current column
        if len(col_missing_idxes) == 0:
            continue
        else:
            X_test[col_missing_idxes, col_id] = col_mean[col_id]
        # predict missing values
        y_test = clf.predict(X_test)
        x[missing_idxes, column] = y_test
        return x