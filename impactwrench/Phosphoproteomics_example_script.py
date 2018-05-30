# Phosphoproteomics MaxQuant Post Processing Script
import pandas as pd
import numpy as np
from padua import *
import csv
import matplotlib.pyplot as plt
from postQuant import *


df = read_proteinGroups('proteinGroups.txt')

phos = read_proteinGroups('Phospho(STY) Sites.txt')

# preliminary file processing

df = remove_potential_contaminants(df)

df = remove_only_identified_by_site(df)

df = remove_reverse(df)

df = localization_probability_filter(df)

df = delta_score_threshold(df)

df = localization_cut_off(df)

df = log_2_transformation(df.loc[:,'Intensity'])

df = guassian_missing_value_imputation(df, width =0.3, downshift = -1.8, prefix = None)

df = anova(df)


#missing value imputation using k-nearest neighbour instead of guassian imputation

df = df.values
clf = fit_knn(df, column)
df = transform (df, column, clf)

# ratio between two groups protein R-site multiplicity at different time points
phos_ratio = phospho_ratio(phos)



#scatter histogram
sh = ScatterHist(df)
sh.plot(kargs_scatter={'c':'r', 's':30, 'alpha':.3},kargs_histy={'color':'g'})

h = hist2d.Hist2D(df)
res = h.plot( contour=True)

c = corrplot.Corrplot(df)
c.plot(colorbar=False, method='square', shrink=.9 ,rotation=45)