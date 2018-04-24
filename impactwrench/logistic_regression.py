import pymc3 as pm
import numpy as np
import pandas as pd
import pandas as pd
import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt
import seaborn
import warnings
warnings.filterwarnings('ignore')
from collections import OrderedDict
from time import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import fmin_powell
from scipy import integrate

import theano as thno
import theano.tensor as T

def run_models(combined_psm, upper_order=5):
 

    models, traces = OrderedDict(), OrderedDict()

    for k in range(1,upper_order+1):

        nm = 'k{}'.format(k)
        fml = create_poly_modelspec(k)

        with pm.Model() as models[nm]:

            print('\nRunning: {}'.format(nm))
            pm.glm.GLM.from_formula(fml, combined_psm, family=pm.glm.families.Normal())

            traces[nm] = pm.sample(2000, chains=1, init=None, tune=1000)

    return models, traces

def plot_traces(traces, retain=1000):

    ax = pm.traceplot(traces[-retain:], figsize=(12,len(traces.varnames)*1.5),
        lines={k: v['mean'] for k, v in pm.summary(traces[-retain:]).iterrows()})

    for i, mn in enumerate(pm.summary(traces[-retain:])['mean']):
        ax[i,0].annotate('{:.2f}'.format(mn), xy=(mn,0), xycoords='data'
                    ,xytext=(5,10), textcoords='offset points', rotation=90
                    ,va='bottom', fontsize='large', color='#AA0022')

def create_poly_modelspec(k=1):
    return ('Identity ~ OMSSA:evalue + Exp m/z OMSSA:pvalue')

def define_model(combined_psm):
    with pm.Model() as logistic_model:
        pm.glm.GLM.from_formula('Identity ~ OMSSA:evalue + Exp m/z + OMSSA:pvalue', combined_psm, family=pm.glm.families.Binomial())
        trace_logistic_model = pm.sample(2000, chains=1, tune=1000)

def plot_distribution_of_factors(x):
	trace = x
	seaborn.jointplot(trace['OMSSA:pvalue'], trace['Exp m/z'], kind="hex", color="#4CB391")
	plt.xlabel("pvalue")
	plt.ylabel("Exp m/z")
	plt.show()
