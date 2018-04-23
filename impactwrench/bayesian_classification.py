import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import theano.tensor as T
import pymc3 as pm
import theano.tensor as tt
import theano
from sklearn.preprocessing import scale
from sklearn.cross_validation import train_test_split


X_train, X_test, Y_train, Y_test = train_test_split(X, y,train_size= 30,random_state=4)

plt.scatter(X[:, 0], X[:,1])
plt.show()

X_unc = np.ones_like(X) * 1.
X_unc[:, 0] = .5

with pm.Model() as logistic_regression:
	inputs = pm.Normal('inputs', mu=X, sd=X_unc, shape=X.shape)
	coeffs = pm.Normal('beta', mu=0, sd=1, shape=2)
	
	linear = pm.math.dot(inputs, coeffs)
	prob = pm.math.sigmoid(linear)
	obs = pm.Bernoulli('obs', p=prob, observed=y)
	
	# Run ADVI which returns posterior means, standard deviations, and the evidence lower bound (ELBO)
	mean_field = pm.fit(method='advi', n=50000,
						callbacks=[pm.callbacks.CheckParametersConvergence(diff='relative')])
	
	out_sample = (mean_field.sample_node(prob, size=100, more_replacements={inputs: grid_2d}) > .5)