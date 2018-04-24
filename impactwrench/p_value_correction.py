import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ursgal
from statsmodels.compat.python import range
from statsmodels.compat.collections import OrderedDict


def _ecdf(x):

	nobs = len(x)
	return np.arange(1,nobs+1)/float(nobs)

	multitest_methods_names = {'b': 'Bonferroni',
							'fdr_bh': 'FDR Benjamini-Hochberg',
							'fdr_by': 'FDR Benjamini-Yekutieli'
							}
	_alias_list = [['b', 'bonf', 'bonferroni'],
				['fdr_bh', 'fdr_i', 'fdr_p', 'fdri', 'fdrp'],
				['fdr_by', 'fdr_n', 'fdr_c', 'fdrn', 'fdrcorr']]

	multitest_alias = OrderedDict()
	for m in _alias_list:
		multitest_alias[m[0]] = m[0]
		for a in m[1:]:
			multitest_alias[a] = m[0]
	
def multipletests(pvals, alpha=0.05, method='b', is_sorted=False,
	  returnsorted=False):

	import gc
	pvals = np.asarray(pvals)
	alphaf = alpha

	if not is_sorted:
			sortind = np.argsort(pvals)
			pvals = np.take(pvals, sortind)

	ntests = len(pvals)
	alphacSidak = 1 - np.power((1. - alphaf), 1./ntests)
	alphacBonf = alphaf / float(ntests)
	if method.lower() in ['b', 'bonf', 'bonferroni']:
		reject = pvals <= alphacBonf
		pvals_corrected = pvals * float(ntests)

	elif method.lower() in ['fdr_bh', 'fdr_i', 'fdr_p', 'fdri', 'fdrp']:
			# delegate, call with sorted pvals
			reject, pvals_corrected = fdrcorrection(pvals, alpha=alpha,
													method='indep',
													is_sorted=True)
		
	elif method.lower() in ['fdr_by', 'fdr_n', 'fdr_c', 'fdrn', 'fdrcorr']:
			# delegate, call with sorted pvals
			reject, pvals_corrected = fdrcorrection(pvals, alpha=alpha,
													method='n',
													is_sorted=True)
			ii = np.arange(1, ntests + 1)
			q = (ntests + 1. - ii)/ii * pvals / (1. - pvals)
			pvals_corrected_raw = np.maximum.accumulate(q)

			pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
			del pvals_corrected_raw
			reject = pvals_corrected <= alpha

	else:
		raise ValueError('method not recognized')

	if not pvals_corrected is None: #not necessary anymore
		pvals_corrected[pvals_corrected>1] = 1
	if is_sorted or returnsorted:
		return reject, pvals_corrected, alphacSidak, alphacBonf
	else:
		pvals_corrected_ = np.empty_like(pvals_corrected)
		pvals_corrected_[sortind] = pvals_corrected
		del pvals_corrected
		reject_ = np.empty_like(reject)
		reject_[sortind] = reject
		return reject_, pvals_corrected_, alphacSidak, alphacBonf
				


def fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False):

	pvals = np.asarray(pvals)

	if not is_sorted:
		pvals_sortind = np.argsort(pvals)
		pvals_sorted = np.take(pvals, pvals_sortind)
	else:
		pvals_sorted = pvals  # alias

	if method in ['i', 'indep', 'p', 'poscorr']:
		ecdffactor = _ecdf(pvals_sorted)
	elif method in ['n', 'negcorr']:
		cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))   #corrected this
		ecdffactor = _ecdf(pvals_sorted) / cm

	else:
		raise ValueError('only indep and negcorr implemented')
	reject = pvals_sorted <= ecdffactor*alpha
	if reject.any():
		rejectmax = max(np.nonzero(reject)[0])
		reject[:rejectmax] = True

	pvals_corrected_raw = pvals_sorted / ecdffactor
	pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
	del pvals_corrected_raw
	pvals_corrected[pvals_corrected>1] = 1
	if not is_sorted:
		pvals_corrected_ = np.empty_like(pvals_corrected)
		pvals_corrected_[pvals_sortind] = pvals_corrected
		del pvals_corrected
		reject_ = np.empty_like(reject)
		reject_[pvals_sortind] = reject
		return reject_, pvals_corrected_
	else:
		return reject, pvals_corrected

