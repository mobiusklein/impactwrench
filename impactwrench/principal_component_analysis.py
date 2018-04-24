import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ursgal
from statsmodels.compat.python import range
from statsmodels.compat.collections import OrderedDict
from numpy import linalg as LA
from sklearn.decomposition import PCA
from sklearn.decomposition import IncrementalPCA
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.decomposition import KernelPCA

def local_linear_embedding(x):
    lle = LocallyLinearEmbedding(n_components = 0.95, n_neighbours = 10)
    return lle.fit_transform(x)


        
def calculate_mean_vector(x):
	
	 _extract_pvalues = x.iloc[:, ['OMSSA:pvalue']]
	 _mean_vector = _extract_pvalues.mean()
	 return  _mean_vector



def acquire_eigen_values(x):
	corrected_values_array = _x.as_matrix(columns = None)
	eigen_values = LA.eigvals(corrected_values_array)
	return LA.eigh(corrected_values_array)

def singular_value_decomposition(x):
    U, s, V = np.linalg.svd(x- x.mean(axis = 0))
    c1 = V.T[:, 0]
    c2 = V.T[:,1]
    return c1, c2

def minimal_number_components(x):
    pca = PCA(n_components= 0.95)
    return pca.fit_transform(x)
    

#increamental PCA
def incremental_pca(x,y):
    n_batches= 100
    inc_pca = IncrementalPCA(n_components = 2)
    for X_batch in np.array_split(x , n_batches):
        inc_pca.partial_fit(X_batch)
    return  inc_pca.transform(x)


def randomized_pca(x):
    rnd_pca = PCA(n_components=0.95,svd_solver= "randomized")
    return rnd_pca.fit_transform(x)
        


# selecting best parameters
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline

def best_parameters(train_set):
    clf= Pipeline([
        ("kpca", KernelPCA(n_components = 2)),
        ("log_reg", LogisticRegression())])
    param_grid = [{"kpca__gamma":np.linspace(0.03, 0.05, 10),"kpca__kernel":["rbf", "sigmoid"]}]
    grid_search = GridSearchCV(clf,param_grid, cv = 3)
    grid_search.fit(train_set)
    return grid_search.best_params_

def plot_pca(x):
    pca = KernelPCA(kernel="rbf", fit_inverse_transform=True, gamma=10)
    kpca = KernelPCA()
    X_kpca = kpca.fit_transform(x)
    X_back = kpca.inverse_transform(x_kpca)
    pca = PCA()
    X_pca = pca.fit_transform(x)

# Plot results

    plt.figure()
    plt.subplot(2, 2, 1, aspect='equal')
    plt.title("Original space")
    reds = y == 0
    blues = y == 1

    plt.scatter(x[reds, 0], x[reds, 1], c="red",
                s=20, edgecolor='k')
    plt.scatter(x[blues, 0], x[blues, 1], c="blue",
                s=20, edgecolor='k')
    plt.xlabel("$x_1$")
    plt.ylabel("$x_2$")

    X1, X2 = np.meshgrid(np.linspace(-1.5, 1.5, 50), np.linspace(-1.5, 1.5, 50))
    X_grid = np.array([np.ravel(x), np.ravel(x)]).T
    # projection on the first principal component (in the phi space)
    Z_grid = kpca.transform(X_grid)[:, 0].reshape(X1.shape)
    plt.contour(X1, X2, Z_grid, colors='grey', linewidths=1, origin='lower')

    plt.subplot(2, 2, 2, aspect='equal')
    plt.scatter(X_pca[reds, 0], X_pca[reds, 1], c="red",
                s=20, edgecolor='k')
    plt.scatter(X_pca[blues, 0], X_pca[blues, 1], c="blue",
                s=20, edgecolor='k')
    plt.title("Projection by PCA")
    plt.xlabel("1st principal component")
    plt.ylabel("2nd component")

    plt.subplot(2, 2, 3, aspect='equal')
    plt.scatter(X_kpca[reds, 0], X_kpca[reds, 1], c="red",
                s=20, edgecolor='k')
    plt.scatter(X_kpca[blues, 0], X_kpca[blues, 1], c="blue",
                s=20, edgecolor='k')
    plt.title("Projection by KPCA")
    plt.xlabel("1st principal component in space induced by $\phi$")
    plt.ylabel("2nd component")

    plt.subplot(2, 2, 4, aspect='equal')
    plt.scatter(X_back[reds, 0], X_back[reds, 1], c="red",
                s=20, edgecolor='k')
    plt.scatter(X_back[blues, 0], X_back[blues, 1], c="blue",
                s=20, edgecolor='k')
    plt.title("Original space after inverse transform")
    plt.xlabel("$x_1$")
    plt.ylabel("$x_2$")

    #plt.subplots_adjust(0.02, 0.10, 0.98, 0.94, 0.04, 0.35)

    plt.show()


