import numpy as np
from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve


#whittaker smoothing function

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


    #using adaptive iterative reweighted partial least squares
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

    import matplotlib.colors as colors
    import numpy as np
    import pymzml

    #normalization
def normalize(y):
    norm = colors.Normalize(0, 100)
    for pix in np.nditer(y, op_flags = ['readwrite']):
        return norm(pix)

    from scipy import signal, misc
    import numpy as np
    from scipy.interpolate import interp1d
    import pymzml

#resampling
def resample(y):
        return signal.resample(y, 1000)

#signal smoothing
def smooth(signal_intensity, box_pts):
    box = np.ones(box_pts)/box_pts
    return  np.convolve(signal_intensity, box, mode='same')

#single object extraction
def load_file(filename):
    for f in filename:
        run = pymzml.run.Reader(filename, MS1_Precision = 5e-6, MSn_Precision = 20e-6)
    for spectrum in run:
        if spectrum['ms level'] ==1:
            y = spectrum.i
            x = spectrum.mz

def mass_to_charge():
    x = spectrum.mz
    return(x)

def signal_intensity():
    y = spectrum.i
    return(y)

def total_ion_count():
    return run['TIC'].peaks

from dtw import dtw

# custom distance distances

def distance(x, y):
    my_custom_norm =  (x *x) +(y*y)
    dist, cost, acc, path = dtw(x, y, dist = my_custom_norm)
    imshow(acc.T, origin='lower', cmap=cm.gray, interpolation='nearest')
    plot(path[0], path[1], 'w')
    xlim((-0.5, acc.shape[0]-0.5))
    ylim((-0.5, acc.shape[1]-0.5))


def dtw_2(x, y):
    dist,cost, acc, path = dtw(x, y, dist = lambda x, y: norm(x-y, ord = 1))
    imshow(acc.T, origin='lower', cmap=cm.gray, interpolation='nearest')
    plot(path[0], path[1], 'w')
    xlim((-0.5, acc.shape[0]-0.5))
    ylim((-0.5, acc.shape[1]-0.5))
