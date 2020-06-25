'''
Contains  PCA class.
'''

import scipy as sp
from scipy import linalg
import futils

class PCA(object):
    def __init__(self, a, center=False, keep=None):
        """
        Copied from matplotlib 1.1.
        
        compute the SVD of a and store data for PCA.  Use project to
        project the data onto a reduced set of dimensions. 

        Inputs:

          *a*: a numobservations x numdims matrix

          If center=True, data will be centered using mean and standard deviation.

        Attrs:

          *mu* : a numdims array of means of a

          *sigma* : a numdims array of standard deviation of a

          *fracs* : the proportion of variance of each of the principal components

          *eofs* : Empirical Orthogonal Functions.

          *pcs* : Principal Components


        """
        n, m = a.shape

        self.mu = a.mean(axis=0)
        self.sigma = a.std(axis=0)

        if center:
            a = self.center(a)

        U, s, Vh = linalg.svd(a, full_matrices=False)

        self.pcs=(U*s[sp.newaxis])[:,:keep]

        evals = s**2/n
        self.evals=evals
        sumvar=(self.sigma**2).sum()
        self.sumvar=sumvar
        self.fracs = (evals/sumvar)[:keep]
        self.eofs = Vh.T[:,:keep]


    def center(self, x):
        'center the data using the mean and sigma from training set a'

        return sp.where(self.sigma==0,(x-self.mu),(x-self.mu)/self.sigma)

    def varimax(self, nr=50):
        '''
        Does rotation of EOFs using varimax criteria.
        
        nr: how many EOFs to keep, default - keep all EOFs 
        '''
        # Need to renormalize eofs and pcs

        N=self.eofs.shape[1] # number of EOFs

        eofs=self.eofs*sp.sqrt(self.evals[sp.newaxis])
        pcs=self.pcs/sp.sqrt(self.evals[sp.newaxis])
        
        eofs, pcs=futils.varimax(eofs, pcs, min(nr, N))

        # Normalize EOFs
        nn=sp.array([sp.linalg.norm(xx) for xx in eofs.transpose()])
        eofs/=nn[sp.newaxis]
        pcs*=nn[sp.newaxis]

        # Sort by PC variance
        lam=sp.array([xx.var() for xx in pcs.transpose()])
        ind=sp.argsort(lam)[::-1]
        
        self.evals=lam[ind]
        self.fracs=lam[ind]/self.sumvar
        self.pcs=pcs[:,ind]
        self.eofs=eofs[:,ind]

def pca(field, keep=None, center=False, weight=True):
    '''
    Performss principal component analysis on data field, and stores
    a PCA object. Please, remove climatology, detrend data etc before
    calling this method. 

    If center=True, PCA object will center data using mean and standard deviation.
    If weight=True, multiply data by area weights.
        '''

    nt,km,jm,im=field.data.shape

    # multiply data by area factor, reshape, return matrix
    if weight:
        factor=sp.cos(sp.deg2rad(field.grid['lat']))
        factor[factor<0.]=0.
        factor=sp.sqrt(factor)
    else:
        factor=sp.ones(field.grid['lat'].shape)

    mask=sp.ma.getmaskarray(field.data).copy()
    field.data[mask]=0.0
    field.data*=factor[sp.newaxis,sp.newaxis]
    X=field.data.reshape((nt,km*jm*im)).view(sp.ndarray) 

    field._pc=pca.PCA(X, center=center, keep=keep)

    field.data/=factor[sp.newaxis,sp.newaxis]
    field.data[mask]=field.data.fill_value
    field.data.mask=mask
