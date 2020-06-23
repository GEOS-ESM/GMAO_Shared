'''
Different routines used by field, dset, grid etc.

Note: g2g functionality was moved to gridtile.py due to dependence on f2py.
'''

import scipy as sp

def scalar2slice(*args):
    val=[]
    for i in args:
        if sp.isscalar(i): i=slice(i,i+1)
        val.append(i)
    return tuple(val)


def find_nearest(array,value):
    idx=sp.absolute(array-value).argmin()
    return array[idx]

def find_nearest_index(array,value):
    idx=sp.absolute(array-value).argmin()
    return idx

def gint(y, x, axis=-1):
    """
    Integrate along the given axis using the composite trapezoidal rule. Copied from numpy.
    Should work with masked arrays
    
    Integrate `y` (`x`) along given axis.
    
    Parameters
    ----------
    y : array
    Input array to integrate.
    x : array
    axis : int, optional
    Specify the axis.
    """
    
    if x.ndim == 1:
        d = sp.diff(x)
        # reshape to correct shape
        shape = [1]*y.ndim
        shape[axis] = d.shape[0]
        d = d.reshape(shape)
    else:
        d = sp.diff(x, axis=axis)
    nd = y.ndim
    slice1 = [slice(None)]*nd
    slice2 = [slice(None)]*nd
    slice1[axis] = slice(1,None)
    slice2[axis] = slice(None,-1)

    return sp.ma.sum(d*(y[slice1]+y[slice2])/2.0,axis)
