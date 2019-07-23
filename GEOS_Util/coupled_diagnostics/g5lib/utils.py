'''
Different routines used by field, dset, grid etc.
'''

import scipy as sp
import futils

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


def g2g(datain,tiledata,hdimsout):
    '''
    Interpolates field from one horizontal grid onto another using exchange 
    grid (tiles).
    
    datain - input array (at least 2d)
    tiledata - array of records of size Nt (Nt - number of tiles).
    tiledata['iin'], tiledata['jin'], tiledata['iout'], 
    tiledata['jout'], tiledata['frac']  - 
    arrays of size Nt of inexes on the input grid, indexes on the output grid, 
    corresponding to a tile N, and fraction of tile on the output grid
    hdimsout - (imout,jmout), dimensions of output grid
    '''
    shin=datain.shape
    shout=datain.shape[:-2]+hdimsout
            
    return sp.ma.masked_values(
        futils.g2g(datain.filled().reshape((-1,)+shin[-2:]),
                   tiledata['iin'].astype(int),
                   tiledata['jin'].astype(int),
                   tiledata['iout'].astype(int),
                   tiledata['jout'].astype(int),
                   tiledata['frac'],
                   datain.fill_value,
                   shout[-1],
                   shout[-2]),datain.fill_value).reshape(shout)
