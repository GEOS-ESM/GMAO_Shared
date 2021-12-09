'''
This module provides different diagnostic utilities.
'''

import numpy as np

def integral(da, dim, weight):
    '''
    Computes weighted integral of DataArray along axis.
    Note: weight should have the same land mask as da.
    '''
    
    return (da*weight).sum(dim)

def average(da, dim, weight):
    '''
    Computes weighted average of DataArray along axis.
    Note: weight should have the same land mask as da.
    '''
    
    return (da*weight).sum(dim)/weight.sum(dim)

def standard_deviation(da, dim, weight):
    '''
    Computes weighted average of DataArray along axis.
    Note: weight should have the same land mask as da.
    '''
    mean=average(da, dim, weight)
    return np.sqrt(average((da-mean)**2,dim,weight))

def print_stat(da, dim, weight):
    '''
    Returns a string with da mean and standard deviation.
    '''
    mean=average(da, dim, weight)
    std=standard_deviation(da, dim, weight)
    return f'mean: {mean.values:.2f}\nstd: {std.values:.2f}'

def shift_lon(da, lon0, lon_name='lon'):
    '''
    Shifts the reference longitude of field da to lon0.
    '''
    idx=np.absolute(da[lon_name]-lon0).argmin().item()

    var=da.roll(**{lon_name:-idx}, roll_coords=True)
    var[lon_name]=var[lon_name].where(var[lon_name]>=var[lon_name][0],var[lon_name]+360)
    return var
