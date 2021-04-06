'''
This module provides different diagnostics utilities.
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

def print_stat(da, dim, weight):
    '''
    Returns a string with da mean and standard deviation.
    '''
    mean=average(da, dim, weight)
    std=np.sqrt(average((da-mean)**2,dim,weight))
    return f'mean: {mean.values:.2f}\nstd: {std.values:.2f}'
