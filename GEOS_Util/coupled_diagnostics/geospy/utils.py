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
