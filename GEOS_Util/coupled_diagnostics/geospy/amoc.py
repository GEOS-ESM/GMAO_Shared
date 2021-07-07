#!/usr/bin/env python3

'''
Calculates AMOC stream function.
'''
import sys
import numpy as np
import matplotlib.pyplot as pl
import geosdset

def mkamoc(exp,dset):
    # Read Atlantic mask
    # xx=dset.ty_trans*atl_mask.values
    # amoc=reversed_cumsum(xx,'st_ocean',skipna=True).sum('xt_ocean',skipna=True)+
    # (dset.ty_trams_gm*atl_mask.values).sum('xt_ocean',skipna=True)
    # check if ty_trans_gm needs to be vertically integrated
    pass

def mkplots(exps,dsets):
    pass

if __name__=='__main__':
    exps=[geosdset.load_exps(sys.argv[1])[0]]
    dsets=geosdset.load_collection(exps,'ocean_month',type='MOM')
    mkplots(exps,dsets)
#    geosdset.close(dsets)
