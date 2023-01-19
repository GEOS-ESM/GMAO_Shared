#!/usr/bin/env python3

'''
Plots sea level.
'''

plotname='SLV'
defaults={'name': 'SLV', 
          'colname': 'geosgcm_ocn2dT', 
          'coltype': 'GEOSTripolar'}

import sys
import numpy as np
import matplotlib.pyplot as pl
import geosdset

def mkclim(exp,dset):
    '''
    Computes climatology for given experiment.
    '''
    vardata=exp['plots'].get(plotname,defaults)
    varname=vardata['name']

    ds=dset[[varname]].sel(time=slice(*exp['dates']))
    ds=ds.groupby('time.season').mean('time')
    ds['weight']=dset['mask']*dset['area']
    return ds

def mkglobal(exp,ds):
    vardata=exp['plots'].get(plotname,defaults)
    varname=vardata['name']
    var=ds[varname]
    wght=ds['mask']*ds['area']
    return utils.average(var,('y','x'),wght)



def mkplots(exps,dsets):
    pass

def main(exps):
    dsets=geosdset.load_data(exps, plotname, defaults)
    mkplots(exps,dsets)
    geosdset.close(dsets)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    main(exps)
