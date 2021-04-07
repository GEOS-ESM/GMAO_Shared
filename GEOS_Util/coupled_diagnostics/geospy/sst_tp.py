#!/usr/bin/env python3

'''
Plots tropical Pacific SST
'''
import sys
import numpy as np
import matplotlib.pyplot as pl
import geosdset, plots, utils

def mkclim(exp,ds):
    varname='TS'
    TFREEZE=273.16
    var=ds[varname].sel(time=slice(*exp['dates']))-TFREEZE
    var=var.groupby('time.month').mean('time')
    mask=1.0-np.isnan(var)
    wght=mask*ds.dy
    var=utils.average(var.sel(lat=slice(-2.1,2.1)),'lat',wght.sel(lat=slice(-2.1,2.1)))
    return utils.shift_lon(var,30)

def plot_equatorial(exps,clims):
    pl.figure(1); pl.clf()
    for exp,clim in zip(exps,clims):
        clim.mean('month').plot()
    ax=pl.gca()
    ax.legend()
    pl.grid()
    pl.show()

def mkplots(exps,dsets):
    # Calculate seasonal equatorial means
    clims=[]
    for exp,dset in zip(exps,dsets):
        clims.append(mkclim(exp,dset))

    # Plots
    plot_equatorial(exps,clims)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    mkplots(exps,dsets)
