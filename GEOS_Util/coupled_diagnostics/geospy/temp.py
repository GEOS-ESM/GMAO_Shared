#!/usr/bin/env python3

import importlib
import matplotlib.pyplot as pl
import cmocean
import geosdset

def plot_clim(exp, ds):
    varid='temp'
    var=ds._ds[var].mean('Time')

    plotopts={'yincrease': False,
              'vmin': 0,
              'vmax': 30,
              'cmap': cmocean.cm.thermal}

    cbaropts={}

    # Replace data below with proper spatially averaged data using xgcm
    pl.figure(1)
    temp.mean('xh').plot(**plotopts)

    pl.figure(2)
    temp.mean('yh').plot(**plotopts)    
    
    pl.show()

def plot_diff(exp, ds1, ds2):
    pass

def plots(exps, dsets):
    obs=[]
    plot_clim(exps[0] ,dsets[0])

    for ds in dsets:
        plot_diff(exps[0], dsets[0], ds)
        
    for ds in obs:
        plot_diff(exps[0], dsets[0], ds)

if __name__=='__main__':
    import sys
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'prog_z',type='MOM')
    plots(exps,dsets)
    
