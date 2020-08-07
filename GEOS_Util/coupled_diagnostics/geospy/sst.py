#!/usr/bin/env python3

'''
Plots SST.
'''

import importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
import geosdset
from plotutils import quickmap, drawstats

def plot_clim(exp, ds):
    varid='TS'
    FREEZE=273.16
    var=ds[varid].sel(time=slice(*exp.dates)); var-=FREEZE
    clim=var.groupby('time.season').mean('time')
    ann=clim.mean('season')

    cbar_kwargs={'orientation': 'horizontal',
                 'shrink': 0.8,
                 'label': '$^0C$'}
    
    fill_opts={'cmap': cmocean.cm.thermal, 
              'levels': np.arange(0.0,31.0,2.0),
              'cbar_kwargs': cbar_kwargs}

    contour_opts={'levels': np.arange(0.0,31.0,4.0),
               'colors': 'black'}

    ff=pl.figure(1) 
    ax=quickmap()
    clim.sel(season='DJF').plot.contourf(ax=ax,**fill_opts)
    cs=clim.sel(season='DJF').plot.contour(ax=ax,**contour_opts)
    cs.clabel(fmt='%1.0f')
    ax.set_title('SST, DJF')
    pl.savefig(exp.plot_path+'/sst_djf.png')

    pl.figure(2)
    ax=quickmap()
    clim.sel(season='JJA').plot.contourf(ax=ax,**fill_opts)
    cs=clim.sel(season='JJA').plot.contour(ax=ax,**contour_opts)
    cs.clabel(fmt='%1.0f')
    ax.set_title('SST, JJA')
    pl.savefig(exp.plot_path+'/sst_jja.png')

    pl.figure(3)
    ax=quickmap()
    ann.plot.contourf(ax=ax,**fill_opts)
    cs=ann.plot.contour(ax=ax,**contour_opts)
    cs.clabel(fmt='%1.0f')
    ax.set_title('SST, Annual Mean')
    pl.savefig(exp.plot_path+'/sst_am.png')
    pl.show()

def plot_diff(exp, ds1, ds2):
    pass

def plots(exps, dsets):
    obs=['OISSTv2']
    plot_clim(exps[0] ,dsets[0])

    for ds in dsets:
        plot_diff(exps[0], dsets[0], ds)
        
    for ds in obs:
        plot_diff(exps[0], dsets[0], ds)

if __name__=='__main__':
    import sys
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    plots(exps,dsets)
    
