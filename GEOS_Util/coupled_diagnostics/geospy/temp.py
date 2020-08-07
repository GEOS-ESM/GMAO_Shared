#!/usr/bin/env python3

'''
Plots temperature depth profiles.
'''

import importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import geosdset

def plot_clim(exp, ds):
    varid='temp'
    var=ds[varid].sel(Time=slice(*exp.dates)).mean('Time')

    cbar_kwargs={'label': '$^0C$'}
    
    fill_opts={'yincrease': False,
               'levels': np.arange(0.0,31.0,2.0),
               'cmap': cmocean.cm.thermal,
               'cbar_kwargs': cbar_kwargs}

    contour_opts={'levels': np.arange(0.0,31.0,4.0),
                  'colors': 'black',
                  'yincrease': False}

    '''
    Note: to do proper spatiall average of 3d fields, 
    3d metrics with 3d land mask are needed. 
    '''
    pl.figure(1)
    plot_var=var.sel(z_l=slice(0.0,3000.0),yh=slice(-85.0,60.0)).mean('xh')
    plot_var.plot.contourf(**fill_opts)
    cs=plot_var.plot.contour(**contour_opts)
    cs.clabel(fmt='%1.0f')
    ax=pl.gca()
    ax.set_title('T, zonal mean')
    ax.set_xlabel('latitude')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.grid()
    pl.savefig(exp.plot_path+'/'+varid+'_lat_depth.png')

    pl.figure(2)
    plot_var=var.sel(yh=slice(-2.1,2.1),z_l=slice(0.0,500.0)).mean('yh')
    plot_var.plot.contourf(**fill_opts)
    cs=plot_var.plot.contour(**contour_opts)
    cs.clabel(fmt='%1.0f')
    ax=pl.gca()
    ax.set_title('T, equatorial')
    ax.set_xlabel('longitude')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.grid()
    pl.savefig(exp.plot_path+'/'+varid+'_eq_depth.png')

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
    
