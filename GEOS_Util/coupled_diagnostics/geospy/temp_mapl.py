#!/usr/bin/env python3

'''
Plots temperature depth profiles.
'''

import importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import geosdset, plots

def plot_clim(exp, da):
    var=da.sel(time=slice(*exp.dates)).mean('time')

    cbar_kwargs={'label': '$^0C$',
                 'extend': 'both'
    }
    
    fill_opts={'yincrease': False,
               'levels': np.arange(-2.0,33.0,2.0),
               'cmap': cmocean.cm.thermal,
               'cbar_kwargs': cbar_kwargs}

    contour_opts={'levels': np.arange(0.0,33.0,4.0),
                  'colors': 'black',
                  'yincrease': False}

    '''
    Note: to do proper spatiall average of 3d fields, 
    3d metrics with 3d land mask are needed. 
    '''
    
    plot=plots.Plot2d(fill_opts=fill_opts, contour_opts=contour_opts)

    pl.figure(1)
    plot_var=var.sel(lev=slice(0.0,3000.0),lat=slice(-85.0,60.0)).mean('lon')
    ax=plot.contour(plot_var)
    ax.set_title('T, zonal mean')
    ax.set_xlabel('latitude')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.grid()
    pl.savefig(exp.plot_path+'/T_lat_depth.png')

    pl.figure(2)
    plot_var=var.sel(lat=slice(-2.1,2.1),lev=slice(0.0,500.0)).mean('lat')
    ax=plot.contour(plot_var)
    ax.set_title('T, equatorial')
    ax.set_xlabel('longitude')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.grid()
    pl.savefig(exp.plot_path+'/T_eq_depth.png')

    pl.show()

def plot_diff(exp, ds1, ds2, type='dif'):
    pass

def mkplots(exps, dsets):
    varname='T'
    TFREEZE=273.16
    da=dsets[0][varname]; da-=TFREEZE
    da=da.assign_coords({'lev': -da.lev})
    plot_clim(exps[0] ,da)

    for ds in dsets[1:]:
        da1=ds[varname]; da1-=TFREEZE
        da1=da1.assign_coords({'lev': -da1.lev})
        plot_diff(exps[0], da, da1, type='dif')

    obs=['woa13']        
    for name in obs:
        obsvarname='t_an'
        da1=importlib.import_module('verification.'+name).ds_t[obsvarname]
        plot_diff(exps[0], da, da1, type='obs')

if __name__=='__main__':
    import sys
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn3d')
    mkplots(exps,dsets)
    
