#!/usr/bin/env python3

'''
Plots temperature depth profiles.
'''

import importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import geosdset, plots, utils

def plot_clim(exp, ds):
    varname='T'
    TFREEZE=273.16
    var=ds[varname].sel(time=slice(*exp.dates)).mean('time')-TFREEZE
    var=var.assign_coords({'lev': -var.lev})
    mask=1.0-np.isnan(var)

    cbar_kwargs={'label': '$^0C$',
    }
    
    fill_opts={'yincrease': False,
               'levels': np.arange(-2.0,33.0,2.0),
               'cmap': cmocean.cm.thermal,
               'cbar_kwargs': cbar_kwargs}

    contour_opts={'levels': np.arange(0.0,33.0,4.0),
                  'colors': 'black',
                  'yincrease': False}

    plot=plots.Plot2d(fill_opts=fill_opts, contour_opts=contour_opts)

    '''
    Plot zonal mean.
    '''

    pl.figure(1)
    wght=mask*ds.dx
    plot_var=utils.average(var,'lon',wght).sel(lat=slice(-85.0,60.0))
    ax=plot.contour(plot_var)
    ax.set_ylim(1500.0,0.0)
    ax.set_title('T, zonal mean')
    ax.set_xlabel('latitude')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.grid()
    pl.savefig(exp.plot_path+'/T_lat_depth.png')

    '''
    Plot equatoral profile (-2S -- 2N).
    '''

    pl.figure(2)
    wght=mask*ds.dy
    plot_var=utils.average(var.sel(lat=slice(-2.1,2.1)),'lat',wght.sel(lat=slice(-2.1,2.1)))
    ax=plot.contour(plot_var)
    ax.set_ylim(500.0,0.0)
    ax.set_title('T, equatorial')
    ax.set_xlabel('longitude')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.grid()
    pl.savefig(exp.plot_path+'/T_eq_depth.png')

    pl.show()
    
def plot_diff(exp, ds1, ds2, type='dif'):
    varname='T'
    TFREEZE=273.16
    var1=ds1[varname].sel(time=slice(*exp.dates)).mean('time')-TFREEZE
    var1=var1.assign_coords({'lev': -var1.lev})
    mask=1.0-np.isnan(var)

    if type=='dif':
        var2=ds2[varname].sel(time=slice(*exp.dates)).mean('time')-TFREEZE
        var2=var1.assign_coords({'lev': -var1.lev})
    elif type=='obs':
        varname='t_an'
        var2=ds2[varname].sel(time=slice(*exp.dates)).mean('time')
    else:
        print('Wrong plot type. Type should be eighter "dif" or "obs"')

    cbar_kwargs={'label': '$^0C$',
    }
    
    fill_opts={'yincrease': False,
               'levels': np.arange(-2.0,33.0,2.0),
               'cmap': cmocean.cm.thermal,
               'cbar_kwargs': cbar_kwargs}

    contour_opts={'levels': np.arange(0.0,33.0,4.0),
                  'colors': 'black',
                  'yincrease': False}

    plot=plots.Plot2d(fill_opts=fill_opts, contour_opts=contour_opts)

    '''
    Plot zonal mean.
    '''

    wght=mask*ds1.dx
    prof1=utils.average(var1,'lon',wght)
    if type=='dif':
        wght=mask*ds2.dx
        prof2=utils.average(var2,'lon',wght)
    elif type=='obs':
        prof2=var2.mean('lon')

def mkplots(exps, dsets):
    plot_clim(exps[0], dsets[0])

    for ds in dsets[1:]:
        plot_diff(exps[0], dsets[0], ds, type='dif')

    obs=['woa13']        
    for name in obs:
        ds=importlib.import_module('verification.'+name).ds_t[obsvarname]
        plot_diff(exps[0], dsets[0], ds, type='obs')

if __name__=='__main__':
    import sys
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn3d')
    mkplots(exps,dsets)
    
