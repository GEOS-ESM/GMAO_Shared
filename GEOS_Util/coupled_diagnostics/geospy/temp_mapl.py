#!/usr/bin/env python3

'''
Plots zonal/meridional depth profiles.
'''

import importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import geosdset, plots, utils

# Globals
varname='T'
TFREEZE=273.16
obs=['woa13']
obsvarname='t_an'

def plot_clim(exp, ds):
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

    pl.figure(1); pl.clf()
    wght=mask*ds.dx
    plot_var=utils.average(var,'lon',wght)
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

    pl.figure(2); pl.clf()
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
    
def plot_diff(exp, ds1, ds2, ftype='dif'):
    var1=ds1[varname].sel(time=slice(*exp.dates)).mean('time')-TFREEZE
    var1=var1.assign_coords({'lev': -var1.lev})
    mask=1.0-np.isnan(var1)
    wght=mask*ds1.dx
    zonal1=utils.average(var1,'lon',wght)
    wght=mask*ds1.dy
    equatorial1=utils.average(var1.sel(lat=slice(-2.1,2.1)),'lat',wght.sel(lat=slice(-2.1,2.1)))

    if ftype=='dif':
        var2=ds2[varname].sel(time=slice(*exp.dates)).mean('time')-TFREEZE
        var2=var1.assign_coords({'lev': -var1.lev})
        mask=1.0-np.isnan(var2)
        wght=mask*ds2.dx
        zonal2=utils.average(var2,'lon',wght)
        wght=mask*ds2.dy
        equatorial2=utils.average(var.sel(lat=slice(-2.1,2.1)),'lat',wght.sel(lat=slice(-2.1,2.1)))
    elif ftype=='obs':
        var2=ds2[obsvarname].mean('time')
        zonal2=var2.mean('lon').rename({'depth':'lev'})
        equatorial2=var2.sel(lat=slice(-2.1,2.1)).mean('lat').rename({'depth':'lev'})
    else:
        print('Wrong plot type. Type should be eighter "dif" or "obs"'); raise

    cbar_kwargs={'label': '$^0C$',
    }
    
    fill_opts={'yincrease': False,
               'levels': np.arange(-5.0,5.1,0.5),
               'cmap': cmocean.cm.diff,
               'cbar_kwargs': cbar_kwargs}

    contour_opts={'levels': np.arange(-5.0,5.1,1.0),
                  'colors': 'black',
                  'yincrease': False}

    plot=plots.Plot2d(fill_opts=fill_opts, contour_opts=contour_opts)

    '''
    Plot zonal mean profile.
    '''
    pl.figure(1); pl.clf()
    dif=zonal1-zonal2.interp_like(zonal1)
    ax=plot.contour(dif)
    ax.set_ylim(1500.0,0.0)
    ax.set_title('T-'+ftype+', zonal mean')
    ax.set_xlabel('latitude')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.grid()
    pl.savefig(exp.plot_path+'/T_'+ftype+'_lat_depth.png')
    pl.show()

    '''
    Plot equatorial profile.
    '''
    pl.figure(2); pl.clf()
    dif=equatorial1-equatorial2.interp_like(equatorial1)
    ax=plot.contour(dif)
    ax=pl.gca()
    ax.set_ylim(500.0,0.0)
    ax.set_title('T-'+ftype+', equatorial')
    ax.set_xlabel('longitude')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.grid()
    pl.savefig(exp.plot_path+'/T_'+ftype+'_eq_depth.png')
    pl.show()

def mkplots(exps, dsets):
#    plot_clim(exps[0], dsets[0])

    for ds in dsets[1:]:
        plot_diff(exps[0], dsets[0], ds, ftype='dif')

    for name in obs:
        ds=importlib.import_module('verification.'+name).ds_t
        plot_diff(exps[0], dsets[0], ds, ftype='obs')

if __name__=='__main__':
    import sys
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn3d')
    mkplots(exps,dsets)
    
