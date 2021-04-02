#!/usr/bin/env python3

'''
Plots SSS.
'''

import sys, importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
import cartopy.crs as ccrs
import xesmf
import geosdset, plots

def plot_clim(exp, da):
    var=da.sel(time=slice(*exp['dates']))
    clim=var.groupby('time.season').mean('time')

    cbar_kwargs={'orientation': 'horizontal',
                 'shrink': 0.8,
                 'label': '$^0C$'}
    
    fill_opts={'cmap': cmocean.cm.haline, 
              'levels': np.arange(32.0,38.1,0.4),
               'cbar_kwargs': cbar_kwargs
    }

    contour_opts={'levels': np.arange(32.0,38.1,0.8),
                  'colors': 'black'
    }

    projection=ccrs.PlateCarree(central_longitude=210.)
    plotmap=plots.PlotMap(projection=projection, fill_opts=fill_opts, contour_opts=contour_opts)

    pl.figure(1); pl.clf() 
    ax=plotmap.contour(clim.sel(season='DJF'))
    ax.set_title('SSS, DJF')
    pl.savefig(f'{exp["plot_path"]}/sss_djf.png')
    
    pl.figure(2); pl.clf()
    ax=plotmap.contour(clim.sel(season='JJA'))
    ax.set_title('SSS, JJA')
    pl.savefig(f'{exp["plot_path"]}/sss_jja.png')

    pl.figure(3); pl.clf()
    ax=plotmap.contour(clim.mean('season'))
    ax.set_title('SSS, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/sss_am.png')
    pl.show()

def plot_diff(exp, da1, da2, ftype='dif'):
    clim1=da1.groupby('time.season').mean('time')
    clim2=da2.groupby('time.season').mean('time')
    rr=xesmf.Regridder(da2,da1,'bilinear',periodic=True)
    dif=clim1-rr(clim2)

    cbar_kwargs={'orientation': 'horizontal',
                 'shrink': 0.8,
                 'label': '$^0C$'}
    
    fill_opts={'cmap': cmocean.cm.diff, 
               'levels': np.arange(-5.,5.1,0.5),
               'cbar_kwargs': cbar_kwargs
    }

    contour_opts={'levels': np.arange(-5.,5.1,1.0),
                  'colors': 'black'
    }

    projection=ccrs.PlateCarree(central_longitude=210.)
    plotmap=plots.PlotMap(projection=projection, fill_opts=fill_opts, contour_opts=contour_opts)

    pl.figure(1); pl.clf() 
    ax=plotmap.contour(dif.sel(season='DJF'))
    ax.set_title(f'SSS-{ftype}, DJF')
    pl.savefig(f'{exp["plot_path"]}/sss-{ftype}_djf.png')
    
    pl.figure(2); pl.clf()
    ax=plotmap.contour(dif.sel(season='JJA'))
    ax.set_title(f'SSS-{ftype}, JJA')
    pl.savefig(f'{exp["plot_path"]}/sss-{ftype}_jja.png')

    pl.figure(3); pl.clf()
    ax=plotmap.contour(dif.mean('season'))
    ax.set_title(f'SSS-{ftype}, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/sss-{ftype}_am.png')
    pl.show()

    rr.clean_weight_file()

def mkplots(exps, dsets):

    varname='SS'
    da=dsets[0][varname]
    plot_clim(exps[0] ,da)

    for ds in dsets[1:]:
        da1=ds[varname]
        plot_diff(exps[0], da, da1, ftype='dif')
        
    obs=['woa13']
    for name in obs:
        obsvarname='s_an'
        da1=importlib.import_module('verification.'+name).ds_s[obsvarname].sel(depth=slice(0.,10.)).mean('depth')
        plot_diff(exps[0], da, da1, ftype='obs')

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    mkplots(exps,dsets)
    
