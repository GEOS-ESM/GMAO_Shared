#!/usr/bin/env python3

'''
Plots SST.
'''

import importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
import cartopy.crs as ccrs
import xesmf
import geosdset, plots

def plot_clim(exp, da):
    var=da.sel(time=slice(*exp.dates))
    clim=var.groupby('time.season').mean('time')

    cbar_kwargs={'orientation': 'horizontal',
                 'shrink': 0.8,
                 'label': '$^0C$'}
    
    fill_opts={'cmap': cmocean.cm.thermal, 
              'levels': np.arange(0.0,31.0,2.0),
               'cbar_kwargs': cbar_kwargs
    }

    contour_opts={'levels': np.arange(0.0,31.0,4.0),
                  'colors': 'black'
    }

    projection=ccrs.PlateCarree(central_longitude=210.)
    plotmap=plots.PlotMap(projection=projection, fill_opts=fill_opts, contour_opts=contour_opts)

    pl.figure(1); pl.clf() 
    ax=plotmap.contour(clim.sel(season='DJF'))
    ax.set_title('SST, DJF')
    pl.savefig(exp.plot_path+'/sst_djf.png')
    
    pl.figure(2); pl.clf()
    ax=plotmap.contour(clim.sel(season='JJA'))
    ax.set_title('SST, JJA')
    pl.savefig(exp.plot_path+'/sst_jja.png')

    pl.figure(3); pl.clf()
    ax=plotmap.contour(clim.mean('season'))
    ax.set_title('SST, Annual Mean')
    pl.savefig(exp.plot_path+'/sst_am.png')
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
              'levels': np.arange(-10.,10.1,1.0),
               'cbar_kwargs': cbar_kwargs
    }

    contour_opts={'levels': np.arange(-10.,10.1,2.0),
                  'colors': 'black'
    }

    projection=ccrs.PlateCarree(central_longitude=210.)
    plotmap=plots.PlotMap(projection=projection, fill_opts=fill_opts, contour_opts=contour_opts)

    pl.figure(1); pl.clf() 
    ax=plotmap.contour(dif.sel(season='DJF'))
    ax.set_title('SST-'+ftype+', DJF')
    pl.savefig(exp.plot_path+'/sst-'+ftype+'_djf.png')
    
    pl.figure(2); pl.clf()
    ax=plotmap.contour(dif.sel(season='JJA'))
    ax.set_title('SST-'+ftype+', JJA')
    pl.savefig(exp.plot_path+'/sst-'+ftype+'_jja.png')

    pl.figure(3); pl.clf()
    ax=plotmap.contour(dif.mean('season'))
    ax.set_title('SST-'+ftype+', Annual Mean')
    pl.savefig(exp.plot_path+'/sst-'+ftype+'_am.png')
    pl.show()

    rr.clean_weight_file()

def mkplots(exps, dsets):

    varname='TS'
    TFREEZE=273.16
    da=dsets[0][varname]; da-=TFREEZE
    plot_clim(exps[0] ,da)

    for ds in dsets[1:]:
        da1=ds[varname]; da1-=TFREEZE
        plot_diff(exps[0], da, da1, ftype='dif')
        
    obs=['OISSTv2']
    for name in obs:
        obsvarname='sst'
        da1=importlib.import_module('verification.'+name).ds[obsvarname]
        plot_diff(exps[0], da, da1, ftype='obs')

if __name__=='__main__':
    import sys
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    mkplots(exps,dsets)
    
