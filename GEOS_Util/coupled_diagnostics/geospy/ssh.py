#!/usr/bin/env python3

'''
Plots SSH
'''

import sys
import numpy as np
import matplotlib.pyplot as pl
import cartopy.crs as ccrs
import cmocean
import geosdset, plots, utils

plotname='SSH'
defaults={'name': 'SSH', 
          'colname': 'geosgcm_ocn2dT', 
          'coltype': 'GEOSTripolar'}

def mkclim(exp,dset):
    '''
    Computes monthly climatology for given experiment.
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

def plot_clim(plotter, exp, clim):
    '''
    Makes climaology plots.
    '''
    vardata=exp['plots'].get(plotname,defaults)
    varname=vardata['name']
    
    pl.figure(1); pl.clf()
    var=clim[varname].sel(season='DJF')    
    ax=plotter.contour(var)
    ax.set_title(f'{exp["expid"]} SSH, DJF')
    pl.savefig(f'{exp["plot_path"]}/ssh_djf.png')
    
    pl.figure(2); pl.clf() 
    var=clim[varname].sel(season='JJA')
    ax=plotter.contour(var)
    ax.set_title(f'{exp["expid"]} SSH, JJA')
    pl.savefig(f'{exp["plot_path"]}/ssh_jja.png')

def plot_gm(exp,da):
    plotter=plots.Plot1d()
    pl.figure(3); pl.clf()
    ax=plotter.line(da)
    ax.set_title('Global SSH')
    ax.set_ylabel('m')
    ax.set_xlabel('')
    pl.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/ssh_gm.png')
    
def mkplots(exps,dsets):
    clim=mkclim(exps[0],dsets[0])
    gm=mkglobal(exps[0],dsets[0])
    
    # Plot parameters
    cbar_kwargs={'orientation': 'horizontal',
                 'shrink': 0.8,
                 'label': '$m$'}
    
    fill_opts={'cmap': cmocean.cm.diff, 
               'levels': np.arange(-2.,2.1,0.5),
               'cbar_kwargs': cbar_kwargs,
               'x': 'lon',
               'y': 'lat'
    }

    projection=ccrs.PlateCarree(central_longitude=210.)
    plotmap=plots.PlotMap(projection=projection, fill_opts=fill_opts)

    # Plots
    plot_clim(plotmap, exps[0], clim)
    plot_gm(exps[0],gm)
    
def main(exps):
    dsets=geosdset.load_data(exps, plotname, defaults)
    mkplots(exps,dsets)
    geosdset.close(dsets)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    main(exps)
