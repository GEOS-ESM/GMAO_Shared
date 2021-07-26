#!/usr/bin/env python3

'''
Plots surface currents
'''
import sys
import numpy as np
import matplotlib.pyplot as pl
import cartopy.crs as ccrs
import cmocean
import geosdset, plots

def mkclim(exp,dset):
    '''
    Computes climatology for given experiment.
    '''
    varname=['US','VS']
    ds=dset[varname].sel(time=slice(*exp['dates']))
    ds=ds.groupby('time.season').mean('time')
    ds['weight']=dset['MASKO'][0]*dset['dx']*dset['dy']
    return ds

def plot_clim(plotter, exp, clim):
    '''
    Makes climaology plots.
    '''
    pl.figure(1); pl.clf() 
    ds=clim.sel(season='DJF')
    ax=plotter.streamplot(ds, x='lon', y='lat', u='US', v='VS')
    ax.set_title(f'{exp["expid"]} Usurf, DJF')
    pl.savefig(f'{exp["plot_path"]}/Usurf_djf.png')
    
    pl.figure(2); pl.clf() 
    ds=clim.sel(season='JJA')
    ax=plotter.streamplot(ds, x='lon', y='lat', u='US', v='VS')
    ax.set_title(f'{exp["expid"]} Usurf, JJA')
    pl.savefig(f'{exp["plot_path"]}/Usurf_jja.png')

    pl.show()

def mkplots(exps,dsets):
    # Calculate climatologies for experiments to comapre
    clim=mkclim(exps[0],dsets[0])

    # Plot parameters
    cbar_kwargs={'orientation': 'horizontal',
                 'shrink': 0.8,
                 'label': '$m/s$'}
    
    fill_opts={'cmap': cmocean.cm.speed, 
               'levels': (0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.8,1.5),
               'cbar_kwargs': cbar_kwargs
    }

    stream_opts={'density': (4,2), 
                 'color': 'black'}

    projection=ccrs.PlateCarree(central_longitude=210.)
    plotmap=plots.PlotMap(projection=projection, fill_opts=fill_opts,
                          stream_opts=stream_opts)

    # Plots
    plot_clim(plotmap, exps[0], clim)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    mkplots(exps,dsets)
    geosdset.close(dsets)
