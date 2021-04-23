#!/usr/bin/env python3

'''
Docstring
'''
import sys
import numpy as np
import matplotlib.pyplot as pl
import cmocean
import cartopy.crs as ccrs
import xesmf
import geosdset, plots, utils

def mkclim(exp,dset):
    '''
    Computes climatology for given experiment.
    '''
    varname=['TAUX','TAUY']
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
    ax=plotter.quiver(ds, x='lon', y='lat', u='TAUX', v='TAUY')
#    ax.set_title(f'{exp["expid"]} SST, DJF')
#    pl.savefig(f'{exp["plot_path"]}/sst_djf.png')
#    
#    pl.figure(2); pl.clf()
#    var=clim.sel(season='JJA')
#    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),exp['weight']))
#    ax.set_title(f'{exp["expid"]} SST, JJA')
#    pl.savefig(f'{exp["plot_path"]}/sst_jja.png')
#
#    pl.figure(3); pl.clf()
#    var=clim.mean('season')
#    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),exp['weight']))
#    ax.set_title(f'{exp["expid"]} SST, Annual Mean')
#    pl.savefig(f'{exp["plot_path"]}/sst_am.png')
    pl.show()

def mkplots(exps,dsets):
    # Calculate climatologies for experiments to comapre
    clims=[]
    for exp,dset in zip(exps,dsets):
        clims.append(mkclim(exp,dset))

    # Plot parameters
    projection=ccrs.PlateCarree(central_longitude=210.)
    plotmap=plots.PlotMap(projection=projection)

    # Plots
    plot_clim(plotmap, exps[0], clims[0])

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    mkplots(exps,dsets)
    geosdset.close(dsets)
