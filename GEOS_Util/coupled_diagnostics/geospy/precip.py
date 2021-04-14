#!/usr/bin/env python3

'''
Plots precipitation
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
    varname='RAIN'
    factor=3600*24 # converts from kg/m^2/sec to mm/day
    var=dset[varname].sel(time=slice(*exp['dates']))*factor
    return var.groupby('time.season').mean('time')

def plot_clim(plotter, exp, clim):
    '''
    Makes climaology plots.
    '''
    pl.figure(1); pl.clf() 
    var=clim.sel(season='DJF')
    ax=plotter.contour(var, mode='filled', stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]} PRECIP, DJF')
    pl.savefig(f'{exp["plot_path"]}/precip_djf.png')
    
    pl.figure(2); pl.clf()
    var=clim.sel(season='JJA')
    ax=plotter.contour(var, mode='filled', stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]} PRECIP, JJA')
    pl.savefig(f'{exp["plot_path"]}/precip_jja.png')

    pl.figure(3); pl.clf()
    var=clim.mean('season')
    ax=plotter.contour(var, mode='filled', stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]} PRECIP, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/precip_am.png')
    pl.show()

def mkplots(exps,dsets):
    # Calculate climatologies for experiments to comapre
    clims=[]
    for exp,dset in zip(exps,dsets):
        clims.append(mkclim(exp,dset))

    mask=1-np.isnan(clims[0][0])
    exps[0]["weight"]=mask*dsets[0]["dx"]*dsets[0]["dy"]

    # Plot parameters
    cbar_kwargs={'orientation': 'horizontal',
                 'shrink': 0.8,
                 'label': 'mm/day'}
    
    fill_opts={'cmap': cmocean.cm.rain, 
               'levels': (.2, .4, .6, .8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.3, 2.6, 3, 
                          3.3, 3.6, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10, 11, 13, 15),
               'cbar_kwargs': cbar_kwargs
    }

    projection=ccrs.PlateCarree(central_longitude=210.)
    plotmap=plots.PlotMap(projection=projection, fill_opts=fill_opts)

    # Plots
    plot_clim(plotmap, exps[0], clims[0]) 

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    mkplots(exps,dsets)
    geosdset.close(dsets)
