#!/usr/bin/env python3

'''
Plots precipitation
'''
import sys, importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
import cartopy.crs as ccrs
import xesmf
import geosdset, plots, utils

varname='TPREC'

def mkclim(exp,dset):
    '''
    Computes climatology for given experiment.
    '''
    factor=3600*24 # converts from kg/m^2/sec to mm/day
    ds=dset[[varname]].sel(time=slice(*exp['dates']))
    ds=ds.groupby('time.season').mean('time')
    ds[varname]*=factor
    ds['weight']=dset['dx']*dset['dy']
    return ds

def plot_clim(plotter, exp, clim):
    '''
    Makes climaology plots.
    '''
    pl.figure(1); pl.clf() 
    var=clim[varname].sel(season='DJF')
    ax=plotter.contour(var, mode='filled', stat=utils.print_stat(var,('lon','lat'),clim['weight']))
    ax.set_title(f'{exp["expid"]} PRECIP, DJF')
    pl.savefig(f'{exp["plot_path"]}/precip_djf.png')
    
    pl.figure(2); pl.clf()
    var=clim[varname].sel(season='JJA')
    ax=plotter.contour(var, mode='filled', stat=utils.print_stat(var,('lon','lat'),clim['weight']))
    ax.set_title(f'{exp["expid"]} PRECIP, JJA')
    pl.savefig(f'{exp["plot_path"]}/precip_jja.png')

    pl.figure(3); pl.clf()
    var=clim[varname].mean('season')
    ax=plotter.contour(var, mode='filled', stat=utils.print_stat(var,('lon','lat'),clim['weight']))
    ax.set_title(f'{exp["expid"]} PRECIP, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/precip_am.png')
    pl.show()

def plot_diff(plotter, exp, cmpexp, clim, cmpclim):
    '''
    Plots climatology difference between two experiments.
    '''

    rr=xesmf.Regridder(cmpclim[varname],clim[varname],'bilinear',periodic=True)
    dif=clim[varname]-rr(cmpclim[varname])

    pl.figure(1); pl.clf() 
    var=dif.sel(season='DJF')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),clim['weight']), mode='filled')
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} PRECIP, DJF')
    pl.savefig(f'{exp["plot_path"]}/precip-{cmpexp["expid"]}_djf.png')
    
    pl.figure(2); pl.clf()
    var=dif.sel(season='JJA')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),clim['weight']), mode='filled')
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} PRECIP, JJA')
    pl.savefig(f'{exp["plot_path"]}/precip-{cmpexp["expid"]}_jja.png')

    pl.figure(3); pl.clf()
    var=dif.mean('season')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),clim['weight']), mode='filled')
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} PRECIP, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/precip-{cmpexp["expid"]}_am.png')
    pl.show()

    rr.clean_weight_file()

def plot_diffobs(plotter, exp, clim, obsclim, obsname):
    '''
    Plots climatology difference against observations
    '''
    rr=xesmf.Regridder(obsclim,clim[varname],'bilinear',periodic=True)
    dif=clim[varname]-rr(obsclim)

    pl.figure(1); pl.clf() 
    var=dif.sel(season='DJF')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),clim['weight']), mode='filled')
    ax.set_title(f'{exp["expid"]}-{obsname} PRECIP, DJF')
    pl.savefig(f'{exp["plot_path"]}/precip-{obsname}_djf.png')
    pl.savefig(f'{exp["plot_path"]}/precip-obs_djf.png')
    
    pl.figure(2); pl.clf()
    var=dif.sel(season='JJA')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),clim['weight']), mode='filled')
    ax.set_title(f'{exp["expid"]}-{obsname} PRECIP, JJA')
    pl.savefig(f'{exp["plot_path"]}/precip-{obsname}_jja.png')
    pl.savefig(f'{exp["plot_path"]}/precip-obs_jja.png')

    pl.figure(3); pl.clf()
    var=dif.mean('season')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),clim['weight']), mode='filled')
    ax.set_title(f'{exp["expid"]}-{obsname} PRECIP, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/precip-{obsname}_am.png')
    pl.savefig(f'{exp["plot_path"]}/precip-obs_am.png')
    pl.show()

    rr.clean_weight_file()

def mkplots(exps,dsets):
    # Calculate climatologies for experiments to comapre
    clims=[]
    for exp,dset in zip(exps,dsets):
        clims.append(mkclim(exp,dset))

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

    plotmap.fill_opts['levels']=(-8, -7, -6, -5, -4, -3, -2, -1, -.5, 0., 
                                 .5, 1, 2, 3, 4, 5, 6, 7, 8)
    plotmap.fill_opts['cmap']=cmocean.cm.diff

    for exp,clim in zip(exps[1:],clims[1:]):
        plot_diff(plotmap, exps[0], exp, clims[0], clim)

    obs={'GPCP': 'precip'} # Names of observational data set and variable in this data set.
    for obsname,obsvarname in obs.items():
        da=importlib.import_module('verification.'+obsname).ds[obsvarname].sel(lev=0.0)
        obsclim=da.groupby('time.season').mean('time')
        plot_diffobs(plotmap, exps[0], clims[0], obsclim, obsname)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_surf')
    mkplots(exps,dsets)
    geosdset.close(dsets)
