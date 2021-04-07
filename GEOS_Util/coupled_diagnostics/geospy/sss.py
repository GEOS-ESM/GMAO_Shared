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
import geosdset, plots, utils

def mkclim(exp,dset):
    '''
    Computes climatology for given experiment.
    '''
    varname='SS'
    var=dset[varname].sel(time=slice(*exp['dates']))
    return var.groupby('time.season').mean('time')
    
def plot_clim(plotter, exp, clim):
    '''
    Makes climaology plots.
    '''
    pl.figure(1); pl.clf() 
    var=clim.sel(season='DJF')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]} SSS, DJF')
    pl.savefig(f'{exp["plot_path"]}/sss_djf.png')
    
    pl.figure(2); pl.clf()
    var=clim.sel(season='JJA')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]} SSS, JJA')
    pl.savefig(f'{exp["plot_path"]}/sss_jja.png')

    pl.figure(3); pl.clf()
    var=clim.mean('season')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]} SSS, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/sss_am.png')
    pl.show()

def plot_diff(plotter, exp, cmpexp, clim, cmpclim):
    '''
    Plots climatology difference between two experiments.
    '''

    rr=xesmf.Regridder(cmpclim,clim,'bilinear',periodic=True)
    dif=clim-rr(cmpclim)

    pl.figure(1); pl.clf() 
    var=dif.sel(season='DJF')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} SSS, DJF')
    pl.savefig(f'{exp["plot_path"]}/sss-{cmpexp["expid"]}_djf.png')
    
    pl.figure(2); pl.clf()
    var=dif.sel(season='JJA')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} SSS, JJA')
    pl.savefig(f'{exp["plot_path"]}/sss-{cmpexp["expid"]}_jja.png')

    pl.figure(3); pl.clf()
    var=dif.mean('season')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} SSS, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/sss-{cmpexp["expid"]}_am.png')
    pl.show()

    rr.clean_weight_file()

def plot_diffobs(plotter, exp, clim, obsclim, obsname):
    '''
    Plots climatology difference against observations
    '''
    rr=xesmf.Regridder(obsclim,clim,'bilinear',periodic=True)
    dif=clim-rr(obsclim)

    pl.figure(1); pl.clf() 
    var=dif.sel(season='DJF')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]}-{obsname} SSS, DJF')
    pl.savefig(f'{exp["plot_path"]}/sss-{obsname}_djf.png')
    pl.savefig(f'{exp["plot_path"]}/sss-obs_djf.png')
    
    pl.figure(2); pl.clf()
    var=dif.sel(season='JJA')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]}-{obsname} SSS, JJA')
    pl.savefig(f'{exp["plot_path"]}/sss-{obsname}_jja.png')
    pl.savefig(f'{exp["plot_path"]}/sss-obs_jja.png')

    pl.figure(3); pl.clf()
    var=dif.mean('season')
    ax=plotter.contour(var, stat=utils.print_stat(var,('lon','lat'),exp['weight']))
    ax.set_title(f'{exp["expid"]}-{obsname} SSS, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/sss-{obsname}_am.png')
    pl.savefig(f'{exp["plot_path"]}/sss-obs_am.png')
    pl.show()

    rr.clean_weight_file()

def mkplots(exps, dsets):
# Calculate climatologies for experiments to comapre
    clims=[]
    for exp,dset in zip(exps,dsets):
        clims.append(mkclim(exp,dset))

    mask=1-np.isnan(clims[0][0])
    exps[0]["weight"]=mask*dsets[0]["dx"]*dsets[0]["dy"]

# Plot parameters
    cbar_kwargs={'orientation': 'horizontal',
                 'shrink': 0.8,
                 'label': 'PSU'}
    
    fill_opts={'cmap': cmocean.cm.haline, 
              'levels': np.arange(32.0,38.1,0.4),
               'cbar_kwargs': cbar_kwargs
    }

    contour_opts={'levels': np.arange(32.0,38.1,0.8),
                  'colors': 'black'
    }

    projection=ccrs.PlateCarree(central_longitude=210.)
    plotmap=plots.PlotMap(projection=projection, fill_opts=fill_opts, 
                          contour_opts=contour_opts)

    # Plots
    plot_clim(plotmap, exps[0], clims[0]) 

    plotmap.fill_opts['levels']=np.arange(-5.,5.1,0.5)
    plotmap.fill_opts['cmap']=cmocean.cm.diff
    plotmap.contour_opts['levels']=np.arange(-5.,5.1,1.0)

    for exp,clim in zip(exps[1:],clims[1:]):
        plot_diff(plotmap, exps[0], exp, clims[0], clim)
        
    obs={'woa13': 's_an'} # Names of observational data set and SST variable in this data set.
    for obsname,obsvarname in obs.items():
        ds=importlib.import_module('verification.'+obsname).ds_s
        da=ds[obsvarname].sel(depth=slice(0.,10.)).mean('depth')
        obsclim=da.groupby('time.season').mean('time')
        plot_diffobs(plotmap, exps[0], clims[0], obsclim, obsname)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    mkplots(exps,dsets)
    
