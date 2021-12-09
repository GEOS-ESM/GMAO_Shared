#!/usr/bin/env python3

'''
Plots SST.
'''

import sys, importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
import cartopy.crs as ccrs
import xesmf
import geosdset, plots, utils

varname='TS'
TFREEZE=273.16    

def mkclim(exp,dset):
    '''
    Computes climatology for given experiment.
    '''
    ds=dset[[varname]].sel(time=slice(*exp['dates']))
    ds=ds.groupby('time.season').mean('time')
    ds[varname]-=TFREEZE
    ds['weight']=dset['mask']*dset['area']
    return ds

def mkglobal(exp,ds):
    var=ds[varname]-TFREEZE
    wght=ds['mask']*ds['area']
    return utils.average(var,('y','x'),wght)
    
def plot_clim(plotter, exp, clim):
    '''
    Makes climaology plots.
    '''
    pl.figure(1); pl.clf() 
    var=clim[varname].sel(season='DJF')
    ax=plotter.contour(var, stat=utils.print_stat(var,('x','y'),clim['weight']))
    ax.set_title(f'{exp["expid"]} SST, DJF')
    pl.savefig(f'{exp["plot_path"]}/sst_djf.png')
    
    pl.figure(2); pl.clf()
    var=clim[varname].sel(season='JJA')
    ax=plotter.contour(var, stat=utils.print_stat(var,('x','y'),clim['weight']))
    ax.set_title(f'{exp["expid"]} SST, JJA')
    pl.savefig(f'{exp["plot_path"]}/sst_jja.png')

    pl.figure(3); pl.clf()
    var=clim[varname].mean('season')
    ax=plotter.contour(var, stat=utils.print_stat(var,('x','y'),clim['weight']))
    ax.set_title(f'{exp["expid"]} SST, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/sst_am.png')
    pl.show()

def plot_diff(plotter, exp, cmpexp, clim, cmpclim):
    '''
    Plots climatology difference between two experiments.
    '''

    rr=xesmf.Regridder(cmpclim[varname],clim[varname],'bilinear',periodic=True)
    dif=clim[varname]-rr(cmpclim[varname])

    pl.figure(1); pl.clf() 
    var=dif.sel(season='DJF')
    ax=plotter.contour(var, stat=utils.print_stat(var,('x','y'),clim['weight']))
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} SST, DJF')
    pl.savefig(f'{exp["plot_path"]}/sst-{cmpexp["expid"]}_djf.png')
    
    pl.figure(2); pl.clf()
    var=dif.sel(season='JJA')
    ax=plotter.contour(var, stat=utils.print_stat(var,('x','y'),clim['weight']))
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} SST, JJA')
    pl.savefig(f'{exp["plot_path"]}/sst-{cmpexp["expid"]}_jja.png')

    pl.figure(3); pl.clf()
    var=dif.mean('season')
    ax=plotter.contour(var, stat=utils.print_stat(var,('x','y'),clim['weight']))
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} SST, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/sst-{cmpexp["expid"]}_am.png')
    pl.show()

    #rr.clean_weight_file()

def plot_diffobs(plotter, exp, clim, obsclim, obsname):
    '''
    Plots climatology difference against observations
    '''
    rr=xesmf.Regridder(obsclim,clim[varname],'bilinear',periodic=True)
    # Just doing clim[varname]-rr(obsclim) does not work, because x,y coords have duplicate values
    dif=clim[varname]
    dif.values-=rr(obsclim).values

    pl.figure(1); pl.clf() 
    var=dif.sel(season='DJF')
    ax=plotter.contour(var, stat=utils.print_stat(var,('x','y'),clim['weight']))
    ax.set_title(f'{exp["expid"]}-{obsname} SST, DJF')
    pl.savefig(f'{exp["plot_path"]}/sst-{obsname}_djf.png')
    pl.savefig(f'{exp["plot_path"]}/sst-obs_djf.png')
    
    pl.figure(2); pl.clf()
    var=dif.sel(season='JJA')
    ax=plotter.contour(var, stat=utils.print_stat(var,('x','y'),clim['weight']))
    ax.set_title(f'{exp["expid"]}-{obsname} SST, JJA')
    pl.savefig(f'{exp["plot_path"]}/sst-{obsname}_jja.png')
    pl.savefig(f'{exp["plot_path"]}/sst-obs_jja.png')

    pl.figure(3); pl.clf()
    var=dif.mean('season')
    ax=plotter.contour(var, stat=utils.print_stat(var,('x','y'),clim['weight']))
    ax.set_title(f'{exp["expid"]}-{obsname} SST, Annual Mean')
    pl.savefig(f'{exp["plot_path"]}/sst-{obsname}_am.png')
    pl.savefig(f'{exp["plot_path"]}/sst-obs_am.png')
    pl.show()

    #rr.clean_weight_file()

def plot_gm(exp,da):
    plotter=plots.Plot1d()
    pl.figure(4); pl.clf()
    ax=plotter.line(da)
    ax.set_title('Global SST')
    ax.set_ylabel('$^0C$')
    ax.set_xlabel('')
    pl.grid()
    pl.tight_layout()
    pl.show()
    pl.savefig(f'{exp["plot_path"]}/sst_gm.png')

def mkplots(exps, dsets):
    # Calculate climatologies for experiments to comapre
    clims=[]
    for exp,dset in zip(exps,dsets):
        clims.append(mkclim(exp,dset))

    gm=mkglobal(exps[0],dsets[0])

    # Plot parameters
    cbar_kwargs={'orientation': 'horizontal',
                 'shrink': 0.8,
                 'label': '$^0C$'}
    
    fill_opts={'cmap': cmocean.cm.thermal, 
              'levels': np.arange(0.0,31.0,2.0),
               'cbar_kwargs': cbar_kwargs,
               'x': 'lon',
               'y': 'lat'
    }

    contour_opts={'levels': np.arange(0.0,31.0,4.0),
                  'colors': 'black',
                  'x': 'lon',
                  'y': 'lat'
    }

    clab_opts={'fmt': '%1.0f'}

    projection=ccrs.PlateCarree(central_longitude=210.)
    plotmap=plots.PlotMap(projection=projection, fill_opts=fill_opts, 
                          contour_opts=contour_opts, clab_opts=clab_opts)

    # Plots
    plot_clim(plotmap, exps[0], clims[0]) 

    plotmap.fill_opts['levels']=np.arange(-10.,10.1,1.0)
    plotmap.fill_opts['cmap']=cmocean.cm.diff
    plotmap.contour_opts['levels']=np.arange(-10.,10.1,2.0)

    for exp,clim in zip(exps[1:],clims[1:]):
        plot_diff(plotmap, exps[0], exp, clims[0], clim)
        
    obs={'OISSTv2': 'sst'} # Names of observational data set and SST variable in this data set.
    for obsname,obsvarname in obs.items():
        da=importlib.import_module('verification.'+obsname).ds[obsvarname]
        obsclim=da.groupby('time.season').mean('time')
        plot_diffobs(plotmap, exps[0], clims[0], obsclim, obsname)

    plot_gm(exps[0],gm)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2dT', type='GEOSTripolar')
    mkplots(exps,dsets)
    geosdset.close(dsets)    
