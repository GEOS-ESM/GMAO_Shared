#!/usr/bin/env python3

'''
Docstring
'''
import sys, importlib
import matplotlib.pyplot as pl
import xarray as xr
import cartopy.crs as ccrs
import xesmf, cmocean
import geosdset, plots

def mkclim(exp,dset):
    '''
    Computes climatology for given experiment.
    '''
    varname=['TAUX','TAUY']
    ds=dset[varname].sel(time=slice(*exp['dates']))
    ds=ds.groupby('time.season').mean('time')
    ds['weight']=dset['mask']*dset['area']
    return ds

def plot_clim(plotter, exp, clim):
    '''
    Makes climaology plots.
    '''
    pl.figure(1); pl.clf() 
    ds=clim.sel(season='DJF')
    ax=plotter.quiver(ds, x='lon', y='lat', u='TAUX', v='TAUY')
    ax.set_title(f'{exp["expid"]} TAU, DJF')
    pl.savefig(f'{exp["plot_path"]}/tau_djf.png')
    
    pl.figure(2); pl.clf() 
    ds=clim.sel(season='JJA')
    ax=plotter.quiver(ds, x='lon', y='lat', u='TAUX', v='TAUY')
    ax.set_title(f'{exp["expid"]} TAU, JJA')
    pl.savefig(f'{exp["plot_path"]}/tau_jja.png')

    pl.figure(3); pl.clf() 
    ds=clim.mean('season')
    ax=plotter.quiver(ds, x='lon', y='lat', u='TAUX', v='TAUY')
    ax.set_title(f'{exp["expid"]} TAU, Annual mean')
    pl.savefig(f'{exp["plot_path"]}/tau_am.png')

    pl.show()

def plot_diff(plotter, exp, cmpexp, clim, cmpclim):
    '''
    Plots climatology difference between two experiments
    '''
    rr=xesmf.Regridder(cmpclim,clim,'bilinear',periodic=True)
    cmp_out=rr(cmpclim)
    dif=xr.Dataset()
    dif['TAUX']=clim['TAUX']-cmp_out['TAUX']
    dif['TAUY']=clim['TAUY']-cmp_out['TAUY']    

    pl.figure(1); pl.clf() 
    ds=dif.sel(season='DJF')
    ax=plotter.quiver(ds, x='lon', y='lat', u='TAUX', v='TAUY')
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} TAU, DJF')
    pl.savefig(f'{exp["plot_path"]}/tau-{cmpexp["expid"]}_djf.png')

    pl.figure(2); pl.clf() 
    ds=dif.sel(season='JJA')
    ax=plotter.quiver(ds, x='lon', y='lat', u='TAUX', v='TAUY')
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} TAU, JJA')
    pl.savefig(f'{exp["plot_path"]}/tau-{cmpexp["expid"]}_jja.png')

    pl.figure(3); pl.clf() 
    ds=dif.mean('season')
    ax=plotter.quiver(ds, x='lon', y='lat', u='TAUX', v='TAUY')
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} TAU, Annual mean')
    pl.savefig(f'{exp["plot_path"]}/tau-{cmpexp["expid"]}_am.png')

    pl.show()

    #rr.clean_weight_file()

def plot_diffobs(plotter, exp, clim, obsclim, obsname,obsvarname):
    '''
    Plots climatology difference against observations
    '''
    rr=xesmf.Regridder(obsclim,clim,'bilinear',periodic=True)
    obs_out=rr(obsclim)
    dif=xr.Dataset()
    dif['TAUX']=clim['TAUX']
    dif['TAUY']=clim['TAUY']
    dif['TAUX'].values-=obs_out[obsvarname[0]].values
    dif['TAUY'].values-=obs_out[obsvarname[1]].values

    pl.figure(1); pl.clf() 
    ds=dif.sel(season='DJF')
    ax=plotter.quiver(ds, x='lon', y='lat', u='TAUX', v='TAUY')
    ax.set_title(f'{exp["expid"]}-{obsname} TAU, DJF')
    pl.savefig(f'{exp["plot_path"]}/tau-{obsname}_djf.png')
    pl.savefig(f'{exp["plot_path"]}/tau-obs_djf.png')

    pl.figure(2); pl.clf() 
    ds=dif.sel(season='JJA')
    ax=plotter.quiver(ds, x='lon', y='lat', u='TAUX', v='TAUY')
    ax.set_title(f'{exp["expid"]}-{obsname} TAU, JJA')
    pl.savefig(f'{exp["plot_path"]}/tau-{obsname}_jja.png')
    pl.savefig(f'{exp["plot_path"]}/tau-obs_jja.png')

    pl.figure(3); pl.clf() 
    ds=dif.mean('season')
    ax=plotter.quiver(ds, x='lon', y='lat', u='TAUX', v='TAUY')
    ax.set_title(f'{exp["expid"]}-{obsname} TAU, Annual mean')
    pl.savefig(f'{exp["plot_path"]}/tau-{obsname}_am.png')
    pl.savefig(f'{exp["plot_path"]}/tau-obs_am.png')

    pl.show()

    #rr.clean_weight_file()

def mkplots(exps,dsets):
    # Calculate climatologies for experiments to comapre
    clims=[]
    for exp,dset in zip(exps,dsets):
        clims.append(mkclim(exp,dset))

    # Plot parameters
    cbar_kwargs={'orientation': 'horizontal',
                 'shrink': 0.8,
                 'label': '$N/m^2$'}
    
    fill_opts={'cmap': cmocean.cm.speed, 
              'levels': (0.01,0.02,0.04,0.06,0.08,0.1,0.15,0.2,0.25,0.3),
               'cbar_kwargs': cbar_kwargs,
               'x': 'lon',
               'y': 'lat'
    }

    projection=ccrs.PlateCarree(central_longitude=210.)
    plotmap=plots.PlotMap(projection=projection, fill_opts=fill_opts)

    # Plots
    plot_clim(plotmap, exps[0], clims[0])

    for exp,clim in zip(exps[1:],clims[1:]):
        plot_diff(plotmap, exps[0], exp, clims[0], clim)

    obs={'QSCAT': ['taux','tauy']} # Names of observational data set and variable in this data set.
    for obsname,obsvarname in obs.items():
        ds=importlib.import_module('verification.'+obsname).ds
        obsclim=ds.groupby('time.season').mean('time')
        plot_diffobs(plotmap, exps[0], clims[0], obsclim, obsname, obsvarname)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2dT',type='GEOSTripolar')
    mkplots(exps,dsets)
    geosdset.close(dsets)
