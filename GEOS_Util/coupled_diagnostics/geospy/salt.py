#!/usr/bin/env python3

'''
Plots zonal/meridional depth profiles.
'''

import sys,importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
import geosdset, plots, utils

# Globals
varname='S'

def mkzonal(exp,ds):
    var=ds[varname].sel(time=slice(*exp['dates'])).mean('time')
    var.coords['lev']=-var['lev']
    mask=1.0-np.isnan(var)
    wght=mask*ds.dx
    return utils.average(var,'lon',wght)

def mkequatorial(exp,ds):
    var=ds[varname].sel(time=slice(*exp['dates'])).mean('time')
    var.coords['lev']=-var['lev']
    mask=1.0-np.isnan(var)
    wght=mask*ds.dy
    var=utils.average(var.sel(lat=slice(-2.1,2.1)),'lat',wght.sel(lat=slice(-2.1,2.1)))
    return utils.shift_lon(var,30)

def plot_zonal(plotter, exp, zonal):
    '''
    Plots zonal mean profile.
    '''
    pl.figure(1); pl.clf()
    ax=plotter.contour(zonal)
    ax.set_ylim(1400.0,0.0)
    ax.set_xlim(-85.0,65.0)
    ax.set_title(f'{exp["expid"]} {varname}, zonal mean')
    ax.set_xlabel('')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(plots.LATITUDE_FORMATTER)
    ax.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/{varname}_lat_depth.png')
    pl.show()

def plot_equatorial(plotter, exp, equatorial):
    '''
    Plots equatoral profile (-2S -- 2N).
    '''
    pl.figure(2); pl.clf()
    ax=plotter.contour(equatorial)
    ax.set_ylim(500.0,0.0)
    ax.set_title(f'{exp["expid"]} {varname}, equatorial')
    ax.set_xlabel('')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(plots.LONGITUDE_FORMATTER)
    ax.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/{varname}_eq_depth.png')
    pl.show()
    
def plot_zonal_diff(plotter, exp, cmpexp, zonal, cmpzonal):
    '''
    Plots zonal mean profile.
    '''
    pl.figure(1); pl.clf()
    dif=zonal-cmpzonal.interp_like(zonal)
    ax=plotter.contour(dif)
    ax.set_ylim(1400.0,0.0)
    ax.set_xlim(-85.0,65.0)
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} {varname}, zonal mean')
    ax.set_xlabel('')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(plots.LATITUDE_FORMATTER)
    ax.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/{varname}-{cmpexp["expid"]}_lat_depth.png')
    pl.show()

def plot_equatorial_diff(plotter, exp, cmpexp, eq, cmpeq):
    '''
    Plots equatorial profile.
    '''
    pl.figure(2); pl.clf()
    dif=eq-cmpeq.interp_like(eq)
    ax=plotter.contour(dif)
    ax=pl.gca()
    ax.set_ylim(500.0,0.0)
    ax.set_title(f'{exp["expid"]}-{cmpexp["expid"]} {varname}, equatorial')
    ax.set_xlabel('')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(plots.LONGITUDE_FORMATTER)
    ax.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/{varname}-{cmpexp["expid"]}_eq_depth.png')
    pl.show()

def plot_zonal_diffobs(plotter, exp, zonal, obszonal, obsname):
    '''
    Plots zonal mean profile.
    '''
    pl.figure(1); pl.clf()
    dif=zonal-obszonal.interp_like(zonal)
    ax=plotter.contour(dif)
    ax.set_ylim(1400.0,0.0)
    ax.set_xlim(-85.0,65.0)
    ax.set_title(f'{exp["expid"]}-{obsname} {varname}, zonal mean')
    ax.set_xlabel('')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(plots.LATITUDE_FORMATTER)
    ax.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/{varname}-{obsname}_lat_depth.png')
    pl.savefig(f'{exp["plot_path"]}/{varname}-obs_lat_depth.png')
    pl.show()

def plot_equatorial_diffobs(plotter, exp, eq, obseq, obsname):
    '''
    Plots equatorial profile.
    '''
    pl.figure(2); pl.clf()
    dif=eq-obseq.interp_like(eq)
    ax=plotter.contour(dif)
    ax=pl.gca()
    ax.set_ylim(500.0,0.0)
    ax.set_title(f'{exp["expid"]}-{obsname} {varname}, equatorial')
    ax.set_xlabel('')
    ax.set_ylabel('depth')
    ax.xaxis.set_major_formatter(plots.LONGITUDE_FORMATTER)
    ax.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/{varname}-{obsname}_eq_depth.png')
    pl.savefig(f'{exp["plot_path"]}/{varname}-obs_eq_depth.png')
    pl.show()

def mkplots(exps, dsets):
    # Calculate zonal and equatorial profiles for experiments to comapre
    zonals=[]
    equatorials=[]
    for exp,dset in zip(exps,dsets):
        zonals.append(mkzonal(exp,dset))
        equatorials.append(mkequatorial(exp,dset))

    # Plots
    cbar_kwargs={'label': 'PSU',
    }
    
    fill_opts={'yincrease': False,
               'levels': np.arange(33.0,36.1,0.2),
               'cmap': cmocean.cm.haline,
               'cbar_kwargs': cbar_kwargs}

    contour_opts={'levels': np.arange(33.0,36.1,0.4),
                  'colors': 'black',
                  'yincrease': False}

    plotter=plots.Plot2d(fill_opts=fill_opts, contour_opts=contour_opts)

    # Plot zonal mean profile
    plot_zonal(plotter, exps[0], zonals[0])

    # Plot equatorial profile
    plot_equatorial(plotter, exps[0], equatorials[0])

    # Validations
    plotter.fill_opts['levels']=np.arange(-2.0,2.1,0.2)
    plotter.fill_opts['cmap']=cmocean.cm.diff
    plotter.contour_opts['levels']=np.arange(-2.0,2.1,0.4)

    for exp,zonal in zip(exps[1:],zonals[1:]):
        plot_zonal_diff(plotter,exps[0],exp,zonals[0],zonal)

    for exp,eq in zip(exps[1:],equatorials[1:]):
        plot_equatorial_diff(plotter,exps[0],exp,equatorials[0],eq)

    obs={'woa13':'s_an'}
    for obsname,obsvarname in obs.items():
        ds=importlib.import_module('verification.'+obsname).ds_s
        var=ds[obsvarname].mean('time').rename({'depth':'lev'})
        obszonal=var.mean('lon')
        obsequatorial=utils.shift_lon(var.sel(lat=slice(-2.1,2.1)).mean('lat'),30)
        plot_zonal_diffobs(plotter,exps[0],zonals[0],obszonal,obsname)
        plot_equatorial_diffobs(plotter,exps[0],equatorials[0],obsequatorial,obsname)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])

    try:
        dsets=geosdset.load_collection(exps,'geosgcm_ocn3d')
    except OSError:
        dsets=geosdset.load_collection(exps,'prog_z',type='MOM')

    mkplots(exps,dsets)
    geosdset.close(dsets)    
