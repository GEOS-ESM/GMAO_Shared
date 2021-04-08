#!/usr/bin/env python3

'''
Plots tropical Pacific SST
'''
import sys, importlib
import numpy as np
import matplotlib.pyplot as pl
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER
import geosdset, plots, utils

def mkequatorial(exp,ds):
    varname='TS'
    TFREEZE=273.16
    var=ds[varname].sel(time=slice(*exp['dates']))-TFREEZE
    mask=1.0-np.isnan(var)
    wght=mask*ds.dy
    var=utils.average(var.sel(lat=slice(-2.1,2.1)),'lat',wght.sel(lat=slice(-2.1,2.1)))
    return utils.shift_lon(var,30).sel(lon=slice(130,280))

def mkequatorial_obs(obsname,obsvarname):
    var=importlib.import_module('verification.'+obsname).ds[obsvarname]
    var=var.sel(lat=slice(2.1,-2.1)).mean('lat')
    return utils.shift_lon(var,30).sel(lon=slice(130,280))

def plot_equatorial(exps,das,obs,obsdas):
    plotter=plots.Plot1d()

    # Plot equatorial means
    pl.figure(1); pl.clf()
    ax=pl.gca()

    for exp,da in zip(exps,das):
        plotter.line_opts['label']=exp['expid']
        plotter.line(da.mean('time'),ax=ax)
    
    for obsname,obsda in zip(obs,obsdas):
        plotter.line_opts['label']=obsname
        plotter.line(obsda.mean('time'),'--',ax=ax)

    ax.legend()
    ax.set_title(f'SST equatorial')
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.set_ylim(20,33)
    ax.set_ylabel('$^0C$')
    ax.set_xlabel('')
    pl.grid()
    pl.tight_layout()
    pl.show()
    pl.savefig(f'{exps[0]["plot_path"]}/sst_eq_am.png')

    # Plot equatorial standard deviations
    pl.figure(2); pl.clf()
    ax=pl.gca()
    for exp,da in zip(exps,das):
        plotter.line_opts['label']=exp['expid']
        mean=da.mean('time')
        std=np.sqrt(((da-mean)**2).mean('time'))
        plotter.line(std,ax=ax)

    for obsname,obsda in zip(obs,obsdas):
        plotter.line_opts['label']=obsname
        mean=obsda.mean('time')
        std=np.sqrt(((obsda-mean)**2).mean('time'))
        plotter.line(std,'--',ax=ax)

    ax.legend()
    ax.set_title(f'SST std, equatorial')
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.set_ylim(0,3)
    ax.set_ylabel('$^0C$')
    ax.set_xlabel('')
    pl.grid()
    pl.tight_layout()
    pl.show()
    pl.savefig(f'{exps[0]["plot_path"]}/sst_eq_std.png')

def mkplots(exps,dsets):
    # Calculate seasonal equatorial means
    equatorials=[]
    for exp,dset in zip(exps,dsets):
        equatorials.append(mkequatorial(exp,dset))

    obs={'OISSTv2': 'sst'}
    obsequatorials=[]
    for obsname,obsvarname in obs.items():
        obsequatorials.append(mkequatorial_obs(obsname,obsvarname))

    # Plots
    plot_equatorial(exps,equatorials,obs,obsequatorials)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    mkplots(exps,dsets)
