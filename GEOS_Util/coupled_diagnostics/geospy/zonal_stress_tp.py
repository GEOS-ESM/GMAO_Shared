#!/usr/bin/env python3

'''
Docstring
'''
import sys, importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
import geosdset, utils, plots

plotname='Zonal_stress'
defaults={'name': 'TAUX', 
          'colname': 'geosgcm_ocn2dT', 
          'coltype': 'GEOSTripolar'}

def mkequatorial(exp,ds):
    vardata=exp['plots'].get(plotname,defaults)
    varname=vardata['name']

    var=ds[varname].sel(time=slice(*exp['dates']))
    wght=ds['mask']*ds['dy']
    eq_region={'y':slice(-2.1,2.1)}
    var=utils.average(var.sel(**eq_region),'y',wght.sel(**eq_region))
    return utils.shift_lon(var,30,'x').sel(x=slice(130,280))

def mkequatorial_obs(obsname,obsvarname):
    var=importlib.import_module('verification.'+obsname).ds[obsvarname]
    var=var.sel(lat=slice(-2.1,2.1)).mean('lat')
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
    ax.set_title(f'TAUX equatorial')
    ax.xaxis.set_major_formatter(plots.LONGITUDE_FORMATTER)
    ax.set_ylim(-.1,.05)
    ax.set_ylabel('$N/m^2$')
    ax.set_xlabel('')
    pl.grid()
    pl.tight_layout()
    pl.savefig(f'{exps[0]["plot_path"]}/taux_eq_am.png')

def plot_equatorial_ac(plotter,exp,da):
    ac=da.groupby('time.month').mean('time')-da.mean('time')

    pl.figure(2); pl.clf()
    ax=plotter.contour(ac)
    ax.set_title(f'{exp["expid"]} SST, eq. annual cycle anom.')
    ax.xaxis.set_major_formatter(plots.LONGITUDE_FORMATTER)
    ax.set_xlabel('')
    ax.yaxis.set_major_locator(plots.MONTH_LOCATOR)        
    ax.yaxis.set_major_formatter(plots.MONTH_FORMATTER)    
    ax.set_ylabel('')
    pl.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/taux_eq_ac.png')

def mkplots(exps,dsets):
    # Calculate equatorial means
    equatorials=[]
    for exp,dset in zip(exps,dsets):
        equatorials.append(mkequatorial(exp,dset))

    obs={'QSCAT': 'taux'}
    obsequatorials=[]
    for obsname,obsvarname in obs.items():
        obsequatorials.append(mkequatorial_obs(obsname,obsvarname))

    # Plots
    plot_equatorial(exps,equatorials,obs,obsequatorials)

    cbar_kwargs={'label': '$N/m^2$'}
    
    fill_opts={'cmap': cmocean.cm.diff, 
               'levels': np.arange(-0.02,0.021,0.002),
               'cbar_kwargs': cbar_kwargs
    }

    contour_opts={'levels': np.arange(-0.02,0.021,0.004),
                  'colors': 'black'
    }

    clab_opts={'fmt': '%1.3f'}

    plotter=plots.Plot2d(fill_opts=fill_opts, contour_opts=contour_opts, clab_opts=clab_opts)

    plot_equatorial_ac(plotter, exps[0],equatorials[0])

def main(exps):
    dsets=geosdset.load_data(exps, plotname, defaults)
    mkplots(exps,dsets)
    geosdset.close(dsets)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    main(exps)
