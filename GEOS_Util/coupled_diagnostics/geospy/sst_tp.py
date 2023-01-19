#!/usr/bin/env python3

'''
Plots tropical Pacific SST
'''
import sys, importlib
import numpy as np
import matplotlib.pyplot as pl
import cmocean
import geosdset, plots, utils

plotname='SST'
defaults={'name': 'TS', 
          'colname': 'geosgcm_ocn2dT', 
          'coltype': 'GEOSTripolar'}

TFREEZE=273.16

def mkequatorial(exp,ds):
    vardata=exp['plots'].get(plotname,defaults)
    varname=vardata['name']

    var=ds[varname].sel(time=slice(*exp['dates']))-TFREEZE
    wght=ds['mask']*ds['dy']
    eq_region={'y':slice(-2.1,2.1)}
    var=utils.average(var.sel(**eq_region),'y',wght.sel(**eq_region))
    return utils.shift_lon(var,30,'x').sel(x=slice(130,280))

def mkequatorial_obs(obsname,obsvarname):
    var=importlib.import_module('verification.'+obsname).ds[obsvarname]
    var=var.sel(lat=slice(2.1,-2.1)).mean('lat')
    return utils.shift_lon(var,30).sel(lon=slice(130,280))

def mknino3(exp,ds):
    vardata=exp['plots'].get(plotname,defaults)
    varname=vardata['name']

    var=ds[varname]-TFREEZE
    wght=ds['mask']*ds['area']
    n3region={'y':slice(-5.1,5.1),
              'x':slice(-150,-90)}
    return utils.average(var.sel(**n3region),('y','x'),
                         wght.sel(**n3region))

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
    ax.xaxis.set_major_formatter(plots.LONGITUDE_FORMATTER)
    ax.set_ylim(20,33)
    ax.set_ylabel('$^0C$')
    ax.set_xlabel('')
    pl.grid()
    pl.tight_layout()
    pl.savefig(f'{exps[0]["plot_path"]}/sst_eq_am.png')

    # Plot equatorial standard deviations
    pl.figure(2); pl.clf()
    ax=pl.gca()
    for exp,da in zip(exps,das):
        plotter.line_opts['label']=exp['expid']
        plotter.line(da.std('time'),ax=ax)

    for obsname,obsda in zip(obs,obsdas):
        plotter.line_opts['label']=obsname
        plotter.line(obsda.std('time'),'--',ax=ax)

    ax.legend()
    ax.set_title(f'SST std, equatorial')
    ax.xaxis.set_major_formatter(plots.LONGITUDE_FORMATTER)
    ax.set_ylim(0,3)
    ax.set_ylabel('$^0C$')
    ax.set_xlabel('')
    pl.grid()
    pl.tight_layout()
    pl.savefig(f'{exps[0]["plot_path"]}/sst_eq_std.png')

def plot_equatorial_ac(plotter,exp,da):
    ac=da.groupby('time.month').mean('time')-da.mean('time')

    pl.figure(3); pl.clf()
    ax=plotter.contour(ac)
    ax.set_title(f'{exp["expid"]} SST, eq. annual cycle anom.')
    ax.xaxis.set_major_formatter(plots.LONGITUDE_FORMATTER)
    ax.set_xlabel('')
    ax.yaxis.set_major_locator(plots.MONTH_LOCATOR)        
    ax.yaxis.set_major_formatter(plots.MONTH_FORMATTER)    
    ax.set_ylabel('')
    pl.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/sst_eq_ac.png')

def plot_hovm(plotter,exp,da):
    anom=da.groupby('time.month')-da.groupby('time.month').mean('time')
    pl.figure(4); pl.clf()
    ax=plotter.contour(anom)
    ax.set_title(f'{exp["expid"]} SST anom.')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.xaxis.set_major_formatter(plots.LONGITUDE_FORMATTER)
    pl.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/hov_tp.png')

def plot_nino3(exp,da):
    anom=da.groupby('time.month')-da.groupby('time.month').mean('time')
    plotter=plots.Plot1d()
    pl.figure(5); pl.clf()
    ax=plotter.line(anom)
    ax.set_title('Nino3 SST anomaly')
    ax.set_ylabel('$^0C$')
    pl.grid()
    pl.tight_layout()
    pl.savefig(f'{exp["plot_path"]}/n3.png')

def mkplots(exps,dsets):
    # Calculate equatorial means
    equatorials=[]
    for exp,dset in zip(exps,dsets):
        equatorials.append(mkequatorial(exp,dset))

    obs={'OISSTv2': 'sst'}
    obsequatorials=[]
    for obsname,obsvarname in obs.items():
        obsequatorials.append(mkequatorial_obs(obsname,obsvarname))

    nino3=mknino3(exps[0],dsets[0])

    # Plots
    plot_equatorial(exps,equatorials,obs,obsequatorials)

    cbar_kwargs={'label': '$^0C$'}
    
    fill_opts={'cmap': cmocean.cm.diff, 
              'levels': np.arange(-2.4,2.5,0.3),
               'cbar_kwargs': cbar_kwargs
    }

    contour_opts={'levels': np.arange(-2.4,2.5,0.6),
                  'colors': 'black'
    }

    plotter=plots.Plot2d(fill_opts=fill_opts, contour_opts=contour_opts)

    plot_equatorial_ac(plotter, exps[0],equatorials[0])

    plot_hovm(plotter, exps[0],equatorials[0]) 

    plot_nino3(exps[0],nino3)

def main(exps):
    dsets=geosdset.load_data(exps, plotname, defaults)
    mkplots(exps,dsets)
    geosdset.close(dsets)

if __name__=='__main__':
    exps=geosdset.load_exps(sys.argv[1])
    main(exps)
