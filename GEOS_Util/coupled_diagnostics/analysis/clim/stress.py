#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
from importlib import import_module
import scipy as sp
import matplotlib.pyplot as pl
from matplotlib import colors
from mpl_toolkits.basemap.cm import s3pcpn_l, sstanom
from matplotlib import dates
import g5lib.field as field
import g5lib.plotters as ptrs
from  g5lib import cmaps as g5cmaps
from g5lib import g5dset

# Read data
exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.Ctl(exp,'geosgcm_surf')

tx=exp.ctl('TAUX',dates=exp.dates).clim(12); tx.data*=10
ty=exp.ctl('TAUY',dates=exp.dates).clim(12); ty.data*=10

exp.stress=field.cmplx(tx,ty); exp.stress.name=exp.ctl.name+' TAU'
exp.stress.shiftgrid(30.);

ind=[0,1,11]; exp.djf=exp.stress.subset(tind=ind).mean(0); exp.djf.name+=', DJF'
ind=[5,6,7]; exp.jja=exp.stress.subset(tind=ind).mean(0); exp.jja.name+=', JJA'
exp.am=exp.stress.mean(0); exp.am.name+=', Annual Mean'

# Calculate equatorial profile
exp.eqprof=exp.am(lons=(130,280), lats=(-2.1,2.1)).ave(2)

# Equatorial Annual Cycle
exp.eqac=exp.stress(lons=(130,280), lats=(-2.1,2.1)).ave(2)
exp.eqac.data-=exp.eqac.mean(0).data
exp.eqac.name=exp.stress.name+', Eq. Annual Cycle'


# Read data
exp1=g5dset.read_exp(exp.cmpexp)
exp1.ctl=g5dset.Ctl(exp1,'geosgcm_surf')

tx=exp1.ctl('TAUX',dates=exp1.dates).clim(12); tx.data*=10
ty=exp1.ctl('TAUY',dates=exp1.dates).clim(12); ty.data*=10

tx.shiftgrid(30.);
ty.shiftgrid(30.);

if tx.dims[2:] != exp.stress.dims[2:]:
    tx.regrid(exp.stress.grid)
    ty.regrid(exp.stress.grid)

exp1.stress=field.cmplx(tx,ty); exp1.stress.name=exp1.ctl.name+' TAU'

ind=[0,1,11]; exp1.djf=exp1.stress.subset(tind=ind).mean(0); exp1.djf.name+=', DJF'
ind=[5,6,7]; exp1.jja=exp1.stress.subset(tind=ind).mean(0); exp1.jja.name+=', JJA'
exp1.am=exp1.stress.mean(0); exp1.am.name+=', Annual Mean'

# Calculate equatorial profile
exp1.eqprof=exp1.am(lons=(130,280), lats=(-2.1,2.1)).ave(2)

# Equatorial Annual Cycle
exp1.eqac=exp1.stress(lons=(130,280), lats=(-2.1,2.1)).ave(2)
exp1.eqac.data-=exp1.eqac.mean(0).data
exp1.eqac.name=exp1.stress.name+', Eq. Annual Cycle'

# Read validation data set
obs=import_module('stress_mon_clim')

tx=obs.ctl.fromfile('taux',kind=0).clim(12); tx.data*=10
ty=obs.ctl.fromfile('tauy',kind=0).clim(12); ty.data*=10

tx.shiftgrid(30.);
ty.shiftgrid(30.);

if tx.dims[2:] != exp.stress.dims[2:]:
    tx.regrid(exp.stress.grid)
    ty.regrid(exp.stress.grid)

obs.stress=field.cmplx(tx,ty); obs.stress.name=obs.ctl.name+' TAU'

ind=[0,1,11]; obs.djf=obs.stress.subset(tind=ind).mean(0); obs.djf.name+=', DJF'
ind=[5,6,7]; obs.jja=obs.stress.subset(tind=ind).mean(0); obs.jja.name+=', JJA'
obs.am=obs.stress.mean(0); obs.am.name+=', Annual Mean'

# Calculate equatorial profile
obs.eqprof=obs.am(lons=(130,280), lats=(-2.1,2.1)).ave(2)

# Equatorial Annual Cycle
obs.eqac=obs.stress(lons=(130,280), lats=(-2.1,2.1)).ave(2)
obs.eqac.data-=obs.eqac.mean(0).data
obs.eqac.name=obs.stress.name+', Eq. Annual Cycle'

# Plots
path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

copts=dict(levels=(0.1,0.2,0.4,0.6,0.8,1.,1.5,2,2.5,3), cmap=s3pcpn_l)

def plot_map(figure,F,copts):
    Nq=int(F.dims[-2]/30)
    p=ptrs.GeoPlotter(copts=copts)
    x=field.absolute(F)
    pl.figure(figure); pl.clf()
    p(x)
    p.method=p.map.quiver
    p.copts={}
    p(F, skip=Nq)
    pl.show()

# DJF
plot_map(1,exp.djf,copts)
pl.savefig(path+'/tau_djf.png')

dif=exp.djf.subset(); dif.data-=exp1.djf.data
dif.name=exp.ctl.name+'-'+exp1.ctl.name+' DJF'
plot_map(2,dif,copts)
pl.savefig(path+'/tau_dif_djf.png')

dif=exp.djf.subset(); dif.data-=obs.djf.data;
dif.name=exp.ctl.name+'-'+obs.ctl.name+' DJF'
plot_map(3,dif,copts)
pl.savefig(path+'/tau-obs_djf.png')

# JJA
plot_map(1,exp.jja,copts)
pl.savefig(path+'/tau_jja.png')

dif=exp.jja.subset(); dif.data-=exp1.jja.data
dif.name=exp.ctl.name+'-'+exp1.ctl.name+' JJA'
plot_map(2,dif,copts)
pl.savefig(path+'/tau_dif_jja.png')

dif=exp.jja.subset(); dif.data-=obs.jja.data;
dif.name=exp.ctl.name+'-'+obs.ctl.name+' JJA'
plot_map(3,dif,copts)
pl.savefig(path+'/tau-obs_jja.png')

# AM
plot_map(1,exp.am,copts)
pl.savefig(path+'/tau_am.png')

dif=exp.am.subset(); dif.data-=exp1.am.data
dif.name=exp.ctl.name+'-'+exp1.ctl.name+' Annual mean'
plot_map(2,dif,copts)
pl.savefig(path+'/tau_dif_am.png')

dif=exp.am.subset(); dif.data-=obs.am.data;
dif.name=exp.ctl.name+'-'+obs.ctl.name+' Annual mean'
plot_map(3,dif,copts)
pl.savefig(path+'/tau-obs_am.png')


# Plot Equatorial Annual Cycle
pl.figure(4);pl.clf()
p=ptrs.Plotter2d(copts=dict(levels=sp.arange(-0.2,0.21,0.02), cmap=sstanom))
p.formatters['time']=dates.DateFormatter('%b')
p(exp.eqac)
del p.copts['cmap']
p.method=pl.contour
p.copts.update(colors='black',levels=sp.arange(-0.2,0.21,0.04))
p(exp.eqac)
ax=p.axis; ax.yaxis.set_major_locator(dates.MonthLocator())
ax.set_title(exp.ctl.name+' Eq. Annual cycle')
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/taux_eq_ac.png')

# Draw equatorial profile
pl.figure(5); pl.clf()
p=ptrs.Plotter1d()
p(exp.eqprof); p(exp1.eqprof); p(obs.eqprof)
ax=p.axis; ax.legend((exp.ctl.name, exp1.ctl.name, obs.ctl.name))
ax=p.axis; ax.set_ylim((-1,0.5)); ax.set_title('Equatorial TAUX')
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/taux_eq_am.png')


