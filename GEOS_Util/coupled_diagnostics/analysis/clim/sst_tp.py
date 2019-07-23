#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
from importlib import import_module
import scipy as sp
import matplotlib.pyplot as pl
from mpl_toolkits.basemap.cm import sstanom
from matplotlib.cm import jet
from matplotlib import dates
import g5lib.plotters as ptrs
from g5lib import g5dset

exp=g5dset.read_exp(sys.argv[1])
exp1=g5dset.read_exp(exp.cmpexp)

varname='TS'

# Read data
exp.ctl=g5dset.Ctl(exp,'geosgcm_ocn2d')
exp.sst=exp.ctl(varname, lats=(-2, 2), levs=(0,), dates=exp.dates).ave(2); exp.sst.data-=273.16
exp.sst.shiftgrid(30.);
exp.clim=exp.sst.clim(12,anom=True)

# Calculate mean equatorial profile
exp.mean=exp.clim(lons=(130, 280)).mean(0)

# Calculate equatorial profile of std
exp.std=exp.sst(lons=(130, 280)).mean(0,ret_std=True)[1]

# Equatorial Annual Cycle
exp.clim.data-=exp.clim.mean(0).data
exp.eqac=exp.clim(lons=(130, 280))
exp.eqac.name=exp.clim.name+', Eq. Annual Cycle'

# Read cmp data 
exp1.ctl=g5dset.Ctl(exp1,'geosgcm_ocn2d')
exp1.sst=exp1.ctl(varname, lats=(-2, 2), levs=(0,), dates=exp1.dates).ave(2); exp1.sst.data-=273.16
exp1.sst.shiftgrid(30.);
exp1.clim=exp1.sst.clim(12,anom=True)

# Calculate mean equatorial profile
exp1.mean=exp1.clim(lons=(130, 280)).mean(0)

# Calculate equatorial profile of std
exp1.std=exp1.sst(lons=(130, 280)).mean(0,ret_std=True)[1]

# Read validation data set
obs=import_module('reynolds')
obs.sst=obs.ctl('sst',lats=(-2, 2), levs=(0,)).ave(2)
obs.clim=obs.sst.clim(12,anom=True)

# Calculate mean equatorial profile
obs.mean=obs.clim(lons=(130, 280)).mean(0)

# Calculate equatorial profile of std
obs.std=obs.sst(lons=(130, 280)).mean(0,ret_std=True)[1]

# Plots

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

# Plot Equatorial Annual Cycle
pl.figure(1);pl.clf()
p2d=ptrs.Plotter2d(copts=dict(levels=sp.arange(-2.4,2.5,0.3),
                            cmap=sstanom))
p2d.formatters['time']=dates.DateFormatter('%b')
p2d(exp.eqac)
del p2d.copts['cmap']
p2d.method=pl.contour
p2d.copts.update(colors='black')
p2d(exp.eqac)
p2d.axis.yaxis.set_major_locator(dates.MonthLocator())
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/sst_eq_ac.png')
    
# Draw equatorial profile
pl.figure(2); pl.clf()
p1d=ptrs.Plotter1d()
p1d(exp.mean); p1d(exp1.mean); p1d(obs.mean); 
ax=p1d.axis; ax.set_ylim((20,33)) 
ax.legend((exp.ctl.name, exp1.ctl.name, obs.ctl.name)); ax.set_title('Equatorial SST')
ax.set_ylabel('$^0$C')
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/sst_eq_am.png')

# Draw equatorial profile of std
pl.figure(3); pl.clf()
p1d(exp.std); p1d(exp1.std); p1d(obs.std); 
ax=p1d.axis; ax.set_ylim((0,2))
ax.legend((exp.ctl.name, exp1.ctl.name, obs.ctl.name), loc=4); ax.set_title('Equatorial SST std.')
ax.set_ylabel('$^0$C')
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/sst_eq_std.png')


