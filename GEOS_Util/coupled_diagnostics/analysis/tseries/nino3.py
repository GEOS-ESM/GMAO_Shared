#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
import matplotlib.pyplot as pl
from matplotlib import dates
import g5lib.plotters as ptrs
from g5lib import g5dset

exp=g5dset.read_exp(sys.argv[1])
varname='TS'
exp.ctl=g5dset.Ctl(exp,'geosgcm_ocn2d')

exp.n3=exp.ctl(varname,lats=(-5,5),lons=(-150,-120)).aave();
exp.n3.clim(12,anom=True)

# Plot

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

pl.clf()
exp.n3.name=exp.ctl.name +' Nino3 SST'
p=ptrs.Plotter1d()
p(exp.n3)
ax=p.axis
if(exp.n3.time.size > 2500):
    myloc=dates.YearLocator((exp.n3.time.size/1000)*10)
    ax.xaxis.set_major_locator(myloc)
ax.set_ylim((-5,5)); ax.set_ylabel('T anomaly, $^0C$')
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/n3.png')


