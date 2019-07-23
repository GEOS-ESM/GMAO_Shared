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

exp.gm=exp.ctl(varname).aave(); exp.gm.data-=273.16

# Plot

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

pl.clf()
exp.gm.name=exp.ctl.name +' Global SST'
p=ptrs.Plotter1d()
p(exp.gm)
ax=p.axis
if(exp.gm.time.size > 2500):
    myloc=dates.YearLocator((exp.gm.time.size/1000)*10)
    ax.xaxis.set_major_locator(myloc)
ax.set_ylim((16,22)); ax.set_ylabel('T, $^0C$')
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/sst_gm.png')


