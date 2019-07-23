#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
import matplotlib.pyplot as pl
from matplotlib import dates
import numpy as np
from numpy import ma
import g5lib.plotters as ptrs
from g5lib import g5dset

exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.Ctl(exp,'geosgcm_ocn3dT')
exp.area=g5dset.Ctl(exp,'geosgcm_ocn2dT')
exp.gvm=exp.ctl.subset(iind=0,jind=0,kind=0)

for tt in xrange(exp.gvm.time.size):
    exp.A  = exp.area.fromfile('AREA',tind=tt)
    exp.T  = exp.ctl.fromfile('T',tind=tt)
    exp.DH = exp.ctl.fromfile('DH',tind=tt)
    xx=(exp.T*exp.DH).aint(area=exp.A).sum(1)  
    yy=exp.DH.aint(area=exp.A).sum(1)
    if xx.data[0,0,0,0] is ma.masked:
        exp.gvm.data[tt] = ma.masked
    else:
        exp.gvm.data[tt] = xx.data[0] / yy.data[0]


exp.gvm.data-=273.16
# Plot

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

pl.clf()
exp.gvm.name=exp.ctl.name +' Global Volume Mean Temp'
p=ptrs.Plotter1d()
p(exp.gvm)
ax=p.axis
if(exp.gvm.time.size > 2500):
    myloc=dates.YearLocator((exp.gm.time.size/1000)*10)
    ax.xaxis.set_major_locator(myloc)
ax.set_ylim((3.6, 4.0)); ax.set_ylabel('T, $^0C$')
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/temp_gvm.png')


