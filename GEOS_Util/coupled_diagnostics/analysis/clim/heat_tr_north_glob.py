#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys, pickle
import netCDF4 as nc
import matplotlib.pyplot as pl
import scipy as sp
from matplotlib import dates, ticker
from mpl_toolkits.basemap.cm import sstanom
import g5lib.bilinear_scale as bscale
import g5lib.plotters as ptrs
from g5lib import g5dset

exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.CtlMOM(exp)

tind=exp.ctl.index(dates=exp.dates)['tind']
ty=exp.ctl.subset(kind=0,tind=tind)
for tt in xrange(tind.start,tind.stop):
    xx=exp.ctl.fromfile('temp_yflux_adv',tind=tt)
    ty.data[tt-tind.start,0,:,:]=xx.data.sum(1)

var=ty.clim(12).mean(0)
var.data/=1e15; var.name=exp.ctl.name+' Northward Heat Transport'

# Global transport
exp.htrans_g=var(lats=(-80,65)).sum(3);

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

pl.figure(1); pl.clf()
p=exp.htrans_g.d()
ax=p.axis; ax.set_ylim(-3,3)
ax.set_ylabel('PW')
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/htrans.png')




