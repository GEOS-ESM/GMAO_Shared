#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
import scipy as sp
import matplotlib.pyplot as pl
from mpl_toolkits.basemap.cm import sstanom
from matplotlib import colors
import g5lib.plotters as ptrs
from  g5lib import cmaps as g5cmaps
from g5lib import g5dset

exp=g5dset.read_exp(sys.argv[1])
varname='TS'
exp.ctl=g5dset.Ctl(exp,'geosgcm_ocn2d')

exp.sst=exp.ctl(varname,lats=(-2, 2), levs=(0,), dates=exp.dates).ave(2); exp.sst.data-=273.16
exp.sst.clim(12,anom=True)

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

pl.clf()
n=colors.Normalize()
clevs=sp.arange(-2.4,2.5,0.3)
n.autoscale(clevs)
m=g5cmaps.FilledCmap(sstanom, fill_range=n((-0.3,0.3)))
p=ptrs.Plotter2d(copts=dict(levels=sp.arange(-2.4,2.5,0.3), cmap=m, norm=n))
p(exp.sst); pl.grid(); pl.tight_layout(); pl.show() 
pl.savefig(path+'/hov_tp.png')

