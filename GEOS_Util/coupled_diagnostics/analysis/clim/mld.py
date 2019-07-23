#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
import scipy as sp
import matplotlib.pyplot as pl
from matplotlib.cm import coolwarm as g5cmap
from g5lib import plotters, g5dset

# Read data
exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.CtlMOM(exp)
exp.mld=exp.ctl('mld',dates=exp.dates,levs=(0,)).clim(12)

ind=[1,2,3]
exp.mld.jfm=exp.mld.subset(tind=ind).mean(0); exp.mld.jfm.name+=', JFM'

ind=[7,8,9]
exp.mld.jas=exp.mld.subset(tind=ind).mean(0); exp.mld.jas.name+=', JAS'

# Plots

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

clevs=sp.array([5, 10, 20, 50, 100, 150, 200, 300])

pp=plotters.GeoPlotter()
pp.copts.update(levels=clevs, cmap=g5cmap)

pl.figure(1); pl.clf()
pp(exp.mld.jfm)
pl.show()
pl.savefig(path+'/mld_jfm.png')

pl.figure(2); pl.clf()
pp(exp.mld.jas)
pl.show()
pl.savefig(path+'/mld_jas.png')

