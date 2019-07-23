#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import os, sys
import scipy as sp
import matplotlib.pyplot as pl
from mpl_toolkits.basemap.cm import sstanom
from g5lib import plotters, g5dset

varname='U'

# Read variable
exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.Ctl(exp,'geosgcm_ocn3d')
if exp.ctl.grid['lev'][-1] < 0.0:
   exp.ctl.grid['lev'][:]*=-1 
exp.am=exp.ctl(varname, dates=exp.dates, lats=(0.,), levs=(0.,500.)).mean(0)
exp.am.shiftgrid(30.)

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass

pl.figure(1); pl.clf()
pp=plotters.Plotter2d()
pp.copts.update(levels=sp.arange(-1.,1.1,0.2), cmap=sstanom)
pp(exp.am)
ax=pl.gca(); ax.set_ylabel('depth, m'); ax.invert_yaxis()
pl.grid(); pl.tight_layout(); pl.show()
pl.savefig(path+'/'+varname+'_eq_depth.png')
