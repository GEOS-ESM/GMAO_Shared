#!/bin/env python

import os
import scipy as sp
import matplotlib.pyplot as pl
from matplotlib import ticker
import my_utils as utl
import my_plots as mpl

print 'Running SSS...'

# Read Levitus annual mean salinity
path=os.environ['NOBACKUP']+'/verification/levitus'
execfile(path+'/ctl.py')

name='Levitus S'

var=ctl.fromfile('salt').ave(0)
var.shiftgrid(30.)
var.grid['lon']=sp.where(var.grid['lon']<29.,var.grid['lon']+360.,var.grid['lon'])

# Latitude-depth
lat_depth=var.ave(3); lat_depth.name=name+', Annual Mean'

# Equatorial depth profile
#lonind=sp.logical_and(var.grid['lon'][0]>=130.0,var.grid['lon'][0]<=280.0)
lonind=slice(None,None)
latind=sp.logical_and(var.grid['lat'][:,0]>=-1.,var.grid['lat'][:,0]<=1.)
eq_depth=var.subset(iind=lonind,jind=latind).ave(2)
eq_depth.name='Eq. '+name+', Annual Mean'

print '...done'
###################### Do plots #######################################################
path=os.environ['NOBACKUP']+'/verification/levitus/pics'

clevs=sp.arange(33.,36.1,0.2)

pl.figure(1)
pl.clf()
lat_depth.copts={'levels' : clevs}
lat_depth.plot2d(); lat_depth.copts.clear()
lat_depth.copts={'levels' : clevs,\
                 'colors' : 'black',\
                 'func'   : pl.contour
                 }
lat_depth.plot2d()
ax=pl.gca(); ax.set_ylim(0.,3000.); ax.invert_yaxis(); ax.set_ylabel('depth, m')
ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
pl.grid(); pl.show()
pl.savefig(path+'/'+var.name+'_lat_depth.png')

clevs=sp.arange(34.,36.1,0.2)

pl.figure(2)
pl.clf()
eq_depth.copts={'levels' : clevs}
eq_depth.plot2d(); eq_depth.copts.clear()
eq_depth.copts={'levels' : clevs,\
                'colors' : 'black',\
                'func'   : pl.contour
                }
eq_depth.plot2d()
ax=pl.gca(); ax.set_ylim(0.,500.); ax.invert_yaxis(); ax.set_ylabel('depth, m')
pl.grid(); pl.show()
pl.savefig(path+'/'+var.name+'_eq_depth.png')
