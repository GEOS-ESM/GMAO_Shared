#!/bin/env python

import os
import scipy as sp
import matplotlib.pyplot as pl
from matplotlib import ticker

# Read variable
execfile('ctl.py')

iind=300
s=ctl.fromfile('salt',iind=iind).ave(0)
s.name='S at 60W'

###################### Do plots #######################################################
clevs=sp.arange(33.,36.1,0.2)

pl.figure(1)
pl.clf()
s.copts={'func': pl.contourf,\
         'levels' : clevs,\
         }
s.plot2d(); s.copts.clear()
s.copts={'levels' : clevs[0::2],\
         'colors' : 'black',\
         'func': pl.contour
         }
s.plot2d()
ax=pl.gca(); ax.set_ylim(0.,3000.); ax.invert_yaxis(); ax.set_ylabel('depth, m')
ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
pl.grid(); pl.show()
pl.savefig('pics/s_profile/s_60W.png')


