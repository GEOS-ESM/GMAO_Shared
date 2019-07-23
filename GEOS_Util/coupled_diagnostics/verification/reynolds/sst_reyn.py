#!/bin/env python

import os
import scipy as sp
import matplotlib.pyplot as pl
from matplotlib import dates
import my_plots as mpl
from mpl_toolkits.basemap.cm import sstanom
from matplotlib.cm import jet

# Read Reynolds SST climatology
path=os.environ['NOBACKUP']+'/verification/reynolds'
execfile(path+'/ctl.py')

obs={}
obs['name']='Reynolds SST'
obs['ctl']=ctl

# Calculate climatology
obs['clim']=obs['ctl'].fromfile('sst',tind=slice(1,None)).clim(12)
obs['clim'].shiftgrid(30.)
obs['clim'].grid['lon']=sp.where(obs['clim'].grid['lon']<29.,obs['clim'].grid['lon']+360.,obs['clim'].grid['lon'])

# Calculate DJF, JJA and annual mean
obs['djf']=obs['clim'].subset(tind=[0,1,11]).ave(0); obs['djf'].name+=', DJF'
obs['jja']=obs['clim'].subset(tind=[5,6,7]).ave(0); obs['jja'].name+=', JJA'
obs['am']=obs['clim'].ave(0)
obs['am'].name+=', Annual Mean'


# Equatorial annual cycle

lonind=sp.logical_and(obs['clim'].grid['lon'][0]>=130.0,obs['clim'].grid['lon'][0]<=280.0)
latind=sp.logical_and(obs['clim'].grid['lat'][:,0]>=-2.1,obs['clim'].grid['lat'][:,0]<=2.0)
obs['eqac']=obs['clim'].subset(iind=lonind,jind=latind).ave(2)
obs['eqac'].data-=obs['eqac'].ave(0).data
obs['eqac'].name=obs['clim'].name+', Eq. Annual Cycle'

######################### Do plots ######################################################
path=os.environ['NOBACKUP']+'/verification/reynolds/pics'

def plot_field(field, fig, clevs, cmap):
        pl.figure(fig); pl.clf()
        field.copts={'levels': clevs,\
                     'cmap': cmap}
        field.plot_map()
        field.copts={'levels': clevs[0::2],\
                     'func': field.map.contour,\
                     'colors': 'black'}
        field.plot_map()
        mean,std=field.aave(ret_std=True)
        mpl.draw_stat((mean.data.squeeze(),std.data.squeeze()))
        pl.show()

clevs1=sp.arange(0.0,32.0,2.0); clevs2=sp.arange(-5.,5.1,.5)
    
# DJF
plot_field(obs['djf'],1,clevs1,jet)
pl.savefig(path+'/sst_lev_djf.png')

# JJA
plot_field(obs['jja'],2,clevs1,jet)
pl.savefig(path+'/sst_lev_jja.png')

# Annual Mean
plot_field(obs['am'],3,clevs1,jet)
pl.savefig(path+'/sst_lev_am.png')

# Plot annual cycle
pl.figure(4);pl.clf()
obs['eqac'].copts={'levels': sp.arange(-2.4,2.5,0.3),\
                   'cmap' : sstanom,\
                   'timefmt': dates.DateFormatter('%b')}
obs['eqac'].plot2d()
obs['eqac'].copts={'func': pl.contour,\
                   'colors': 'black',\
                   'levels': sp.arange(-2.4,2.5,0.3),\
                   'timefmt': dates.DateFormatter('%b')}
obs['eqac'].plot2d()
ax=pl.gca(); ax.yaxis.set_major_locator(dates.MonthLocator())
pl.grid(); pl.show()
    
pl.savefig(path+'/sst_lev_eq_ac.png')


