#!/bin/env python

import matplotlib.pyplot as pl
from matplotlib import dates, ticker
import my_utils as utl
import scipy as sp

def moc_stream(field):
    km=field.grid['lev'].size
    moc=field.subset(iind=0)
    var=field.gint(3)
    for k in range(km)[-2::-1]:
        kind=slice(km-1,k-1,-1)
        moc.data[:,k,:,:]=var.subset(kind=kind).gint(1).data/1e6

    return moc

execfile('ctl.py')
tind=sp.arange(ctl.time.size)
#tind=sp.arange(10)
# Atlantic MOC (two parts, north and south)
lon=sp.where(ctl.grid['lon']>180,ctl.grid['lon']-360.,ctl.grid['lon'])
iind_north=sp.where(sp.logical_and(lon[0]>=-100.,lon[0]<=10))[0]
jind_north=sp.where(ctl.grid['lat'][:,0]>=20.)[0]
moc_north=ctl.subset(iind=0,tind=tind,jind=jind_north)

iind_mid=sp.where(sp.logical_and(lon[0]>=-90.,lon[0]<=0.))[0]
jind_mid=sp.where(sp.logical_and(ctl.grid['lat'][:,0]>10.,ctl.grid['lat'][:,0]<20.))[0]
moc_mid=ctl.subset(iind=0,tind=tind,jind=jind_mid)

iind_south=sp.where(sp.logical_and(lon[0]>=-70.,lon[0]<=20))[0]
jind_south=sp.where(sp.logical_and(ctl.grid['lat'][:,0]<=10.,ctl.grid['lat'][:,0]>=-30.))[0]
moc_south=ctl.subset(iind=0,tind=tind,jind=jind_south)

for t in xrange(tind.size):
    print t
    x=ctl.fromfile('vo',tind=tind[t])
    moc_north.data[t]=moc_stream(x.subset(iind=iind_north,jind=jind_north)).data
    moc_mid.data[t]=moc_stream(x.subset(iind=iind_mid,jind=jind_mid)).data
    moc_south.data[t]=moc_stream(x.subset(iind=iind_south,jind=jind_south)).data


moc=ctl.subset(tind=0, iind=0)
moc.data[:,:,jind_north,:]=moc_north.ave(0).data
moc.data[:,:,jind_mid,:]=moc_mid.ave(0).data
moc.data[:,:,jind_south,:]=moc_south.ave(0).data
moc=moc.subset(jind=sp.concatenate((jind_south,jind_mid,jind_north)))
moc.name+=' Atl. MOC Stream Function'

x=moc
pl.figure(1); pl.clf();
x.copts={'levels': sp.arange(-30,31,3),\
         'func': pl.contourf}
x.cbar_opts={'orientation' : 'horizontal'}
x.plot2d()
x.copts={'levels': sp.arange(-30,31,3),\
         'func': pl.contour,\
         'colors': 'black'}
x.plot2d()
ax=pl.gca(); ax.invert_yaxis(); ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
pl.grid(); pl.show()
pl.savefig('pics/atl_moc.png')

