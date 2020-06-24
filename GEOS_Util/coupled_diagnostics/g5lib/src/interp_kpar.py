#!/usr/bin/env python

import scipy as sp
import netCDF4 as nc
from scipy import interpolate

#gi=nc.Dataset('/discover/nobackup/yvikhlia/coupled/Forcings/a288x181_o720x410/INPUT/grid_spec.nc')
#loni=gi.variables['x_T'][:]; lati=gi.variables['y_T'][:]
#jm,im=loni.shape

loni,lati=sp.meshgrid(sp.arange(-179.875,180,0.25),sp.arange(-89.875,90,0.25))
loni[loni>79.875]-=360.
jm,im=720,1440

vlist=[]
vlist.append(('head','f4',16))
vlist.append(('h1','i4',1))
vlist.append(('kpar','f4',(jm,im)))
vlist.append(('h2','i4',1))

#a=sp.memmap('/discover/nobackup/yvikhlia/coupled/Forcings/a288x181_o720x410/SEAWIFS_KPAR_mon_clim.720x410',dtype=vlist,mode='r')
a=sp.memmap('/discover/nobackup/projects/gmao/share/dao_ops/fvInput/g5gcm/bcs/realtime/SST/1440x720/SEAWIFS_KPAR_mon_clim.1440x720',dtype=vlist,mode='r')

#go=nc.Dataset('/home/yvikhlya/nobackup/coupled/4_0_beta10/GEOSagcm/src/GEOSgcs_GridComp/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/Shared/Raster/data/MOM/1440x1080/grid_spec.nc')
go=nc.Dataset('/home/yvikhlia/nobackup/coupled/Forcings/MOM6/CF0090x6C_TM0360xTM0210/INPUT/ocean_hgrid.nc')
#lono=go.variables['x_T'][:]; lato=go.variables['y_T'][:]
lono=go.variables['x'][1::2,1::2]; lato=go.variables['y'][1::2,1::2]
jm,im=lono.shape

vlist=[]
vlist.append(('head','f4',16))
vlist.append(('h1','i4',1))
vlist.append(('kpar','f4',(jm,im)))
vlist.append(('h2','i4',1))
b=sp.zeros(14,dtype=vlist)

for i,x in enumerate(a):
    print i
    xx=x['kpar'][:]
    b[i]['head'][:]=x['head'][:]
    b[i]['h1']=im*jm*4
    b[i]['kpar'][:]=interpolate.griddata(zip(loni.flatten(),lati.flatten()),xx.flatten(),(lono,lato))
    b[i]['h2']=im*jm*4

b['head'][:,13]=im; b['head'][:,14]=jm
b['kpar'][sp.isnan(b['kpar'])]=1.0
b['kpar'][:,1060:,:]=1.0
b.tofile('SEAWIFS_KPAR_mon_clim.360x210')
