#!/usr/bin/env python

import numpy as np
import netCDF4 as nc
from scipy import interpolate

loni,lati=np.meshgrid(np.arange(-179.875,180,0.25),np.arange(-89.875,90,0.25))
loni[loni>79.875]-=360.
jm,im=720,1440

vlist=[]
vlist.append(('head','f4',16))
vlist.append(('h1','i4'))
vlist.append(('kpar','f4',(jm,im)))
vlist.append(('h2','i4'))

a=np.memmap('/discover/nobackup/projects/gmao/share/dao_ops/fvInput/g5gcm/bcs/realtime/SST/1440x720/SEAWIFS_KPAR_mon_clim.1440x720',dtype=vlist,mode='r')

go=nc.Dataset('/home/yvikhlia/nobackup/coupled/Forcings/MOM6/CF0090x6C_TM0360xTM0320/MAPL_Tripolar.nc')
lono=go.variables['lon_centers'][:]; lato=go.variables['lat_centers'][:]
jm,im=lono.shape

vlist=[]
vlist.append(('head','f4',16))
vlist.append(('h1','i4'))
vlist.append(('kpar','f4',(jm,im)))
vlist.append(('h2','i4'))
b=np.zeros(14,dtype=vlist)

for i,x in enumerate(a):
    print(i)
    xx=x['kpar'][:]
    b[i]['head'][:]=x['head'][:]
    b[i]['h1']=im*jm*4
    b[i]['kpar'][:]=interpolate.griddata(list(zip(loni.flatten(),lati.flatten())),xx.flatten(),(lono,lato), 
                                         method='nearest', fill_value=1.0)
    b[i]['h2']=im*jm*4

b['head'][:,13]=im; b['head'][:,14]=jm
b.tofile('SEAWIFS_KPAR_mon_clim.360x320')
