import os
import scipy as sp
import netCDF4 as nc
from g5lib import dset

__all__=['ctl']

oceanval=os.environ.get('OCEANVAL',
                        '/discover/nobackup/projects/gmao/oceanval/verification')

class Ctl(dset.NCDset):
    def __init__(self):
        name='Reynolds'
        
        flist=[oceanval+'/reynolds_sst/sst_NOAA_OI_v2.nc']

        f=nc.Dataset(flist[0])
        time=f.variables['time']
        time=sp.array(nc.num2date(time[:],time.units),dtype='|O')
        f.close()

        super(Ctl,self).__init__(flist,time=time,name=name)

    def fromfile(self,varname,iind=slice(None),jind=slice(None),kind=slice(None),tind=slice(None), dtype=sp.float32):

        ii,jj,kk,tt=dset.scalar2slice(iind,jind,kind,tind)
        
        var=super(Ctl,self).fromfile(varname,iind=iind,jind=jind,kind=kind,tind=tind,dtype=dtype,maskandscale=False)
        
        # Applly land mask
        data=var.data*0.01
        fmask=nc.Dataset(oceanval+'/reynolds_sst/lsmask.nc')
        x=fmask.variables['mask'][0][jj][:,ii]; mask=sp.zeros(data.shape); mask[:]=x
        var.data=sp.ma.masked_where(mask==0.0,data)
        del fmask

        # Flip north-south
        var.grid['lat']=var.grid['lat'].copy()[-1::-1]
        var.data=var.data.copy()[:,:,-1::-1,:]
        
        return var


ctl=Ctl()
