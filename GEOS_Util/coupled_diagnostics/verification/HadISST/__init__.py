from g5lib import dset
import netCDF4 as nc
import scipy as sp

__all__=['ctl']

class Ctl(dset.NCDset):
    def __init__(self):
        name='HadISST'

        flist=['/gpfsm/dnb42/projects/p16/ssd/ocean/kovach/odas-2/obs/HADISST/HadISST_sst.nc']
        f=nc.Dataset(flist[0])
        tt=f.variables['time']
        time=nc.num2date(tt[:],'hours since 0001-01-01 00:00:00')
        undef=f.variables['sst'].missing_value
        f.close()

        super(Ctl,self).__init__(flist,time=time,undef=undef,name=name)
        
    def fromfile(self,varname,iind=slice(None),jind=slice(None),kind=slice(None),tind=slice(None), dtype=sp.float32):

        var=super(Ctl,self).fromfile(varname,iind=iind,jind=jind,kind=kind,tind=tind,dtype=dtype)
        
        # Flip north-south
        var.grid['lat']=var.grid['lat'].copy()[-1::-1]
        var.data=var.data.copy()[:,:,-1::-1,:]
        
        return var
        
ctl=Ctl()
