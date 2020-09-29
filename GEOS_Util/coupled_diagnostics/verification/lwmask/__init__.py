import scipy as sp
from datetime import date
import scipy as sp
from g5lib import dset

__all__=['mask']

class Ctl(dset.GADset):
    def __init__(self):
        name='lwmask'
        undef=1e15
        dir='/discover/nobackup/projects/gmao/share/dao_ops/verification/lwmask'

        flist=[dir+'/lwmask1440721.data']
        
        lon=sp.arange(-180,180,0.25); lat=sp.arange(-90,90.1,0.25); lev=sp.zeros(1)
        grid=dset.Grid(lon,lat,lev)
        vlist=[('mask','>f4',grid.dims)]
        time=sp.array([date(2000,1,1)],dtype='|O')
        super(Ctl,self).__init__(flist,vlist,grid,time,undef,name)

mask=Ctl()('mask')
