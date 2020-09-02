import os
import scipy as sp
from g5lib import dset, grid
import datetime
from dateutil import rrule

__all__=['ctl']

class Ctl(dset.NCDset):
    def __init__(self):
        name='Levitus'
        
        oceanval=os.environ.get('OCEANVAL',
                                '/discover/nobackup/projects/gmao/oceanval/verification')
        flist=[oceanval+'/levitus/levitus_grd.nc']

        t=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(0001,1,1),count=12)
        time=sp.array(t[:],dtype='|O')
        
#        super(Ctl,self).__init__(flist,levname='depth',\
#                                time=time,name=name,undef=10e11)

        super(Ctl,self).__init__(flist,\
                                     time=time,name=name,undef=10e11)

# Something wrong with levels in a datafile
        lev=sp.array((0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600,\
                          700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500,\
                          3000, 3500, 4000, 4500, 5000, 5500))
        lon=self.grid['lon']; lat=self.grid['lat']

        self.grid=grid.Grid(lon=lon,lat=lat,lev=lev)

    def fromfile(self,varname,iind=slice(None),jind=slice(None),kind=slice(None),\
                 tind=slice(None), dtype=sp.float32):
        
        var=super(Ctl,self).fromfile(varname,iind=iind,jind=jind,kind=kind,\
                                     tind=tind,maskandscale=False,dtype=dtype)
        
        # Applly land mask
        var.data=sp.ma.masked_invalid(var.data)
        
        return var

ctl=Ctl()
