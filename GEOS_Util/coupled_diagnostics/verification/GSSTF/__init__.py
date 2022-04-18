import scipy as sp
import datetime
import dateutil.rrule as rrule
from g5lib import dset

class Ctl(dset.GADset):
    def __init__(self):
        name='GSSTF'
        undef=-999
        
        dir='/discover/nobackup/projects/gmao/share/dao_ops/verification/GSSTF_Surf_Flux'
        rr=rrule.rrule(rrule.MONTHLY, dtstart=datetime.date(2000,1,1),count=12)
        flist=[dir+'/GSSTF.Jul1987_Dec2000_Clim.'
               +str(date.month).zfill(2)+'.data' 
               for date in rr]
        
        lon=sp.arange(-179.5,180.,1); lat=sp.arange(-89.5,90,1.); lev=sp.zeros(1)
        grid=dset.Grid(lon,lat,lev)
        vlist=[]
        vlist.append(('','i4',1))
        vlist.append(('eflux','>f4',grid.dims))
        vlist.append(('','i4',1))
        vlist.append(('','i4',1))
        vlist.append(('uflux','>f4',grid.dims))
        vlist.append(('','i4',1))
        vlist.append(('','i4',1))
        vlist.append(('vflux','>f4',grid.dims))
        vlist.append(('','i4',1))
        vlist.append(('','i4',1))
        vlist.append(('hflux','>f4',grid.dims))
        vlist.append(('','i4',1))
        vlist.append(('','i4',1))
        vlist.append(('q10m','>f4',grid.dims))
        vlist.append(('','i4',1))
        vlist.append(('','i4',1))
        vlist.append(('tpw','>f4',grid.dims))
        vlist.append(('','i4',1))
        vlist.append(('','i4',1))
        vlist.append(('s10','>f4',grid.dims))
        vlist.append(('','i4',1))
        vlist.append(('','i4',1))
        vlist.append(('dq','>f4',grid.dims))
        vlist.append(('','i4',1))
        
        time=sp.array(rr[:],dtype='|O')
        
        super(Ctl,self).__init__(flist,vlist,grid,time,undef,name)


ctl=Ctl()
