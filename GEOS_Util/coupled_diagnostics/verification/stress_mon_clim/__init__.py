import scipy as sp
import os
from g5lib import dset
import datetime
import dateutil.rrule as rrule

class Ctl(dset.GADset):
    def __init__(self):
        
        name='QSCAT'
        undef=-9999.

        oceanval=os.environ.get('OCEANVAL',
                                '/discover/nobackup/projects/gmao/oceanval/verification')
        flist=[oceanval+'/QSCAT/qscat_clim.dat']

        lon=sp.arange(-179.5,180,1); lat=sp.arange(-89.5,90,1); lev=sp.zeros(1)
        grid=dset.Grid(lon,lat,lev)
        
        vlist=[]
        vlist.append(('taux','f4',grid.dims))
        vlist.append(('tauy','f4',grid.dims))

        t=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(0001,1,1),count=12)
        time=sp.array(t[:],dtype='|O')

        super(Ctl,self).__init__(flist,vlist,grid,time,undef,name)


ctl=Ctl()
