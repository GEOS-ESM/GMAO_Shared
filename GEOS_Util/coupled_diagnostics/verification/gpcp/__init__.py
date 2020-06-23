import os
from g5lib import dset
import scipy as sp
import datetime
import dateutil.rrule as rrule

class Ctl(dset.GADset):
    def __init__(self):
        name='GPCP'
        undef=-9.99990000e+04
        
        syr=1982; eyr=2006
        yrs=[str(i) for i in range(syr,eyr)]
        mon=[str(i).zfill(2) for i in range(1,13)]
        dates=[i+j for i in yrs for j in mon]
        oceanval=os.environ.get('OCEANVAL',
                                '/discover/nobackup/projects/gmao/oceanval/verification')
        flist=sp.array([oceanval+'/gpcp/gpcp_v2_psg.'+i+'.data' for i in dates])

        lon=sp.arange(1.25,361,2.5); lat=sp.arange(-88.75,91,2.5); lev=sp.zeros(1)
        grid=dset.Grid(lon,lat,lev)

        vlist=[('precip','>f4',grid.dims)]

        r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(syr,1,1),count=len(flist))
        time=sp.array(r[:],dtype='|O')

        super(Ctl,self).__init__(flist,vlist,grid,time,undef,name)

ctl=Ctl()
