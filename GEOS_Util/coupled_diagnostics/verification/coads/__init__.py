import os
import scipy as sp
import datetime
import dateutil.rrule as rrule
from g5lib import dset 
from vlist import vlist

__all__=['ctl']
class Ctl(dset.GADset):
    def __init__(self):
        name='COADS'
        undef=-1e10

        flist=[os.environ['NOBACKUP']+'/verification/coads/coads2x25_climate.dat']

        lon=sp.arange(0.0,360.,2.5); lat=sp.arange(-90.,91,2.); lev=sp.zeros(1)
        grid=dset.Grid(lon,lat,lev)
        rr=rrule.rrule(rrule.MONTHLY, dtstart=datetime.date(1980,1,15),count=12)
        time=sp.array(rr[:],dtype='|O')
        
        super(Ctl,self).__init__(flist,vlist,grid,time,undef,name)

ctl=Ctl()
