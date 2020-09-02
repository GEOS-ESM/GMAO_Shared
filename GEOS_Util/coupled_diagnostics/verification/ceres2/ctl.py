from g5lib import dset
import scipy as sp
import datetime
import dateutil.rrule as rrule

class Ctl(dset.NCDset):
    def __init__(self):
        name='CERES2'
        undef=1e15
        
        dir='/discover/nobackup/projects/gmao/share/dao_ops/verification/JPOTTER'
        flist=sp.array([dir+'/CERES2_'+str(mo).zfill(2)+'_climo.nc' for mo in xrange(1,13)])
        
        t=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(0001,1,1),count=12)
        time=sp.array(t[:],dtype='|O')
        
        super(Ctl,self).__init__(flist, time=time,name=name, undef=1e20)

ctl=Ctl()
