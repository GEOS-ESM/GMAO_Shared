import datetime
import dateutil.rrule as rrule
import scipy as sp
from g5lib import dset

class Ctl(dset.NCDset):
    def __init__(self,collection='tavg1_2d_rad_Nx'):
        '''
        See directory for possible collections

        '''
        stream='d5124_rpit_jan00'
        flist=['/home/dao_ops/'+stream+'/run/.../scratch/diag/Y'
            +str(yr)+'/M'+str(mo).zfill(2)+'/'+stream+'.'+collection+'.monthly.'+str(yr)+str(mo).zfill(2)+'.nc4'
            for yr in xrange(2000,2004) for mo in xrange(1,13)]
        stream='d5124_rpit_jan04'
        flist+=['/home/dao_ops/'+stream+'/run/.../scratch/diag/Y'
            +str(yr)+'/M'+str(mo).zfill(2)+'/'+stream+'.'+collection+'.monthly.'+str(yr)+str(mo).zfill(2)+'.nc4'
            for yr in xrange(2004,2012) for mo in xrange(1,13)]
        stream='d5124_rpit_jan12'
        flist+=['/home/dao_ops/'+stream+'/run/.../scratch/diag/Y'
            +str(yr)+'/M'+str(mo).zfill(2)+'/'+stream+'.'+collection+'.monthly.'+str(yr)+str(mo).zfill(2)+'.nc4'
            for yr in xrange(2012,2019) for mo in xrange(1,13)]
        r=rrule.rrule(rrule.MONTHLY,
                      dtstart=datetime.date(2000,1,1),
                      until=datetime.date(2018,12,1))
        time=sp.array(r[:],dtype='|O')

        super(Ctl,self).__init__(sp.array(flist),time=time,name='FPIT',nrecs=sp.ones(sp.size(flist)))
