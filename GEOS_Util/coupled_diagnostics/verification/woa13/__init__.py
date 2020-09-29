import os,datetime
from dateutil import rrule
import scipy as sp
from g5lib import dset

__all__=['tctl','sctl','rhoctl']

class Ctl(dset.NCDset):
    def __init__(self,var='t'):
        name='WOA13 '+var.upper()
        
        oceanval=os.environ.get('OCEANVAL',
                                '/discover/nobackup/projects/gmao/oceanval/verification')
        data_dir=oceanval+'/WOA13'
        if var=='I':
            ver='01'
        else:
            ver='01v2'
        flist=sp.array([data_dir+'/woa13_decav_'+var+str(mm).zfill(2)+'_'+ver+'.nc' for mm in range(1,13)],dtype='|O')
        
        rr=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(2013,1,1),count=12)
        time=sp.array(rr[:],dtype='|O')

        super(Ctl,self).__init__(flist,time=time,name=name,levname='depth')
        
tctl=Ctl('t')
sctl=Ctl('s')
rhoctl=Ctl('I')
