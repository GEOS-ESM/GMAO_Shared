import os
import scipy as sp
from g5lib import dset
import datetime
import dateutil.rrule as rrule

class Ctl(dset.NCDset):
    def __init__(self,**opts):

        expid=opts.get('name','MERRA')
        collection=opts.get('collection','inst3_3d_asm_Cp')
        expdir=opts.get('data_path',os.environ['NOBACKUP']+\
                            '/verification/merra/'+collection)
        syr=opts.get('start_year',1979)
        eyr=opts.get('end_year',2010)
        
        # Create meta-data
        yrs=[str(i).zfill(4) for i in range(syr,eyr+1)]
        mon=[str(i).zfill(2) for i in range(1,13)]
        dates=[i+j for i in yrs for j in mon]
        flist=sp.array([expdir+'/d5_merra.inst3_3d_asm_Cp.monthly.'+date+'.nc' \
                        for date in dates])

        r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(syr,1,1),count=len(flist))
        time=sp.array(r[:],dtype='|O')

        super(Ctl,self).__init__(flist,time=time,name=expid)


