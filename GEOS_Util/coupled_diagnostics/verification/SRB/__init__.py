import scipy as sp
import datetime
import dateutil.rrule as rrule
from g5lib import dset

class Ctl(dset.GADset):
    def __init__(self, col='SW',subcol='FTOA'):
        '''
        col is 'LW' or 'SW'
        subcol
           for SW: 'CWV', 'FABS', 'FALL', 'FCLR', 'FTOA', 'SALB', 'TCLD'
           for LW: 'sfc_down', 'sfc_up', 'toa_up',
                  'clr_sfc_down', 'clr_sfc_up', 'clr_toa_up'
        '''
        name='SRB_'+col
        undef=-999
        
        dir='/discover/nobackup/projects/gmao/share/dao_ops/verification/SRB'
        rr=rrule.rrule(rrule.MONTHLY, dtstart=datetime.date(1983,7,1),until=datetime.date(2000,6,1))
        flist=[dir+'/'+col+'/'+subcol+'_monthly_'
               +str(date.year)+str(date.month).zfill(2)+'.binary' 
               for date in rr]
        
        lon=sp.arange(0.5,360.1,1); lat=sp.arange(-89.5,90,1.); lev=sp.zeros(1)
        grid=dset.Grid(lon,lat,lev)
        vlist=[(subcol,'>f4',grid.dims)]
        time=sp.array(rr[:],dtype='|O')
        
        super(Ctl,self).__init__(flist,vlist,grid,time,undef,name)

ctl=Ctl()
