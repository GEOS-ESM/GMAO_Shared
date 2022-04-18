from g5lib import dset
import scipy as sp
import datetime
import dateutil.rrule as rrule

class Ctl(dset.GADset):
    def __init__(self):
        name='ECMWF_Interim'
        undef=-1e-15
        
        r=rrule.rrule(rrule.YEARLY,dtstart=datetime.date(1980,1,1),until=datetime.date(2011,1,1))        
        dir='/discover/nobackup/projects/gmao/share/dao_ops/verification/ECMWF_Interim_FullRes/surface'
        flist=sp.array([dir+'/Interim_monthly_FLX_'+str(date.year)+'.bin'
                        for date in r[:]])
        
        lon=sp.arange(0,360.5,0.75); lat=sp.arange(-90,90.5,0.75); lev=sp.zeros(1)
        grid=dset.Grid(lon,lat,lev)
        
        vlist=[('cp','f4',grid.dims),
               ('e','f4',grid.dims),
               ('slhf','f4',grid.dims),
               ('sshf','f4',grid.dims),
               ('ssr','f4',grid.dims),
               ('ssrc','f4',grid.dims),
               ('ssrd','f4',grid.dims),
               ('str','f4',grid.dims),
               ('strc','f4',grid.dims),
               ('strd','f4',grid.dims),
               ('tisr','f4',grid.dims),
               ('tp','f4',grid.dims),
               ('tsr','f4',grid.dims),
               ('tsrc','f4',grid.dims),
               ('ttr','f4',grid.dims),
               ('ttrc','f4',grid.dims),
               ('csf','f4',grid.dims),
               ('es','f4',grid.dims),
               ('par','f4',grid.dims),
               ('ro','f4',grid.dims),
               ('sf','f4',grid.dims),
               ('smlt','f4',grid.dims),
               ('lsf','f4',grid.dims),
               ('lsp','f4',grid.dims),
               ('parcs','f4',grid.dims),
               ('uvb','f4',grid.dims),
               ]
        r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(1980,1,1),until=datetime.date(2011,12,1))
        time=sp.array(r[:],dtype='|O')

        super(Ctl,self).__init__(flist,vlist,grid,time,undef,name)

ctl=Ctl()
