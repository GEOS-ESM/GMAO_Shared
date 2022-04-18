import scipy as sp
import datetime
import dateutil.rrule as rrule
from g5lib  import dset

class Ctl(dset.GADset):
    def __init__(self):
        name='Reynolds MERRA'
        undef=-999.0
        
        flist=['/discover/nobackup/projects/gmao/share/dao_ops/fvInput/g5gcm/bcs/realtime/SST/360x180/dataoceanfile_MERRA_sst_1971-current.360x180.LE']

        lon=sp.arange(-179.5,180,1); lat=sp.arange(-89.5,90,1); lev=sp.zeros(1)
        grid=dset.Grid(lon,lat,lev)

        vlist=[('','i4',1),
               ('head','f4',14),
               ('','i4',1),
               ('','i4',1),
               ('sst','f4',grid.dims),
               ('','i4',1)]

        time=sp.repeat(datetime.date(1,1,1),1732)
        a=sp.memmap(flist[0],dtype=vlist,mode='r')
        
        for i,hh in enumerate(a['head']):
            time[i]=datetime.date(hh[0],hh[1],hh[2])

        a.close()

        super(Ctl,self).__init__(flist,vlist,grid,time,undef,name)
