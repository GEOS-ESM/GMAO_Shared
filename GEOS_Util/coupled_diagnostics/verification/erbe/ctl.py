from g5lib import dset
import scipy as sp
import datetime
import dateutil.rrule as rrule

class Ctl(dset.GADset):
    def __init__(self):
        name='ERBE'
        undef=1e15
        
        dir='/discover/nobackup/projects/gmao/share/dao_ops/verification/Clouds_radiation/erbe'
        flist=sp.array([dir+'/erbe2x25.data'])

        lon=sp.arange(0.,360,2.5); lat=sp.arange(-90.,91,2.); lev=sp.zeros(1)
        grid=dset.Grid(lon,lat,lev)

        print grid.dims

        vlist=[('OLR','>f4',grid.dims)]
        vlist.append(('OLRCLR','>f4',grid.dims))
        vlist.append(('NSR','>f4',grid.dims))
        vlist.append(('NSRCLR','>f4',grid.dims))
        vlist.append(('OSR','>f4',grid.dims))
        vlist.append(('OSRCLR','>f4',grid.dims))

        r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(1985,1,1),count=60)
        time=sp.array(r[:],dtype='|O')

        super(Ctl,self).__init__(flist,vlist,grid,time,undef,name)

    def fromfile(self,varname,iind=slice(None),jind=slice(None),kind=slice(None),\
                 tind=slice(None)):
        
        var=super(Ctl,self).fromfile(varname,iind=iind,jind=jind,kind=kind,\
                                     tind=tind)
        
        # Applly land mask
        var.data=sp.ma.masked_values(var.data, self.undef)
        
        return var


ctl=Ctl()
