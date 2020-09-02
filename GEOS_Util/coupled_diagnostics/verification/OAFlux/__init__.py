import os
from g5lib import dset
import scipy as sp
import datetime
import dateutil.rrule as rrule

oceanval=os.environ.get('OCEANVAL',
                        '/discover/nobackup/projects/gmao/oceanval/verification')

class CtlQnet(dset.GADset):
    def __init__(self):
        expid='OAFlux'
        undef=10e15
        
        expdir=oceanval+'/OAfluxes/netheat_1983-2004/grads'
        r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(1983,7,1),count=258)
        flist=sp.array([expdir+'/a.qnet.'+str(year)+'.grads' for year in xrange(1983,2005)])
        lon=sp.arange(0.5,360,1); lat=sp.arange(-89.5,90,1); lev=sp.zeros(1)
        grid=dset.Grid(lon,lat,lev)
        vlist=[('qnet','<f4',grid.dims)]
        time=sp.array(r[:],dtype='|O')
        super(CtlQnet,self).__init__(flist,vlist,grid,time,undef,name=expid)

class CtlTurb(dset.GADset):
    def __init__(self,collection='lh'):
        expid='OAFlux'
        undef=10e15
        
        expdir=oceanval+'/OAfluxes/turbulence_1958-2006/'+collection+'/grads'
        r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(1958,1,1),count=588)
        flist=sp.array([expdir+'/a.'
                        +collection+'_oaflux.'
                        +str(year)
                        +'.grads' 
                        for year in xrange(1958,2007)])
        lon=sp.arange(0.5,360,1); lat=sp.arange(-89.5,90,1); lev=sp.zeros(1)
        grid=dset.Grid(lon,lat,lev)
        vlist=[(collection,'<f4',grid.dims),
               ('err','<f4',grid.dims)]
        time=sp.array(r[:],dtype='|O')
        super(CtlTurb,self).__init__(flist,vlist,grid,time,undef,name=expid)

class CtlLHSH(dset.NCDset):
    def __init__(self):
        expid='OAFlux'
        expdir=oceanval+'/OAFlux'
        r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(1981,1,1),until=datetime.date(2001,12,1))
        flist=sp.array([expdir+'/lhsh_oaflux_'+str(yr)+'.nc'  for yr in xrange(1981,2002)])
        time=sp.array(r[:],dtype='|O')
        nrecs=sp.ones(flist.size)*12
        super(CtlLHSH,self).__init__(flist,time=time,name=expid,nrecs=nrecs)

class CtlLWSW(dset.NCDset):
    def __init__(self):
        expid='OAFlux'
        expdir=oceanval+'/OAFlux'
        r=rrule.rrule(rrule.MONTHLY,dtstart=datetime.date(1983,7,1),until=datetime.date(2001,12,1))
        flist=sp.array([expdir+'/lwsw_isccp_'+str(yr)+'.nc'  for yr in xrange(1983,2002)])
        time=sp.array(r[:],dtype='|O')
        nrecs=sp.ones(flist.size)*12; nrecs[0]=6
        super(CtlLWSW,self).__init__(flist,time=time,name=expid,nrecs=nrecs)

ctlqnet=CtlQnet()
ctllwsw=CtlLWSW()
ctllhsh=CtlLHSH()
ctllh=CtlTurb()
ctlqa=CtlTurb('qa')
ctlsh=CtlTurb('sh')
ctlta=CtlTurb('ta')
ctlts=CtlTurb('ts')
ctlws=CtlTurb('ws')
