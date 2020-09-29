import os, pickle
from g5lib import domain
import scipy as sp

path=os.environ['NOBACKUP']+'/verification/merra'
execfile(path+'/ctl.py')

ctl=Ctl()
dates=(datetime.date(1979,1,1),datetime.date(2009,12,1))
dom=domain.Domain(lons=(-180,180),lats=(-90,90),levs=(1000,),dates=dates)

var='SLP'

slp=ctl.fromfile(var,**dom(ctl.grid,ctl.time)); slp.data/=100.
slp.clim(12,anom=True); slp.detrend()

# NAO index
d=domain.Domain(lons=(-27,),lats=(34,),levs=dom.levs,dates=dom.dates)
nao=slp.subset(**d(slp.grid,slp.time))
d=domain.Domain(lons=(-32,),lats=(60,),levs=dom.levs,dates=dom.dates)
nao.data-=slp.subset(**d(slp.grid,slp.time)).data

# Make composites
pos,neg=slp.composites(nao)

# Save data
try:
    os.makedirs(path+'/data')
except OSError:
    pass 

f=open(path+'/data/'+var+'_nao.dat','w')
pickle.dump(nao,f)
f.close()

f=open(path+'/data/'+var+'_nao_plus.dat','w')
pickle.dump(pos,f)
f.close()

f=open(path+'/data/'+var+'_nao_minus.dat','w')
pickle.dump(neg,f)
f.close()

