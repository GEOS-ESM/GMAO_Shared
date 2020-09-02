import os, pickle
from g5lib import domain
from datetime import datetime
import scipy as sp

path=os.environ['NOBACKUP']+'/verification/HadISST'
execfile(path+'/ctl.py')

ctl=Ctl()
dates=(datetime(1870,1,1),datetime(2008,12,1))
dom=domain.Domain(lons=(-180,180),lats=(90,-90),dates=dates)

var='sst'

sst=ctl.fromfile(var,**dom(ctl.grid,ctl.time))
sst.clim(12,anom=True); sst.detrend()

sst.shiftgrid(30.)
sst.grid['lon']=sp.where(sst.grid['lon']>=29.5,sst.grid['lon'],sst.grid['lon']+360.)

# AMO index
amodom=domain.Domain(lons=(270,360),lats=(0,80),dates=dates)
amo=sst.subset(**amodom(sst.grid,sst.time)).aave()
amo.running_mean(121)

# Make composites
pos,neg=sst.composites(amo)

# Save data
try:
    os.makedirs(path+'/data')
except OSError:
    pass 

f=open(path+'/data/'+var+'_amo.dat','w')
pickle.dump(amo,f)
f.close()

f=open(path+'/data/'+var+'_amo_plus.dat','w')
pickle.dump(pos,f)
f.close()

f=open(path+'/data/'+var+'_amo_minus.dat','w')
pickle.dump(neg,f)
f.close()

