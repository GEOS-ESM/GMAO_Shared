import os, pickle
from g5lib import domain
#from datetime import datetime
import scipy as sp

path=os.environ['NOBACKUP']+'/verification/merra'
execfile(path+'/ctl.py')

ctl=Ctl()
dates=(datetime.date(1979,1,1),datetime.date(2009,12,1))
dom=domain.Domain(lons=(-180,180),lats=(-90,90),levs=(500,),dates=dates)

var='H'

h500=ctl.fromfile(var,**dom(ctl.grid,ctl.time)); 
h500.clim(12,anom=True); h500.detrend()

h500.shiftgrid(30.)
h500.grid['lon']=sp.where(h500.grid['lon']>=29.,h500.grid['lon'],h500.grid['lon']+360.)

# PNA index
d=domain.Domain(lons=(200,),lats=(20,),levs=dom.levs,dates=dom.dates)
pna=h500.subset(**d(h500.grid,h500.time))
d=domain.Domain(lons=(195,),lats=(45,),levs=dom.levs,dates=dom.dates)
pna.data-=h500.subset(**d(h500.grid,h500.time)).data
d=domain.Domain(lons=(245,),lats=(55,),levs=dom.levs,dates=dom.dates)
pna.data+=h500.subset(**d(h500.grid,h500.time)).data
d=domain.Domain(lons=(275,),lats=(30,),levs=dom.levs,dates=dom.dates)
pna.data-=h500.subset(**d(h500.grid,h500.time)).data

# Make composites
pos,neg=h500.composites(pna)

# Save data
try:
    os.makedirs(path+'/data')
except OSError:
    pass 

f=open(path+'/data/'+var+'_pna.dat','w')
pickle.dump(pna,f)
f.close()

f=open(path+'/data/'+var+'_pna_plus.dat','w')
pickle.dump(pos,f)
f.close()

f=open(path+'/data/'+var+'_pna_minus.dat','w')
pickle.dump(neg,f)
f.close()

