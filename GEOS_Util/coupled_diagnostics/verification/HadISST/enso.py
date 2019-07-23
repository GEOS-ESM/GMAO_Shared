import os, pickle
from datetime import datetime
import scipy as sp
import HadISST

ctl=HadISST.Ctl()
var='sst'
sst=ctl(var,dates=(datetime(1950,1,1),datetime(2008,12,1)))
sst.clim(12,anom=True); sst.detrend()

sst.shiftgrid(30.)

# First PC
pc,eof=sst.get_pc(1)

# Save data
path=os.environ['NOBACKUP']+'/verification/HadISST'
try:
    os.makedirs(path+'/data')
except OSError:
    pass 

f=open(path+'/data/'+var+'_pc1.dat','w')
pickle.dump(pc,f)
f.close()

f=open(path+'/data/'+var+'_eof1.dat','w')
pickle.dump(eof,f)
f.close()

