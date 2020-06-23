import os, pickle
from  g5lib import domain
from datetime import datetime
import scipy as sp
import HadISST

path=os.environ['NOBACKUP']+'/verification/HadISST'

ctl=HadISST.Ctl()
dates=(datetime(1950,1,1),datetime(2008,12,1))

var='sst'

sst=ctl(var,dates=dates)
sst.clim(12,anom=True); sst.detrend()

sst.shiftgrid(30.)

# PDO index and pattern
pc1,eof1=sst(lons=(120,260),lats=(20,60)).get_pc(1); 

pos,neg=sst.composites(pc1)

# Save data
try:
    os.makedirs(path+'/data')
except OSError:
    pass 

f=open(path+'/data/'+var+'_pc1_np.dat','w')
pickle.dump(pc1,f)
f.close()

f=open(path+'/data/'+var+'_eof1_np.dat','w')
pickle.dump(eof1,f)
f.close()

f=open(path+'/data/'+var+'_pdo_plus.dat','w')
pickle.dump(pos,f)
f.close()

f=open(path+'/data/'+var+'_pdo_minus.dat','w')
pickle.dump(neg,f)
f.close()

