import cPickle, os
from datetime import datetime
import merra

varname='H'

ctl=merra.Ctl()
var=ctl(varname,levs=(350,))
var.clim(12,anom=True)
var.pca()

path=os.environ['HOME']+'/verification/merra'

try:
    os.makedirs(path+'/data')
except OSError:
    pass 

f=open(path+'/data/'+varname+'_eofs.dat','w')
cPickle.dump(var,f,-1)
f.close()
