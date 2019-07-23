import cPickle, os
from datetime import datetime

path=os.environ['NOBACKUP']+'/verification/HadISST'
execfile(path+'/ctl.py')

ctl=Ctl()
dates=(datetime(1951,1,1),datetime(2008,12,1))

varname='sst'

var=ctl(varname,dates=dates)
var.clim(12,anom=True)
var.pca()

try:
    os.makedirs(path+'/data')
except OSError:
    pass 

f=open(path+'/data/'+varname+'_eofs.dat','w')
cPickle.dump(var,f,-1)
f.close()
