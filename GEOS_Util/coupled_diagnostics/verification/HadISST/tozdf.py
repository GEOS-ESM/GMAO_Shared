import os, cPickle
import scipy as sp
from datetime import datetime

execfile(os.environ['HOME']+'/verification/HadISST/ctl.py')

ctl=Ctl()
dates=(datetime(1951,1,1),datetime(2008,12,1))
x=ctl('sst',levs=(0,),dates=dates)


nt,km,jm,im=x.dims
x.clim(12,anom=True)

rec=[]
rec.append(('','i4'))
rec.append(('head','f4',100))
rec.append(('','i4'))
rec.append(('','i4'))
rec.append(('field','f4',(jm,im)))
rec.append(('','i4'))

xx=sp.zeros(nt, dtype=rec)

i1=100*4
i2=jm*im*4

xx['f0']=xx['f2']=i1
xx['f3']=xx['f5']=i2

xx['field']=x.data[:,0,:,:].view(sp.ndarray)

for i,tt in enumerate(x.time):
    xx['head'][i,0]=tt.year
    xx['head'][i,1]=tt.month
    xx['head'][i,2]=tt.day
    xx['head'][i,3]=tt.hour
    xx['head'][i,4]=tt.minute
    xx['head'][i,5]=tt.second

xx['head'][:,6]=1
xx['head'][:,7]=1
xx['head'][:,8]=im
xx['head'][:,9]=jm
xx['head'][:,10]=1
xx['head'][:,11]=-179.5
xx['head'][:,12]=-89.5
xx['head'][:,13]=1
xx['head'][:,14]=1
xx['head'][:,15]=1
xx['head'][:,16]=2
xx['head'][:,17]=4
xx['head'][:,18]=-1e30
xx['head'][:,19]=1

xx.tofile(os.environ['HOME']+'/verification/HadISST/data/temp.zdf')
    
