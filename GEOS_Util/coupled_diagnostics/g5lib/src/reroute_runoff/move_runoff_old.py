import os
import scipy as sp
import netCDF4 as nc

# Read ocean tiles which has discharge
print 'Reading discharge tile indexes'
dir=os.environ['NOBACKUP']+'/workdir/move_runoff'
runoff_table=dir+'/runoff.bin'

Nt=sp.fromfile(runoff_table,dtype='i4',count=3)[1]

vlist=[('','i4'),('Nt','i4'),('','i4'),
       ('','i4'),('sind','i4',Nt),('','i4'),
       ('','i4'),('dind','i4',Nt),('','i4'),
       ('','i4'),('weight','f4',Nt),('','i4')]

ind=sp.unique(sp.fromfile(runoff_table,dtype=vlist)['dind']); ind-=1 # indises start with 0 

# Read ocean i,j corresponding to tmp from tile file
print 'Reading discharge grid indexes'
tfile=dir+'/tile.dat'
sind=sp.loadtxt(tfile, skiprows=8,dtype='f4')[ind,9:7:-1].astype(int); sind-=1;
sp.save(dir+'/sind.npy',sind)
#sind=sp.load(dir+'/sind.npy')

# Read ocean land mask
print 'Reading ocean land mask'
gfile=dir+'/grid_spec.nc'
ff=nc.Dataset(gfile)
lon=ff.variables['x_T'][:]
lat=ff.variables['y_T'][:]
wet=ff.variables['wet'][:]
sh=wet.shape
tmp=sp.zeros((sh[0]+2,sh[1]+2))
tmp[1:-1,1:-1]=wet
nwet=tmp[1:-1,1:-1]+tmp[:-2,:-2]+tmp[1:-1,:-2]\
    +tmp[2:,:-2]+tmp[2:,1:-1]+tmp[2:,2:]\
    +tmp[1:-1,2:]+tmp[:-2,2:]+tmp[:-2,1:-1]


# Move discharge.
print 'Moving discharge'
def move_discharge(ji,nwet):
  sj,ej=max(ji[0]-2,0),min(ji[0]+3,sh[0])
  si,ei=max(ji[1]-2,0),min(ji[1]+3,sh[1])
  tmp=nwet[sj:ej,si:ei]
  if tmp.max() >= 6.0:
    j,i=sp.unravel_index(tmp.argmax(),tmp.shape)
    ji[0],ji[1]=sj+j,si+i
  return ji
  

dind=sind.copy()
for i,ji in enumerate(dind):
  if nwet[ji[0],ji[1]] < 6.0 or wet[ji[0],ji[1]]==0.0:
    # move discharge
    dind[i]=move_discharge(ji,nwet)

# Move discharge out of Uruguay Bay
ind=sp.argmax(sp.all(sind==(350,891),axis=1)); dind[ind]=sp.array((345,896))
ind=sp.argmax(sp.all(sind==(349,891),axis=1)); dind[ind]=sp.array((344,896))
ind=sp.argmax(sp.all(sind==(350,890),axis=1)); dind[ind]=sp.array((345,895))
ind=sp.argmax(sp.all(sind==(349,889),axis=1)); dind[ind]=sp.array((344,894))
ind=sp.argmax(sp.all(sind==(349,887),axis=1)); dind[ind]=sp.array((344,892))
ind=sp.argmax(sp.all(sind==(351,888),axis=1)); dind[ind]=sp.array((346,893))
ind=sp.argmax(sp.all(sind==(352,887),axis=1)); dind[ind]=sp.array((345,897))
ind=sp.argmax(sp.all(sind==(347,891),axis=1)); dind[ind]=sp.array((342,896))
ind=sp.argmax(sp.all(sind==(349,893),axis=1)); dind[ind]=sp.array((344,898))

# Move discharge out of Amazon Bay
ind=sp.argmax(sp.all(sind==(502,919),axis=1)); dind[ind]=sp.array((505,924))
ind=sp.argmax(sp.all(sind==(500,918),axis=1)); dind[ind]=sp.array((500,924))
ind=sp.argmax(sp.all(sind==(497,922),axis=1)); dind[ind]=sp.array((502,927))
ind=sp.argmax(sp.all(sind==(498,921),axis=1)); dind[ind]=sp.array((498,926))
ind=sp.argmax(sp.all(sind==(501,921),axis=1)); dind[ind]=sp.array((501,926))
ind=sp.argmax(sp.all(sind==(499,923),axis=1)); dind[ind]=sp.array((499,928))

disch_old=wet.copy(); disch_old[sind[:,0,],sind[:,1]]=2
disch_new=wet.copy(); disch_new[dind[:,0,],dind[:,1]]=2

ind=sp.unique(sp.where(sind != dind)[0])
tmp=sp.concatenate((sp.roll(sind[ind]+1,1,1),sp.roll(dind[ind]+1,1,1)),1)
sp. savetxt('/home/yvikhlia/nobackup/workdir/move_runoff/ocean_routing_1440x1080.txt',tmp,fmt='%4.1i')
