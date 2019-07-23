import os, sys
import scipy as sp
from matplotlib import pyplot as pl
from g5lib import g5dset

exp=g5dset.read_exp(sys.argv[1])
exp.ctl=g5dset.Ctl('geosgcm_surf')
rad=exp.ctl('RADSRF',dates=exp.dates).clim(12).mean(0)
lhf=exp.ctl('LHFX',dates=exp.dates).clim(12).mean(0)
shf=exp.ctl('SHFX',dates=exp.dates).clim(12).mean(0)
fro=exp.ctl('FROCEAN',dates=exp.dates).clim(12).mean(0)
fro.data=sp.ma.masked_values(fro.data,0.0)
netheat=rad(); netheat.data+=(-lhf.data-shf.data)
netheat.data*=fro.data
netheat.data-=netheat.aave().data

zonal=netheat.gint(3)
imptrans=zonal()
jm=imptrans.dims[2]-1
for j in xrange(jm-1,-1,-1):
    imptrans.data[:,:,j,:]=-zonal.subset(jind=slice(j,jm)).gint(2).data

imptrans.data/=1e15
GWlats=[47, 24, -19, -30];
GWerr=[0.1, 0.3, 0.6, 0.3];
GWoht=[0.6, 1.8, -0.8, -0.6];

trenb=sp.ma.masked_values(sp.genfromtxt('/home/yvikhlia/verification/implied/ANNUAL_TRANSPORTS_1985_1989.ascii',skiprows=1),-999.0)/100.

path=exp.plot_path
try:
    os.makedirs(path)
except OSError:
    pass


pl.figure(1); pl.clf()
imptrans.d();
pl.plot(trenb[:,0],trenb[:,6], linewidth=2,color='red')
pl.plot(GWlats,GWoht,'*',color='green')
pl.errorbar(GWlats,GWoht,yerr=GWerr,fmt='*',color='green')
pl.plot(trenb[:,0],trenb[:,6]+trenb[:,13],color='red')
pl.plot(trenb[:,0],trenb[:,6]-trenb[:,13],color='red')
ax=pl.gca()
ax.set_title('Implied Ocean Heat Transport')
ax.set_ylabel('PW')
ax.set_xlim(-90.,90)
ax.set_ylim(-3,3)
ax.legend((exp.__name__,'Trenberth-Caron','Ganachaud-Wunsch'),loc=4)
pl.grid()
pl.tight_layout()
pl.show()
pl.savefig(path+'/implied.png')
