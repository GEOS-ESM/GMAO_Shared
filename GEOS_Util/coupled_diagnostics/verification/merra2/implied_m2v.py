import os
from matplotlib import pyplot as pl
import scipy as sp
import netCDF4 as nc
import merra2, merra, ecmwf_int,lwmask, OAFlux

# Read MERRA2 surface fluxes
lw_m2=merra2.CtlClim('tavgM_2d_rad_Nx')('lwgnt').mean(0)
sw_m2=merra2.CtlClim('tavgM_2d_rad_Nx')('swgnt').mean(0)
lhf_m2=merra2.CtlClim('tavgM_2d_flx_Nx')('eflux').mean(0); lhf_m2.data*=-1
shf_m2=merra2.CtlClim('tavgM_2d_flx_Nx')('hflux').mean(0); shf_m2.data*=-1
netheat_m2=sw_m2()
netheat_m2.data+=(lw_m2.data+lhf_m2.data+shf_m2.data)

# Read ocean mask
f=nc.Dataset('/gpfsm/dnb02/projects/p53/merra2/scratch/MERRA2_100/stage/MERRA2_100.const_2d_asm_Nx.00000000.nc4')
frocean=sp.ma.masked_values(f.variables['FROCEAN'][0],0)
f.close()
netheat_m2.data*=frocean
aave_m2=netheat_m2.aave()
netheat_m2.data-=aave_m2.data

zonal=netheat_m2.gint(3)
imptrans_m2=zonal()
jm=imptrans_m2.dims[2]-1
for j in xrange(jm-1,-1,-1):
    imptrans_m2.data[:,:,j,:]=-zonal.subset(jind=slice(j,jm)).gint(2).data
imptrans_m2.data/=1e15

# Read MERRA surface fluxes
lw_m=merra.Ctl('tavgM_2d_rad_Nx')('LWGNT').clim(12).mean(0)
sw_m=merra.Ctl('tavgM_2d_rad_Nx')('SWGNT').clim(12).mean(0)
lhf_m=merra.Ctl('tavgM_2d_flx_Nx')('EFLUX').clim(12).mean(0); lhf_m.data*=-1
shf_m=merra.Ctl('tavgM_2d_flx_Nx')('HFLUX').clim(12).mean(0); shf_m.data*=-1
netheat_m=sw_m(); netheat_m.data+=(lw_m.data+lhf_m.data+shf_m.data)

# Read ocean mask
f=nc.Dataset('/discover/nobackup/dao_ops/scratch/d5_merra_jan79/stage/d5_merra_jan79.const_2d_asm_Nx.hdf')
frocean=sp.ma.masked_values(f.variables['FROCEAN'][0],0)
f.close()
netheat_m.data*=frocean
aave_m=netheat_m.aave()
netheat_m.data-=aave_m.data

zonal=netheat_m.gint(3)
imptrans_m=zonal()
jm=imptrans_m.dims[2]-1
for j in xrange(jm-1,-1,-1):
    imptrans_m.data[:,:,j,:]=-zonal.subset(jind=slice(j,jm)).gint(2).data
imptrans_m.data/=1e15

# Read Interim surface fluxes
lw_int=ecmwf_int.Ctl()('str').clim(12).mean(0); lw_int.data/=3600*24
sw_int=ecmwf_int.Ctl()('ssr').clim(12).mean(0); sw_int.data/=3600*24
lhf_int=ecmwf_int.Ctl()('slhf').clim(12).mean(0); lhf_int.data/=3600*24
shf_int=ecmwf_int.Ctl()('sshf').clim(12).mean(0); shf_int.data/=3600*24
netheat_int=sw_int()
netheat_int.data+=(lw_int.data+lhf_int.data+shf_int.data)

# Read ocean mask
mask=lwmask.mask
mask.shiftgrid(0.)
mask.regrid(netheat_int.grid)
netheat_int.data*=sp.ma.masked_less(mask.data,0.5)
aave_int=netheat_int.aave()
netheat_int.data-=aave_int.data

zonal=netheat_int.gint(3)
imptrans_int=zonal()
jm=imptrans_int.dims[2]-1
for j in xrange(jm-1,-1,-1):
    imptrans_int.data[:,:,j,:]=-zonal.subset(jind=slice(j,jm)).gint(2).data
imptrans_int.data/=1e15

# Read WHOI OAFlux
lw_whoi=OAFlux.ctllwsw('nlwrs').clim(12).mean(0); lw_whoi.data*=0.1
sw_whoi=OAFlux.ctllwsw('nswrs').clim(12).mean(0); sw_whoi.data*=0.1
lhf_whoi=OAFlux.ctllhsh('lhtfl').clim(12).mean(0); lhf_whoi.data*=0.1
shf_whoi=OAFlux.ctllhsh('shtfl').clim(12).mean(0); shf_whoi.data*=0.1
netheat_whoi=OAFlux.ctlqnet('qnet').clim(12).mean(0)
aave_whoi=netheat_whoi.aave()
netheat_whoi.data-=aave_whoi.data

zonal=netheat_whoi.gint(3)
imptrans_whoi=zonal()
jm=imptrans_whoi.dims[2]-1
for j in xrange(jm-1,-1,-1):
    imptrans_whoi.data[:,:,j,:]=-zonal.subset(jind=slice(j,jm)).gint(2).data
imptrans_whoi.data/=1e15

GWlats=[47, 24, -19, -30];
GWerr=[0.1, 0.3, 0.6, 0.3];
GWoht=[0.6, 1.8, -0.8, -0.6];
trenb=sp.ma.masked_values(sp.genfromtxt('/home/yvikhlia/verification/implied/ANNUAL_TRANSPORTS_1985_1989.ascii',skiprows=1),-999.0)/100.

path=os.environ['NOBACKUP']+'/verification/merra2/plots'
try:
    os.makedirs(path)
except OSError:
    pass

sw_m2.data/=1e9
sw_m.data/=1e9
sw_int.data/=1e9
sw_whoi.data/=1e9

lw_m2.data/=1e9
lw_m.data/=1e9
lw_int.data/=1e9
lw_whoi.data/=1e9 

lhf_m2.data/=1e9
lhf_m.data/=1e9
lhf_int.data/=1e9
lhf_whoi.data/=1e9

shf_m2.data/=1e9
shf_m.data/=1e9
shf_int.data/=1e9
shf_whoi.data/=1e9

netheat_m2.data/=1e9
netheat_m.data/=1e9
netheat_int.data/=1e9
netheat_whoi.data/=1e9

pl.figure(1,figsize=(12,15)); pl.clf()
ax=pl.subplot(321)
ax.text(0,1.1,'a)',transform=ax.transAxes)
sw_m2.gint(3).d()
sw_m.gint(3).d()
sw_int.gint(3).d()
sw_whoi.gint(3).d()
ax.set_xlim(-90,90)
ax.set_title('Net Surface Short Wave Radiation')
ax.set_ylabel('GW/m')
pl.grid()

ax=pl.subplot(322)
ax.text(0,1.1,'b)',transform=ax.transAxes)
lw_m2.gint(3).d()
lw_m.gint(3).d()
lw_int.gint(3).d()
lw_whoi.gint(3).d()
ax.set_xlim(-90,90)
ax.set_title('Net Surface Long Wave Radiation')
ax.set_ylabel('GW/m')
pl.grid()

ax=pl.subplot(323)
ax.text(0,1.1,'c)',transform=ax.transAxes)
lhf_m2.gint(3).d()
lhf_m.gint(3).d()
lhf_int.gint(3).d()
lhf_whoi.gint(3).d()
ax.set_xlim(-90,90)
ax.set_title('Latent Heat Flux')
ax.set_ylabel('GW/m')
pl.grid()

ax=pl.subplot(324)
ax.text(0,1.1,'d)',transform=ax.transAxes)
shf_m2.gint(3).d()
shf_m.gint(3).d()
shf_int.gint(3).d()
shf_whoi.gint(3).d()
ax.set_xlim(-90,90)
ax.set_title('Sensible Heat Flux')
ax.set_ylabel('GW/m')
pl.grid()

ax=pl.subplot(325)
ax.text(0,1.1,'e)',transform=ax.transAxes)
netheat_m2.gint(3).d()
netheat_m.gint(3).d()
netheat_int.gint(3).d()
netheat_whoi.gint(3).d()
ax.set_xlim(-90,90)
ax.set_title('Net Surface Heat Flux Anomaly')
ax.set_ylabel('GW/m')
pl.grid()

ax=pl.subplot(326)
ax.text(0,1.1,'f)',transform=ax.transAxes)
imptrans_m2.d();
imptrans_m.d();
imptrans_int.d();
imptrans_whoi.d()
pl.plot(trenb[:,0],trenb[:,6], linewidth=2,color='black')
pl.plot(GWlats,GWoht,'*',color='black')
pl.errorbar(GWlats,GWoht,yerr=GWerr,fmt='*',color='black')
pl.plot(trenb[:,0],trenb[:,6]+trenb[:,13],color='black')
pl.plot(trenb[:,0],trenb[:,6]-trenb[:,13],color='black')
ax.set_title('Implied Ocean Heat Transport')
ax.set_ylabel('PW')
ax.set_xlim(-90.,90)
ax.set_ylim(-3,3)
ax.legend(('MERRA-2','MERRA','ERA-interim','WHOI OAFlux','Trenberth-Caron','Ganachaud-Wunsch'),loc=4, fontsize='x-small')
pl.grid()
pl.show()
pl.savefig(path+'/implied.png')
