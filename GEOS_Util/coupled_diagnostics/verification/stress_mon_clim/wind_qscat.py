import scipy as sp
import os
import matplotlib.pyplot as pl
from mpl_toolkits.basemap.cm import sstanom
import my_plots as mpl
from matplotlib import dates
path=os.environ['NOBACKUP']+'/verification/stress_mon_clim'
execfile(path+'/ctl.py')
ctl=Ctl()

taux={}
taux['clim']=ctl.fromfile('taux'); taux['clim'].data*=10.
taux['clim'].shiftgrid(.1)
taux['clim'].grid['lon']=sp.where(taux['clim'].grid['lon']<0.,taux['clim'].grid['lon']+360, \
                                  taux['clim'].grid['lon'])
ind=[0,1,11]; taux['djf']=taux['clim'].subset(tind=ind).ave(0); taux['djf'].name+=', DJF'
ind=[5,6,7]; taux['jja']=taux['clim'].subset(tind=ind).ave(0); taux['jja'].name+=', JJA'
taux['am']=taux['clim'].subset(tind=ind).ave(0); taux['am'].name+=', Annual Mean'

tauy={}
tauy['clim']=ctl.fromfile('tauy'); tauy['clim'].data*=10.
tauy['clim'].shiftgrid(.1)
tauy['clim'].grid['lon']=sp.where(tauy['clim'].grid['lon']<0.,tauy['clim'].grid['lon']+360, \
                                  tauy['clim'].grid['lon'])
ind=[0,1,11]; tauy['djf']=tauy['clim'].subset(tind=ind).ave(0); tauy['djf'].name+=', DJF'
ind=[5,6,7]; tauy['jja']=tauy['clim'].subset(tind=ind).ave(0); tauy['jja'].name+=', JJA'
tauy['am']=tauy['clim'].subset(tind=ind).ave(0); tauy['am'].name+=', Annual Mean'

# Equatorial annual cycle
lonind=sp.logical_and(taux['clim'].grid['lon'][0]>=130.0,taux['clim'].grid['lon'][0]<=280.0)
latind=sp.logical_and(taux['clim'].grid['lat'][:,0]>=-2.1,taux['clim'].grid['lat'][:,0]<=2.0)
taux['eqac']=taux['clim'].subset(iind=lonind,jind=latind).ave(2)
taux['eqac'].data-=taux['eqac'].ave(0).data
taux['eqac'].name=taux['clim'].name+', Eq. Annual Cycle'

# Plots

path+='/pics'

def plot_field(field, fig):
    clevs=sp.arange(-2,2.1,0.2)
    pl.figure(fig); pl.clf()
    field.copts={'levels': clevs,\
                 'cmap': sstanom}
    field.plot_map()
    mean,std=field.aave(ret_std=True)
    mpl.draw_stat((mean.data.squeeze(),std.data.squeeze()))
    pl.show()


#DJF
plot_field(taux['djf'],1)
pl.savefig(path+'/taux_djf.png')
plot_field(tauy['djf'],2)
pl.savefig(path+'/tauy_djf.png')

#JJA
plot_field(taux['jja'],1)
pl.savefig(path+'/taux_jja.png')
plot_field(tauy['jja'],2)
pl.savefig(path+'/tauy_jja.png')

#AM
plot_field(taux['am'],1)
pl.savefig(path+'/taux_am.png')
plot_field(tauy['am'],2)
pl.savefig(path+'/tauy_am.png')

# Plot annual cycle
pl.figure(3);pl.clf()
taux['eqac'].copts={'levels': sp.arange(-0.2,0.21,0.05),\
                   'cmap' : sstanom,\
                   'timefmt': dates.DateFormatter('%b')}
taux['eqac'].plot2d()
taux['eqac'].copts={'func': pl.contour,\
                   'colors': 'black',\
                   'levels': sp.arange(-0.2,0.21,0.05),\
                   'timefmt': dates.DateFormatter('%b')}
taux['eqac'].plot2d()
ax=pl.gca(); ax.yaxis.set_major_locator(dates.MonthLocator())
pl.grid(); pl.show()
    
pl.savefig(path+'/taux_eq_ac.png')
