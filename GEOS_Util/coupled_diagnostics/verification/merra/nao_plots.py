import os, pickle
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import scipy as sp
import mpl_toolkits.basemap as bm

dsetname='merra'
varname='SLP'
indname='nao'
path=os.environ['NOBACKUP']+'/verification/'+dsetname
indfile=path+'/data/'+varname+'_'+indname+'.dat'
indpic=path+'/pics/'+varname+'_'+indname+'.png'
indtitle='Merra NAO Index (DJF)';xlab='years'
indylim=(-20,20); tint=5
posfile=path+'/data/'+varname+'_'+indname+'_plus.dat'
pospic=path+'/pics/'+varname+'_'+indname+'_plus.png'
postitle='Merra SLP > NAO_std'
negfile=path+'/data/'+varname+'_'+indname+'_minus.dat'
negpic=path+'/pics/'+varname+'_'+indname+'_minus.png'
negtitle='Merra SLP < -NAO_std'
units='mbar'
copts={'levels': sp.arange(-7,7.1,1),\
           'cmap': bm.cm.sstanom}
cbar_opts={'orientation': 'vertical'}
fmap=bm.Basemap(projection='npstere',boundinglat=0.0,lon_0=-60)

try:
    os.makedirs(path+'/pics')
except OSError:
    pass

# Plot index, DJF means
f=open(indfile); x=pickle.load(f); f.close()
tind=range(12,x.time.size,12)
ind=x.subset(tind=tind); 

for i,tt in enumerate(ind.time):
    ind.data[i]=x.data[tind[i]-1:tind[i]+2].mean(0)

ind.name=indtitle

pl.figure(1,figsize=(12,4)); pl.clf()
ind.plot1d(); ax=pl.gca()
ax.set_xlabel(xlab); ax.set_ylabel(units); ax.set_ylim(indylim)
ax.xaxis.set_major_locator(mdates.YearLocator(tint))
pl.grid(); pl.show()
pl.savefig(indpic)

# Positive composite
f=open(posfile); x=pickle.load(f); f.close()
x.name=postitle; x.units=units; x.copts=copts; x.cbar_opts=cbar_opts; x.map=fmap
pl.figure(2); pl.clf()
x.plot_mapfc()
pl.grid(); pl.show()
pl.savefig(pospic)

# Negative composite
f=open(negfile); x=pickle.load(f); f.close()
x.name=negtitle; x.units=units; x.copts=copts; x.cbar_opts=cbar_opts; x.map=fmap
pl.figure(3); pl.clf()
x.plot_mapfc()
pl.grid(); pl.show()
pl.savefig(negpic)
