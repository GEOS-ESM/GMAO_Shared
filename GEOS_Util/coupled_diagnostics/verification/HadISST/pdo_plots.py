import os, pickle
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import scipy as sp
import mpl_toolkits.basemap as bm
from mpl_toolkits.basemap.cm import sstanom
from g5lib import plotters

dsetname='HadISST'
varname='sst'
indname='pdo'
path=os.environ['NOBACKUP']+'/verification/'+dsetname
indfile=path+'/data/'+varname+'_pc1_np.dat'
indpic=path+'/pics/'+varname+'_pc1_np.png'
indtitle='North Pacific HadlSST PC1 (DJF)';xlab='years'
indylim=(-3,3); tint=10
posfile=path+'/data/'+varname+'_'+indname+'_plus.dat'
pospic=path+'/pics/'+varname+'_'+indname+'_plus.png'
postitle='HadlSST > PDO_std'
negfile=path+'/data/'+varname+'_'+indname+'_minus.dat'
negpic=path+'/pics/'+varname+'_'+indname+'_minus.png'
negtitle='HadlSST < -PDO_std'
units='$^0$C'
copts={'levels': sp.arange(-1,1.1,0.1),\
           'cmap': sstanom}
cbar_opts={'orientation': 'vertical'}

try:
    os.makedirs(path+'/pics')
except OSError:
    pass

# Plot index
f=open(indfile); ind=pickle.load(f); f.close()

# Normalize
s=ind.ave(0,ret_std=True)[1].data
ind.data/=s;

ind.name=indtitle

pl.figure(1,figsize=(12,4)); pl.clf()
ind.d(); ax=pl.gca()
ax.set_xlabel(xlab); ax.set_ylabel(units); ax.set_ylim(indylim)
ax.xaxis.set_major_locator(mdates.YearLocator(tint))
pl.grid(); pl.show()
pl.savefig(indpic)

# Positive composite
f=open(posfile); pos=pickle.load(f); f.close()
pp=plotters.GeoPlotter()
pos.name=postitle; pos.units=units; pp.copts.update(copts); pp.cbar_opts.update(cbar_opts)
pl.figure(2); pl.clf()
pp(pos)
pl.grid(); pl.show()
pl.savefig(pospic)

# Negative composite
f=open(negfile); neg=pickle.load(f); f.close()
neg.name=negtitle; neg.units=units; 
pl.figure(3); pl.clf()
pp(neg)
pl.grid(); pl.show()
pl.savefig(negpic)
