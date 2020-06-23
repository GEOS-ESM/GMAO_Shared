import os, pickle
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import scipy as sp
from mpl_toolkits.basemap.cm import sstanom
from g5lib import plotters

dsetname='HadISST'
varname='sst'
path=os.environ['NOBACKUP']+'/verification/'+dsetname
indfile=path+'/data/'+varname+'_pc1.dat'
indpic=path+'/pics/'+varname+'_pc1.png'
indtitle='HadISST PC1'; xlab='years'
indylim=(-3,3); tint=10
eoffile=path+'/data/'+varname+'_eof1.dat'
eofpic=path+'/pics/'+varname+'_eof1.png'
units='$^0$C'
copts={'levels': sp.arange(-1,1.1,0.1),\
           'cmap': sstanom}
cbar_opts={'orientation': 'vertical'}

try:
    os.makedirs(path+'/pics')
except OSError:
    pass

# read data
f=open(indfile); pc=pickle.load(f); f.close()
f=open(eoffile); eof=pickle.load(f); f.close()

# Normalize pc
s=pc.ave(0,ret_std=True)[1].data
pc.data/=s
eof.data*=s

pc.name=indtitle

pl.figure(1,figsize=(12,4)); pl.clf()
pc.d(); ax=pl.gca()
ax.set_xlabel(xlab); ax.set_ylim(indylim)
ax.xaxis.set_major_locator(mdates.YearLocator(tint))
pl.grid(); pl.show()
pl.savefig(indpic)

# Plot EOF1
pp=plotters.GeoPlotter()
eof.units=units; pp.copts.update(copts); pp.cbar_opts.update(cbar_opts)
pl.figure(2); pl.clf()
pp(eof)
pl.grid(); pl.show()
pl.savefig(eofpic)

