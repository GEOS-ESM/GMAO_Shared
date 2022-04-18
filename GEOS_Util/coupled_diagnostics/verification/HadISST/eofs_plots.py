import os, pickle,sys
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import scipy as sp
from mpl_toolkits.basemap.cm import sstanom
from g5lib import plotters as pltrs
from datetime import datetime

dsetname='HadISST'
varname='sst'
path=os.environ['NOBACKUP']+'/verification/'+dsetname

varname='sst'

try:
    os.makedirs(path+'/pics')
except OSError:
    pass

f=open(path+'/data/sst_eofs.dat'); var=pickle.load(f); f.close()
var._pc.varimax(nr=76)


def eof_plot(var,N):
    copts={'levels': sp.arange(-1,1.1,0.2),
           'cmap': sstanom}
    cbar_opts={'orientation': 'vertical'}
    xlabel='time, yrs'
    ylabel='SST $^0$C'
    indylim=(-5,5)
    tint=10

    pc,eof=var.get_pc(N)
    
    pc.data/=sp.sqrt(var._pc.evals[N-1])
    eof.data*=sp.sqrt(var._pc.evals[N-1])

    eof.name=eof.name.replace('EOF','REOF')
    pc.name=pc.name.replace('PC','RPC')
    
    pl.figure(1,figsize=(12,4)); pl.clf()
    p=pltrs.Plotter1d(); p(pc)
    ax=p.axis
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_ylim(indylim)
    ax.xaxis.set_major_locator(mdates.YearLocator(tint))
    pl.grid(); pl.show()
    pl.savefig(path+'/pics/'+varname+'_rpc'+str(N)+'.png')
    
    # Plot EOF
    eof.shiftgrid(30.)
    p=pltrs.GeoPlotter(copts=copts, cbar_opts=cbar_opts)
    pl.figure(2); pl.clf()
    p(eof)
    pl.grid(); pl.show()
    pl.savefig(path+'/pics/'+varname+'_reof'+str(N)+'.png')
    
for n in xrange(1,11):
    print 'plotting EOF '+str(n)+'\n'
    eof_plot(var,n)
    raw_input('Press ENTER')

