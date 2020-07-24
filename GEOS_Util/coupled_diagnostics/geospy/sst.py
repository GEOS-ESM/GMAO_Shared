import importlib
import matplotlib.pyplot as pl
import geosdset

def plot_clim(exp, ds):
    var='TS'
    sst=ds[var]
    sst.mean('time').plot()
    pl.show()

def plot_diff(exp, ds1, ds2):
    pas    

def plots(exps, dsets):
    obs=[]
    plot_clim(exps[0] ,dsets[0])

    for ds in dsets:
        plot_diff(exps[0], dsets[0], ds)
        
    for ds in obs:
        plot_diff(exps[0], dsets[0], ds)

if __name__=='__main__':
    import sys
    exps=geosdset.load_exps(sys.argv[1])
    dsets=geosdset.load_collection(exps,'geosgcm_ocn2d')
    plost(exps,dsets)
    
