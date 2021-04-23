'''
Some classes for standard plots.
'''

import numpy as np
import matplotlib.pyplot as pl
import matplotlib.ticker as mticker
import cartopy.crs as ccrs

class Plot(object):
    '''
    Abstract plot class.
    '''
    pass

class Plot1d(Plot):
    '''
    Class for 1d plots.
    '''
    def __init__(self, **kwargs):
        super(Plot1d,self).__init__()
        self.line_opts=kwagrs.get('line_opts',{})

    def line(self, da, fmt='-',ax=None):
        if ax is None:
            ax=pl.gca()

        da.plot.line(fmt,ax=ax,**self.line_opts)
        return ax

class Plot2d(Plot):
    '''
    Class for general 2d plots.
    '''

    def __init__(self, **kwargs):
        super(Plot2d,self).__init__()
        self.fill_opts=kwargs.get('fill_opts',{})
        self.contour_opts=kwargs.get('contour_opts',{})
        self.quiver_opts=kwargs.get('quiver_opts',{})
        self.clab_opts={'fmt': '%1.1f'}
        self.clab_opts.update(kwargs.get('clab_opts',{}))
        
    def contour(self, da, ax=None, mode='both',stat=None):
        '''
        Makes contour plot of data array.
        
        Parameters
        ----------
        da: DataArray
        ax: axes
        mode: 'filled', 'contour' or 'both'
        '''
        
        if ax is None:
            ax=pl.gca()

        if mode in ('filled', 'both'):
            cs=da.plot.contourf(ax=ax,**self.fill_opts)
        
        if mode in ('contour', 'both'):
            cs=da.plot.contour(ax=ax,**self.contour_opts)
            cs.clabel(**self.clab_opts)

        if stat is not None:
            ylim=ax.get_ylim(); dely=(max(ylim)-min(ylim))/20
            xlim=ax.get_xlim(); delx=(max(xlim)-min(xlim))/10
            ax.text(min(xlim)-delx,max(ylim)+dely,stat)
        
        return ax

    def quiver(self, ds, x, y, u, v, ax=None):
        '''
        Makes quiver plot of vector data.
        
        Parameters
        ----------
        ds: Dataset with vector data components
        x: name of x axis
        y: name of y axis
        u: name of U component
        v: name of V component
        ax: axes
        '''
        
        if ax is None:
            ax=pl.gca()

        skip=int(ds[u].shape[0]/20)
        if ds[x].ndim==1:
            x,y=np.meshgrid(ds[x].values,ds[y].values)
        else:
            x,y=ds[x].values,ds[y].values

        args=[x[::skip,::skip],
              y[::skip,::skip],
              ds[u].values[::skip,::skip],
              ds[v].values[::skip,::skip]]
        cs=pl.quiver(*args,**self.quiver_opts)

        return ax
        
class PlotMap(Plot2d):
    '''
    Class for map plots.
    '''
    
    def __init__(self, projection=ccrs.PlateCarree(), **kwargs):
        super(PlotMap,self).__init__(**kwargs)
        self.projection=projection
        self.fill_opts['transform']=ccrs.PlateCarree()
        self.contour_opts['transform']=ccrs.PlateCarree()
        self.quiver_opts['transform']=ccrs.PlateCarree()

    def plot_map(self):
        ax=pl.axes(projection=self.projection)
        ax.set_global()
        ax.coastlines()
        gl=ax.gridlines(draw_labels=True)
        gl.top_labels=False
        gl.right_labels=False
        
        return ax

    def contour(self, da, ax=None, mode='both', stat=None):
        '''
        Makes contour plot of data array.
        
        Parameters
        ----------
        da: DataArray
        ax: axes
        mode: 'filled', 'contour' or 'both'
        '''
        if ax is None:
            ax=self.plot_map()

        super(PlotMap,self).contour(da, ax=ax, mode=mode, stat=stat)
        
        return ax

    def quiver(self, ds, x, y, u, v, ax=None):
        '''
        Makes quiver plot of vector data.
        
        Parameters
        ----------
        ds: Dataset with vector data components
        x: name of x axis
        y: name of y axis
        u: name of U component
        v: name of V component
        ax: axes
        '''
        
        if ax is None:
            ax=self.plot_map()

        super(PlotMap,self).quiver(ds, x, y, u, v, ax=ax)

        return ax

# Longitude and Latitude formatters 
def format_lon(lon,pos=None):
    lon = lon+360 if lon<-180 else lon-360 if lon>180 else lon

    if lon == -180.:
        return f'{-lon:g}'
    elif (lon > -180.) and (lon < 0.):
        return f'{-lon:g}W'
    elif (lon == 0.) or (lon == 180.):
        return f'{lon:g}'
    else:
        return f'{lon:g}E'
    
def format_lat(lat,pos=None):
    return f'{-lat:g}S' if lat<0 else f'{lat:g}N' if lat>0 else 'EQ'

LONGITUDE_FORMATTER=mticker.FuncFormatter(format_lon)
LATITUDE_FORMATTER=mticker.FuncFormatter(format_lat)

MONTH_LOCATOR=mticker.FixedLocator(range(1,13))
MONTH_FORMATTER=mticker.FixedFormatter(('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
