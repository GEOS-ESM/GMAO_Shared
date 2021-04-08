'''
Some classes for standard plots.
'''

import matplotlib.pyplot as pl
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

class Plot(object):
    '''
    Abstract plot class.
    '''
    pass

class Plot1d(Plot):
    '''
    Class for 1d plots.
    '''
    def __init__(self, line_opts={}):
        super(Plot1d,self).__init__()
        self.line_opts=line_opts

    def line(self, da, fmt='-',ax=None):
        if ax is None:
            ax=pl.gca()

        da.plot.line(fmt,ax=ax,**self.line_opts)
        return ax

class Plot2d(Plot):
    '''
    Class for general 2d plots.
    '''

    def __init__(self, fill_opts={}, contour_opts={}, clab_opts={}):
        super(Plot2d,self).__init__()
        self.fill_opts=fill_opts
        self.contour_opts=contour_opts
        self.clab_opts={'fmt': '%1.1f'}
        self.clab_opts.update(clab_opts)
        
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
        
class PlotMap(Plot2d):
    '''
    Class for map plots.
    '''
    
    def __init__(self, projection=ccrs.PlateCarree(), fill_opts={}, contour_opts={}, clab_opts={}):
        super(PlotMap,self).__init__(fill_opts, contour_opts, clab_opts)
        self.projection=projection
        self.fill_opts['transform']=ccrs.PlateCarree()
        self.contour_opts['transform']=ccrs.PlateCarree()

    def plot_map(self):
        ax=pl.axes(projection=self.projection)
        ax.set_global()
        ax.coastlines()
        gl=ax.gridlines(draw_labels=True)
        gl.top_labels=False
        gl.right_labels=False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        
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


# Longitude and Latitude formatters 
def format_lon(lon,pos=None):
    if (lon >= -360.) and (lon < -180.):
        lon+=360.
    elif (lon > 180.) and (lon <= 360.):
        lon-=360.
    if lon == -180.:
        return '%i' %-lon
    elif (lon > -180.) and (lon < 0.):
        return '%iW' %-lon
    elif (lon == 0.) or (lon == 180.):
        return '%i' %lon
    else:
        return '%iE' %lon
    
def format_lat(lat,pos=None):
    if lat < 0:
        return f'{-lat:g}S'
    elif lat > 0:
        return f'{lat:g}N'
    else:
        return 'EQ'

LONGITUDE_FORMATTER=mticker.FuncFormatter(format_lon)
LATITUDE_FORMATTER=mticker.FuncFormatter(format_lat)
