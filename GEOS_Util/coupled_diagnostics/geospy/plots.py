'''
Some classes for standard plots.
'''

import matplotlib.pyplot as pl
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

class Plot(object):
    '''
    Abstract plot class.
    '''
    pass

class Plot1d(Plot):
    '''
    Class for 1d plots.
    '''
    pass

class Plot2d(Plot):
    '''
    Class for general 2d plots.
    '''

    def __init__(self, fill_opts={}, contour_opts={}):
        super(Plot2d,self).__init__()
        self.fill_opts=fill_opts
        self.contour_opts=contour_opts
        
    def contour(self, da, ax=None, mode='both', clab_fmt='%1.0f'):
        '''
        Makes contour plot of data array.
        
        Parameters
        ----------
        da: DataArray
        ax: axes
        mode: 'filled', 'contour' or 'both'
        clab_fmt: format for contour labels
        '''
        
        if ax is None:
            ax=pl.gca()

        if mode in ('filled', 'both'):
            da.plot.contourf(ax=ax,**self.fill_opts)
        
        if mode in ('contour', 'both'):
            cs=da.plot.contour(ax=ax,**self.contour_opts)
            cs.clabel(fmt=clab_fmt)
        
        return ax
        
class PlotMap(Plot2d):
    '''
    Class for map plots.
    '''
    
    def __init__(self, projection=ccrs.PlateCarree(), fill_opts={}, contour_opts={}):
        super(PlotMap,self).__init__(fill_opts, contour_opts)
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

    def contour(self, da, ax=None, mode='both', clab_fmt='%1.0f'):
        '''
        Makes contour plot of data array.
        
        Parameters
        ----------
        da: DataArray
        ax: axes
        mode: 'filled', 'contour' or 'both'
        clab_fmt: format for contour labels
        '''
        if ax is None:
            ax=self.plot_map()

        super(PlotMap,self).contour(da, ax=ax, mode=mode, clab_fmt=clab_fmt)
        
        return ax
