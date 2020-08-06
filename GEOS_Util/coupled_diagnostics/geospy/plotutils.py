'''
Frequently used plotting utilities.
'''

import matplotlib.pyplot as pl
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def quickmap(projection=ccrs.PlateCarree()):
    ax=pl.axes(projection=projection)
    ax.set_global()
    ax.coastlines()
    gl=ax.gridlines(draw_labels=True)
    gl.top_labels=False
    gl.right_labels=False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    return ax

def drawstats(ax, mean='', std=''):
    pass
