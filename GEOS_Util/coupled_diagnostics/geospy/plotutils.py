import matplotlib.pyplot as pl
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def quickmap(projection=ccrs.PlateCarree()):
    ax=pl.axes(projection=projection)
    ax.set_global()
    ax.coastlines()
    gl=ax.gridlines(draw_labels=True)
    gl.xlabels_top=False
    gl.ylabels_right=False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    return ax

def drawstats(ax, mean='', std=''):
    pass
