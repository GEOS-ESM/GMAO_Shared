'''
This module defines classes derived from basemap.Basemap for drawing maps.
'''
import mpl_toolkits.basemap as bm
import scipy as sp

class BaseMapper(bm.Basemap):
    '''
    A class for drawing general maps.
    '''

    def __init__(self, mopts={}, popts={}, fill=False,**bopts):
        '''
        bopts - options for bm.Basemap
        mopts - options for bm.Basemap.drawmeridians
        popts - options for bm.Basemap.drawparrallels
        '''
    
        super(BaseMapper, self).__init__(**bopts)
        self.mopts=dict(meridians=sp.arange(-360,361,60), labels=(0,0,0,1))
        self.mopts.update(mopts)
        self.popts=dict(circles=sp.arange(-90,91,30), labels=(1,0,0,0))
        self.popts.update(popts)
        self.fill=fill # Fills continetns, default=False

    def draw(self):
        '''
        Draw map
        '''
        if self.fill:
            self.fillcontinents()
        self.drawcoastlines()
        self.drawparallels(**self.popts)
        self.drawmeridians(**self.mopts)

def mapper(grid,**bopts):
    '''
    Returns a BaseMapper object based on grid domain. 
    '''

    opts=dict(llcrnrlat=grid['lat'].min(),

              llcrnrlon=grid['lon'].min(),
              urcrnrlat=grid['lat'].max(),
              urcrnrlon=grid['lon'].max())
    opts.update(bopts)

    return BaseMapper(**opts)
