'''
This module define a converter from world-time coordinates to grid indexes.
'''

import scipy as sp
import my_lib.utils as utl
import matplotlib.dates as mdt
from dateutil import parser

class Domain(object):

    def __init__(self,lons=(0.,360.), lats=(-90.,90.), levs=(0,), \
                     dates=(mdt.datetime.datetime.now(),)):
        '''
        Provide the domain extension
        lons - tuple of longitudes (eastern boundary,  western boundary)
        lats - tuple of latitudes (southern boundary, northern boundary)
        levs - tuple of levels (min level, max level)
        dates -  tuple (start time, end time) of datetime objects or strings in format 'yyyy-mo-dy hr:mn:sc' 

        If any of tuples contain only one entry, the index of nearest coordinate will be
        returned by __call__.

        Default is global domain.
        '''
        self.lons=lons
        self.lats=lats
        self.levs=levs
        dd=[]
        for tt in dates:
            if isinstance(tt,str): 
                dd.append(parser.parse(tt,yearfirst=True))
            else:
                dd.append(tt)
        self.dates=tuple(dd) 

    def __call__(self,grid,time):
        '''
        Returns a dictionary, containing indices for a given grid, which correspond to self domain.
        Result is suitable for dset.fromfile
        '''

        iind1=utl.find_nearest_index(grid['lon'][0],self.lons[0])
        jind1=utl.find_nearest_index(grid['lat'][:,0],self.lats[0])
        kind1=utl.find_nearest_index(grid['lev'],self.levs[0])
        tind1=utl.find_nearest_index(mdt.date2num(time),mdt.date2num(self.dates[0]))

        iind2=utl.find_nearest_index(grid['lon'][0],self.lons[1]) if  len(self.lons)>1 else iind1
        jind2=utl.find_nearest_index(grid['lat'][:,0],self.lats[1]) if len(self.lats)>1 else jind1
        kind2=utl.find_nearest_index(grid['lev'],self.levs[1]) if len(self.levs)>1 else kind1
        tind2=utl.find_nearest_index(mdt.date2num(time),mdt.date2num(self.dates[1])) \
            if len(self.dates)>1 else tind1

        iind=slice(iind1,iind2+1,1) if iind2>=iind1 else slice(iind2,iind1+1,1)
        jind=slice(jind1,jind2+1,1) if jind2>=jind1 else slice(jind2,jind1+1,1)
        kind=slice(kind1,kind2+1,1) if kind2>=kind1 else slice(kind2,kind1+1,1)
        tind=slice(tind1,tind2+1,1)
        return {'iind': iind,\
                'jind': jind,\
                'kind': kind,\
                'tind': tind}

    @property
    def dict(self):
        '''
        Returns a dictionary with domain span for convenience
        '''
        return {'lons': self.lons,\
                    'lats': self.lats,\
                    'levs': self.levs,\
                    'dates': self.dates}

    def __repr__(self):
        return repr(self.dict)

# Put some preset domains here (i.e. Nino3,3.4,4, etc)


