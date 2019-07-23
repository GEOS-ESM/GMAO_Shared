'''
Contains plotter classes for plotting 1-d, 2-d and geographical data.
'''

import matplotlib.pylab as pl
import matplotlib as mlib
from matplotlib import animation
import scipy as sp
import mappers

class Plotter(object):
    '''
    General plotter class.
    '''
    def __init__(self, method=None, axes=None,):

        self.fig=None
        self.axis=None
        self.axes=axes
        self.method=method

        # Set formatters only for lons and lats, others are handled by matplotlib
        self.formatters=dict(lon=mlib.ticker.FuncFormatter(self.format_lon),
                             lat=mlib.ticker.FuncFormatter(self.format_lat))
        

    def __call__(self, field):

        if self.axes is None: self._find_axes(field)

        self.fig=pl.gcf()
        self.axis=pl.gca()        
        self.axis.set_title(field.name)

        # Set formatters for lons and lats
        try:
            self.axis.xaxis.set_major_formatter(self.formatters[self.axes[0]])
            self.axis.yaxis.set_major_formatter(self.formatters[self.axes[1]])
        except KeyError:
            pass

    def _find_axes(self,field):
        '''
        Finds varying dimensions which will be plotted along the x and y axes.
        Should be overloaded in the subclass.
        '''
        raise NotImplementedError("subclass must implement this method")

    def draw_stat(self,ss):
        ax=self.axis
        ylim=ax.get_ylim(); dely=(max(ylim)-min(ylim))/20
        xlim=ax.get_xlim(); delx=(max(xlim)-min(xlim))/10
        ax.text(min(xlim)-delx,max(ylim)+dely,ss)
        
    def format_lon(self,lon,pos=None):
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
    
    def format_lat(self,lat,pos=None):
        if lat < 0:
            return '%iS' %-lat
        elif lat > 0:
            return '%iN' %lat
        else:
            return 'EQ'



class Plotter1d(Plotter):
    '''
    For 1-d plots
    '''

    def __init__(self, method=pl.plot, axes=None, style='-',lopts={}):
        '''
        method - a plotting function
        axes   - dimension environment, by default inferred from field shape
        style  - line style
        lopts  - keyword arguments for plotting function
        '''
        super(Plotter1d,self).__init__(method, axes)
        self.method=method
        self.style=style
        self.lopts=dict(linewidth=2)
        self.lopts.update(lopts)

    def __call__(self, field):
        '''
        Plot field.
        '''        
        super(Plotter1d,self).__call__(field)

        #Find ploting dimension
        dims=[0]*4
        dimnames=('time', 'lev', 'lat', 'lon')
        dimname=self.axes[0] if self.axes[1]=='data' else self.axes[1]
        ind=dimnames.index(dimname)
        dims[ind]=slice(None)

        # Find plot arguments
        dd=dict(lon=field.grid['lon'][0,:],
                lat=field.grid['lat'][:,0],
                lev=field.grid['lev'],
                time=field.time,
                data=field.data[dims])
        
        args=[]
        args.append(dd[self.axes[0]])
        args.append(dd[self.axes[1]])
        args.append(self.style)

        # Plot
        self.method(*args,**self.lopts)

        # If x-axis is time, apply default time formatting
        if self.axes[0]=='time':
            fig=pl.gcf()
            fig.autofmt_xdate()

    def _find_axes(self,field):
        '''Finds varying dimension which will be plotted along the x or y axis.'''

        dimnames=('time','lev','lat','lon')
        dims=dict(zip(dimnames,field.dims))

        for name in dimnames:
            if dims[name] > 1:
                self.axes=(name,'data')
                return
        

class Plotter2d(Plotter):
    '''
    For 2-d plots.
    '''

    def __init__(self, method=pl.contourf, axes=None, copts={}, cbar_opts={}, clab_opts={}):
        '''
        method - a plotting function
        axes   - dimension environment, by default inferred from field shape
        copts  - keyword arguments for plotting function
        cbar_opts - keyword arguments for colorbar()
        clab_opts - keyword arguments for clabel()
        '''
        
        super(Plotter2d,self).__init__(method, axes)
        
        self.method=method
        if method is pl.contourf:
            self.copts=dict(extend='both')
            self.copts.update(copts)
            self.cbar_opts=None
            if cbar_opts is not None:
                self.cbar_opts=dict(shrink=0.8, extend='both')
                self.cbar_opts.update(cbar_opts)
        else:
            self.copts=dict()
            self.copts.update(copts)
            self.cbar_opts=dict()
            self.cbar_opts.update(cbar_opts)

        self.clab_opts=dict(fmt='%1.1f', colors='black', fontsize=12)
        self.clab_opts.update(clab_opts)
        
    def __call__(self, field, skip=1):
        '''
        Plot field.
        '''
        super(Plotter2d,self).__call__(field)

        dims=[0]*4
        dimnames=('time', 'lev', 'lat', 'lon')
        xind=dimnames.index(self.axes[0])
        dims[xind]=slice(None)
        yind=dimnames.index(self.axes[1])
        dims[yind]=slice(None)
    
        # Find plot arguments
        dd=dict(lon=field.grid['lon'][0,:],
                lat=field.grid['lat'][:,0],
                lev=field.grid['lev'],
                time=field.time)

        args=[]
        args.append(dd[self.axes[0]][::skip])
        args.append(dd[self.axes[1]][::skip])
        
        if xind>yind:
            z=field.data[dims][::skip,::skip]
        else:
            z=sp.ma.transpose(field.data[dims][::skip,::skip])

        args.append(sp.real(z))
        if sp.any(sp.iscomplex(z.view(sp.ndarray))): args.append(sp.imag(z))
        
        cs=self.method(*args, **self.copts)

        if isinstance(cs, mlib.contour.ContourSet):
            if cs.filled:
                if self.cbar_opts is not None:
                    cb=pl.colorbar(**self.cbar_opts)
                    units=getattr(field,'units',"")
                    cb.set_label(units)
            else:
                cs.clabel(**self.clab_opts)

        # If x-axis is time, apply default time formatting
        if self.axes[0]=='time':
            fig=pl.gcf()
            fig.autofmt_xdate()

        return cs
        
    def _find_axes(self,field):
        '''Finds varying dimensions which will be plotted along the x and y axes.'''

        dimnames=('time','lev','lat','lon')
        dims=dict(zip(dimnames,field.dims))
        axes=[]
        for name in reversed(dimnames):
            if dims[name]>1 and len(axes)<2: 
                axes.append(name)
        self.axes=axes

class GeoPlotter(Plotter2d):
    '''
    For map plots.
    '''
    
    def __init__(self, map=None, method=None, copts={}, cbar_opts={}, clab_opts={}):
        
        super(GeoPlotter,self).__init__(method, ('lon','lat'), copts, cbar_opts, clab_opts)
        self.map=map
        
        if type(cbar_opts) is dict:
            self.cbar_opts.update(orientation='horizontal')
        
    def __call__(self, field, skip=1, stat=False):
        
        self.fig=pl.gcf()
        self.axis=pl.gca()        
        self.axis.set_title(field.name)

        if self.map is None: self.map=mappers.mapper(field.grid)
        if self.method is None: 
            self.method=self.map.contourf
            self.copts.update(extend='both')

        self.map.draw()

        args=[]
        X,Y=self.map(field.grid['lon'], field.grid['lat'])
        args.append(X[::skip,::skip])
        args.append(Y[::skip,::skip])
        z=field.data[0,0,::skip,::skip]
        args.append(sp.real(z))
        if sp.any(sp.iscomplexobj(z.view(sp.ndarray))): args.append(sp.imag(z))

        cs=self.method(*args, **self.copts)

        if isinstance(cs, mlib.contour.ContourSet):
            if cs.filled:
                if self.cbar_opts is not None:
                    cb=pl.colorbar(**self.cbar_opts)
                    units=getattr(field,'units',"")
                    cb.set_label(units)
            else:
                cs.clabel(**self.clab_opts)
        
        if stat:
            mean,std=field.aave(ret_std=True)
            ss='mean: %1.2f\nstd: %1.2f' %(mean.data.squeeze(), std.data.squeeze())
            self.draw_stat(ss)
        
        return cs

    def set_map(self,newmap,method):
        self.map=newmap
        self.method=method

class Animator(object):
    '''
    This should animate the plot. To be changed when nccs upgrade matplotlib to 1.2.
    '''
    def __init__(self, plotter,loopdim='time',ani_opts={}):
        """
        plotter is a Plotter object
        loopdim is one of 'time', 'lev', 'lat', 'lon'. Default - 'time'
        ani_opts - different options for animation.FuncAnimation
        """
        self.plotter=plotter
        self.loopdim=loopdim
        self.ani_opts=ani_opts
    
    def __call__(self, field):
        dims={'time':(0,'tind'), 
              'lev': (1,'kind'), 
              'lat': (2,'jind'), 
              'lon': (3,'iind')}
        loopdimind=dims[self.loopdim][0]
        loopindname=dims[self.loopdim][1]

        arg={loopindname: 0}
        img=self.plotter(field.subset(**arg))
        fargs=(field, loopindname) 
        nframes=field.dims[loopdimind]
        ani = animation.FuncAnimation(self.plotter.fig, 
                                      self.__update__, 
                                      fargs=fargs,
                                      frames=nframes,
                                      **self.ani_opts)
        return ani


    def __update__(self,ind,field,loopindname):
        self.plotter.fig.clf()
        arg={loopindname: ind}
        img=self.plotter(field.subset(**arg))
        ss=str(field.time[ind])
        self.plotter.draw_stat(ss) 
        return img

def plotter(dims):
    '''
    Based on field shape
    '''
    
    ndims=sp.sum(sp.array(dims) > 1)

    if ndims<1: raise Exception('Field has no varying dimensions.')
    if ndims>2: raise Exception('Can not plot {0}-d fields.'.format(ndims))
    if ndims==1: return Plotter1d()
    if ndims==2: 
        if dims[2]>1 and dims[3]>1:
            return GeoPlotter()
        else:
            return Plotter2d()

        

    
