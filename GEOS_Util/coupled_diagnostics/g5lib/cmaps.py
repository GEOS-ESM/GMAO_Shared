from matplotlib import colors
import scipy as sp

class FilledCmap(colors.Colormap):
    '''
    This colormap fills the values in fill_range with fill_color.
    Range should be normalized. 
    '''

    def __init__(self,cmap,name='filled_map',fill_color='white',fill_range=(0.45,0.55)):
        self.cmap = cmap
        self.fill_color=fill_color
        self.fill_range = fill_range
        self.monochrome=cmap.monochrome
        colors.Colormap.__init__(self,name,cmap.N)

    def __call__(self, X, alpha=1.0, bytes=False):

        rgba=sp.array(self.cmap( X, alpha=1.0, bytes=False))
        rgba[...,0] = sp.where(sp.logical_and(X>=self.fill_range[0],X<=self.fill_range[1]), 1.0, rgba[...,0])
        rgba[...,1] = sp.where(sp.logical_and(X>=self.fill_range[0],X<=self.fill_range[1]), 1.0, rgba[...,1])
        rgba[...,2] = sp.where(sp.logical_and(X>=self.fill_range[0],X<=self.fill_range[1]), 1.0, rgba[...,2])
        return rgba

