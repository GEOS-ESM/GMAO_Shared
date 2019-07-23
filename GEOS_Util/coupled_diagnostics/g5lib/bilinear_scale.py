import scipy as sp
from matplotlib import scale as mscale
from matplotlib import transforms as mtransforms

class BiLinearScale(mscale.ScaleBase):
    """
    Allows to put depth profiles on irregular vertical scale
    """
    
    name = 'bilinear'


    def __init__(self, axis, threshold, multiplier):
        """
        threshold - threshold depth
        multiplier - scale factor for upper part
        """

        super(BiLinearScale,self).__init__()

        self.thresh = threshold
        self.mult   = multiplier

    def get_transform(self):
        """
        Override this method to return a new instance that does the
        actual transformation of the data.
        """
        return self.Transform(self.thresh, self.mult)

    def set_default_locators_and_formatters(self, axis):
        """
        Override to set up the locators and formatters to use with the
        scale.
        """
        pass


    class Transform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, thresh, mult):
            super(BiLinearScale.Transform,self).__init__()
            self.thresh = thresh
            self.mult   = mult

        def transform(self, a):
            """
            This transform takes an Nx1 ``numpy`` array and returns a
            transformed copy.
            """

#            aa=a.copy()
            a[a<self.thresh]*=self.mult
            a[a>=self.thresh]+=self.thresh*(self.mult-1)
            return a

        def inverted(self):
            """
            Override this method so matplotlib knows how to get the
            inverse transform for this transform.
            """
            return BiLinearScale.InvertedTransform(self.thresh, self.mult)

    class InvertedTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, thresh, mult):
            super(BiLinearScale.InvertedTransform,self).__init__()
            self.thresh = thresh
            self.mult   = mult

        def transform(self, a):
#            aa=a.copy()
            a[a<self.thresh*self.mult]/=self.mult
            a[a>=self.thresh*self.mult]-=self.thresh*(self.mult-1)
            return a

        def inverted(self):
            return BiLinearScale.Transform(self.thresh)

