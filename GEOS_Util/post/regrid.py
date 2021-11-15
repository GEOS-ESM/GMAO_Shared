#!/usr/bin/env python
#
# source g5_modules
#

import yaml
from regrider_base  import *
from regrider_upper import *
from regrider_surf  import *
from regrider_ana  import *

if __name__ == '__main__' :

   # load input yaml
   stream = open("regrid.yaml", 'r')
   config = yaml.full_load(stream)

   # upper air
   upper = upperair(config)
   upper.regrid()
   
   # surface
   surf  = surface(config)
   surf.regrid()

   # analysis
   ana = analysis(config)
   ana.regrid()
   

