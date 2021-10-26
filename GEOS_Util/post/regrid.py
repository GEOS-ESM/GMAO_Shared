#!/usr/bin/env python
#
# module load python/GEOSpyD/Min4.9.2_py3.9 
#

import yaml
from regrider_base  import *
from regrider_upper import *
from regrider_surf  import *

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

   # what ever
#   forcing = force(config)
#   forcing.regrid()
   

