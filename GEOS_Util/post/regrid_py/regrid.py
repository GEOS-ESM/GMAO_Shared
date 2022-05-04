#!/usr/bin/env python3
#
# source install/bin/g5_modules
#
# Newer GEOS code should load a module with GEOSpyD Python3 if not run:
#   module load python/GEOSpyD/Min4.10.3_py3.9
#

import sys, getopt
from regrid_questions import get_config_from_questionary 
from regrid_params import *
from regrid_upper import *
from regrid_lake_landice_saltwater import *
from regrid_analysis  import *
from regrid_catchment  import *

def main(argv):
  config_yaml = ''
  try:
    opts, args = getopt.getopt(argv,"hc:", ['config_file='])
  except getopt.GetoptError:
    print('Usage: regrid.py -c regrid.yaml or ./regrid.py ')
    sys.exit('command line error')
  for opt, arg in opts:
    if opt == '-h':
      print('''\nThere are two ways to use this script to regrid restarts. \n 
              1) Use a config file to regrid: \n
                ./regrid.py -c my_config.yaml \n \n
              2) Use questionary to convert template regrid_params.tpl to \n
                 regrid_params.yaml and regrid. \n
                ./regrid.py \n
              \nHelp message: \n
              1) The rst_dir directory should have three sub-directories: \n
                upperair, surface and analysis which contain restart files respectively. \n''')
      sys.exit(0)
    if opt in("-c", "--config_file"):
      config_yaml = arg

  if config_yaml == '':
    config = get_config_from_questionary()
    params = regrid_params(config)
    params.convert_to_yaml()
    config_yaml = 'regrid_params.yaml'

  # upper air
  upper = upperair(config_yaml)
  upper.regrid()

  # lake, landice and saltwater
  lls  = lake_landice_saltwater(config_yaml)
  lls.regrid()

  # catchment
  catch  = catchment(config_yaml)
  catch.regrid()

  # analysis
  ana = analysis(config_yaml)
  ana.regrid()

if __name__ == '__main__' :
  main(sys.argv[1:])

