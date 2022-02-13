#!/usr/bin/env python3
#
# source g5_modules
# module load python/GEOSpyD/Min4.9.2_py3.9 
#

import sys, getopt
import yaml
from regrid_questions import get_config_from_questionary 
from regrider_base  import *
from regrider_upper import *
from regrider_surf  import *
from regrider_ana  import *

def main(argv):
  config_yaml = ''
  try:
    opts, args = getopt.getopt(argv,"hc:", ['config_file='])
  except getopt.GetoptError:
    print('Usage: regrid.py -c regrid.yaml or ./regrid.py ')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print('\nUse a config file to regrid: \n'
            './regrid.py -c my_config.yaml \n \n'
            'Use questionary to generate regrid.yaml and regrid \n'
            './regrid.py \n'
            '\nHelp message: \n'
            '1) .... \n'
            '2) .... ')
      sys.exit()
    if opt in("-c", "--config_file"):
      config_yaml = arg

  if config_yaml == '':
    get_config_from_questionary()
    config_yaml = 'regrid.yaml'

  # load input yaml
  stream = open(config_yaml, 'r')
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

if __name__ == '__main__' :
  main(sys.argv[1:])

