#!/usr/bin/env python3
#
# source install/bin/g5_modules
#
# Newer GEOS code should load a module with GEOSpyD Python3 if not run:
#   module load python/GEOSpyD/Min4.10.3_py3.9
#

import sys, getopt
import ruamel.yaml
import questionary
from regrid_questions import get_config_from_questionary 
from regrid_params import *
from regrid_upper import *
from regrid_lake_landice_saltwater import *
from regrid_analysis  import *
from regrid_catchANDcn  import *

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
                 regrid_params.yaml and then regrid. \n
                ./regrid.py \n
              \nHelp message: \n
              1) Each individual script can be executed independently
              2) regrid_questions.py generates raw_answer.yaml
              3) regrid_params.py uses raw_answer.yaml and regrid_params.tpl as inputs and generates regrid_params.yaml
              4) regrid_upper.py uses regrid_params.yaml as input for regriding
              5) regrid_lake_landice_saltwater.py uses regrid_params.yaml as input for regriding
              6) regrid_catchANDcn.py uses regrid_params.yaml as input for regriding
              7) regrid_analysis.py uses regrid_params.yaml as input for regriding ''')
      sys.exit(0)
    if opt in("-c", "--config_file"):
      config_yaml = arg

  if config_yaml == '':
    config = get_config_from_questionary()
    params = regrid_params(config)
    params.convert_to_yaml()
    config_yaml = 'regrid_params.yaml'

  with  open(config_yaml, 'r') as f:
    for line in f.readlines():
      trimmed_line = line.rstrip()
      if trimmed_line: # Don't print blank lines
          print(trimmed_line)

  print('\n')
  questions = [
        {
            "type": "confirm",
            "name": "Continue",
            "message": "Above is the YAML config file, would you like to continue?",
            "default": True
        },]  
  answer = questionary.prompt(questions)

  if not answer['Continue'] :
     print("\nYou answered not to continue, exiting.\n")
     sys.exit(0)
  # upper air
  upper = upperair(config_yaml)
  upper.regrid()

  # lake, landice and saltwater
  lls  = lake_landice_saltwater(config_yaml)
  lls.regrid()

  # catchANDcn
  catch  = catchANDcn(config_yaml)
  catch.regrid()

  # analysis
  ana = analysis(config_yaml)
  ana.regrid()

if __name__ == '__main__' :
  main(sys.argv[1:])

