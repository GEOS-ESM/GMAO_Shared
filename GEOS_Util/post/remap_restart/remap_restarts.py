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
from remap_questions import get_config_from_questionary 
from remap_params import *
from remap_upper import *
from remap_lake_landice_saltwater import *
from remap_analysis  import *
from remap_catchANDcn  import *

def main(argv):
  config_yaml = ''
  try:
    opts, args = getopt.getopt(argv,"hc:", ['config_file='])
  except getopt.GetoptError:
    print('Usage: remap_restarts.py -c remap_params.yaml or ./remap_restarts.py ')
    sys.exit('command line error')
  for opt, arg in opts:
    if opt == '-h':
      print('''\nThere are two ways to use this script to remap restarts. \n 
              1) Use an exsiting config file to remap: \n
                ./remap_restarts.py -c my_config.yaml \n \n
              2) Use questionary to convert template remap_params.tpl to \n
                 remap_params.yaml and then remap. \n
                ./remap_restarts.py \n
              \nHelp message: \n
              1) Each individual script can be executed independently
              2) remap_questions.py generates raw_answer.yaml
              3) remap_params.py uses raw_answer.yaml and remap_params.tpl as inputs and generates remap_params.yaml
              4) remap_upper.py uses remap_params.yaml as input for remapping
              5) remap_lake_landice_saltwater.py uses remap_params.yaml as input for remapping
              6) remap_catchANDcn.py uses remap_params.yaml as input for remapping
              7) remap_analysis.py uses remap_params.yaml as input for remapping ''')
      sys.exit(0)
    if opt in("-c", "--config_file"):
      config_yaml = arg

  params = ''
  if config_yaml == '':
    config = get_config_from_questionary()
    params = remap_params(config)
    params.convert_to_yaml()
    config_yaml = 'remap_params.yaml'

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
  
  # copy merra2 files from archives
  if params:
     params.copy_merra2()

  # upper air
  upper = upperair(params_file=config_yaml)
  upper.remap()

  # lake, landice and saltwater
  lls  = lake_landice_saltwater(params_file=config_yaml)
  lls.remap()

  # catchANDcn
  catch  = catchANDcn(params_file=config_yaml)
  catch.remap()

  # analysis
  ana = analysis(params_file=config_yaml)
  ana.remap()

if __name__ == '__main__' :
  main(sys.argv[1:])

