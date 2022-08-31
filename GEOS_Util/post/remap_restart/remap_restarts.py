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
from remap_utils import *
from remap_questions import get_config_from_questionary 
from remap_params import *
from remap_upper import *
from remap_lake_landice_saltwater import *
from remap_analysis  import *
from remap_catchANDcn  import *

def main(argv):

  question_flag = False
  cmd2yaml_flag = False
  yaml_flag     = False
  config        = ''

  try:
    opts, args = getopt.getopt(argv,"hc:o", ['config_file='])
  except getopt.GetoptError:
    print('Usage: remap_restarts.py -c remap_params.yaml or ./remap_restarts.py or ./remap_restarts.py -o key=value ....')
    sys.exit('command line error')
  for opt, arg in opts:
    if opt == '-h':
      print('''\nThere are two ways to use this script to remap restarts. \n 
              1) Use an exsiting config file to remap: \n
                ./remap_restarts.py -c my_config.yaml \n \n
              2) Use questionary to convert template remap_params.tpl to \n
                 remap_params.yaml and then remap. \n
                ./remap_restarts.py \n
              3) Use command line to input a flattened yaml file: \n
                ./remap_restarts.py -o input:air:drymass=1 input:air:hydrostatic=0 ... 
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
      yaml_flag = True
      
    if opt == '-o':
      print("\n Use command line to input a flattened yaml config file \n")
      cmd2yaml_flag = True
      config_yaml = 'remap_params.yaml'

  if yaml_flag :
     config  = yaml_to_config(config_yaml)

  if cmd2yaml_flag :
     config  = args_to_config(args)

  if not ( yaml_flag or cmd2yaml_flag) :
     raw_config = get_config_from_questionary()
     params = remap_params(raw_config)
     config = params.config
     question_flag = True  
     config_yaml = 'remap_params.yaml'
 
  print('\n')
  print_config(config)
  
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
 
  if yaml_flag     or question_flag :  write_cmd(config)
  if cmd2yaml_flag or question_flag :  config_to_yaml(config, config_yaml)

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

