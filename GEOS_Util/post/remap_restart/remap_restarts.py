#!/usr/bin/env python3
#
# source install/bin/g5_modules
#
# Newer GEOS code should load a module with GEOSpyD Python3 if not run:
#   module load python/GEOSpyD/Min4.11.0_py3.9
#

import sys
import argparse
import textwrap
import ruamel.yaml
import questionary
from remap_utils import *
from remap_questions import get_config_from_questionary
from remap_params import *
from remap_upper import *
from remap_lake_landice_saltwater import *
from remap_analysis  import *
from remap_catchANDcn  import *

# Define the argument parser
def parse_args():

    program_description = textwrap.dedent(f'''
      USAGE:

      There are three ways to use this script to remap restarts.

      1. Use an existing config file to remap:
           ./remap_restarts.py -c my_config.yaml

      2. Use questionary to convert template remap_params.tpl to
         remap_params.yaml and then remap:
           ./remap_restarts.py

      3. Use command line to input a flattened yaml file:
           ./remap_restarts.py -o input:air:drymass=1 input:air:hydrostatic=0 ...

      NOTE: Each individual script can be executed independently
        1. remap_questions.py generates raw_answer.yaml
        2. remap_params.py uses raw_answer.yaml and remap_params.tpl as inputs and generates remap_params.yaml
        3. remap_upper.py uses remap_params.yaml as input for remapping
        4. remap_lake_landice_saltwater.py uses remap_params.yaml as input for remapping
        5. remap_catchANDcn.py uses remap_params.yaml as input for remapping
        6. remap_analysis.py uses remap_params.yaml as input for remapping
    ''')

    parser = argparse.ArgumentParser(description='Remap restarts',epilog=program_description,formatter_class=argparse.RawDescriptionHelpFormatter)
    # define a mutually exclusive group of arguments
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-c', '--config_file', help='YAML config file')
    group.add_argument('-o', '--flattened_yaml', help='Flattened YAML config', metavar='input:air:drymass=1 input:air:hydrostatic=0 ...')

    # Parse using parse_known_args so we can pass the rest to the remap scripts
    # If config_file is used, then extra_args will be empty
    # If flattened_yaml is used, then extra_args will be populated
    args, extra_args = parser.parse_known_args()
    return args, extra_args

def main():

  question_flag = False
  config        = ''

  # Parse the command line arguments from parse_args() capturing the arguments and the rest
  command_line_args, extra_args = parse_args()
  config_yaml = command_line_args.config_file
  flattened_yaml = command_line_args.flattened_yaml

  if config_yaml:
      config = yaml_to_config(config_yaml)
  elif flattened_yaml:
      config_yaml = 'remap_params.yaml'
      config  = args_to_config(extra_args)
  else:
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

  if config_yaml or question_flag: write_cmd(config)
  if flattened_yaml or question_flag: config_to_yaml(config, config_yaml)

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
  main()

