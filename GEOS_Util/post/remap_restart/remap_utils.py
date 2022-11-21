#!/usr/bin/env python3
#
import os
import subprocess
import ruamel.yaml
from collections import OrderedDict
import shutil
import questionary

def config_to_yaml(config, yaml_file):
   if os.path.exists(yaml_file) :
      overwrite = questionary.confirm("Do you want to overwrite " + yaml_file + "?" , default=False).ask()
      if not overwrite :
         while True:
           new_name = questionary.text("What's the backup name?  ", default=yaml_file +'.1').ask()
           if os.path.exists(new_name):
              print('\n'+ new_name + ' exists, please enter a new one. \n')
           else:
              shutil.move(yaml_file, new_name)
              break
   yaml = ruamel.yaml.YAML()
   with open(yaml_file, "w") as f:
      yaml.dump(config, f)

def yaml_to_config(yaml_file):
   yaml = ruamel.yaml.YAML()
   stream =''
   with  open(yaml_file, 'r') as f:
     stream = f.read()
   config = yaml.load(stream)
   return config

def write_cmd(config) :
   def flatten_nested(nested_dict, result=None, prefix=''):
     if result is None:
       result = dict()
     for k, v in nested_dict.items():
       new_k = ':'.join((prefix, k)) if prefix else k
       if not (isinstance(v, dict) or isinstance(v, OrderedDict)):
         result.update({new_k: v})
       else:
         flatten_nested(v, result, new_k)
     return result

   out_dir = config['output']['shared']['out_dir']
   if not os.path.exists(out_dir) : os.makedirs(out_dir)
   bin_path = os.path.dirname(os.path.realpath(__file__))

   cmd = '#!/usr/local/bin/csh \n'
   cmd = cmd + 'set BINPATH=' + bin_path +'\n'
   cmd = cmd + 'source $BINPATH/g5_modules \n'

   flat_dict = flatten_nested(config)

   k = 1
   for key, value in flat_dict.items():
     if isinstance(value, int) or isinstance(value, float) or isinstance(value, bool) or isinstance(value, type(None)): value = str(value)
     if k == 1:
       cmd = cmd + 'set FLAT_YAML="' + key+"="+ value+ '"\n'
     else:
       cmd = cmd + 'set FLAT_YAML="$FLAT_YAML '+ key+"="+ value+ '"\n'
     k = k+1

   cmd = cmd + '$BINPATH/remap_restarts.py -o $FLAT_YAML'

   with open(out_dir + '/remap_restarts.CMD', 'w') as f:
     f.write(cmd)
   subprocess.call(['chmod', '+x',out_dir + '/remap_restarts.CMD'])

def args_to_config(args):
   # template file should be with this util file
   remap_tpl = os.path.dirname(os.path.realpath(__file__)) + '/remap_params.tpl'
   config = yaml_to_config(remap_tpl)
   # fill in the config with args
   for values in args:
     [keys, value] = values.split("=")
     key = keys.split(':')
     if value.lower() in ['false', 'null', 'none'] :
        value = False
     elif value.lower() in ['true'] :
        value = True
     if len(key) == 2:
        config[key[0]][key[1]] = value
     if len(key) == 3:
        config[key[0]][key[1]][key[2]] = value

   return config

def print_config( config, indent = 0 ):
   for k, v in config.items():
     if isinstance(v, dict):
        print("   " * indent, f"{k}:")
        print_config(v, indent+1)
     else:
        print("   " * indent, f"{k}: {v}")

def merra2_expid(config):
   if not config['MERRA-2']:
      return config
   yyyymm = int(config.get('yyyymmddhh')[0:6])
   if yyyymm < 197901 :
      exit("Error. MERRA-2 data < 1979 not available\n")
   elif (yyyymm < 199201):
      expid = "d5124_m2_jan79"
   elif (yyyymm < 200106):
      expid = "d5124_m2_jan91"
   elif (yyyymm < 202106):
      expid = "d5124_m2_jan00"
   # There was a rewind in MERRA2 from Jun 2021 to Sept 2021
   elif (yyyymm < 202110):
      expid = "d5124_m2_jun21"
   else:
      expid = "d5124_m2_jan10"
   config['expid'] = expid

   return config

if __name__ == '__main__' :
   config = yaml_to_config('c24Toc12.yaml')
   print_config(config)
