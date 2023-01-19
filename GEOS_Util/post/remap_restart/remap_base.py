#!/usr/bin/env python3
#
import os
import ruamel.yaml
import shutil
import subprocess

class remap_base(object):
  def __init__(self, **configs):
     for key, value in configs.items():
        if (key == 'params_file'):
          print( "use Config yaml file: " + value)
          yaml = ruamel.yaml.YAML()
          stream =''
          with  open(value, 'r') as f:
            stream = f.read()
          self.config = yaml.load(stream)
          out_dir    = self.config['output']['shared']['out_dir']
          if not os.path.exists(out_dir) : os.makedirs(out_dir)
          f = os.path.basename(value)
          dest = out_dir+'/'+f
          try:
            shutil.copy(value, dest)
          except shutil.SameFileError:
            pass
        if (key == 'config_obj'):
          print( "use Config obj")
          self.config = value
          out_dir    = self.config['output']['shared']['out_dir']
          if not os.path.exists(out_dir) : os.makedirs(out_dir)
        break
  def remove_merra2(self):
    if self.config['input']['shared']['MERRA-2']:
      print(" remove temporary folder that contains MERRA-2 archived files ... \n")
      subprocess.call(['/bin/rm', '-rf', self.config['input']['shared']['rst_dir']])
