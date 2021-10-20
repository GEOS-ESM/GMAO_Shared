#!/usr/bin/env python
#
import os
import subprocess
import shutil
import glob

tags_to_bcs ={"Icarus_Reynolds": "/discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/bcs/Icarus_Updated/Icarus_Reynolds/", 
              "Icarus-NLv3_MERRA-2": "/discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/bcs/Icarus-NLv3/Icarus-NLv3_MERRA-2/"}

def grid_to_string(grid):
  if(grid[0] == 'C'):
    n = int(grid[1:])
    return f'CF{n:04d}x6C' 
  else:
    xy = grid.upper().split('X')
    x = int(xy[0])
    y = int(xy[1])
    if (y/x == 6) :
      return f'CF{x:04d}x6C'
    else:
      return f'DE{x:04d}xPE{y:04d}'

class regrider:
  def __init__(self, config):
     self.common_in  = config['options_in']['COMMON']
     self.common_out = config['options_out']['COMMON']
     self.restarts_in = config['restarts_in']
     self.slurm_options = config['slurm_options']

     # deduce the options
     agrid_str = grid_to_string(self.common_in['agrid'])
     ogrid_str = grid_to_string(self.common_in['ogrid'])
     ao_str = agrid_str + '_' + ogrid_str
     self.in_bcsdir = tags_to_bcs[self.common_in['tag']] + ao_str
     self.in_til    = self.in_bcsdir + '/'+ao_str + '-Pfafstetter.til' 
     
     agrid_str = grid_to_string(self.common_out['agrid'])
     ogrid_str = grid_to_string(self.common_out['ogrid'])
     ao_str = agrid_str + '_' + ogrid_str
     self.out_bcsdir = tags_to_bcs[self.common_out['tag']] + ao_str
     self.out_til    = self.out_bcsdir + '/'+ao_str + '-Pfafstetter.til' 

