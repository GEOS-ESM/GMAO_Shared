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
import glob
import remap_restarts 
from remap_questions import get_config_from_questionary
from remap_params import *
from remap_upper import *
from remap_lake_landice_saltwater import *
from remap_analysis  import *
from remap_catchANDcn  import *

def compare(base, result):
  #1) comparing nc4
  bases   =  glob.glob(base + '/*_rst*.nc4')
  results =  glob.glob( result + '/*_rst*.nc4')
  assert ( len(bases) == len(results), " number of restart out should be the same")

  for f in bases:
     bname = os.basename(f)
     core = banme.split('_rst')[0].split('.')[-1]+'_rst'
     for r in results:
        if core in r :
          cmd = 'nccmp -dmgfs '+ f + ' ' + r
          p = subprocess.Popen(shlex.split(cmd), universal_newlines=True, stdout=sp.PIPE, stderr=sp.PIPE)
          rc = p.wait()
          if not rc : 
            print ( f + ' is different from ' + r)
            return False
  return True

if __name__ == '__main__' :
  # test c24_to_c12
  base_c12 = '/gpfsm/dnb32/wjiang/pl_c24Tc12'
  yaml = ruamel.yaml.YAML()
  stream =''
  c24toc12 = 'c24Toc12.yaml'
  with  open(c24toc12, 'r') as f:
     stream = f.read()
  results_c12 = yaml.load(stream)['output']['shared']['out_dir']
  remap_restarts.main -c c24toc12
  
  rc = compare(base_c12, results_c12)
  assert( rc, 'failed to compare ' + base_c12)

