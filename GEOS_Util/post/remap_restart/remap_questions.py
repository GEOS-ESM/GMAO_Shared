#!/usr/bin/env python3
#
# source install/bin/g5_modules
#
# Newer GEOS code should load a module with GEOSpyD Python3 if not run:
#   module load python/GEOSpyD/Min4.10.3_py3.9
#

import os
import subprocess
import shlex
import ruamel.yaml
import shutil
import questionary
import glob

def has_fvcore(x):
  ymdh = x['yyyymmddhh']
  time = ymdh[0:8] + '_'+ymdh[8:10]
  files = glob.glob(x['rst_dir']+'/*fvcore_*'+time+'*')
  if len(files) ==1 :
    return True
  else:
    print('\nDo not have fvcore for this date\n try "fvcore_internal_rst"')
    if os.path.exists(x['rst_dir']+'/fvcore_internal_rst'):
       return True 
    return False


def ask_common_in():
   questions = [
        {
            "type": "confirm",
            "name": "MERRA-2",
            "message": "Would you like to get restarts from MERRA-2 archive?",
            "default": False
        },
        {
            "type": "path",
            "name": "rst_dir",
            "message": "Enter the directory containing restart files to be remapped: \n  (If it is MERRA-2, files from achives will be copied to this directory automatically). \n",
        },
        {
            "type": "text",
            "name": "yyyymmddhh",
            "message": "What time would you like to remap from?(must be 10 digits: yyyymmddhh)",
            "validate": lambda text: len(text)==10 ,
        },
        {
            "type": "text",
            "name": "agrid",
            "message": "Enter input atmospheric grid: \n C12   C180  C1000  \n C24   C360  C1440 \n C48   C500  C2880 \n C90   C720  C5760 \n ",
            "default": 'C360',
            # if it is merra-2 or has_fvcore, agrid is deduced
            "when": lambda x: not x['MERRA-2'] and not has_fvcore(x),
        },
        {
            "type": "select",
            "name": "model",
            "message": "Select ocean model:",
            "choices": ["data", "MOM5", "MOM6"],
            "default": "data",
            "when": lambda x: not x['MERRA-2']
        },

        {
            "type": "select",
            "name": "ogrid",
            "message": "Input Ocean grid: \n \
             Data Ocean Grids \n \
             ------------------- \n \
             360X180  (Reynolds) \n \
             1440X720 (MERRA-2) \n \
             2880X1440  (OSTIA) \n \
             CS = same as atmospere grid (OSTIA cubed-sphere) \n",
            "choices": ['360X180','1440X720','2880X1440','CS'],
            "when": lambda x: x.get('model') == 'data' and not x['MERRA-2'],
        },
        {
            "type": "select",
            "name": "ogrid",
            "message": "Input Ocean grid: \n \
             Coupled Ocean Grids \n \
             ------------------- \n \
             72X36 \n \
             360X200 \n \
             720X410 \n \
             1440X1080 \n ",
            "choices": ['72X36','360X200','720X410','1440X1080'],
            "when": lambda x: x.get('model') == 'MOM5' or x.get('model')== 'MOM6'
        },
        {
            "type": "text",
            "name": "tag",
            "message": "Enter GCM or DAS tag for input: \n \
Sample GCM tags \n \
--------------- \n \
G30   : Ganymed-3_0  .........  Ganymed-3_0_p1 \n \
G40   : Ganymed-4_0  .........  Heracles-5_4_p3 \n \
ICA   : Icarus  ..............  Jason \n \
GITOL : 10.3  ................  10.18 \n \
INL   : Icarus-NL  ...........  Jason-NL \n \
GITNL : 10.19  ...............  10.23 \n \
\n \
Sample DAS tags \n \
--------------- \n \
5B0  : GEOSadas-5_10_0_p2  ..  GEOSadas-5_11_0 \n \
512  : GEOSadas-5_12_2  .....  GEOSadas-5_16_5\n \
517  : GEOSadas-5_17_0  .....  GEOSadas-5_24_0_p1\n \
525  : GEOSadas-5_25_1  .....  GEOSadas-5_29_4\n",
            "default": "INL",
            "when": lambda x: not x["MERRA-2"],
        },
        {
            "type": "select",
            "name": "bc_base",
            "message": "Select bcs base \n \
             discover_ops: /discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/bcs \n \
             discover_lt: /discover/nobackup/ltakacs/bcs \n \
             discover_couple: /discover/nobackup/projects/gmao/ssd/aogcm/atmosphere_bcs \n",
            "choices": ["discover_ops", "discover_lt", "discover_couple", "other"],
            "when": lambda x: not x['MERRA-2'],
        },
        {
            "type": "path",
            "name": "alt_bcs",
            "message": "Specify your own bcs absolute path (do not contain grid info) for restarts: \n ",
            "when": lambda x: x.get("bc_base")=="other",
        },
   ]
   common_in = questionary.prompt(questions)
   if not common_in.get('model') :
      common_in['model'] = 'data'
   common_in['rst_dir'] = os.path.abspath(common_in['rst_dir'])
   return common_in

def ask_common_out():
   questions = [
        {
            "type": "path",
            "name": "out_dir",
            "message": "Enter the directory for new restarts:\n"
        },
        {
            "type": "text",
            "name": "expid",
            "message": "Enter new restarts expid:",
            "default": "",
        },
        {
            "type": "text",
            "name": "agrid",
            "message": "Enter new atmospheric grid: \n C12   C180  C1000  \n C24   C360  C1440 \n C48   C500  C2880 \n C90   C720  C5760 \n ",
            "default": 'C360',
        },
        {
            "type": "text",
            "name": "tag",
            "message": "Enter GCM or DAS tag for new restarts: \n \
Sample GCM tags \n \
--------------- \n \
F14   : Fortuna-1_4  .........  Fortuna-1_4_p1 \n \
F20   : Fortuna-2_0  .........  Fortuna-2_0 \n \
F21   : Fortuna-2_1  .........  Fortuna-2_5_pp2 \n \
G10   : Ganymed-1_0  .........  Ganymed-1_0_BETA4 \n \
G10p  : Ganymed-1_0_p1  ......  Ganymed-1_0_p6 \n \
G20   : Ganymed-2_0  .........  Ganymed-2_1_p6 \n \
G30   : Ganymed-3_0  .........  Ganymed-3_0_p1 \n \
G40   : Ganymed-4_0  .........  Heracles-5_4_p3 \n \
ICA   : Icarus  ..............  Jason \n \
GITOL : 10.3  ................  10.18 \n \
INL   : Icarus-NL  ...........  Jason-NL \n \
GITNL : 10.19  ...............  10.23 \n \
\n \
Sample DAS tags \n \
--------------- \n \
214  : GEOSdas-2_1_4  .......  GEOSdas-2_1_4-m4 \n \
540  : GEOSadas-5_4_0  ......  GEOSadas-5_5_3 \n \
561  : GEOSadas-5_6_1  ......  GEOSadas-5_7_3_p2 \n \
580  : GEOSadas-5_8_0  ......  GEOSadas-5_9_1 \n \
591p : GEOSadas-5_9_1_p1  ...  GEOSadas-5_9_1_p9 \n \
5A0  : GEOSadas-5_10_0  .....  GEOSadas-5_10_0_p1 \n \
5B0  : GEOSadas-5_10_0_p2  ..  GEOSadas-5_11_0 \n \
512  : GEOSadas-5_12_2  .....  GEOSadas-5_16_5\n \
517  : GEOSadas-5_17_0  .....  GEOSadas-5_24_0_p1\n \
525  : GEOSadas-5_25_1  .....  GEOSadas-5_29_4\n",
            "default": "INL",
        },
        {
            "type": "select",
            "name": "model",
            "message": "Select ocean model for new restarts:",
            "choices": ["data", "MOM5", "MOM6"],
            "default": "data",
        },
        {
            "type": "select",
            "name": "ogrid",
            "message": "Select new Ocean grid: \n \
             Data Ocean Grids \n \
             ------------------- \n \
             360X180  (Reynolds) \n \
             1440X720 (MERRA-2) \n \
             2880X1440  (OSTIA) \n \
             CS = same as atmospere grid (OSTIA cubed-sphere) \n",
            "choices": ['360X180','1440X720','2880X1440','CS'],
            "when": lambda x: x['model'] == 'data',
        },
        {
            "type": "select",
            "name": "ogrid",
            "message": "Select new Ocean grid: \n \
             Coupled Ocean Grids \n \
             ------------------- \n \
             72X36 \n \
             360X200 \n \
             720X410 \n \
             1440X1080 \n ",
            "choices": ['72X36','360X200','720X410','1440X1080'],
            "when": lambda x: x['model'] != 'data',
        },
        {
            "type": "select",
            "name": "bc_base",
            "message": "Select bcs base for new restarts: \n \
             discover_ops: /discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/bcs \n \
             discover_lt: /discover/nobackup/ltakacs/bcs \n \
             discover_couple: /discover/nobackup/projects/gmao/ssd/aogcm/atmosphere_bcs \n",
            "choices": ["discover_ops", "discover_lt", "discover_couple", "other"],
        },
        {
            "type": "path",
            "name": "alt_bcs",
            "message": "Specify your own bcs path (do not contain grid info) for new restarts: \n ",
            "when": lambda x: x.get("bc_base")=="other",
        },

   ]
   common_out = questionary.prompt(questions)
   if common_out.get('ogrid') == 'CS':
      common_out['ogrid'] = common_out['agrid']
   common_out['out_dir'] = os.path.abspath(common_out['out_dir'])

   return common_out
def ask_upper_out():
   questions = [
        {
            "type": "text",
            "name": "nlevel",
            "message": "Enter new atmospheric levels: (71 72 91 127 132 137 144 181)",
            "default": "72",
        },
   ]
   return questionary.prompt(questions)
def ask_surface_in(common_in):
   wemin_default = '26'
   tag = common_in.get('tag')
   if tag in ['INL','GITNL', '525'] : wemin_default = '13'
   questions = [
        {
            "type": "text",
            "name": "wemin",
            "message": "What is value of Wemin?",
            "default": wemin_default
        },
        {
            "type": "text",
            "name": "zoom",
            "message": "What is value of zoom [1-8]?",
            "default": '8'
        },
   ]
   return questionary.prompt(questions)

def ask_surface_out(common_out):
   wemout_default = '26'
   tag = common_out.get('tag')
   if tag in ['INL','GITNL', '525'] : wemout_default = '13'
   questions = [
        {
            "type": "text",
            "name": "wemout",
            "message": "What is value of Wemout?",
            "default": wemout_default
        },
   ]
   return questionary.prompt(questions)

def ask_analysis_in():
   questions = [
        {
            "type": "confirm",
            "name": "bkg",
            "message": "Regrid bkg files?",
            "default": False,
        },
        {
            "type": "confirm",
            "name": "lcv",
            "message": "Write lcv?",
            "default": False,
        },
   ]
   return questionary.prompt(questions)

def ask_slurm_options():
   cmd = 'id -gn'
   p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
   (accounts, err) = p.communicate()
   p_status = p.wait()
   accounts = accounts.decode().split()

   questions = [
        {
            "type": "text",
            "name": "qos",
            "message": "qos?",
            "default": "debug",
        },
        {
            "type": "text",
            "name": "account",
            "message": "account?",
            "default": accounts[0],
        },
        {
            "type": "select",
            "name": "partition",
            "message": "partition?",
            "choices": ['hasw', 'sky'],
        },
   ]
   return questionary.prompt(questions)

def get_config_from_questionary():
    print("\nCommon parameters INPUT: \n")
    common_in = ask_common_in()
    print("\nSurface parameters INPUT: \n")
    surface_in = ask_surface_in(common_in)
    print("\nAnalysis parameters INPUT: \n")
    analysis_in = ask_analysis_in()
    print("\nCommon parameters OUTPUT: \n")
    common_out = ask_common_out()
    print("\nUPPER parameters OUTPUT: \n")
    upper_out = ask_upper_out()
    print("\nSurface parameters OUTPUT: \n")
    surface_out = ask_surface_out(common_out)
    print("\nSlurm parameters OUTPUT: \n")
    slurm = ask_slurm_options()

    config={}
    inputs={}
    outputs={}
    parameters_in={}
    parameters_in['COMMON']  = common_in
    parameters_in['SURFACE'] = surface_in
    inputs['parameters'] = parameters_in

    config['input']= inputs

    parameters_out={}
    parameters_out['COMMON']  = common_out
    parameters_out['UPPERAIR'] = upper_out
    parameters_out['SURFACE'] = surface_out
    parameters_out['ANALYSIS'] = analysis_in
    outputs['parameters'] = parameters_out

    config['output']= outputs

    config['slurm_options'] = slurm

    return config

if __name__ == "__main__":
  config = get_config_from_questionary()
  yaml = ruamel.yaml.YAML()
  with open("raw_answers.yaml", "w") as f:
    yaml.dump(config, f)

