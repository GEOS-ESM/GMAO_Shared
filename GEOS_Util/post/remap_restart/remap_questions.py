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

def fvcore_name(x):
  ymdh = x['input:shared:yyyymmddhh']
  time = ymdh[0:8] + '_'+ymdh[8:10]
  files = glob.glob(x['input:shared:rst_dir']+'/*fvcore_*'+time+'*')
  if len(files) ==1 :
    fname = files[0]
    print('\nFound ' + fname) 
    return fname
  else:
    fname = x['input:shared:rst_dir']+'/fvcore_internal_rst'
    if os.path.exists(fname):
       print('\nFound ' + fname) 
       return fname
    return False

def tmp_merra2_dir(x):
   cmd = 'whoami'
   p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
   (user, err) = p.communicate()
   p_status = p.wait()
   print(user)
   user = user.decode().split()
   tmp_merra2 = '/discover/nobackup/'+user[0]+'/merra2_tmp'+x['input:shared:yyyymmddhh']+'/'
   return tmp_merra2

def we_default(tag):
   default_ = '26'
   if tag in ['INL','GITNL', '525'] : default_ = '13'
   return default_

def zoom_default(x):
   zoom_ = '8'   
   fvcore = fvcore_name(x)
   if fvcore :
      fvrst = os.path.dirname(os.path.realpath(__file__)) + '/fvrst.x -h '
      cmd = fvrst + fvcore
      print(cmd +'\n')
      p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
      (output, err) = p.communicate()
      p_status = p.wait()
      ss = output.decode().split()
      x['input:shared:agrid'] = "C"+ss[0] # save for air parameter
      lat = int(ss[0])
      lon = int(ss[1])
      if (lon != lat*6) :
         sys.exit('This is not a cubed-sphere grid fvcore restart. Please contact SI team')
      ymdh  = x.get('input:shared:yyyymmddhh')
      ymdh_ = str(ss[3]) + str(ss[4])[0:2]
      if (ymdh_ != ymdh) :
         print("Warning: The date in fvcore is different from the date you input\n")
      zoom = lat /90.0
      zoom_ = str(int(zoom))
      if zoom < 1 : zoom_ = '1'
      if zoom > 8 : zoom_ = '8'
   if x['input:shared:MERRA-2'] :
      zoom_ = '2'
   return zoom_

def get_account():
   cmd = 'id -gn'
   p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
   (accounts, err) = p.communicate()
   p_status = p.wait()
   accounts = accounts.decode().split()
   return accounts[0]

def ask_questions():

   questions = [
        {
            "type": "confirm",
            "name": "input:shared:MERRA-2",
            "message": "Would you like to remap archived MERRA-2 restarts?",
            "default": False,
        },
        {
            "type": "path",
            "name": "input:shared:rst_dir",
            "message": "Enter the directory containing restart files to be remapped:",
            "when": lambda x: not x['input:shared:MERRA-2'],
        },
        {
            "type": "text",
            "name": "input:shared:yyyymmddhh",
            "message": "From what restart date/time would you like to remap? (must be 10 digits: yyyymmddhh)",
            "validate": lambda text: len(text)==10 ,
        },
        {
            "type": "path",
            "name": "input:shared:rst_dir",
            "message": "Enter a directory to which the archived MERRA-2 archive files can be copied: ",
            "default": lambda x: tmp_merra2_dir(x),
            "when": lambda x: x['input:shared:MERRA-2'],
        },

        {
            "type": "path",
            "name": "output:shared:out_dir",
            "message": "Enter the directory for new restarts:\n"
        },

        {
            "type": "text",
            "name": "input:shared:agrid",
            "message": "Enter input atmospheric grid: \n C12   C180  C1000  \n C24   C360  C1440 \n C48   C500  C2880 \n C90   C720  C5760 \n ",
            "default": 'C360',
            # if it is merra-2 or has_fvcore, agrid is deduced
            "when": lambda x: not x['input:shared:MERRA-2'] and not fvcore_name(x),
        },

        {
            "type": "text",
            "name": "output:shared:agrid",
            "message": "Enter new atmospheric grid: \n C12   C180  C1000  \n C24   C360  C1440 \n C48   C500  C2880 \n C90   C720  C5760 \n ",
            "default": 'C360',
        },

        {
            "type": "text",
            "name": "output:air:nlevel",
            "message": "Enter new atmospheric levels: (71 72 91 127 132 137 144 181)",
            "default": "72",
        },

        {
            "type": "select",
            "name": "input:shared:model",
            "message": "Select input ocean model:",
            "choices": ["data", "MOM5", "MOM6"],
            "default": "data",
            "when": lambda x: not x['input:shared:MERRA-2']
        },

        {
            "type": "select",
            "name": "input:shared:ogrid",
            "message": "Input Ocean grid: \n \
             Data Ocean Grids \n \
             ------------------- \n \
             360X180  (Reynolds) \n \
             1440X720 (MERRA-2) \n \
             2880X1440  (OSTIA) \n \
             CS = same as atmospere grid (OSTIA cubed-sphere) \n",
            "choices": ['360X180','1440X720','2880X1440','CS'],
            "when": lambda x: x.get('input:shared:model') == 'data' and not x['input:shared:MERRA-2'],
        },

        {
            "type": "select",
            "name": "output:shared:model",
            "message": "Select ocean model for new restarts:",
            "choices": ["data", "MOM5", "MOM6"],
            "default": "data",
        },
        {
            "type": "select",
            "name": "output:shared:ogrid",
            "message": "Select new ocean grid:", 
            "choices": ['360X180','1440X720','2880X1440','CS'],
            "when": lambda x: x['output:shared:model'] == 'data',
        },

        {
            "type": "select",
            "name": "input:shared:ogrid",
            "message": "Input ocean grid: \n \
             Coupled Ocean Grids \n \
             ------------------- \n \
             72X36 \n \
             360X200 \n \
             720X410 \n \
             1440X1080 \n ",
            "choices": ['72X36','360X200','720X410','1440X1080'],
            "when": lambda x: x.get('input:shared:model') == 'MOM5' or x.get('input:shared:model')== 'MOM6'
        },
        {
            "type": "select",
            "name": "output:shared:ogrid",
            "message": "Select new ocean grid: \n \
             Coupled Ocean Grids \n \
             ------------------- \n \
             72X36 \n \
             360X200 \n \
             720X410 \n \
             1440X1080 \n ",
            "choices": ['72X36','360X200','720X410','1440X1080'],
            "when": lambda x: x['output:shared:model'] != 'data',
        },

        {
            "type": "text",
            "name": "input:shared:tag",
            "message": "Enter GCM or DAS tag for input: \n \
Sample GCM tags \n \
--------------- \n \
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
            "when": lambda x: not x["input:shared:MERRA-2"],
        },
        {
            "type": "text",
            "name": "output:shared:tag",
            "message": "Enter GCM or DAS tag for new restarts:",
            "default": "INL",
        },

        {
            "type": "select",
            "name": "input:shared:bc_base",
            "message": "Select bcs base \n \
             discover_ops: /discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/bcs \n \
             discover_lt: /discover/nobackup/ltakacs/bcs \n \
             discover_couple: /discover/nobackup/projects/gmao/ssd/aogcm/atmosphere_bcs \n",
            "choices": ["discover_ops", "discover_lt", "discover_couple", "other"],
            "when": lambda x: not x['input:shared:MERRA-2'],
        },
        {
            "type": "path",
            "name": "input:shared:alt_bcs",
            "message": "Specify your own bcs absolute path (do not contain grid info) for restarts: \n ",
            "when": lambda x: x.get("input:shared:bc_base")=="other",
        },

        {
            "type": "select",
            "name": "output:shared:bc_base",
            "message": "Select bcs base for new restarts:",
            "choices": ["discover_ops", "discover_lt", "discover_couple", "other"],
        },
        {
            "type": "path",
            "name": "output:shared:alt_bcs",
            "message": "Specify your own bcs path (do not contain grid info) for new restarts: \n ",
            "when": lambda x: x.get("output:shared:bc_base")=="other",
        },

        {
            "type": "confirm",
            "name": "output:air:remap",
            "message": "Would you like to remap upper air?",
            "default": True,
        },
        {
            "type": "confirm",
            "name": "output:surface:remap",
            "message": "Would you like to remap surface?",
            "default": True,
        },
        {
            "type": "confirm",
            "name": "output:analysis:bkg",
            "message": "Regrid bkg files?",
            "default": False,
        },
        {
            "type": "confirm",
            "name": "output:analysis:lcv",
            "message": "Write lcv?",
            "default": False,
        },
        {
            "type": "text",
            "name": "input:surface:wemin",
            "message": "What is value of Wemin?",
            "default": lambda x: we_default(x.get('input:shared:tag'))
        },
        {
            "type": "text",
            "name": "output:surface:wemout",
            "message": "What is value of Wemout?",
            "default": lambda x: we_default(x.get('output:shared:tag'))
        },
        {
            "type": "text",
            "name": "input:surface:zoom",
            "message": "What is value of zoom [1-8]?",
            "default": lambda x: zoom_default(x)
        },
        {
            "type": "text",
            "name": "output:shared:expid",
            "message": "Enter new restarts expid:",
            "default": "",
        },

        {
            "type": "text",
            "name": "slurm:qos",
            "message": "qos?",
            "default": "debug",
        },

        {
            "type": "text",
            "name": "slurm:account",
            "message": "account?",
            "default": get_account(),
        },
        {
            "type": "select",
            "name": "slurm:constraint",
            "message": "constraint?",
            "choices": ['hasw', 'sky', 'cas'],
        },
   ]
   answers = questionary.prompt(questions)
   if not answers.get('input:shared:model') :
      answers['input:shared:model'] = 'data'
   answers['input:shared:rst_dir'] = os.path.abspath(answers['input:shared:rst_dir'])
   if answers.get('output:shared:ogrid') == 'CS':
      answers['output:shared:ogrid'] = answers['output:shared:agrid']
   answers['output:shared:out_dir']  = os.path.abspath(answers['output:shared:out_dir'])

   return answers

def get_config_from_questionary():
   answers = ask_questions()
   config  = {}
   config['input'] = {}
   config['input']['shared'] = {}
   config['input']['surface'] = {}
   config['output'] = {}
   config['output']['shared'] = {}
   config['output']['air'] = {}
   config['output']['surface'] = {}
   config['output']['analysis'] = {}
   config['slurm'] = {} 
   for key, value in answers.items():
     keys = key.split(":")
     if len(keys) == 2:
       config[keys[0]][keys[1]] = value 
     if len(keys) == 3:
       config[keys[0]][keys[1]][keys[2]] = value 

   return config

if __name__ == "__main__":
   config = get_config_from_questionary()
   yaml = ruamel.yaml.YAML()
   with open("raw_answers.yaml", "w") as f:
     yaml.dump(config, f)

