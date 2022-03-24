#!/usr/bin/env python3
#
import os
import subprocess
import shutil
import glob
import yaml
import shlex

class lakelandice(object):
  def __init__(self, params_file):
     stream = open(params_file, 'r')
     self.config = yaml.full_load(stream)
     stream.close()

  def regrid(self):
     print("\nRegridding land, landice, saltwater.....\n")
     config = self.config
     bindir  = os.getcwd()
     in_bcsdir  = config['input']['shared']['bcs_dir']
     out_bcsdir = config['output']['shared']['bcs_dir']
     out_dir    = config['output']['shared']['out_dir']
     if not os.path.exists(out_dir) : os.makedirs(out_dir)
     print( "cd " + out_dir)
     os.chdir(out_dir)

     InData_dir = out_dir+'/InData/'
     if os.path.exists(InData_dir) : subprocess.call(['rm', '-rf',InData_dir])
     print ("mkdir " + InData_dir)
     os.makedirs(InData_dir)

     OutData_dir = out_dir+'/OutData/'
     if os.path.exists(OutData_dir) : subprocess.call(['rm', '-rf',OutData_dir])
     print ("mkdir " + OutData_dir)
     os.makedirs(OutData_dir)

     rst_dir = config['input']['surface']['rst_dir']
     if (not rst_dir): rst_dir = config['input']['shared']['rst_dir']
     restarts_in = glob.glob(rst_dir +'surface/*') 
     yyyymmddhh_ = str(config['input']['shared']['yyyymmddhh'])
     suffix = yyyymmddhh_[0:8]+'_'+yyyymmddhh_[8:10]+'z.nc4'
     if (not restarts_in) :
         expid = config['input']['shared']['expid']
         suffix = yyyymmddhh_[0:8]+'_'+yyyymmddhh_[8:10]+'z.bin'
         # from merra-2
         restarts_in = [rst_dir + expid + '.lake_internal_rst.'     + suffix,
                        rst_dir + expid + '.landice_internal_rst.'  + suffix,
                        rst_dir + expid + '.saltwater_internal_rst.'+ suffix]

     def find_rst(rsts,rst_name):
        for rst in rsts:
           f = os.path.basename(rst)
           if (f.find(rst_name) != -1):
              return f
        return None  
     
     saltwater = find_rst(restarts_in,'saltwater')
     openwater = find_rst(restarts_in,'openwater')
     seaice  = find_rst(restarts_in,'seaice') 
     landice = find_rst(restarts_in,'landice')
     lake  = find_rst(restarts_in,'lake')
     route = find_rst(restarts_in,'route')

     for rst in restarts_in:
        f = os.path.basename(rst)
        dest = InData_dir+'/'+f
        if os.path.exists(dest) : shutil.remove(dest)
        print('\nCopy ' + rst + ' to ' +dest)
        shutil.copy(rst,dest)

     in_tile_file  = glob.glob(in_bcsdir+ '/*-Pfafstetter.til')[0]
     out_tile_file = glob.glob(out_bcsdir+ '/*-Pfafstetter.til')[0]

     in_til = InData_dir+'/' + in_tile_file.split('/')[-1]
     out_til = OutData_dir+'/'+ out_tile_file.split('/')[-1]

     if os.path.exists(in_til)  : shutil.remove(in_til)
     if os.path.exists(out_til) : shutil.remove(out_til)
     print('\n Copy ' + in_tile_file + ' to ' + in_til)
     shutil.copy(in_tile_file, in_til)
     print('\n Copy ' + out_tile_file + ' to ' + out_til)
     shutil.copy(out_tile_file, out_til)

     exe = bindir + '/mk_LakeLandiceSaltRestarts '
     zoom = config['output']['surface']['zoom']

     if (saltwater):
       cmd = exe + out_til + ' ' + in_til + ' InData/'+ saltwater + ' 0 ' + str(zoom)
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))
  
       # split Saltwater
       if  config['output']['surface']['split_saltwater']:
         print("\nSplitting Saltwater...\n")
         cmd = bindir+'/SaltIntSplitter ' + out_til + ' ' + 'OutData/' + saltwater
         print('\n'+cmd)
         subprocess.call(shlex.split(cmd))
         openwater = None
         seaice  = None

     if (openwater):
       cmd = exe + out_til + ' ' + in_til + ' InData/' + openwater + ' 0 ' + str(zoom)
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))

     if (seaice):
       cmd = exe + out_til + ' ' + in_til + ' InData/' + seaice + ' 0 ' + str(zoom)
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))

     if (lake):
       cmd = exe + out_til + ' ' + in_til + ' InData/' + lake + ' 19 ' + str(zoom)
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))

     if (landice):
       cmd = exe + out_til + ' ' + in_til + ' InData/' + landice + ' 20 ' + str(zoom)
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))

     if (route):
       route = bindir + '/mk_RouteRestarts '
       cmd = route + out_til + ' ' + yyyymmddhh_[0:6]
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))

     expid = config['output']['shared']['expid']
     if (expid) :
        expid = expid + '.'
     else:
        expid = ''
     suffix = '_rst.' + suffix
     for out_rst in glob.glob("OutData/*_rst*"):
       filename = expid + os.path.basename(out_rst).split('_rst')[0].split('.')[-1]+suffix
       print('\n Move ' + out_rst + ' to ' + out_dir+"/"+filename)
       shutil.move(out_rst, out_dir+"/"+filename)
     os.chdir(bindir)

if __name__ == '__main__' :
   lli = lakelandice('regrid_params.yaml')
   lli.regrid()

