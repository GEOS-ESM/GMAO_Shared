#!/usr/bin/env python3
#
import os
import subprocess
import shutil
import glob
import ruamel.yaml
import shlex
from remap_base import remap_base

class lake_landice_saltwater(remap_base):
  def __init__(self, **configs):
     super().__init__(**configs)
     self.copy_merra2()

  def remap(self):
     if not self.config['output']['surface']['remap_water']:
        return

     restarts_in = self.find_rst()
     if len(restarts_in) == 0:
       return

     print("\nRemapping land, landice, saltwater.....\n")
     config = self.config
     cwdir  = os.getcwd()
     bindir  = os.path.dirname(os.path.realpath(__file__)) 
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

     types = 'z.bin'
     type_str = subprocess.check_output(['file','-b', restarts_in[0]])
     type_str = str(type_str)
     if 'Hierarchical' in type_str:
        types = 'z.nc4'
     yyyymmddhh_ = str(config['input']['shared']['yyyymmddhh'])
     suffix = yyyymmddhh_[0:8]+'_'+yyyymmddhh_[8:10]+ types

     saltwater = ''
     seaice    = ''
     landice   = ''
     lake      = ''
     route     = ''
     openwater = ''
     for rst in restarts_in:
        f = os.path.basename(rst)
        dest = InData_dir+'/'+f
        if os.path.exists(dest) : shutil.remove(dest)
        print('\nCopy ' + rst + ' to ' +dest)
        shutil.copy(rst,dest)
        if 'saltwater' in f : saltwater = f
        if 'seaice'    in f : seaice    = f
        if 'landice'   in f : landice   = f
        if 'lake'      in f : lake      = f
        if 'roue'      in f : route     = f
        if 'openwater' in f : openwater = f

     in_tile_file  = glob.glob(in_bcsdir+ '/*-Pfafstetter.til')[0]
     out_tile_file = glob.glob(out_bcsdir+ '/*-Pfafstetter.til')[0]

     in_til = InData_dir+'/' + os.path.basename(in_tile_file)
     out_til = OutData_dir+'/'+ os.path.basename(out_tile_file)

     if os.path.exists(in_til)  : shutil.remove(in_til)
     if os.path.exists(out_til) : shutil.remove(out_til)
     print('\n Copy ' + in_tile_file + ' to ' + in_til)
     shutil.copy(in_tile_file, in_til)
     print('\n Copy ' + out_tile_file + ' to ' + out_til)
     shutil.copy(out_tile_file, out_til)

     exe = bindir + '/mk_LakeLandiceSaltRestarts.x '
     zoom = config['input']['surface']['zoom']

     if (saltwater):
       cmd = exe + out_til + ' ' + in_til + ' InData/'+ saltwater + ' 0 ' + str(zoom)
       print('\n'+cmd)
       subprocess.call(shlex.split(cmd))
  
       # split Saltwater
       if  config['output']['surface']['split_saltwater']:
         print("\nSplitting Saltwater...\n")
         cmd = bindir+'/SaltIntSplitter.x ' + out_til + ' ' + 'OutData/' + saltwater
         print('\n'+cmd)
         subprocess.call(shlex.split(cmd))
         openwater = ''
         seaice  = ''

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
       route = bindir + '/mk_RouteRestarts.x '
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
     print('cd ' + cwdir)
     os.chdir(cwdir)

     self.remove_merra2()

  def find_rst(self):
     surf_restarts =[
                 "route_internal_rst"       ,
                 "lake_internal_rst"        ,
                 "landice_internal_rst"     ,
                 "openwater_internal_rst"   ,
                 "saltwater_internal_rst"   ,
                 "seaicethermo_internal_rst"]

     rst_dir = self.config['input']['shared']['rst_dir']
     yyyymmddhh_ = str(self.config['input']['shared']['yyyymmddhh'])
     time = yyyymmddhh_[0:8]+'_'+yyyymmddhh_[8:10]
     restarts_in=[]
     for f in surf_restarts :
        files = glob.glob(rst_dir+ '/*'+f+'*'+time+'*')
        if len(files) >0:
          restarts_in.append(files[0])
     if (len(restarts_in) == 0) :
        print("\n try restart file names without time stamp\n")
        for f in surf_restarts :
           fname = rst_dir+ '/'+f
           if os.path.exists(fname):
             restarts_in.append(fname)

     return restarts_in

  def copy_merra2(self):
    if not self.config['input']['shared']['MERRA-2']:
      return

    expid = self.config['input']['shared']['expid']
    yyyymmddhh_ = str(self.config['input']['shared']['yyyymmddhh'])
    yyyy_ = yyyymmddhh_[0:4]
    mm_   = yyyymmddhh_[4:6]
    dd_   = yyyymmddhh_[6:8]
    hh_   = yyyymmddhh_[8:10]

    suffix = yyyymmddhh_[0:8]+'_'+ hh_ + 'z.bin'
    merra_2_rst_dir = '/archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/'+expid +'/rs/Y'+yyyy_ +'/M'+mm_+'/'
    rst_dir = self.config['input']['shared']['rst_dir'] + '/'
    os.makedirs(rst_dir, exist_ok = True)
    print(' Copy MERRA-2 surface restarts \n from \n    ' + merra_2_rst_dir + '\n to\n    '+ rst_dir +'\n')

    surfin = [ merra_2_rst_dir +  expid+'.lake_internal_rst.'     + suffix,
               merra_2_rst_dir +  expid+'.landice_internal_rst.'  + suffix,
               merra_2_rst_dir +  expid+'.saltwater_internal_rst.'+ suffix]

    for f in surfin :
       fname = os.path.basename(f)
       dest = rst_dir + '/'+fname
       print("Copy file "+f +" to " + rst_dir)
       shutil.copy(f, dest)

if __name__ == '__main__' :
   lls = lake_landice_saltwater(params_file='remap_params.yaml')
   lls.remap()
