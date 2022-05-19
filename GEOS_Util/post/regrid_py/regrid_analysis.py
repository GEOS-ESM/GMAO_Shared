#!/usr/bin/env python3
#
# source install/bin/g5_modules
#
# Newer GEOS code should load a module with GEOSpyD Python3 if not run:
#   module load python/GEOSpyD/Min4.10.3_py3.9
#

import os
from datetime import datetime, timedelta
import subprocess
import shlex
import shutil
import glob
import fileinput
import ruamel.yaml

class analysis(object):
  def __init__(self, params_file):
     yaml = ruamel.yaml.YAML()
     stream =''
     with  open(params_file, 'r') as f:
        stream = f.read()
     self.config = yaml.load(stream)

  def regrid(self):
     config = self.config
     bkg = config['output']['analysis']['bkg']
     if ( not bkg ): return

     agrid_in  = config['input']['shared']['agrid']
     agrid_out = config['output']['shared']['agrid']

     if (self.get_grid_kind(agrid_in.upper()) == self.get_grid_kind(agrid_out.upper())):
       print(" No need to regrid anaylysis file")
       return

     bindir  = os.getcwd()
     in_bcsdir  = config['input']['shared']['bcs_dir']
     out_bcsdir = config['output']['shared']['bcs_dir']
     out_dir    = config['output']['shared']['out_dir']
     if not os.path.exists(out_dir) : os.makedirs(out_dir)
     print( "cd " + out_dir)
     os.chdir(out_dir)

     tmpdir = out_dir+'/ana_data/'
     if os.path.exists(tmpdir) : subprocess.call(['rm', '-rf',tmpdir])
     print ("mkdir " + tmpdir)
     os.makedirs(tmpdir)

     print( "cd " + tmpdir)
     os.chdir(tmpdir)

     rst_dir_orig = config['input']['analysis']['rst_dir']
     if (not rst_dir_orig): rst_dir_orig = config['input']['shared']['rst_dir']
     anafiles = glob.glob(rst_dir_orig +'/analysis/*')
     yyyymmddhh_ = str(config['input']['shared']['yyyymmddhh'])
     yyyy_ = yyyymmddhh_[0:4]
     mm_   = yyyymmddhh_[4:6]
     dd_   = yyyymmddhh_[6:8]
     hh_   = yyyymmddhh_[8:10]
     rst_time = datetime(year=int(yyyy_), month=int(mm_), day=int(dd_), hour = int(hh_))
     expid_in  = config['input']['shared']['expid']
     expid_out = config['output']['shared']['expid']
     if (expid_out) :
        expid_out = expid_out + '.'
     else:
        expid_out = ''

     if len(anafiles) == 0 :
       anafiles=[]
       for h in [3,4,5,6,7,8,9]:
          delt = timedelta(hours = h-3)
          new_time = rst_time + delt
          yyyy = "Y"+str(new_time.year)
          mm   = 'M%02d'%new_time.month
          ymd  = '%04d%02d%02d'%(new_time.year,new_time.month, new_time.day)
          hh   = '%02d'%h
          newhh= '%02d'%new_time.hour
          rst_dir = rst_dir_orig.replace('Y'+yyyy_,yyyy).replace('M'+mm_,mm)
          # bkg files
          for ftype in ['sfc', 'eta']:
             fname = expid_in+'.bkg'+hh+'_'+ftype+'_rst.'+ymd+'_'+newhh+'z.nc4'
             f = rst_dir+'/'+fname
             if(os.path.isfile(f)):
               anafiles.append(f)
             else:
               print('Warning: Cannot find '+f)

          # cbkg file
          fname = expid_in + '.cbkg' + hh + '_eta_rst.' + ymd + '_' + newhh + 'z.nc4'
          f = rst_dir+'/'+fname
          if(os.path.isfile(f)):
            anafiles.append(f)
          else:
            print('Warning: Cannot find '+f)
          # gaas_bkg_sfc files
          if (h==6 or h==9):
             fname = expid_in+'.gaas_bkg_sfc_rst.'+ymd+'_'+newhh+'z.nc4'
             f = rst_dir+'/'+fname
             if (os.path.isfile(f)):
               anafiles.append(f)
             else:
               print('Warning: Cannot find '+f)
       # trak.GDA.rst file
       delt = timedelta(hours = 3)
       new_time = rst_time - delt
       yyyy = "Y"+str(new_time.year)
       mm   = 'M%02d'%new_time.month
       ymdh = '%04d%02d%02d%02d'%(new_time.year, new_time.month, new_time.day, new_time.hour)
       rst_dir = rst_dir_orig.replace('Y'+yyyy_,yyyy).replace('M'+mm_,mm)
       fname = expid_in+'.trak.GDA.rst.'+ymdh+'z.txt'
       f = rst_dir+'/'+fname
       if (os.path.isfile(f)): anafiles.append(f)

     aqua   = config['output']['analysis']['aqua']
     local_fs=[]
     for f in anafiles:
       fname    = os.path.basename(f)
       out_name = fname.replace(expid_in + '.', expid_out)
       f_tmp = tmpdir+'/'+out_name
       local_fs.append(f_tmp) 
       shutil.copy(f,f_tmp)
       if out_name.find('satbias') != -1 :
         if (aqua):
            f_orig = f_tmp+'.orig'
            shutil.move(f_tmp, f_orig)
            f_ = open(f_tmp, 'w')
            for line in fileinput.input(f_orig):
               f_.write(line.replace('airs281SUBSET_aqua', 'airs281_aqua      '))
            f_.close() 

     
     nlevel = config['output']['air']['nlevel']
     flags = "-g5 -res " + self.get_grid_kind(agrid_out.upper()) + " -nlevs " + str(nlevel)
     bkg_files = glob.glob(tmpdir+'/*.bkg??_eta_rst*')
     for f in bkg_files:
        f_orig = f + ".orig"
        shutil.move(f,f_orig)
        cmd = bindir + '/dyn2dyn.x ' + flags + ' -o ' + f + ' ' + f_orig
        print(cmd)
        subprocess.call(shlex.split(cmd))
     for f in local_fs:
       fname = os.path.basename(f)
       shutil.move(f, out_dir+'/'+fname)
    # write lcv 
     lcv = config['output']['analysis']['lcv']
     if lcv :
       ymd_ = yyyymmddhh_[0:8]
       hh_  = yyyymmddhh_[8:10]
       hms_ = hh_+'0000'
       rstlcvOut = out_dir+'/'+expid_out+'rst.lcv.'+ymd_+'_'+hh_+'z.bin'
       cmd = bindir+'/mkdrstdate.x ' + ymd_ + ' ' + hms_ +' ' + rstlcvOut
       print(cmd)
       subprocess.call(shlex.split(cmd))
     os.chdir(bindir)

  def get_grid_kind(this, grid):
     hgrd = {}
     hgrd['C12']   = 'a'
     hgrd['C24']   = 'a'
     hgrd['C48']   = 'b'
     hgrd['C90']   = 'c'
     hgrd['C180']  = 'd'
     hgrd['C360']  = 'd'
     hgrd['C500']  = 'd'
     hgrd['C720']  = 'e'
     hgrd['C1000'] = 'e'
     hgrd['C1440'] = 'e'
     hgrd['C2000'] = 'e'
     hgrd['C2880'] = 'e'
     hgrd['C5760'] = 'e'
     return hgrd[grid]

if __name__ == '__main__' :
   ana = analysis('regrid_params.yaml')
   ana.regrid()
