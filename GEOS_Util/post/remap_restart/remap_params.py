#!/usr/bin/env python3
#
import os,sys
import ruamel.yaml
import shutil
import glob
import time
import shlex
import subprocess
import questionary
from datetime import datetime
from datetime import timedelta

class remap_params(object):
  def __init__(self, config_from_question):
     self.common_in     = config_from_question['input']['parameters']['COMMON']
     self.common_out    = config_from_question['output']['parameters']['COMMON']
     self.upper_out     = config_from_question['output']['parameters']['UPPERAIR']
     self.slurm_options = config_from_question['slurm_options']
     self.surf_in  = config_from_question['input']['parameters']['SURFACE']
     self.surf_out = config_from_question['output']['parameters']['SURFACE']
     #self.ana_in   = config_from_question['input']['parameters']['ANALYSIS']
     self.ana_out  = config_from_question['output']['parameters']['ANALYSIS']

     self.init_time()
     self.init_tags()
     self.init_merra2()


     # load input yaml
     yaml = ruamel.yaml.YAML() 
     stream = ''
     remap_tpl = os.path.dirname(os.path.realpath(__file__)) + '/remap_params.tpl'
     with  open(remap_tpl, 'r') as f:
       stream = f.read()
     config_tpl = yaml.load(stream)

     # params for shared
     config_tpl['input']['shared']['agrid']    = self.common_in.get('agrid')
     config_tpl['input']['shared']['ogrid']    = self.common_in.get('ogrid')
     config_tpl['input']['shared']['rst_dir']  = self.common_in['rst_dir']+'/'
     config_tpl['input']['shared']['expid']    = self.common_in.get('expid')
     config_tpl['input']['shared']['yyyymmddhh'] = self.common_in['yyyymmddhh']

     config_tpl['output']['air']['nlevel']     = self.upper_out.get('nlevel')
     config_tpl['output']['shared']['agrid']   = self.common_out['agrid']
     config_tpl['output']['shared']['ogrid']   = self.common_out['ogrid']
     config_tpl['output']['shared']['out_dir'] = self.common_out['out_dir'] + '/'
     config_tpl['output']['shared']['expid']   = self.common_out['expid']

     # params for upper air
     config_tpl = self.params_for_air(config_tpl)
     config_tpl = self.params_for_surface(config_tpl)
     config_tpl = self.params_for_analysis(config_tpl)
     config_tpl = self.options_for_slurm(config_tpl)

     # get bc directory and tile file
     in_bcsdir  = self.get_bcdir("IN")
     out_bcsdir = self.get_bcdir("OUT")
     config_tpl['input']['shared']['bcs_dir']    = in_bcsdir+ '/'
     config_tpl['output']['shared']['bcs_dir']   = out_bcsdir + '/'

     self.config = config_tpl

  def convert_to_yaml(self) :
     if os.path.exists('remap_params.yaml') :
       overwrite = questionary.confirm("Do you want to overwrite remap_params.yaml file?", default=False).ask()
       if not overwrite :
         while True:
           new_name = questionary.text("What's the backup name for existing remap_params.yaml?", default='remap_params.yaml.1').ask()
           if os.path.exists(new_name):
              print('\n'+ new_name + ' exists, please enter a new one. \n')
           else:
              shutil.move('remap_params.yaml', new_name)
              break
     yaml = ruamel.yaml.YAML()
     with open("remap_params.yaml", "w") as f:
        yaml.dump(self.config, f)

  def init_tags(self):
     # copy and paste from remap.pl
     # minor change. Add "D" to the number for each group
     # BCS Tag: Fortuna-1_4
     F14  = ( 'F14',              'Fortuna-1_4',            'Fortuna-1_4_p1' )
     D214 = ( 'D214',              'GEOSdas-2_1_4',          'GEOSdas-2_1_4-m1', 
              'GEOSdas-2_1_4-m2', 'GEOSdas-2_1_4-m3',       'GEOSdas-2_1_4-m4' )
     D540 = ( 'D540',              'GEOSadas-5_4_0',         'GEOSadas-5_4_0_p1',
              'GEOSadas-5_4_0_p2',  'GEOSadas-5_4_0_p3',    'GEOSadas-5_4_0_p4',
              'GEOSadas-5_4_1',     'GEOSadas-5_4_1_p1',    'GEOSadas-5_4_2',
              'GEOSadas-5_4_3',     'GEOSadas-5_4_4',       'GEOSadas-5_5_0',
              'GEOSadas-5_5_1',     'GEOSadas-5_5_2',       'GEOSadas-5_5_3' )

    # BCS Tag: Fortuna-2_0
    #---------------------
     F20  = ( 'F20',                    'Fortuna-2_0')

    # BCS Tag: Fortuna-2_1
    #---------------------
     F21  = ( 'F21',                  'Fortuna-2_1',         'Fortuna-2_1_p1',
              'Fortuna-2_1_p2',       'Fortuna-2_1_p3',      'Fortuna-2_2',
              'Fortuna-2_2_p1',       'Fortuna-2_2_p2',      'Fortuna-2_3',
              'Fortuna-2_3_p1',       'Fortuna-2_4',         'Fortuna-2_4_p1',
              'Fortuna-2_4_p2',       'Fortuna-2_5',         'Fortuna-2_5_BETA0',
              'Fortuna-2_5_p1',       'Fortuna-2_5_p2',      'Fortuna-2_5_p3',
              'Fortuna-2_5_p4',       'Fortuna-2_5_p5',      'Fortuna-2_5_p6',
              'Fortuna-2_5_pp2' )
     D561 = ( 'D561',                  'GEOSadas-5_6_1',      'GEOSadas-5_6_1_p1',
              'GEOSadas-5_6_1_p2',    'GEOSadas-5_6_1_p3',   'GEOSadas-5_6_1_p4',
              'GEOSadas-5_6_2',       'GEOSadas-5_6_2_p1',   'GEOSadas-5_6_2_p2',
              'GEOSadas-5_6_2_p3',    'GEOSadas-5_6_2_p4',   'GEOSadas-5_6_2_p5',
              'GEOSadas-5_6_2_p6',    'GEOSadas-5_7_1',      'GEOSadas-5_7_1_p1',
              'GEOSadas-5_7_1_p2',    'GEOSadas-5_7_2',      'GEOSadas-5_7_2_p1',
              'GEOSadas-5_7_2_p2',    'GEOSadas-5_7_2_p2_m1','GEOSadas-5_7_2_p3',
              'GEOSadas-5_7_2_p3_m1', 'GEOSadas-5_7_2_p3_m2','GEOSadas-5_7_2_p4',
              'GEOSadas-5_7_2_p5',    'GEOSadas-5_7_2_p5_m1','GEOSadas-5_7_3',
              'GEOSadas-5_7_3_p1',    'GEOSadas-5_7_3_p2',   'GEOSadas-5_7_3_p2' )

   # BCS Tag: Ganymed-1_0
    #---------------------
     G10 =  ( 'G10',                  'Ganymed-1_0',          'Ganymed-1_0_BETA',
              'Ganymed-1_0_BETA1',    'Ganymed-1_0_BETA2',    'Ganymed-1_0_BETA3',
              'Ganymed-1_0_BETA4' )

     D580 = ( 'D580',                  'GEOSadas-5_8_0',       'GEOSadas-5_9_0',
              'GEOSadas-5_9_1' )

    # BCS Tags: Ganymed-1_0_M and Ganymed-1_0_D
    #------------------------------------------
     G10p = ( 'G10p',                 'Ganymed-1_0_p1',       'Ganymed-1_0_p2',
              'Ganymed-1_0_p3',       'Ganymed-1_0_p4',       'Ganymed-1_0_p5',
              'Ganymed-1_0_p6' )

     D591p= ( 'D591p',                 'GEOSadas-5_9_1_p1',    'GEOSadas-5_9_1_p2',
              'GEOSadas-5_9_1_p3',    'GEOSadas-5_9_1_p4',    'GEOSadas-5_9_1_p5',
              'GEOSadas-5_9_1_p6',    'GEOSadas-5_9_1_p7',    'GEOSadas-5_9_1_p8',
              'GEOSadas-5_9_1_p9' )

    # BCS Tags: Ganymed-1_0_M and Ganymed-1_0_D w/ new landice rst
    #------------------------------------------------------------------------
     G20  = ( 'G20',                  'Ganymed-2_0',          'Ganymed-2_1',
              'Ganymed-2_1_p1',       'Ganymed-2_1_p2',       'Ganymed-2_1_p3',
              'Ganymed-2_1_p4',       'Ganymed-2_1_p5',       'Ganymed-2_1_p6' )
     D5A0 = ( 'D5A0',                  'GEOSadas-5_10_0',      'GEOSadas-5_10_0_p1' )


    # BCS Tags: Ganymed-1_0_Reynolds and Ganymed-1_0_Ostia
    #-----------------------------------------------------
     G30  = ( 'G30',                  'Ganymed-3_0',         'Ganymed-3_0_p1' )
     D5B0 = ( 'D5B0',                  'GEOSadas-5_10_0_p2',  'GEOSadas-5_11_0' )

    # BCS Tags: Ganymed-4_0_Reynolds, Ganymed-4_0_MERRA-2, and Ganymed-4_0_Ostia
    #---------------------------------------------------------------------------
     G40  = ( 'G40',                  'Ganymed-4_0',         'Ganymed-4_0_p1',
              'Ganymed-4_1',          'Heracles-1_0',        'Heracles-1_1',
              'Heracles-2_0',         'Heracles-2_1',        'Heracles-3_0',
              'Heracles-4_0',         'Heracles-5_4_p3' )
     D512 = ( 'D512',                  'GEOSadas-5_12_2',     'GEOSadas-5_12_4',
              'GEOSadas-5_12_4_p1',   'GEOSadas-5_12_4_p2',  'GEOSadas-5_12_4_p3',
              'GEOSadas-5_12_5',      'GEOSadas-5_13_0_p1',  'GEOSadas-5_13_0_p2',
              'GEOSadas-5_13_1',      'GEOSadas-5_16_5' )

    # BCS Tags: Icarus (New Land Parameters, New Topography)
    #---------------------------------------------------------------------------
     ICA  = ( 'ICA',                  'Icarus',              'Jason' )
     D517 = ( 'D517', 'GEOSadas-5_17_0',      'GEOSadas-5_17_1',     'GEOSadas-5_18_0',
              'GEOSadas-5_18_1',      'GEOSadas-5_18_2',     'GEOSadas-5_18_3',
              'GEOSadas-5_18_3_p1',   'GEOSadas-5_19_0',     'GEOSadas-5_20_0',
              'GEOSadas-5_20_0_p1',   'GEOSadas-5_20_0_p2',  'GEOSadas-5_21_0',
              'GEOSadas-5_21_2',      'GEOSadas-5_21_3_p1',  'GEOSadas-5_22_0',
              'GEOSadas-5_22_0_p1',   'GEOSadas-5_22_0_p2',  'GEOSadas-5_23_0',
              'GEOSadas-5_23_0_p1',   'GEOSadas-5_24_0',     'GEOSadas-5_24_0_p1' )
     GITOL = ( 'GITOL', '10.3',  '10.4',  '10.5',
               '10.6',  '10.7',  '10.8',
               '10.9',  '10.10', '10.11',
               '10.12', '10.13', '10.14',
               '10.15', '10.16', '10.17',
               '10.18'  )

    # BCS Tags: Icarus-NLv3 (New Land Parameters)
    #---------------------------------------------------------------------------
     INL  = ( 'INL', 'Icarus-NL', 'Icarus-NLv3', 'Jason-NL' )
     GITNL = ( 'GITNL', '10.19', '10.20', '10.21', '10.22', '10.23' )
     D525  = ( '525', 'GEOSadas-5_25_1', 'GEOSadas-5_25_1_p5', 'GEOSadas-5_25_p7',
                      'GEOSadas-5_27_1', 'GEOSadas-5_29_3',    'GEOSadas-5_29_4' )

     self.bcsTag={}
     for tag in F14:   self.bcsTag[tag]= "Fortuna-1_4"
     for tag in F20:   self.bcsTag[tag]= "Fortuna-2_0"
     for tag in F21:   self.bcsTag[tag]= "Fortuna-2_1"
     for tag in G10:   self.bcsTag[tag]= "Ganymed-1_0"
     for tag in G10p:  self.bcsTag[tag]= "Ganymed-1_0_M"
     for tag in G20:   self.bcsTag[tag]= "Ganymed-1_0_M"
     for tag in G30:   self.bcsTag[tag]= "Ganymed-1_0_Reynolds"
     for tag in G40:   self.bcsTag[tag]= "Ganymed-4_0_Reynolds"
     for tag in ICA:   self.bcsTag[tag]= "Icarus_Reynolds"
     for tag in GITOL: self.bcsTag[tag]= "Icarus_Reynolds"
     for tag in INL:   self.bcsTag[tag]= "Icarus-NLv3_Reynolds"
     for tag in GITNL: self.bcsTag[tag]= "Icarus-NLv3_Reynolds"

     
     for tag in D214:  self.bcsTag[tag]= "Fortuna-1_4"
     for tag in D540:  self.bcsTag[tag]= "Fortuna-1_4"
     for tag in D561:  self.bcsTag[tag]= "Fortuna-2_1"
     for tag in D580:  self.bcsTag[tag]= "Ganymed-1_0"
     for tag in D591p: self.bcsTag[tag]= "Ganymed-1_0_M"
     for tag in D5A0:  self.bcsTag[tag]= "Ganymed-1_0_M"
     for tag in D5B0:  self.bcsTag[tag]= "Ganymed-1_0_Reynolds"
     for tag in D512:  self.bcsTag[tag]= "Ganymed-4_0_Reynolds"
     for tag in D517:  self.bcsTag[tag]= "Icarus_Reynolds"
     for tag in D525:  self.bcsTag[tag]= "Icarus-NLv3_Reynolds"
    
     self.tagsRank ={}
     self.tagsRank['Fortuna-1_4'] = 1 
     self.tagsRank['Fortuna-1_5'] = 2 
     self.tagsRank['Fortuna-2_0'] = 3 
     self.tagsRank['Fortuna-2_1'] = 4 
     self.tagsRank['Ganymed-1_0'] = 5 
     self.tagsRank['Ganymed-1_0_m1'] = 6 
     self.tagsRank['Ganymed-1_0_m2'] = 7 
     self.tagsRank['Ganymed-1_0_M']  = 8 
     self.tagsRank['Ganymed-1_0_D']  = 9 
     self.tagsRank['Ganymed-1_0_Reynolds'] = 10
     self.tagsRank['Ganymed-1_0_Ostia']    = 11
     self.tagsRank['Ganymed-4_0_Reynolds'] = 12
     self.tagsRank['Ganymed-4_0_Ostia']    = 13
     self.tagsRank['Ganymed-4_0_MERRA-2']  = 14
     self.tagsRank['Icarus_Reynolds']      = 15
     self.tagsRank['Icarus_MERRA-2']       = 16
     self.tagsRank['Icarus_Ostia']         = 17
     self.tagsRank['Icarus-NLv3_Reynolds'] = 18
     self.tagsRank['Icarus-NLv3_MERRA-2']  = 19
     self.tagsRank['Icarus-NLv3_Ostia']    = 20

     self.bcbase={}
     self.bcbase['discover_ops'] = "/discover/nobackup/projects/gmao/share/gmao_ops/fvInput/g5gcm/bcs"
     self.bcbase['discover_lt']  = "/discover/nobackup/ltakacs/bcs"
     self.bcbase['discover_couple']  = "/discover/nobackup/projects/gmao/ssd/aogcm/atmosphere_bcs"

  def init_time(self):
     ymdh = self.common_in.get('yyyymmddhh')
     self.yyyymm = ymdh[0:6]
     self.yyyy = ymdh[0:4]  
     self.mm   = ymdh[4:6]  
     self.dd   = ymdh[6:8]  
     self.hh   = ymdh[8:10]  
     self.ymd  = ymdh[0:8]  

  def init_merra2(self):
    def get_grid_kind(grid):
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

    if not self.common_in['MERRA-2']:
      return
    print("\n MERRA-2 sources:\n")
    yyyymm = int(self.yyyymm)
    if yyyymm < 197901 :
      exit("Error. MERRA-2 data < 1979 not available\n")
    elif (yyyymm < 199201):
      self.common_in['expid'] = "d5124_m2_jan79"     
    elif (yyyymm < 200106):
      self.common_in['expid'] = "d5124_m2_jan91"
    elif (yyyymm < 201101):
      self.common_in['expid'] = "d5124_m2_jan00"
    else:
      self.common_in['expid'] = "d5124_m2_jan10"

    self.common_in['agrid'] = 'C180'
    self.common_in['ogrid'] = '1440x720'
    self.common_in['bc_base']= 'discover_ops'
    self.common_in['tag']= 'Ganymed-4_0'
    expid = self.common_in['expid']
    yyyymmddhh_ = str(self.common_in['yyyymmddhh'])
    surfix = yyyymmddhh_[0:8]+'_'+self.hh+'z.bin'
    merra_2_rst_dir = '/archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/'+expid +'/rs/Y'+self.yyyy +'/M'+self.mm+'/'
    rst_dir = self.common_in['rst_dir'] + '/'
    os.makedirs(rst_dir, exist_ok = True)
    print(' Copy MERRA-2 Restart \n from \n    ' + merra_2_rst_dir + '\n to\n    '+ rst_dir +'\n') 

    upperin =[merra_2_rst_dir +  expid+'.fvcore_internal_rst.' + surfix,
              merra_2_rst_dir +  expid+'.moist_internal_rst.'  + surfix,
              merra_2_rst_dir +  expid+'.agcm_import_rst.'     + surfix,
              merra_2_rst_dir +  expid+'.gocart_internal_rst.' + surfix,
              merra_2_rst_dir +  expid+'.pchem_internal_rst.'  + surfix ]

    surfin = [ merra_2_rst_dir +  expid+'.catch_internal_rst.'    + surfix, 
               merra_2_rst_dir +  expid+'.lake_internal_rst.'     + surfix,
               merra_2_rst_dir +  expid+'.landice_internal_rst.'  + surfix,
               merra_2_rst_dir +  expid+'.saltwater_internal_rst.'+ surfix]

    for f in upperin :
       fname = os.path.basename(f)
       dest = rst_dir + '/'+fname
       print("Copy file "+f +" to " + rst_dir)
       shutil.copy(f, dest)

    for f in surfin :
       fname = os.path.basename(f)
       dest = rst_dir + '/'+fname
       print("Copy file "+f +" to " + rst_dir)
       shutil.copy(f, dest)

    # prepare analysis files
    bkg = self.ana_out['bkg']
    if ( not bkg ): return
    yyyy_ = yyyymmddhh_[0:4]
    mm_   = yyyymmddhh_[4:6]
    dd_   = yyyymmddhh_[6:8]
    hh_   = yyyymmddhh_[8:10]
    rst_time = datetime(year=int(yyyy_), month=int(mm_), day=int(dd_), hour = int(hh_))
    expid_in  = self.common_in['expid']
    expid_out = self.common_out['expid']
    if (expid_out) :
       expid_out = expid_out + '.'
    else:
       expid_out = ''

    agrid_in  = self.common_in['agrid']
    agrid_out = self.common_out['agrid']

    if (get_grid_kind(agrid_in.upper()) == get_grid_kind(agrid_out.upper())):
      print(" No need to remap anaylysis file according to air grid in and out")
      return

    anafiles=[]
    for h in [3,4,5,6,7,8,9]:
       delt = timedelta(hours = h-3)
       new_time = rst_time + delt
       yyyy = "Y"+str(new_time.year)
       mm   = 'M%02d'%new_time.month
       ymd  = '%04d%02d%02d'%(new_time.year,new_time.month, new_time.day)
       hh   = '%02d'%h
       newhh= '%02d'%new_time.hour
       m2_rst_dir = merra_2_rst_dir.replace('Y'+yyyy_,yyyy).replace('M'+mm_,mm)
       # bkg files
       for ftype in ['sfc', 'eta']:
          fname = expid_in+'.bkg'+hh+'_'+ftype+'_rst.'+ymd+'_'+newhh+'z.nc4'
          f = m2_rst_dir+'/'+fname
          if(os.path.isfile(f)):
             anafiles.append(f)
          else:
             print('Warning: Cannot find '+f)
       # cbkg file
       fname = expid_in + '.cbkg' + hh + '_eta_rst.' + ymd + '_' + newhh + 'z.nc4'
       f = m2_rst_dir+'/'+fname
       if(os.path.isfile(f)):
         anafiles.append(f)
       else:
         print('Warning: Cannot find '+f)
       # gaas_bkg_sfc files
       if (h==6 or h==9):
          fname = expid_in+'.gaas_bkg_sfc_rst.'+ymd+'_'+newhh+'z.nc4'
          f = m2_rst_dir+'/'+fname
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
    m2_rst_dir = merra_2_rst_dir.replace('Y'+yyyy_,yyyy).replace('M'+mm_,mm)
    fname = expid_in+'.trak.GDA.rst.'+ymdh+'z.txt'
    f = m2_rst_dir+'/'+fname
    if (os.path.isfile(f)): anafiles.append(f)

    for f in anafiles:
      fname    = os.path.basename(f)
      f_tmp = rst_dir+'/'+fname
      print("Copy file "+f +" to " + rst_dir)
      shutil.copy(f,f_tmp)

  def get_bcbase(self, opt):
     base = ''
     model = ''
     
     if opt.upper() == 'IN':
        model = self.common_in.get('model')
        if model == 'MOM6' or model == 'MOM5':
          base = 'discover_couple'
        else:
          base  = self.common_in.get('bc_base')

     if opt.upper() == 'OUT':
        model = self.common_out.get('model')
        if model == 'MOM6' or model == 'MOM5':
          base = 'discover_couple'
        else:
          base = self.common_out.get('bc_base')
     assert base, 'please specify bc_base: discover_ops, discover_lt, discover_couple or an absolute path'
     if base == 'discover_ops' or base == 'discover_lt' or base=='discover_couple':
        return self.bcbase[base]
     else:
        return base 

  def get_bcdir(self, opt):
    tag = self.common_in['tag']
    ogrid = self.common_in['ogrid']
    model = self.common_in['model']
    bcdir = self.common_in.get('alt_bcs', None)
    if opt.upper() == "OUT":
       tag = self.common_out['tag']
       ogrid = self.common_out['ogrid']
       model = self.common_out['model']
       bcdir = self.common_out.get('alt_bcs', None)

    if bcdir is None :
       bc_base = self.get_bcbase(opt)
       bctag = self.get_bcTag(tag,ogrid)
       tagrank = self.tagsRank[bctag]
       if (tagrank >= self.tagsRank['Icarus-NLv3_Reynolds']) :
          bcdir = bc_base+'/Icarus-NLv3/'+bctag+'/'
          if model == 'MOM6' or model == 'MOM5':
             bcdir = bc_base+'/Icarus-NLv3/'+model+'/'
       elif (tagrank >= self.tagsRank['Icarus_Reynolds']):
          if bc_base == self.bcbase['discover_ops']:
             bcdir = bc_base+'/Icarus_Updated/'+bctag+'/'
          else:
             bcdir = bc_base+'/Icarus/'+bctag+'/'
          if model == 'MOM6' or model == 'MOM5':
             bcdir = bc_base+'/Icarus/'+model+'/'
       elif(tagrank >= self.tagsRank["Ganymed-4_0_Reynolds"]):
          bcdir = bc_base + '/Ganymed-4_0/'+bctag+'/'
          if model == 'MOM6' or model == 'MOM5':
             bcdir = bc_base+'/Ganymed/'+model+'/'
       else:
          bcdir = bc_base + '/' + bctag + '/'
          if model == 'MOM6' or model == 'MOM5':
             bcdir = bc_base+'/Ganymed/'+model+'/'

    if not os.path.exists(bcdir):
       exit("Cannot find bc dir " +  bcdir)

    gridStr = self.get_grid_subdir(bcdir,opt)
    bcdir =  bcdir + '/' + gridStr

    return bcdir
        
  def get_grid_subdir(self, bcdir, opt):

     def get_name_with_grid( grid, names, a_o):
       if not grid :
         return names
       namex = []
       if (grid[0].upper() == 'C'):
         n = int(grid[1:])
         s1 ='{n}x6C'.format(n=n)
         j=n*6
         s2 =str(n)
         s3 =str(j)
         # first try
         for aoname in names:
           name = ''
           if(a_o == 'a'):
             name = aoname.split('_')[0]
           else:
             name = aoname.split('_')[1]
           if (name.find(s1) != -1 or (name.find(s2) != -1 and name.find(s3) != -1 )):
              namex.append(aoname)
       else:
         xy = grid.upper().split('X')
         s2 = xy[0]
         s3 = xy[1]
         for aoname in names:
           name = ''
           if(a_o == 'a'):
             name = aoname.split('_')[0]
           else:
             name = aoname.split('_')[1]
           if (name.find(s2) != -1 and name.find(s3) != -1): namex.append(aoname)
       return namex
     #v3.5
     #dirnames = [ f.name for f in os.scandir(bcdir) if f.is_dir()]
     #v2.7
     dirnames = [f for f in os.listdir(bcdir) if os.path.isdir(os.path.join(bcdir,f))]
     agrid_ = self.common_in['agrid']   
     ogrid_ = self.common_in['ogrid'] 
     if opt.upper() == "OUT" :
       agrid_ = self.common_out['agrid']   
       ogrid_ = self.common_out['ogrid'] 
          
     anames = get_name_with_grid(agrid_, dirnames, 'a')
     gridID = get_name_with_grid(ogrid_, anames, 'o')
     if len(gridID) == 0 :
       exit("cannot find the grid string: " + bcdir)
     g = ''
     if len(gridID) == 1 : g = gridID[0]
     if len(gridID) >=2 :
       print("find too many grid strings in " + bcdir)
       print(" gridIDs found", gridID)
       for g_ in gridId:
         if g_.count('_') == 1 :
           g = g_ 
           #WY note, found many string in couple model
           print(" pick the first directory with only one '_' " + g)
           break
     return g

  def get_bcTag(self, tag, ogrid):
    bctag = self.bcsTag[tag]
    if ogrid[0].upper() == "C": 
       bctag=bctag.replace('_Reynolds','_Ostia')
    else:
       xy = ogrid.upper().split('X')
       x = int(xy[0])
       if x == 1440:   bctag=bctag.replace('_Reynolds','_MERRA-2')
       if x == 2880:  
          bctag=bctag.replace('_Reynolds','_Ostia')
          bctag=bctag.replace('_M','_D')
    return bctag

  def params_for_air(self, config_tpl):
     if self.common_in['MERRA-2']:
       return config_tpl 
     # verify agrid
     rst_dir = self.common_in['rst_dir'] + '/'
     time = self.ymd + '_'+ self.hh
     files = glob.glob(rst_dir +'/*fvcore_*'+time+'*')
     if len(files) == 0 :
       fname_ = rst_dir +'/fvcore_internal_rst'
       if os.path.exists(fname_) :
         files.append(fname_)

     # get expid
     if (len(files) >0) : 
        fname = os.path.basename(files[0])
        expid = fname.split('fvcore')[0]
        config_tpl['input']['shared']['expid'] = expid[0:-1] #remove the last '.'

     agrid_ = self.common_in['agrid']
     if self.common_in['ogrid'] == 'CS' :
        config_tpl['input']['shared']['ogrid']  = agrid_ 
        self.common_in['ogrid'] = agrid_
        
     ogrid = config_tpl['input']['shared']['ogrid']
     tagout = self.common_out['tag']
     bctag  = self.get_bcTag(tagout, ogrid)
     tagrank = self.tagsRank[bctag]
     if ( not config_tpl['input']['air']['drymass']) :
        config_tpl['input']['air']['drymass'] = 0
        if tagrank >=12 :
          config_tpl['input']['air']['drymass'] = 1

     return config_tpl

  def options_for_slurm(self, config_tpl):
    config_tpl['slurm']['account'] = self.slurm_options['account']
    config_tpl['slurm']['qos'] = self.slurm_options['qos']
    config_tpl['slurm']['constraint'] = self.slurm_options['constraint']
    return config_tpl

  def params_for_surface(self, config_tpl):
    config_tpl['output']['surface']['surflay'] = 20.
    tagout = self.common_out['tag']
    ogrid = self.common_out['ogrid']
    bctag = self.get_bcTag(tagout, ogrid)
    tagrank = self.tagsRank[bctag]
    if tagrank >=12 :
       config_tpl['output']['surface']['surflay'] = 50.
    if tagrank >= self.tagsRank["Icarus_Reynolds"]:
       config_tpl['output']['surface']['split_saltwater'] = True
    config_tpl['input']['surface']['zoom']= self.surf_in['zoom']
    config_tpl['input']['surface']['wemin']= self.surf_in['wemin']
    config_tpl['output']['surface']['wemin']= self.surf_out['wemout']

    rst_dir = self.common_in['rst_dir'] + '/'
    time = self.ymd + '_'+ self.hh
    files = glob.glob(rst_dir +'/*catch_*'+time+'*')
    if (len(files)== 0) :
       files = glob.glob(rst_dir +'/*catch_*')

    if (len(files) > 0) :
        config_tpl['input']['surface']['catch_model'] = 'catch'

    files = glob.glob(rst_dir +'/*catchcnclm40_*'+time+'*')
    if (len(files)== 0) :
       files = glob.glob(rst_dir +'/*catchcnclm40_*')

    if (len(files) > 0) :
        config_tpl['input']['surface']['catch_model'] = 'catchcnclm40'

    files = glob.glob(rst_dir +'/*catchcnclm45_*'+time+'*')
    if (len(files)== 0) :
       files = glob.glob(rst_dir +'/*catchcnclm45_*')

    if (len(files) > 0) :
        config_tpl['input']['surface']['catch_model'] = 'catchcnclm45'

    return config_tpl

  def params_for_analysis(self, config_tpl):
    config_tpl['output']['analysis']['lcv'] = self.ana_out.get('lcv')
    config_tpl['output']['analysis']['bkg'] = self.ana_out.get('bkg')
    
    ogrid = self.common_out['ogrid']
    tagout = self.common_out['tag']
    bctag  = self.get_bcTag(tagout, ogrid)
    tagrank = self.tagsRank[bctag]
    if tagrank >= self.tagsRank["Ganymed-4_0_Reynolds"] :
      config_tpl['output']['analysis']['aqua'] = True
    return config_tpl

if __name__ == "__main__":
  yaml = ruamel.yaml.YAML()
  stream =''
  with open("raw_answers.yaml", "r") as f:
     stream = f.read()
  config = yaml.load(stream)
  param = remap_params(config) 
  param.convert_to_yaml() 
