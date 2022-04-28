#!/usr/bin/env python3
#
import os
import yaml
import shutil
import glob
import time
import questionary
from datetime import datetime
from datetime import timedelta

class regrid_params(object):
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
     self.init_restarts_in()

     # get bc directory and tile file
     self.in_bcsdir  = self.get_bcdir("IN")
     self.in_til     = glob.glob(self.in_bcsdir+ '/*-Pfafstetter.til')[0] 
     print("\ninput tile file: " + self.in_til)
     self.out_bcsdir = self.get_bcdir("OUT")
     self.out_til     = glob.glob(self.out_bcsdir+ '/*-Pfafstetter.til')[0] 
     print("\noutput tile file: " + self.out_til)

     # load input yaml
     stream     = open('regrid_params.tpl', 'r')
     config_tpl = yaml.full_load(stream)
     
     # params for shared

     config_tpl['input']['shared']['agrid']    = self.common_in['agrid']
     config_tpl['input']['shared']['ogrid']    = self.common_in['ogrid']
     config_tpl['input']['shared']['bcs_dir']  = self.in_bcsdir
     config_tpl['input']['shared']['rst_dir']  = self.common_in['rst_dir']
     config_tpl['input']['shared']['expid']    = self.common_in['expid']
     config_tpl['input']['shared']['yyyymmddhh'] = self.common_in['yyyymmddhh']

     config_tpl['output']['shared']['agrid']   = self.common_out['agrid']
     config_tpl['output']['shared']['ogrid']   = self.common_out['ogrid']
     config_tpl['output']['shared']['bcs_dir'] = self.out_bcsdir
     config_tpl['output']['shared']['out_dir'] = self.common_out['out_dir']
     config_tpl['output']['shared']['expid']   = self.common_out['expid']

     # params for upper air
     config_tpl = self.params_for_air(config_tpl)
     config_tpl = self.params_for_surface(config_tpl)
     config_tpl = self.params_for_analysis(config_tpl)
     config_tpl = self.options_for_slurm(config_tpl)

     self.config = config_tpl

  def convert_to_yaml(self) :
     if os.path.exists('regrid_params.yaml') :
       overwrite = questionary.confirm("Do you want to overwrite regrid_params.yaml file?", default=False).ask()
       if not overwrite :
         while True:
           new_name = questionary.text("What's the backup name for existing regrid_params.yaml?", default='regrid_params.yaml.1').ask()
           if os.path.exists(new_name):
              print('\n'+ new_name + ' exists, please enter a new one. \n')
           else:
              shutil.move('regrid_params.yaml', new_name)
              break
     f = open("regrid_params.yaml", "w")
     yaml.dump(self.config, f, allow_unicode=True, default_flow_style=False)
     f.close()

  def init_tags(self):
     # copy and paste from regrid.pl
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
     yyyymmddhh = str(self.common_in['yyyymmddhh'])
     self.yyyymm = yyyymmddhh[0:6]
     self.yyyy = yyyymmddhh[0:4]  
     self.mm   = yyyymmddhh[4:6]  
     self.dd   = yyyymmddhh[6:8]  
     self.hh   = yyyymmddhh[8:10]  
     self.ymd  = yyyymmddhh[0:8]  

  def init_restarts_in(self):
    if self.common_in['MERRA-2']:
      return
      
    rst_dir = self.common_in.get('rst_dir')+'/'
    self.restarts_in={}
    self.restarts_in['UPPERAIR'] = glob.glob(rst_dir +'upperair/*')
    self.restarts_in['SURFACE'] = glob.glob(rst_dir +'surface/*')
    self.restarts_in['ANALYSIS'] = glob.glob(rst_dir +'analysis/*')
    
  def init_merra2(self):
    if not self.common_in['MERRA-2']:
      return
    print("\nMERRA-2 sources:\n")
    yyyymm = int(self.yyyymm)
    if yyyymm < 197901 :
      print("Error. MERRA-2 data < 1979 not available\n")
      exit()
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
    yyyymmddhh = str(self.common_in['yyyymmddhh'])
    surfix = yyyymmddhh[0:8]+'_'+self.hh+'z.bin'
    self.common_in['rst_dir'] = '/archive/users/gmao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/'+expid +'/rs/Y'+self.yyyy +'/M'+self.mm+'/'
    print('\nMERRA-2 Restart dir: ' + self.common_in['rst_dir'] +'\n') 

    self.restarts_in = {}
    upperin =[self.common_in['rst_dir']+  expid+'.fvcore_internal_rst.' + surfix,
              self.common_in['rst_dir']+  expid+'.moist_internal_rst.'  + surfix,
              self.common_in['rst_dir']+  expid+'.agcm_import_rst.'     + surfix,
              self.common_in['rst_dir']+  expid+'.gocart_internal_rst.' + surfix,
              self.common_in['rst_dir']+  expid+'.pchem_internal_rst.'  + surfix ]
    self.restarts_in['UPPERAIR'] = upperin 

    surfin = [self.common_in['rst_dir']+  expid+'.catch_internal_rst.'    + surfix, 
              self.common_in['rst_dir']+  expid+'.lake_internal_rst.'     + surfix,
              self.common_in['rst_dir']+  expid+'.landice_internal_rst.'  + surfix,
              self.common_in['rst_dir']+  expid+'.saltwater_internal_rst.'+ surfix]
    self.restarts_in['SURFACE'] = surfin 
    self.restarts_in['ANALYSIS'] = []

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
       bcdir = self.common_in.get('alt_bcs', None)

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
       print( "Cannot find bc dir " +  bcdir)
       exit()

    gridStr = self.get_grid_subdir(bcdir,opt)
    bcdir =  bcdir + '/' + gridStr

    return bcdir
        
  def get_grid_subdir(self, bcdir, opt):

     def get_name_with_grid( grid, names):
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
         for name in names:
           if (name.find(s1) != -1):
              namex.append(name)
         if len(namex) ==0:
           for name in names:
             if (name.find(s2) != -1 and name.find(s3) != -1): namex.append(name)
       else:
         xy = grid.upper().split('X')
         s2 = xy[0]
         s3 = xy[1]
         for name in names:
           if (name.find(s2) != -1 and name.find(s3) != -1): namex.append(name)
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
          
     anames = get_name_with_grid(agrid_, dirnames)
     gridID = get_name_with_grid(ogrid_, anames)
     if len(gridID) == 0 :
       print("cannot find the grid string: " + bcdir)
       exit()
     if len(gridID) >=2 :
       print("find too may grid strings in " + bcdir)
       print(" gridIDs found", gridID)
       print(" pick the first one " + gridID[0])
     return gridID[0]

  def get_bcTag(self, tag, ogrid):
    bctag = self.bcsTag[tag]
    if ogrid[0].upper() == "C": 
       bctag=bctag.replace('_Reynolds','_Ostia')
    else:
       xy = ogrid.upper().split('X')
       x = int(xy[0])
       if x == 1440:   bctag=bctag.replace('_Reynolds','_MERRA-2')
       if x == 2800:  
          bctag=bctag.replace('_Reynolds','_Ostia')
          bctag=bctag.replace('_M','_D')
    return bctag

  def params_for_air(self, config_tpl):
     # verify agrid
     agrid = config_tpl['input']['shared']['agrid']
     fvrst = os.path.dirname(os.path.realpath(__file__)) + '/fvrst.x -h '
     if not self.common_in.get('MERRA-2'):
       fvcore = ''
       for f in self.restarts_in:
           if 'fvcore' in f:
             fvcore = f
             break
       cmd = fvrst + fvcore
       print(cmd +'\n')
       p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
       (output, err) = p.communicate()
       p_status = p.wait()
       ss = output.decode().split()
       if (agrid):
         if agrid[0].upper() == "C":
            n=int(agrid[1:])
            o=int(ss[0])
            assert n==o, "input agrid is not consistent with fvcore restart"
       else:
         config_tpl['input']['shared']['agrid'] = "C"+ss[0]

     ogrid = config_tpl['input']['shared']['ogrid']
     tagout = self.common_out['tag']
     bctag  = self.get_bcTag(tagout, ogrid)
     tagrank = self.tagsRank[bctag]
     if ( not config_tpl['input']['air']['drymass']) :
        config_tpl['input']['air']['drymass'] = 0
        if tagrank >=12 :
          config_tpl['input']['air']['drymass'] = 1

     config_tpl['output']['air']['nlevel'] = self.upper_out['nlevel']

     return config_tpl

  def options_for_slurm(self, config_tpl):
    config_tpl['slurm']['account'] = self.slurm_options['account']
    config_tpl['slurm']['debug'] = self.slurm_options['debug']
    config_tpl['slurm']['partition'] = self.slurm_options['partition']
    return config_tpl

  def params_for_surface(self, config_tpl):
    config_tpl['output']['surface']['surflay'] = 20.
    config_tpl['output']['surface']['rescale'] = False
    tagout = self.common_out['tag']
    ogrid = self.common_out['ogrid']
    bctag = self.get_bcTag(tagout, ogrid)
    tagrank = self.tagsRank[bctag]
    if tagrank >=12 :
       config_tpl['output']['surface']['surflay'] = 50.
    if tagrank > self.tagsRank["Fortuna-2_0"]:
       config_tpl['output']['surface']['rescale'] = True
    if tagrank >= self.tagsRank["Icarus_Reynolds"]:
       config_tpl['output']['surface']['split_saltwater'] = True
    config_tpl['output']['surface']['zoom']= self.surf_in['zoom']
    config_tpl['output']['surface']['wemin']= self.surf_in['wemin']
    config_tpl['output']['surface']['wemout']= self.surf_out['wemout']
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
    
