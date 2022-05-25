#
# This template file can be filled with questionary or manually
# 
# By default, the values in child dictionary overwrites the values
# in parent dictionary. For example, rst_dir in input is the location with 
# restart files to be regrrided. The surface restart files can be in
# surfface:rst_dir if it has value 
# 

input:
  air:
    drymass: 1
    hydrostatic: 0
    rst_dir: null
  analysis:
    rst_dir: null
  shared:
    agrid: 
    bcs_dir:
    expid:
    ogrid:
    rst_dir:
    yyyymmddhh:
  surface:
    wemin:
    rst_dir: null
    tile_file: null
    catchcnclm40:
      regrid: false
      rst_file: null
    catchcnclm45:
      regrid: false
      rst_file: null
    catchment:
      regrid: false
      rst_file: null
    lakelandice:
      regrid: false
      rst_file: null

output:
  shared:
    agrid:
    bcs_dir: 
    expid: 
    ogrid: 
    out_dir:
  air:
    nlevel:
  surface:
    zoom:
    split_saltwater: false
    surflay: 20.
    wemin:
    tile_file: none
  analysis:
    bkg: true
    aqua: False
    lcv: false

slurm:
  account: 
  debug: 
  partition: 
