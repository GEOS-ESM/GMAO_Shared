input:
  air:
    drymass: 1
    hydrostatic: 0
    nlevel:
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
    catchCNclm40:
      regrid: false
      rst_dir: null
    catchCNclm45:
      regrid: false
      rst_dir: null
    catchment:
      regrid: false
      rst_dir: null
    lakelandice:
      regrid: false
      rst_dir: null

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
    rescale: true
    surflay: 20.
    wemout:
  analysis:
    bkg: true
    aqua: False
    lcv: false

slurm:
  account: 
  debug: 
  partition: 