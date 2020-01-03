esma_set_this ()

set (srcs
  Chem_MieMod2G.F90
  Chem_MieTableMod.F90
  )

esma_add_library(${this} 
  SRCS ${srcs} 
  DEPENDENCIES MAPL_Base)

target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

#include_directories (${esma_include}/${this})

