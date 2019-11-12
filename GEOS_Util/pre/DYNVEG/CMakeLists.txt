# No library
set (srcs regrid_DYNVEG.F90)

include_directories(${INC_ESMF})
include_directories(${INC_NETCDF})

esma_add_library(${this} SRCS ${srcs} DEPENDENCIES MAPL_Base)

ecbuild_add_executable (TARGET regrid_DYNVEG.x SOURCES regrid_DYNVEG.F90  LIBS MAPL_Base ${this})

install(PROGRAMS regrid_DYNVEG.csh process_modis.py DESTINATION bin)
