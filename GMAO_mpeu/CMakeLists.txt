esma_set_this()

add_subdirectory(mpi0)

set(srcs
  get_zeits.c
  m_mpif.F m_random.F
  m_FileResolv.f90
  m_date.F90
  m_datetime.F90
  m_IndexBin_char.F90
  m_IndexBin_integer.F90
  m_IndexBin_logical.F90
  m_List.F90
  m_MergeSorts.F90
  m_Filename.F90
  m_Group.F90
  m_Permuter.F90
  m_SortingTools.F90
  m_String.F90
  m_StrTemplate.F90
  m_chars.F90
  m_die.F90
  m_dropdead.F90
  m_flow.F90
  m_inpak90.F90
  m_intkinds.F90
  m_ioutil.F90
  m_mall.F90
  m_mergedList.F90
  m_mpif90.F90
  m_mpout.F90
  m_parDOT.F90
  m_rankMerge.F90
  m_realkinds.F90
  m_showDistrib.F90
  m_stdio.F90
  m_zeit.F90
  m_sysclocks.F90
  m_ior4.F90
  m_rGrADS.F90
  m_wGrADS.F90
  m_GrADS.F90
  m_fpe.F90
  )

esma_add_library(${this} SRCS ${srcs} DEPENDENCIES MPI::MPI_Fortran)

esma_add_library(GMAO_eu SRCS ${srcs} DEPENDENCIES mpi0)
target_include_directories (GMAO_eu PRIVATE $<BUILD_INTERFACE:${MPI_Fortran_INCLUDE_PATH}>)
file (MAKE_DIRECTORY ${esma_include}/GMAO_eu)
file (MAKE_DIRECTORY ${esma_etc}/GMAO_eu)

#############################################################
# Exceptional cases:
# One source file uses the preprocessor but is has the .f90 suffix (lower case).
set_source_files_properties (m_FileResolv.f90 PROPERTIES COMPILE_FLAGS ${PP})

# NAG notices if the same procedure is called with varying types in the absence
# of an explicit interface.
target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${MISMATCH}>)
target_compile_options (GMAO_eu PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${MISMATCH}>)

file (COPY assert.H DESTINATION ${include_${this}})
file (COPY assert.H DESTINATION ${esma_include}/GMAO_eu)


