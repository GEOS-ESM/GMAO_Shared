#
# Earth System Modeling Applications (ESMA) PILGRIM makefile fragment.
# This fragment customize GNUmakefile for each architecture. 
#
# REVISION HISTORY:
#
# 01May2006  Todling   Cannot run w/ option -fpe0 (ana5sfc crashes)
#
#--------------------------------------------------------------------------

_STATIC =
	# _STATIC:	Tell the fortran compiler to treate all local
	#		variables as with SAVE attribute (required
	#		by elliptic.f).

_ELLIPTIC = 

#   -----
#   LINUX
#   -----

ifeq ($(ARCH),Linux)

FPE =
ifeq ($(FC),ifort)
  _STATIC = -save
  _ELLIPTIC = -vec-report0   -ftz   -align all   -fno-alias   -fpe0   -fp-model source -fp-model source
endif
ifeq ($(FC),lf95)
  _STATIC = -sav
endif

endif  #    Linux

ifeq ($(ARCH),IRIX64)
  _STATIC = -static
endif	# IRIX64
