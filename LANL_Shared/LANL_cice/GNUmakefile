#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 09Jun2003  da Silva  First crack.
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # Generic stuff

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"

THIS := $(shell basename `pwd`)
LIB   = lib$(THIS).a

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          install clean distclean doc 

export ESMADIR BASEDIR ARCH SITE

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done
	$(MAKE) local_$@

local_esma_install local_install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)
#	$(CP) -p GEOS_ErrLog.h  $(ESMAINC)/$(THIS)

local_esma_clean local_clean:
	-$(RM) *~ *.[aox] *.mod *.x *.ps

local_esma_distclean local_distclean:
	-$(RM) *~ *.[aoxd] *.mod *.x *.ps

local_esma_doc local_doc:
	echo ""; echo Making *GridComp*.[fF]* in `pwd`          ;\
	$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex

CPPDEFS    := -DLINUX
# if add -Dcoupled, cice will use true sea level tilt;
# otherwise, tilt terms are computed based on surface currents
CPPDEFS :=  $(CPPDEFS) -Dcoupled -DGEOS -DMODIFY_TOPOGRAPHY -DUSE_R8

ifeq ($(findstring dataatm,$(options)),dataatm)
    CPPDEFS += -DDATAATM
endif

ifneq ($(PGI_MAJOR),'')
    CPPDEFS += -D__PGI
endif

USER_FDEFS := $(CPPDEFS)

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)
#

#-------------------------------------------------------------------------------
# CVS $Id$
# CVS $Source$
# CVS $Name$
#-------------------------------------------------------------------------------
# Common Makefile: a framework for building all CCSM components and more
#
# Command-line variables
#   MACFILE=<file> ~ the macros definition file to use/include
#   EXEC=<name>    ~ name given to executable, default is a.out
#   VPATH=<vpath>  ~ VPATH            , default is . (cwd only)
#   SRCS=<files>   ~ list of src files, default is all .c .F .F90 files in VPATH
#   VPFILE=<file>  ~ file with list of dirs, used to create VPATH
#   SRCFILE=<file> ~ file with list of src files, used to create SRCS
#   DEPGEN=<exec>  ~ dependency generator utility, default is makdep
#
#   <macro defns>  ~ any macro definitions found in this file or the included 
#                    MACFILE will be over-riden by cmd-line macro definitions
#   MODEL=<model>  ~ a standard macro definition, often found in the included 
#                    MACFILE, used to trigger special compilation flags
#
# Usage examples:
#   % gmake MACFILE=Macros.AIX VPFILE=Filepath MODEL=ccm3 EXEC=atm
#   % gmake MACFILE=Macros.AIX VPFILE=Filepath SRCFILE=Srclist EXEC=pop
#   % gmake MACFILE=Macros.C90 VPATH="dir1 dir2" SRCS="file1.c file2.F90"
#   % gmake MACFILE=Macros.SUN SRCS="test.F"
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# parse cmd-line and establish values for EXEC, VPATH, SRCS, OBJS, etc
#-------------------------------------------------------------------------------

EXEC     := a.out
MACFILE  := NONE 
MODEL    := NONE 
VPFILE   := NONE
#VPATH    := .
VPATH    := drivers/cice4 csm_share mpi source
SRCFILE  := NONE
#SRCS     := NONE
SRCS     := shr_orb_mod.F90 ice_constants.F90 ice_prescribed_mod.F90 \
	ice_boundary.F90 ice_broadcast.F90 ice_communicate.F90 \
	ice_exit.F90 ice_gather_scatter.F90 ice_global_reductions.F90 ice_timers.F90 \
	ice_age.F90 ice_atmo.F90 ice_blocks.F90 ice_calendar.F90 ice_diagnostics.F90 \
	ice_distribution.F90 ice_domain.F90 ice_domain_size.F90 ice_dyn_evp.F90 \
	ice_fileunits.F90 ice_flux.F90 ice_forcing.F90 ice_grid.F90 ice_history.F90 \
	ice_init.F90 ice_itd.F90 ice_kinds_mod.F90 ice_mechred.F90 ice_meltpond.F90 \
	ice_ocean.F90 ice_orbital.F90 ice_read_write.F90 ice_restart.F90 ice_restoring.F90 \
	ice_shortwave.F90 ice_spacecurve.F90 ice_state.F90 ice_step_mod.F90 ice_therm_itd.F90 \
	ice_therm_vertical.F90 ice_transport_driver.F90 ice_transport_remap.F90 ice_work.F90

DEPGEN   := $(ESMABIN)/makdep  # an externally provided dependency generator


ifneq ($(VPATH),.)
  # this variable was specified on cmd line or in an env var
else
  ifneq ($(VPFILE),NONE)
    # explicit list of VPATH dirs is provided
    VPATH := $(wildcard . $(shell cat $(VPFILE) ) )
  endif
endif

ifneq ($(SRCS),NONE)
  # this variable was specified on cmd line or in an env var
else
  ifneq ($(SRCFILE),NONE)
    # explicit list of src files is provided
    SRCS := $(shell cat $(SRCFILE) )
  else
    # list of src files is all .F90 .F .c files in VPATH
    SRCS := $(wildcard $(addsuffix /*.F90 , $(VPATH)) \
		       $(addsuffix /*.[cF], $(VPATH)) )
  endif
endif

OBJS  := $(addsuffix .o, $(sort $(basename $(notdir $(SRCS)))))
DEPS  := $(addsuffix .d, $(sort $(basename $(notdir $(SRCS)))))
INCS  := $(patsubst %,-I%, $(VPATH) $(INC_MAPL_BASE) )

INC_DIRS = . $(INC_MAPL_BASE)
MOD_DIRS = $(INC_DIRS)

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir))

.SUFFIXES:
.SUFFIXES: .F90 .F .c .o

all: $(LIB)

#-------------------------------------------------------------------------------
# include the file that provides macro definitions required by build rules
# note: the MACFILE may not be needed for certain goals
#-------------------------------------------------------------------------------
ifneq ($(MAKECMDGOALS), db_files)
  -include $(MACFILE)
endif

#-------------------------------------------------------------------------------
# echo file names, paths, compile flags, etc. used during build
#-------------------------------------------------------------------------------

db_files:
	@echo " "
	@echo "* EXEC    := $(EXEC)"
	@echo "* MACFILE := $(MACFILE)"
	@echo "* VPFILE  := $(VPFILE)"
	@echo "* VPATH   := $(VPATH)"
	@echo "* SRCFILE := $(SRCFILE)"
	@echo "* INCS    := $(INCS)"
	@echo "* SRCS    := $(SRCS)"
	@echo "* OBJS    := $(OBJS)"
	@echo "* DEPS    := $(DEPS)"
db_flags:
	@echo " "
	@echo "* cpp     := $(CPP) $(CPPFLAGS) $(CPPDEFS) $(INCS) $(INCLDIR)"
	@echo "* cc      := cc -c $(CFLAGS) $(INCS) $(INCLDIR)"
	@echo "* .F.o    := $(FC) -c $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR)"
	@echo "* .F90.o  := $(FC) -c $(FFLAGS) $(FREEFLAGS) $(INCS) $(INCLDIR)"
	@echo "* INCS  := $(INCS)"

#-------------------------------------------------------------------------------
# build rules: MACFILE, cmd-line, or env vars must provide the needed macros
#-------------------------------------------------------------------------------

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(RANLIB_FLAGS) $(LIB)

$(EXEC): $(OBJS)
	$(LD) -o $(EXEC) $(LDFLAGS) $(OBJS) $(ULIBS) $(SLIBS)

mostlyclean:
	$(RM) -f *.f *.f90 

clean:
	$(RM) -f *.f *.f90 *.d *.mod *.o $(LIB) $(DEPGEN)
#	$(RM) -f *.f *.f90 *.d *.$(MOD_SUFFIX) $(OBJS) 

distclean:
	$(RM) -f *.f *.f90 *.d *.mod *.o $(LIB) $(DEPGEN)

#realclean:
#	$(RM) -f *.f *.f90 *.d *.$(MOD_SUFFIX) $(OBJS) 

#-------------------------------------------------------------------------------
# Build & include dependency files
#-------------------------------------------------------------------------------
# ASSUMPTIONS:
# o an externally provided dependency generator, $(DEPGEN), is available, 
#   its cmd line syntax is compatible with the build rules below.   Eg, for 
#   each .o file, there is a corresponding .d (dependency) file, and both
#   will be dependent on the same src file, eg.    foo.o foo.d : foo.F90
#   Also, the dependancy genorator's capabilities, limitations, and assumptions
#   are understood & accepted.
#-------------------------------------------------------------------------------

%.d : %.c
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $< | head -3  > $@
%.d : %.F
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@
%.d : %.F90
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@
%.d : %.H
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@

$(DEPS): $(DEPGEN)
$(DEPGEN): bld/makdep.c
	@ $(CC) -o $@ $<

# the if-tests prevent DEPS files from being created when they're not needed
ifneq ($(MAKECMDGOALS), db_files)
ifneq ($(MAKECMDGOALS), db_flags)
ifneq ($(MAKECMDGOALS), mostlyclean)
ifneq ($(MAKECMDGOALS), clean)
ifneq ($(MAKECMDGOALS), help)
    -include $(DEPS)
endif
endif
endif
endif
endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
