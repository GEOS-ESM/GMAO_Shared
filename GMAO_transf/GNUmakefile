#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 3mar2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva  Standardization
# 03Jun2008  Todling   - Merged Merra w/ FDDA; add AMS/Ravi changes
#                      - temporarily removed fv2ss/ss2fv (need dyn update)
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies
  include ./GMAO_transf_arch.mk           # Costumized configuration

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------


THIS = $(shell basename `pwd`)
LIB  = lib$(THIS).a

BIN_SRCS = fv2ss.f90  ss2fv.f90  ss2gg.f90  ec2fv.F90 g4tog5.f90 dyn_divo.F90 # fv2gg.F  gg2fv.F
BINS = $(addsuffix .x, $(basename $(BIN_SRCS)))


esma_install install: $(LIB) $(BINS)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.a         $(ESMALIB)
	$(CP) -p $(BINS)     $(ESMABIN)
	$(CP) -p ec2fv.csh   $(ESMABIN)
	$(CP) -p *.mod       $(ESMAINC)/$(THIS)
	chmod 755 $(ESMABIN)/ec2fv.csh

esma_clean clean:
	$(RM) *~ *.[aox] *.[Mm][Oo][Dd] m_gmap.F

esma_distclean distclean:
	$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] m_gmap.F

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"

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
	@echo "        FREAL = $(LIB_EU)"

macros: help
	@echo ""
	@echo "Additional macros:"
	@echo "           FC = $(FC)"
	@echo "      _STATIC = $(_STATIC)"

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS =	atod.f				\
	blend_uvsh.F90                  \
	dtoa.F				\
	dynamics_lattice_module.F	\
	m_elliptic.f			\
	fft.F				\
	fv2fv.F90			\
	fv2ll.F				\
	gauss_lat_nmc.f			\
	get_ncep_phis.f			\
	get_ncep_tsx.f			\
	m_ec2fv.F			\
	m_hinterp.F			\
	interp_h.F			\
	ll2fv.F				\
	m_InterleaveScatterer.F90	\
	m_InterleaveScattererComm.F90	\
	m_Interleaved.F90		\
	m_InterleavedComm.F90		\
	m_SubdomainDistributor.F90	\
	m_SubdomainDistributorComm.F90	\
	m_checksums.F90			\
	m_daInterp.F90			\
	m_fgInterp.F90			\
	m_fvGrid.F90			\
	m_fvGridHeader.F90		\
	m_fvGridThickness.F90		\
	m_fvGriddedState.F90		\
	m_fvSurface.F90			\
	m_geosapi.F90			\
	m_gg.F90			\
	m_ggGradient.F90		\
	m_ggGradientSP.F90		\
	m_gmap.F90			\
	m_gsiGriddedState.F90		\
	m_interleavedObject.F90		\
	m_llInterp.F90			\
	m_massadj.F90			\
	m_ncSurface.F90			\
	m_ncepphis.F90			\
	m_parGrADS.F90			\
	m_ppInterp.F90			\
	m_simpleGrid.F90		\
	m_ss.f90			\
	m_swapij.F90			\
	m_undef.F90			\
	m_xform.f90			\
	make_surf.f90			\
	mappm.F				\
	mpi_util.F			\
	mymalloc.F			\
	qsat.f				\
	set_regrid.f90			\
	tick.f				\
	timer.F				\
	util.f90

OBJS := $(addsuffix .o, $(basename $(SRCS)))

ALL_SRCS = $(SRCS) $(BIN_SRCS)
DEPS := $(addsuffix .d, $(basename $(ALL_SRCS)))

ifneq ("$(BOPT)","g")
   FOPT = $(FOPT3) 
endif
#FREAL = 

OMPFLAG = # disable OpenMP until it has been tested.

THIS_SP    = NCEP_sp_r8i4
THIS_W3    = NCEP_w3_r8i4

MOD_DIRS = . $(INC_MPEU) $(INC_SP) $(INC_CRTM) $(INC_IRSSE) \
             $(INC_SFCIO) $(INC_HERMES) $(INC_SIGIO)
USER_FINCS = $(I)$(INC_MPEU) 
USER_FFLAGS = $(BIG_ENDIAN) $(OMPFLAG) $(EXTEND_SOURCE)
USER_FDEFS  = $(D)openmp
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

LDFLAGS += $(OMPFLAG) 

vpath % $(MOD_DIRS)

$(LIB) lib : $(MODS) $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

%.x : $(LIB) %.o 
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB) $(LIB_SCI) $(LIB_SYS)

# Special case: ifort has a hardtime optmizing this on the Altix,
# maybe elsewhere. To be safe, lower optimization everywhere.

m_xform.o : m_xform.f90
	$(FC) -c $(patsubst $(FOPT3),$(FOPT2),$(f90FLAGS)) $<

# Specific targets

m_gmap.F : gmap.F
	$(RM) m_gmap.F
	$(LN) -s gmap.F m_gmap.F

fv2gg.x : $(LIB) m_gmap.o fv2gg.o
	$(LD) $(LDFLAGS) -o fv2gg.x fv2gg.o m_gmap.o $(LIB) \
	$(LIB_GFIO) $(LIB_SP) $(LIB_SDF) $(LIB_SCI) $(LIB_SYS)

gg2fv.x : $(LIB) m_gmap.o gg2fv.o 
	$(LD) $(LDFLAGS) -o gg2fv.x gg2fv.o m_gmap.o $(LIB) \
	$(LIB_SP) $(LIB_GFIO) $(LIB_SDF) $(LIB_SCI) $(LIB_SYS)

fv2ss.x : $(MODS) fv2ss.o $(LIB) $(LIB_SP)
	$(LD) $(LDFLAGS) -o fv2ss.x $(MODS) fv2ss.o \
	$(LIB) $(LIB_W3) $(LIB_HERMES) $(LIB_SFCIO) $(LIB_EU) \
	$(LIB_GFIO) $(LIB_SP) $(LIB_MPEU)  $(LIB_SIGIO) $(LIB_SDF) \
	$(LIB_SCI) $(LIB_SYS)

dyn_divo.x : $(MODS) dyn_divo.o $(LIB) $(LIB_SP)
	$(LD) $(LDFLAGS) -o dyn_divo.x $(MODS) dyn_divo.o \
	$(LIB) $(LIB_W3) $(LIB_HERMES) $(LIB_SFCIO) $(LIB_EU) \
	$(LIB_GFIO) $(LIB_SP) $(LIB_MPEU)  $(LIB_SIGIO) $(LIB_SDF) \
	$(LIB_SCI) $(LIB_SYS)

ec2fv.x : $(MODS) ec2fv.o $(LIB) $(LIB_SP)
	$(LD) $(LDFLAGS) -o ec2fv.x $(MODS) ec2fv.o \
	$(LIB) $(LIB_HERMES) $(LIB_EU) \
	$(LIB_GFIO) $(LIB_SDF) $(LIB_SCI) $(LIB_SYS)

ss2fv.x : $(MODS) ss2fv.o $(LIB) $(LIB_SP)
	$(LD) $(LDFLAGS) -o ss2fv.x $(MODS) ss2fv.o \
	$(LIB) $(LIB_W3) $(LIB_HERMES) $(LIB_SFCIO) $(LIB_EU) \
	$(LIB_GFIO) $(LIB_SP) $(LIB_EU) $(LIB_SIGIO) $(LIB_SDF) \
	$(LIB_SCI) $(LIB_SYS)

# BEWARE: There is much to question the link of SP lib and scslib (precision-related issues)
ss2gg.o: ss2gg.f90
	$(FC) $(patsubst $(FREAL8),,$(FFLAGS)) -c $<
ss2gg.x : ss2gg.o $(LIB_W3) $(LIB_SP) $(LIB_SIGIO)
	$(LD) $(LDFLAGS) -o ss2gg.x ss2gg.o \
	-lNCEP_w3_r4i4 -lNCEP_sp_r4i4 -lNCEP_bacio_r4i4 $(LIB_SIGIO) $(LIB_SCI) $(LIB_SYS)

ssinfo.x : $(MODS) ssinfo.o $(LIB) $(LIBsp)
	$(LD) $(LDFLAGSMP) -o ssinfo.x $(MODS) ssinfo.o \
	$(LIB) $(LIB_W3) $(LIB_HERMES) $(LIB_SIGIO) $(LIB_SFCIO) \
	$(LIB_GFIO) $(LIB_SP) $(LIB_EU) $(LIB_SDF) $(LIB_SCI) $(LIB_SYS)

ssstats.x : $(MODS) ssstats.o $(LIB)  $(LIB_SP)
	$(LD) $(LDFLAGS) -o ssstats.x $(MODS) ssstats.o \
	$(LIB) $(LIB_W3) $(LIB_HERMES) \
	$(LIB_GFIO) $(LIB_SP) $(LIB_EU) $(LIB_SDF) $(LIB_SCI) $(LIB_SYS)

surf_bkg.x: surf_bkg.o $(LIB)
	$(LD) $(LDFLAGS) -o surf_bkg.x surf_bkg.o $(LIB) \
	$(LIB_GFIO) $(LIB_SDF) $(LIB_EU) $(LIB_SFCIO) $(LIB_HERMES)

ssi_surf.x: ssi_surf.o $(LIB)
	$(LD) $(LDFLAGS) -o ssi_surf.x ssi_surf.o $(LIB) 

sfc_regrid.x: sfc_regrid.o $(OBJ)
	$(LD) $(LDFLAGS) -o sfc_regrid.x sfc_regrid.o $(LIB)

fv2ll.o: fv2ll.F
	$(FC) -c $(patsubst $(D)openmp,$(null),$(FFLAGS)) $<

ll2fv.o: ll2fv.F
	$(FC) -c $(patsubst $(D)openmp,$(null),$(FFLAGS)) $<

ut_massadj.x: m_massadj.o ut_massadj.o
	$(FC) $(LDFLAGS) -o ut_massadj.x ut_massadj.o m_massadj.o \
	      $(LIB_HERMES)  $(LIB_GFIO) $(LIB_MPEU) \
	      $(LIB_SDF) $(LIB_SCI) $(LIB_SYS)

#m_elliptic.o: m_elliptic.f
#	$(FC) -c $(_STATIC) $(patsubst $(D)openmp,$(null),$(FFLAGS)) $<
m_elliptic.o: m_elliptic.f
	$(FC) -c $(FFLAGS) $(FREAL8) -I. $(_ELLIPTIC) $<


g4tog5.x : $(LIB) g4tog5.o
	$(LD) $(LDFLAGS) -o g4tog5.x g4tog5.o $(LIB) \
	$(LIB_HERMES) $(LIB_GFIO) $(LIB_SDF) $(LIB_SYS)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
