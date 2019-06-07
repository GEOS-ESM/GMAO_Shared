#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 3mar2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva  Standardization
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

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------


THIS = arpack
LIB  = lib$(THIS).a
BIN  = 

esma_install install: $(LIB) $(BIN)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.a         $(ESMALIB)
	$(CP) -p *.h         $(ESMAINC)/$(THIS)

esma_clean clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

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
	@echo "        FREAL = $(FREAL)"

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS = \
cgetv0.f  cstatn.f            dsconv.f  dstqrb.f  snconv.f  sseigt.f sngets.f \
cnaitr.f  dnconv.f  dseigt.f  sneigh.f  ssesrt.f  zsortc.f \
cnapps.f  dneigh.f  dsesrt.f  sgetv0.f  sneupd.f  sseupd.f  zgetv0.f   zstatn.f \
cnaup2.f  dgetv0.f  dneupd.f  dseupd.f  slaqrb.f  sngets.f  ssgets.f  znaitr.f \
cnaupd.f  dlaqrb.f  dngets.f  dsgets.f  snaitr.f  ssaitr.f  ssortc.f  znapps.f \
cneigh.f  dnaitr.f  dsaitr.f  dsortc.f  snapps.f  ssapps.f  ssortr.f  znaup2.f \
cneupd.f  dnapps.f  dsapps.f  dsortr.f  snaup2.f  ssaup2.f  sstatn.f  znaupd.f \
cngets.f  dnaup2.f  dsaup2.f  dstatn.f  snaupd.f  ssaupd.f  sstats.f  zneigh.f \
csortc.f  dnaupd.f  dsaupd.f  dstats.f            ssconv.f  sstqrb.f  zneupd.f \
sspg.f    dspg.f
 
OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

FREAL = 
FOPT = $(FOPT2)

USER_FFLAGS = $(BIG_ENDIAN) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

%.x : $(LIB) %.o 
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB) $(LIB_EU) $(LIB_GFIO) \
                                    $(LIB_MFHDF3) $(LIB_SDF) $(LIB_SYS)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
