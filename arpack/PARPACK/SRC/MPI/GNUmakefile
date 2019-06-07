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
       ESMADIR := $(PWD)/../../../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------


THIS = parpack
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
pcgetv0.f  pcnaupd.f  pdlamch.f  pdnaupd.f  pdsaitr.f  pdseupd.f   pslamch.f  psnaupd.f  pssaitr.f  psseupd.f  pznapps.f  pzngets.f\
pclarnv.f  pcneigh.f  pdlarnv.f  pdneigh.f  pdsapps.f  pdsgets.f   pslarnv.f  psneigh.f  pssapps.f  pssgets.f  pznaup2.f\
pcnaitr.f  pcneupd.f  pdnaitr.f  pdneupd.f  pdsaup2.f  pdznorm2.f  psnaitr.f  psneupd.f  pssaup2.f  pzgetv0.f  pznaupd.f\
pcnapps.f  pcngets.f  pdnapps.f  pdngets.f  pdsaupd.f  pscnorm2.f  psnapps.f  psngets.f  pssaupd.f  pzlarnv.f  pzneigh.f\
pcnaup2.f  pdgetv0.f  pdnaup2.f  pdnorm2.f  pdseigt.f  psgetv0.f   psnaup2.f  psnorm2.f  psseigt.f  pznaitr.f  pzneupd.f
 
OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

FREAL = 
FOPT = $(FOPT2)
FDP  =

USER_FFLAGS = $(BIG_ENDIAN) -I$(INC_MPI)
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

vpath % $(MOD_DIRS)

$(LIB) lib : $(DEPS) $(OBJS)
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
