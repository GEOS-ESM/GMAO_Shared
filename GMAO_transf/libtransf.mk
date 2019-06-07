#________________________________________
# This make module contains targets for building libtransf.a

#all: libtransf.a

LIBSRCS_F90	=	\
	m_fvGridHeader.F90		\
	m_fvGridThickness.F90		\
	m_fvGriddedState.F90		\
	m_fvGrid.F90			\
	m_fvSurface.F90			\
	m_checksums.F90			\
	m_daInterp.F90			\
	m_fgInterp.F90			\
	m_ggGradient.F90		\
	m_ggGradientSP.F90		\
	m_llInterp.F90			\
	m_ncepphis.F90			\
	m_ncSurface.F90			\
	m_parGrADS.F90			\
	m_ppInterp.F90			\
	m_interleavedObject.F90		\
	m_Interleaved.F90		\
	m_InterleavedComm.F90		\
	m_InterleaveScattererComm.F90	\
	m_InterleaveScatterer.F90	\
	m_SubdomainDistributorComm.F90	\
	m_SubdomainDistributor.F90	\
	m_swapij.F90			\
	m_undef.F90			\
	m_geosapi.F90

LIBSRCS_F	=	\
	fft.F		\
	interp_h.F	\
	mappm.F		\
	getcon.F

LIBSRCS_f	=	\
	atod.f		\
	dtoa.f		\
	elliptic.f	\
	gauss_lat_nmc.f	\
	get_ncep_phis.f

LIBSRCS	= $(LIBSRCS_F90)	\
	  $(LIBSRCS_F)		\
	  $(LIBSRCS_f)

LIBOBJS	= $(LIBSRCS_F90:.F90=.o)	\
	  $(LIBSRCS_F:.F=.o)		\
	  $(LIBSRCS_f:.f=.o)

libtransf.a: $(LIBOBJS)
	$(RM) $@
	$(AR) $@ $(LIBOBJS)

list:
	@ echo "-- sources of libtransf.a --"
	@ for f in $(LIBSRCS); do echo $$f; done

# --- Local module dependencies ---

m_checksums.o: m_undef.o
m_daInterp.o:
m_fgInterp.o: m_daInterp.o m_interleavedObject.o m_llInterp.o m_ncepphis.o m_ppInterp.o m_SubdomainDistributorComm.o m_SubdomainDistributor.o m_swapij.o mytrace.H
m_fvGridThickness.o: m_fvGrid.o m_checksums.o
m_fvGridHeader.o: m_checksums.o m_fvGrid.o m_fvGridThickness.o mytrace.H
m_fvGriddedState.o: m_checksums.o m_fvGrid.o m_fvGridThickness.o m_interleavedObject.o m_InterleaveScattererComm.o m_InterleaveScatterer.o mytrace.H m_daInterp.o
m_fvGrid.o:
m_geosapi.o:
m_ggGradientSP.o: m_geosapi.o
m_ggGradient.o: m_SubdomainDistributor.o m_SubdomainDistributorComm.o m_ggGradientSP.o m_interleavedObject.o m_swapij.o mytrace.H
m_interleavedObject.o: mytrace.H
m_InterleavedComm.o: m_interleavedObject.o mytrace.H
m_Interleaved.o: mytrace.H
m_InterleaveScattererComm.o: m_InterleaveScatterer.o
m_InterleaveScatterer.o: m_interleavedObject.o mytrace.H
m_llInterp.o: m_undef.o
m_ncepphis.o: m_interleavedObject.o
m_parGrADS.o: m_SubdomainDistributor.o m_SubdomainDistributorComm.o m_interleavedObject.o m_InterleavedComm.o m_swapij.o
m_ppInterp.o: m_checksums.o m_geosapi.o m_interleavedObject.o m_SubdomainDistributorComm.o m_SubdomainDistributor.o mytrace.H
m_SubdomainDistributorComm.o: m_interleavedObject.o m_SubdomainDistributor.o mytrace.H
m_SubdomainDistributor.o: m_interleavedObject.o mytrace.H
m_swapij.o:
m_fvSurface.o: m_InterleaveScattererComm.o m_InterleaveScatterer.o m_undef.o m_checksums.o
m_ncSurface.o: m_InterleaveScattererComm.o m_InterleaveScatterer.o m_checksums.o
m_undef.o: m_geosapi.o

interp_h.o:
dtoa.o:
atod.o:
mappm.o:
gauss_lat_nmc.o:
get_ncep_phis.o:
getcon.o:

u_intp.o: m_agInterp.o m_checksums.o
u_intp: u_intp.o m_agInterp.o libtransf.a
	$(FC) -o $@ u_intp.o m_agInterp.o -L. -ltransf $(_L)
