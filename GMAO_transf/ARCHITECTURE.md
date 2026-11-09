# GMAO_transf: architecture

**Audience:** coding agents and developers working on the GMAO_transf library (grid transforms) in GMAO_Shared.
**Describes:** the MAPL2 line, branch `main`.
**Verified against:** `GMAO_Shared` @ `b0ebdca` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`; this directory is `GMAO_transf/`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Library-level document: purpose, contents and users. Read the source for routine detail.

## 1. Purpose

Grid transforms and interpolation for GEOS data assimilation and
pre/post-processing: FV (cubed/lat-lon dynamics) grids to and from
lat-lon, Gaussian, spectral (NCEP/ECMWF) grids, horizontal and vertical
remapping, surface-field regridding, mass adjustment, and parallel
data-layout helpers. It also builds a few stand-alone conversion programs.

## 2. Build

`CMakeLists.txt` (single file, `esma_set_this()` then):

- Library `GMAO_transf`: one live (the second, later) `esma_add_library(${this} SRCS ${srcs}
  DEPENDENCIES GMAO_mpeu NCEP_sp_r4i4 NCEP_bufr_r4i4 NCEP_w3_r8i4
  NCEP_sfcio GMAO_hermes NCEP_sigio ${MKL_LIBRARIES})` . The
  earlier `esma_add_library` (near the top) is commented out (part of a note on
  r4i4 vs r8i4 NCEP libraries); the commented-out form lacks NCEP_bufr.
- Sources (`srcs`, 55 files) are listed explicitly. Every file is compiled
  with the `${FREAL8}` flags (real*8) and `-Dopenmp`; module directory is
  `include_${this}`. Fixed-form `.f`/`.F` and free-form `.f90`/`.F90` mix.
- Executables (`ecbuild_add_executable`): `dyn_divo.x`, `ec2fv.x`, `fv2ss.x`,
  `g4tog5.x`, `ss2fv.x` link `GMAO_transf GMAO_eu GMAO_hermes NCEP_w3_r8i4
  NCEP_sp_r8i4` with `${FREAL8}`; `ss2gg.x` links only NCEP r4i4 libraries
  (no `GMAO_transf`), marked "needs to be worked on". `ec2fv.csh` is
  installed to `bin`.
- Precision: no `_r4`/`_r8` library variants; one real*8 library.
- Present but not compiled by CMake: `m_cwf.F90`, `fv2gg.F`, `gg2fv.F`,
  `geos1fft.f`, `make_date.f90`, `sfc_regrid.f90`, `ssi_surf.f90`,
  `ssinfo.f90`, `ssstats.f90`, `ut_massadj.f90`, `mytrace.H`, `physcons.h`,
  `timer.h`. Legacy GNU-make files: `GMAO_transf_arch.mk`, `libtransf.mk`,
  `Externals.depend` (`m_ioutil`, `m_die`, `m_stdio` need mpeu; `m_dyn`,
  `m_dynp` need hermes). `AnIntro` is an introduction text.

## 3. Contents

Module names are those declared in the files (file name differs in some
cases).

| Group | Module (file) | Role (from names/headers) |
| --- | --- | --- |
| FV <-> other grids | `m_fv2fv` (`fv2fv.F90`), `fv2ll.F`, `ll2fv.F`, `m_ec2fv`, `m_gg`, `m_ggGradient`, `m_ggGradientSP`, `m_ss`, `m_specgrid`, `m_xform`, `m_simpleGrid` | Transforms between FV, lat-lon, Gaussian, spectral grids |
| Interpolation | `m_hinterp`, `m_llInterp`, `m_ppInterp`, `m_fgInterp`, `m_daInterp`, `interp_h.F`, `set_regrid.f90`, `m_compact_diff`, `m_elliptic` | Horizontal/vertical interpolators; first-guess and analysis interpolators |
| FV state | `m_fvGrid`, `m_fvGridHeader`, `m_fvGridThickness`, `m_fvGriddedState`, `m_gsiGriddedState`, `m_massadj`, `m_geosapi` | FV grid/state containers, GSI-gridded state, mass adjustment, API to the GEOS dynamics state |
| Surface | `m_fvSurface`, `m_ncSurface`, `m_NCEPphis` (`m_ncepphis.F90`), `make_surf.f90`, `get_ncep_phis.f`, `get_ncep_tsx.f` | Surface-field handling and NCEP surface/topography input |
| Parallel layout | `m_Interleaved`, `m_InterleavedComm`, `m_interleavedObject`, `m_InterleaveScatterer`, `m_InterleaveScattererComm`, `m_SubdomainDistributor`, `m_SubdomainDistributorComm`, `m_parGrADS`, `m_swapij`, `dynamics_lattice_module` | Distribution/scatter of fields across PEs |
| Utilities | `m_checksums`, `m_undef`, `util.f90`, `mpi_util.F`, `mymalloc.F`, `timer.F`, `tick.f`, `qsat.f` | Misc helpers |
| F77 kernels | `fft.F`, `gmap.F`, `mappm.F`, `atod.f`, `dtoa.F`, `gauss_lat_nmc.f`, `blend_uvsh.F90` | FFT, vertical remap (PPM), Gaussian latitudes, conversions |
| Programs | `dyn_divo.F90`, `ec2fv.F90` (+`ec2fv.csh`), `fv2ss.f90`, `ss2fv.f90`, `g4tog5.f90`, `ss2gg.f90` | Stand-alone converters (ECMWF->FV, FV<->spectral, GEOS-4 -> GEOS-5, spectral -> Gaussian, divergence/vorticity from dynamics file) |

## 4. Who uses it

- `GMAO_stoch`: its `CMakeLists.txt` lists `GMAO_transf` in `DEPENDENCIES`
  and adds `include_GMAO_transf`, but `stoch_module.F90` `use`s none of its
  modules (it uses ESMF, MAPL, mod_param, stoch_data, initialize, m_die).
- `@GMAO_Shared/CMakeLists.txt` lists it as a subdirectory.
- `GMAO_hermes`: no `use` of transf modules found.
- `GEOS_Util` (nested repo): no `GMAO_transf` link in its CMake files. Its
  `post/` programs compile their own `dynamics_lattice_module.F`,
  `mpi_util.F`, `mymalloc.F`, `timer.F` (`post/CMakeLists.txt` SRCS).
- MOM (`@mom`): no reference to `GMAO_transf`; its `gotm-4.0` has its own
  `util` module (see section 5).
- No other user found in the GEOSgcm fixture. GEOSadas (GSI, mkiau) may use
  it: *Unverified.* The CMake note mentions `GEOSmkiau_GridComp` using
  NCEP_sp r4i4.

## 5. Known traps

- Module names differ from file names: `m_fv2fv` in `fv2fv.F90`, `m_NCEPphis`
  in `m_ncepphis.F90`; no `m_fv2ll`/`m_ll2fv` modules exist (`fv2ll.F`,
  `ll2fv.F` are plain routine files).
- `GMAO_eu` is linked by the executables, while the library depends on
  `GMAO_mpeu`.
- Mixed NCEP precision: library uses `NCEP_sp_r4i4` with `NCEP_w3_r8i4`; the
  executables use `NCEP_sp_r8i4`; `ss2gg.x` uses r4i4 only.
- Duplicate module names: `dynamics_lattice_module` exists here and as
  `@GEOS_Util/post/dynamics_lattice_module.F`; `util` exists here
  (`util.f90`, used by 11 GMAO_transf files) and in MOM
  (`gotm-4.0/util/util.F90`, used by gotm-4.0 `turbulence` and `util`). The
  GEOS_Util and MOM users compile their own copies, not GMAO_transf's.
- Source files set `-Dopenmp` globally; `${FREAL8}` promotes reals in all
  library sources.

## 6. Observations (noted, not fixed)

- Two `esma_add_library` lines exist; only the later one is active.
  Comments describe a conflict between mkiau (r4i4 needed) and the
  executables here (r8 needed).
- Roughly a dozen source files in the directory are not in `srcs` (section 2).
- `m_fv2ll`, `m_ll2fv` are not modules; `fv2ll.F`, `ll2fv.F` are plain
  routine files. Only `m_fv2fv` and `m_NCEPphis` differ from their file names.
