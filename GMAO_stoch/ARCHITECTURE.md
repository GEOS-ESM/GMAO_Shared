# GMAO_stoch: architecture

**Audience:** coding agents and developers working on the GMAO_stoch library (stochastic physics patterns) in GMAO_Shared.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GMAO_Shared` @ `595999b` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`; this directory is `GMAO_stoch/`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Library-level document: purpose, contents and users. Read the source for routine detail.

> **Build status on this line:** not built, not ported. On `release/MAPL-v3`, `CMakeLists.txt` (repository root) calls `esma_add_subdirectories` only for `GMAO_mpeu`, `GMAO_hermes` and `GEOS_Shared` (comment: "not yet ported to MAPL3 - only build GMAO_mpeu, GMAO_hermes (light), and GEOS_Shared for FVdycoreCubed_GridComp"); `GMAO_stoch` appears only in the commented-out `if ( NOT HERMES_LIGHT )` block. The code in this directory is the MAPL2 code, kept for a later port.

## 1. Purpose
Library that generates stochastic-physics random patterns on the sphere for
the atmospheric physics: SPPT (stochastically perturbed physics tendencies)
and SKEB (stochastic kinetic energy backscatter). Patterns are built in
spherical-harmonic space with prescribed temporal and spatial correlation
(`patterngenerator.f90`), then transformed to a Gaussian grid and
interpolated to the model grid. The code derives from the NCEP spectral
(GFS) stochastic-physics code and its spectral-transform helpers.

## 2. Build
`GMAO_stoch/CMakeLists.txt` builds one library, `GMAO_stoch`
(`esma_add_library`).
- Dependencies: `GMAO_mpeu GMAO_gfio_r8 GMAO_transf MAPL NCEP_sfcio
  NCEP_sp_r4i4 ESMF::ESMF NetCDF::NetCDF_Fortran`.
- Include directories also add `NCEP_sigio`, `GMAO_hermes`.
- Compiled with the `FREAL8` flags (real*8 default), `-O1` release flags and
  `EXTENDED_SOURCE`.
- `NCEP_sp_r4i4` is linked: a comment in the file says the "correct" r8i4
  variant would break `NCEP_sp` use in `GEOSmkiau_GridComp`, which needs
  the r4i4 interface.
- Present but not compiled: `stochprog.F90` (main program calling
  `stochmain`), `stoch_ll2cc_stubs.F90` (empty `cube2latlon`/`latlon2cube`
  stubs for standalone linking), `sppt.rc`, `stoch.rc`, `pref.txt` (sample
  resources). `stochmain.f90` (the test driver subroutine) is compiled into
  the library.

## 3. Contents
| File | Module / unit | Role |
|---|---|---|
| `stoch_module.F90` | `stoch_module` | Entry point: `setup_pattern`, `clear_pattern`, `sppt_pattern`, `skeb_pattern`, `writout8/4`; uses ESMF/MAPL |
| `stoch_data.f90` | `stoch_data` | Sets up and initializes the random patterns and their state |
| `get_pattern.f90` | `get_pattern` | `get_pattern_sppt/skeb`; scalar/vort-div spectral <-> grid transforms |
| `patterngenerator.f90` | `patterngenerator` | Random patterns with given time/space autocorrelation in spectral space |
| `mersenne_twister.f` | `mersenne_twister` | Random number generator (`random_setseed`, `random_gauss`, `random_stat`) |
| `mod_param.F90` | `mod_param` | Kinds, constants, truncation and namelist parameters (`sppt*`, `skeb*`); includes `mpif.h` |
| `initialize.F90` | `initialize` | Latitude/node arrays and spherical-harmonic setup arrays |
| `mod_four.F90` | `mod_four` | Fourier transform routines |
| `mod_uvdz.F90` | `mod_uvdz` | Wind <-> vorticity/divergence spectral conversion |
| `grd_xform.F90` | `grd_xform` | Grid transform helpers (uses `m_const`) |
| `num_parthd.f90` | `num_parthd` | `num_parthds` returns a hard-coded 8 (OpenMP query commented out) |
| `epslon.f90`, `gauss2grd.f90`, `get_ls_node.f`, `glats.f`, `gozrineo.f`, `pln2eo.f`, `setlats.f`, `wrap_fftpack.f`, `mpi_more.f` | external procedures | Spectral-transform helpers (Legendre, Gaussian latitudes, FFT wrapper, MPI helpers) |
| `function_indlsev.h`, `function_indlsod.h` | include files | Index statement functions |
| `stoch.rc` | namelists `&spptparm`, `&skebparm` | Sample amplitudes, time scales, length scales, seeds |
| `sppt.rc` | resource file | Sample driver config (grid, `SPPT_RUN_DT`) |

## 4. Who uses it
- `GEOSphysics_GridComp` links it: `GEOSagcm_GridComp/GEOSphysics_GridComp/CMakeLists.txt`
  (`DEPENDENCIES ... GMAO_stoch`) and `GEOS_PhysicsGridComp.F90` does
  `use stoch_module`. The component reads `SPPT:` (and a SKEB flag) from its
  config and updates the SPPT/SKEB patterns and perturbed tendencies when enabled
  (`DO_SPPT`, `DO_SKEB`; exports such as `DUDTSTOCH`).
- `Shared/@GMAO_Shared/CMakeLists.txt` lists `GMAO_stoch` as a subdirectory.
- No other user found in the GEOSgcm fixture. *Unverified:* use from GEOSadas.

## 5. Known traps
- **Duplicate module name `mersenne_twister`.** The MOM6 plug's
  `@stochastic_physics` (`MOM6_GEOSPlug/@stochastic_physics`) has its own
  `mersenne_twister.F90` (listed in its `CMakeLists.txt`) and its own
  `stochy_*` modules. It does not use any GMAO_stoch module: its sources
  `use` only `stochy_*`, `spectral_transforms`, `kinddef`, `mpi_wrapper` and
  `mersenne_twister`. *Inferred:* it is a separate UFS/CCPP-derived package
  and a different code base, despite similar function and the shared
  `mersenne_twister` module name (the two are not linked into one target in
  the fixture).
- A third `mersenne_twister.f` is built by `@NCEP_Shared/NCEP_w3`
  (`NCEP_w3/CMakeLists.txt`), so three libraries in the fixture provide a
  file of that name.
- *Inferred:* source is mixed fixed/free form and `.f`/`.F90`;
  `EXTENDED_SOURCE` is applied to the whole target.
- Libraries are r8 (`FREAL8`) but `NCEP_sp_r4i4` is linked (see Build).

## 6. Observations (noted, not fixed)
- `num_parthds` returns a constant 8 irrespective of the actual thread count.
- `stochprog.F90` and `stoch_ll2cc_stubs.F90` are not in the source list.
- `stochmain.f90` reads `teststoch.rc`, while the sample resource files are
  named `sppt.rc` and `stoch.rc`.
- `CMakeLists.txt` includes `NCEP_sigio` and `GMAO_hermes` include
  directories but does not list them as dependencies.
