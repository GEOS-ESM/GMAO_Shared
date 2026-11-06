# LANL_Shared: architecture

**Audience:** coding agents and developers working on LANL_Shared (the CICE4 sea-ice library) in GMAO_Shared.
**Describes:** the MAPL2 line, branch `main`.
**Verified against:** `GMAO_Shared` @ `b0ebdca` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`; this directory is `LANL_Shared/`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Library-level document: purpose, contents and users. Read the source for routine detail.

## 1. Purpose
`LANL_Shared` holds one vendored third-party model: CICE version 4.0, the Los
Alamos (LANL) sea-ice model (`CICE4/README_v4.0`, dated August 2008). GEOS builds
a subset of it as the shared library `CICE4`, which the sea-ice components and
the saltwater surface component link. It is *not* CICE6; CICE6 is a separate
library used by `CICE_GEOSPlug` (see Known traps).

## 2. Build
- `LANL_Shared/CMakeLists.txt`: only `add_subdirectory (CICE4)`.
- `CICE4/CMakeLists.txt`: `esma_set_this()` is called, then `get_filename_component(this ... NAME)` resets `this` to the directory name, so the target is `CICE4`;
  `esma_add_library(... DEPENDENCIES MAPL TYPE SHARED)`.
- Sources (explicit list, 43 files): `csm_share/shr_orb_mod.F90`,
  `drivers/cice4/ice_constants.F90`, `drivers/cice4/ice_prescribed_mod.F90`,
  `mpi/*.F90` (7 files, MPI variant of the communication layer) and
  `source/*.F90` (physics and infrastructure).
- Compile definitions (PRIVATE): `coupled GEOS MODIFY_TOPOGRAPHY USE_R8`; `LINUX`
  on Linux. Option `BUILD_MIT_OCEAN` (default OFF) adds a `BUILD_MIT_OCEAN`
  definition.
- Intel Release builds reset `CMAKE_Fortran_FLAGS_RELEASE` to
  `GEOS_Fortran_FLAGS_NOVECT` (comment: vectorised Release flags segfault CICE4).
- Executable: `makdep` from `bld/makdep.c` (dependency generator of the original
  CICE build).
- Present but not compiled: `drivers/esmf/` (including its own `ice_constants.F90`
  and `CICE_*Mod.F90`), the rest of `drivers/cice4/` (`CICE.F90`,
  `CICE_InitMod/RunMod/FinalMod.F90`, `CICE_RunMod.F90_debug`), `serial/`,
  `cice/` (`comp_ice`, `clean_ice`, `README_v4.0`, an `ice.log.Linux.LANL.coyote`
  run log), `doc/cicedoc.pdf` and top-level copies `cicedoc.pdf`, `comp_ice`,
  `clean_ice`, `ice.log.Linux.LANL.coyote`.

## 3. Contents
Roles are taken from directory and file names.

| Directory / modules | Role |
|---|---|
| `source/` | CICE4 physics and infrastructure, `ice_*` modules: grid and domain (`ice_grid`, `ice_domain`, `ice_blocks`, `ice_distribution`, `ice_spacecurve`), state and fluxes (`ice_state`, `ice_flux`), dynamics (`ice_dyn_evp`, `ice_mechred`, `ice_transport_driver`, `ice_transport_remap`), thermodynamics (`ice_therm_vertical`, `ice_therm_itd`, `ice_itd`, `ice_meltpond`, `ice_shortwave`, `ice_atmo`, `ice_ocean`), I/O and control (`ice_init`, `ice_restart`, `ice_history`, `ice_read_write`, `ice_forcing`, `ice_calendar`, `ice_step_mod`, ...) |
| `mpi/` | `ice_boundary`, `ice_broadcast`, `ice_communicate`, `ice_exit`, `ice_gather_scatter`, `ice_global_reductions`, `ice_timers` |
| `drivers/cice4/` (compiled part) | `ice_constants`, `ice_prescribed_mod` (prescribed-ice and nudging) |
| `csm_share/` | `shr_orb_mod` (orbital parameters) |
| `bld/makdep.c` | `makdep` tool |

## 4. Who uses it
Link target `CICE4` (CMake `DEPENDENCIES`), all under `GEOSgcm_GridComp`:
- `GEOSogcm_GridComp/CMakeLists.txt` (line 10) and
  `GEOSseaice_GridComp/CMakeLists.txt` (line 9, which also lists subcomponents
  `GEOSdataseaice_GridComp GEOSCICEDyna_GridComp CICE_GEOSPlug`).
- `GEOSdataseaice_GridComp`, `GEOSCICEDyna_GridComp`, `GEOSMITDyna_GridComp`
  CMakeLists.
- `GEOSsurface_GridComp/GEOSsaltwater_GridComp/CMakeLists.txt` (line 11).

Fortran `use` citations (verified):
- `GEOSCICEDyna_GridComp/GEOS_CICEDynaGridComp.F90`: `ice_constants`, `ice_state`,
  `ice_communicate`, `ice_domain`, `ice_grid`, `ice_dyn_evp`, `ice_step_mod`,
  `ice_itd`, `ice_flux`, `ice_init`, ... (init and dynamics entry points).
- `GEOSdataseaice_GridComp/GEOS_DataSeaIceGridComp.F90`: `ice_state`,
  `ice_prescribed_mod`; `GEOSseaice_GridComp/GEOS_SeaIceGridComp.F90`:
  `ice_prescribed_mod` (`ice_nudging`).
- `GEOSMITDyna_GridComp/GEOS_MITDynaGridComp.F90`: `ice_domain_size`, `ice_init`,
  `ice_work`.
- `GEOSogcm_GridComp/GEOS_OgcmGridComp.F90` (line 1324): `ice_constants`.
- `GEOSsaltwater_GridComp/GEOS_CICE4ColumnPhysGridComp.F90`: `ice_kinds_mod`.

## 5. Known traps
- **Duplicate module names with CICE6.** `CICE_GEOSPlug/CMakeLists.txt` builds
  its own target `cice6` (from `cice6_cmake/`, sources taken from the separate
  `@cice6` tree) and depends on `MAPL cice6`, not on `CICE4`. CICE6 has its own
  `ice_*` modules (and `CICE_InitMod`, `CICE_RunMod`, `CICE_FinalMod`, used by
  `CICE_GEOSPlug.F90`), so name-based greps for `ice_*`/`CICE_*Mod` return about
  about 120 hits in `@cice6` and about 25 under `@sis2` (the latter are SIS2 sources, not
  CICE4 users). Those are not LANL_Shared users. The `ice_*` module names
  collide, so *Inferred:* the two libraries cannot be assumed to coexist in one link line.
  *Inferred:* `GEOSseaice_GridComp` lists both `CICE4` (dependency) and
  `CICE_GEOSPlug` (subcomponent, linking `cice6`); whether both are in the same
  executable in a given configuration was not checked.
- `ice_constants.F90` exists in both `drivers/cice4/` (compiled) and
  `drivers/esmf/` (not compiled); only the `drivers/cice4` copy is built.
- Compiled sources come from an explicit list: adding a file to a directory
  does not build it.
- Not compiled: `drivers/cice4/CICE.F90`, `CICE_*Mod.F90`; *Inferred:* the GEOS
  components call CICE4 routines directly instead of running the CICE4 driver.

## 6. Observations (noted, not fixed)
- Duplicate copies of docs, `comp_ice`, `clean_ice` and the `coyote` run log sit
  both at `CICE4/` and in `CICE4/cice/`.
- The target and the directory (`CICE4`) differ from the repository-level name
  (`LANL_Shared`); the `LANL_Shared` name appears only in the parent
  `GMAO_Shared/CMakeLists.txt` subdirectory list.
- `GEOSseaice_GridComp` depends on `CICE4` while its subcomponent
  `CICE_GEOSPlug` depends on `cice6`; how the shared `ice_*` names are kept
  apart at link/run time is not visible from the CMake files read.
- Users in `GEOSgcm` are identified; users in GEOSadas or other fixtures:
  *Unverified:*.
