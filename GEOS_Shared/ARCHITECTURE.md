# GEOS_Shared: architecture

**Audience:** coding agents and developers working on the GEOS_Shared library (shared physics utilities) in GMAO_Shared.
**Describes:** the MAPL2 line, branch `main`.
**Verified against:** `GMAO_Shared` @ `b0ebdca` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`; this directory is `GEOS_Shared/`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Library-level document: purpose, contents and users. Read the source for routine detail.

## 1. Purpose

Small Fortran library of utilities shared by the GEOS gridded components:
saturation vapor pressure / humidity functions, tridiagonal solvers, surface
drag coefficients, the atmosphere-ocean interface layer, radar reflectivity,
vertical remapping and a few small tools. Not a gridded component.

## 2. Build

`GEOS_Shared/CMakeLists.txt`: one library target (`esma_set_this()`, name
`GEOS_Shared`) via `esma_add_library`, dependencies `MAPL`, `ESMF::ESMF`,
`NetCDF::NetCDF_Fortran`. No executables, no precision variants.

Compiled sources: `GEOS_Utilities.F90`, `DragCoefficients.F90`, `GEOS.F90`,
`tropovars.F90`, `surfacelayer.F90`, `GEOS_TopoGet.F90`, `GEOS_GmapMod.F90`,
`GEOS_RemapMod.F90`, `G3_MPI_Util_Mod.F90`, `windfix.F90`, `PWSSSP.F`,
`sphere.F`, `calcdbz.F`, `compute_radar_reflectivity.F90`, `OVP.F90`,
`getco2.F90`, `atmOceanIntLayer.F90`. `sphere.F` gets `${FREAL8}`; `PWSSSP.F`
gets `${FREAL8}` plus `${EXTENDED_SOURCE}` when set.

Present but not in `srcs`: `condense.F90` (external subroutine `condensation`),
`tstqsat.F90` (a test
program using `GEOS_Mod`), `qsat.h`, and `qsat_new.code`, `qsat_fast.code`,
`dqsat.code`, `dqsat_new.code`, `dqsat_fast.code`, `qsat.code`.

## 3. Contents

| Module / file | Role |
|---|---|
| `GEOS_UtilsMod` (`GEOS_Utilities.F90`) | Computational utilities: generic interfaces `QSAT`, `QSATLQU`, `QSATICE`, `DQSAT` (public as `GEOS_DQsat`; rank 0-3 variants, as for `QSAT*`), `GEOS_TRILU`, `GEOS_TRISOLVE` (rank 1-3 variants only), plus other helpers |
| `GEOS_Mod` (`GEOS.F90`) | Convenience module: `use MAPL` and `use GEOS_UtilsMod` only |
| `DragCoefficientsMod` | Surface drag/transfer coefficients (momentum, heat, scalar) |
| `sfclayer` (`surfacelayer.F90`) | Surface-layer routines; uses `DragCoefficientsMod` |
| `atmOcnIntlayer` (`atmOceanIntLayer.F90`) | Atmosphere-ocean interface layer; uses `GEOS_UtilsMod` |
| `GEOS_RadarMod` (`compute_radar_reflectivity.F90`) | Radar reflectivity; `calcdbz.F` is a companion file |
| `GEOS_RemapMod`, `GEOS_GmapMod` | Vertical remapping; `GEOS_RemapMod` uses `GEOS_GmapMod` |
| `GEOS_TopoGetMod` | Topography reading |
| `G3_MPI_Util_Mod` | MPI helpers; used by `windfix.F90` |
| `OVP` | Exports for local overpass times, accumulated over 24 h; used by chemistry |
| non-module files | `tropovars.F90` (tropopause variables, external subroutine), `getco2.F90` (`GetCO2`, RCP45 CO2 by year/day), `windfix.F90`, `PWSSSP.F`, `sphere.F` |

## 4. Include (`.code`) files

Verified from source with grep for `#include`:
- `GEOS_Utilities.F90` includes `qsatlqu.code` and `qsatice.code` (in the
  bodies of the `QSATLQU0-3`/`QSATICE0-3` procedures) and `trilu.code` and
  `trisolve.code` (in the `GEOS_TRILU1-3`/`GEOS_TRISOLVE1-3` procedures).
  The comment block at its top lists the same names plus `esatlqu.code` and
  `esatice.code`.
- `qsatlqu.code` includes `esatlqu.code`; `qsatice.code` includes
  `esatice.code` (twice each).
- `DragCoefficients.F90` includes `getcds.code`, `getcdm.code`, `getcdh.code`
  (each twice); `getcdh.code` includes `unused_dummy.H`.
- No Fortran source in this directory includes `qsat.h`, `qsat.code`,
  `qsat_new/fast.code` or `dqsat*.code`.

## 5. Who uses it

Counts are `use <module>` files under `src` (GEOSgcm fixture), by grep.

- `GEOS_UtilsMod` (saturation functions, tridiagonal solvers): Moist 13 files, Surface 13, Turbulence 2, Radiation 2
  (Irrad, Satsim), Chemistry 1.
  `GEOSturbulence_GridComp/LockEntrain.F90`, `edmf.F90`;
  `GEOSmoist_GridComp/GEOS_GF_InterfaceMod.F90` (13 Moist files);
  13 Surface files; Radiation `GEOS_IrradGridComp.F90`,
  `GEOS_SatsimGridComp.F90`; Chemistry
  `Chem_Shared/Lightning_mod.F90`; `GEOS_GcmGridComp.F90`,
  `GEOSmkiau_GridComp/GEOS_mkiauGridComp.F90`; also GMAO_hermes
  (`lcv2prs.F90`, `dyndiff.f90`).
- `GEOS_Mod`: 13 files, e.g. `GEOSdatmodyn_GridComp/NeuralNet.F90`.
- `DragCoefficientsMod` (10 files) and `sfclayer` (9): Surface children,
  e.g. `GEOSsurface_GridComp/GEOSsaltwater_GridComp/GEOS_OpenWaterGridComp.F90`,
  `GEOSlandice_GridComp/GEOS_LandIceGridComp.F90`.
- `atmOcnIntlayer` (3): `GEOS_OpenWaterGridComp.F90`,
  `GEOS_SaltWaterGridComp.F90`, `GEOS_SimpleSeaiceGridComp.F90`.
- `GEOS_RadarMod` (2): Moist `GEOS_GFDL_1M_InterfaceMod.F90`,
  `Process_Library.F90`.
- `GEOS_RemapMod`, `GEOS_TopoGetMod`: `GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`.
  `GEOS_GmapMod`: also `GEOS_mkiauGridComp.F90`.
- `OVP` (4): Chemistry `GEOS_ChemEnvGridComp.F90`, `@TR/TR_GridCompMod.F90`,
  `@GMI/GMIchem_GridCompMod.F90`, `GmiEmiss_GridCompClassMod.F90`.
- `G3_MPI_Util_Mod`: only `windfix.F90` in this library.
- Link dependency `GEOS_Shared` in CMake: `FVdycoreCubed_GridComp`,
  `GEOSsuperdyn`, `GEOSdatmodyn`, `GEOSturbulence`, `GEOSmkiau`, Chemistry
  `@TR`, `@GMI`, `Chem_Shared`, Radiation (`GEOSirrad`, `GEOSsolar`,
  `GEOSsatsim`, RRTMG module dirs), `GMAO_hermes`, `@GEOS_Util/post`.
  FVdycoreCubed has no `use` of `GEOS_UtilsMod`; it only links.
- Not checked for users: `tropovars`,
  `getco2`, `PWSSSP`, `sphere`. *Unverified:* use by GEOSadas.

## 6. Known traps

- `GEOS_Utilities.F90` declares `module GEOS_UtilsMod`: file, library and
  module names differ (also `GEOS.F90` -> `GEOS_Mod`, `surfacelayer.F90` ->
  `sfclayer`, `atmOceanIntLayer.F90` -> `atmOcnIntlayer`).
- `.code` files are textually included into procedure bodies; they are not
  standalone sources and several variants (`_new`, `_fast`) are not built.
- Several procedures are external subroutines (`tropovars`, `windfix`),
  not in a module.

## 7. Observations (noted, not fixed)

- `condense.F90`, `qsat.h`, `tstqsat.F90` and the `*_new`/`*_fast`/`dqsat*`
  `.code` files are not referenced by the build or by any include found.
- `GEOS_DQsat` is public in `GEOS_Utilities.F90` (`DQSAT0-3` procedures,
  not built from `.code` includes) and is in use (`LockEntrain.F90`:
  `use GEOS_UtilsMod, only: DQSAT=>GEOS_DQSAT`); only the `dqsat*.code`
  files are unused.
