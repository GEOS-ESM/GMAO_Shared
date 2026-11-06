# GMAO_gfio: architecture

**Audience:** coding agents and developers working on the GMAO_gfio library (gridded file I/O) in GMAO_Shared.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GMAO_Shared` @ `595999b` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`; this directory is `GMAO_gfio/`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Library-level document: purpose, contents and users. Read the source for routine detail.

> **Build status on this line:** not built, not ported. On `release/MAPL-v3`, `CMakeLists.txt` (repository root) calls `esma_add_subdirectories` only for `GMAO_mpeu`, `GMAO_hermes` and `GEOS_Shared` (comment: "not yet ported to MAPL3 - only build GMAO_mpeu, GMAO_hermes (light), and GEOS_Shared for FVdycoreCubed_GridComp"); `GMAO_gfio` appears only in the commented-out `if ( NOT HERMES_LIGHT )` block (as the `GMAO_gfio_r4`/`GMAO_gfio_r8` special case); `HERMES_LIGHT` is ON, so the targets `GMAO_gfio_r4` and `GMAO_gfio_r8` do not exist. `GMAO_hermes` in its light form does not depend on `GMAO_gfio_r8`. The code in this directory is the MAPL2 code, kept for a later port.

## 1. Purpose
GFIO ("Grided File I/O") is the legacy Data Assimilation Office library for
reading and writing gridded, time-stamped data in self-describing files
(HDF/netCDF), via a Fortran 77 style API (`GFIO_Open`, `GFIO_GetVar`,
`GFIO_Close`, ...). `README` dates it to GFIO 1.0.8 (1999). Most of the API
is plain external routines, not Fortran modules.

## 2. Build
On this line `GMAO_Shared/CMakeLists.txt` sets `HERMES_LIGHT ON` and the block below is commented out, so neither target is built.
`GMAO_gfio/CMakeLists.txt` is added **twice** by `GMAO_Shared/CMakeLists.txt`
(the special case is commented out on this line; live on the MAPL2 line), only when `HERMES_LIGHT` is OFF (the default option value; block
`if ( NOT HERMES_LIGHT )`, which holds the commented-out `GMAO_gfio_r4`/`GMAO_gfio_r8` special case here) and the `GMAO_gfio` directory exists: as `GMAO_gfio_r4` and `GMAO_gfio_r8`, with
`add_dependencies(GMAO_gfio_r4 GMAO_gfio_r8)`. The precision is taken from
the target name.
- Compiled sources: `gfio.F`, `gfioCF.F`, `diffdate.f`, `getdate.f`,
  `timesubs.f`, `gfioutil.F`, `m_GFIO_GetFld.f`, `m_GFIO_PutFld.F`.
- Dependencies: `GMAO_mpeu`, `esmf`, `NetCDF::NetCDF_Fortran`.
- `r8`: adds the `FREAL8` flags and a `link_target` that symlinks the include
  directory `GMAO_gfio` to `GMAO_gfio_r8` (comment: "_r8 is default for gfio").
- `r4` only, if `USE_F2PY` and F2PY3 are found: Python module `GFIO_` from
  `GFIO_py.F90` (exports `gfioopen`, `gfiocreate`, `gfiogetvar`, ...),
  installed to `lib/Python`.
- Executable per precision: `GFIO_mean_${precision}.x` from `GFIO_mean.f90`.
- `gfio.py` is installed to `lib/Python`.
- Present, not in the source list: `m_gfio.f90`, `gfio_init.f`,
  `gfio_flush.F`, `gfio_test.F`, `ut_cyclic.f90`, C files (`gfio_idl.c`,
  `testcall.c`, `cfortran.h`, `gfio_c.h`), `gfio.h`, `unused_dummy.H`, legacy
  `Arch.mk`, `make.dist.csh`.

## 3. Contents
| File | Role (from header comments / names) |
|---|---|
| `gfio.F`, `gfioCF.F` | Core F77 API: open/create/inquire/get/put/close. `gfioCF.F` is the CF-convention variant (by name; *Inferred*) |
| `gfioutil.F` | Helpers used by the core routines |
| `m_GFIO_GetFld.f`, `m_GFIO_PutFld.F` | Modules `m_GFIO_GetFld`, `m_GFIO_PutFld`: higher-level get/put of whole fields |
| `diffdate.f`, `getdate.f`, `timesubs.f` | Date/time arithmetic (`nymd`/`nhms` integers) |
| `m_gfio.f90` | Module `m_gfio` "Grided File I/O Class"; not in the CMake source list |
| `gfio_init.f`, `gfio_flush.F` | Init/flush routines; not compiled |
| `GFIO_mean.f90` | Program: monthly means (and linear combinations) of files |
| `GFIO_py.F90`, `gfio.py` | Python (f2py) wrapper and script |
| `gfio_idl.c`, `cfortran.h`, `gfio_c.h`, `testcall.c` | C/IDL bindings and C test |
| `gfio_test.F`, `ut_cyclic.f90` | Tests; not built by this CMake |

## 4. Who uses it
Callers use the external routines, so a `use` grep finds few.
- `GMAO_hermes` links `GMAO_gfio_r8` (`GMAO_hermes/CMakeLists.txt`, the `DEPENDENCIES` of `esma_add_library` and the include directories, inside the `if ( NOT HERMES_LIGHT )` block); about 27 files there mention `GFIO_*` (e.g. `m_dyn.f90`,
  `GFIO_remap.F90`, `dyn2prs.f90`, `getgfiodim.f90`). `m_dyn_util.F90`,
  `recalcsfc.F90`, `m_vdc2vdc.F90` use the modules `m_GFIO_PutFld` /
  `m_GFIO_GetFld` (not `m_gfio`).
- `GMAO_stoch` links `GMAO_gfio_r8` (`GMAO_stoch/CMakeLists.txt`, `esma_add_library` `DEPENDENCIES`);
  no GMAO_stoch source file mentions gfio.
- `GMAO_transf`: 8 files mention gfio (`fv2gg.F`, `m_cwf.F90`, `ec2fv.F90`,
  `m_fvSurface.F90`, `ss2fv.f90`, `ec2fv.csh`, `gg2fv.F`, `make_surf.f90`);
  not opened.
- GEOS-Chem interface: `GEOSCHEMchem_GridComp/@geos-chem/Interfaces/GEOS/geos_interface.F90`
  calls `GFIO_Open` / `GFIO_Close`.
- Chemistry `Shared/Chem_Base` links `GMAO_gfio_r4`
  (`Chem_Base/CMakeLists.txt`, `esma_add_library` `DEPENDENCIES`). `Chem_Shared`
  (`Chem_UtilMod_geos4.F90`, `Chem_SettlingCalculator.F90`, `Chem_UtilMod.F90`)
  and `Chem_Base` (`Chem_Aod.F90`, `Chem_BundleMod.F90`, ...) call
  `GFIO_Open`, `GFIO_GetVar`, `GFIO_GetVarT1`, `GFIO_Close`.
- `@GEOS_Util/post` scripts mention it (count not verified; nested repo, not read).
- MAPL's own CFIO is separate: a commented line in `GMAO_stoch/CMakeLists.txt`
  lists `GMAO_gfio_r8` together with `MAPL_cfio_r4`.
- Other fixtures (GEOSadas): *Unverified:*.

## 5. Known traps
- Two precision builds; consumers choose `_r4` or `_r8` (Chemistry `r4`,
  hermes/stoch `r8`). The unsuffixed include directory is a symlink to `_r8`.
- Calls are external routines without explicit interfaces for most of the API,
  so precision mismatches may not be caught at compile time (*Inferred*).
- `GMAO_gfioeos` (named in the `esma_add_subdirectories` list of `GMAO_Shared/CMakeLists.txt`) is a different name.

## 6. Observations (noted, not fixed)
- Several files (`m_gfio.f90`, `gfio_init.f`, `gfio_flush.F`, tests, C
  bindings) are not in the CMake source list, yet other libraries mention
  `m_gfio`. Which of them are still intended to build is unclear.
- `README` points to an intranet page and a Postscript manual not in the
  directory.
