# GMAO_Shared: decisions

**Audience:** coding agents and developers working on the GMAO_Shared repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GMAO_Shared` @ `595999b` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Only what is visible in the code now. Likely reasons are inferred.

### Light build option (`HERMES_LIGHT`)
- Pattern: `GEOS_Shared` and `GMAO_hermes` are always added; all other libraries only when `HERMES_LIGHT` is `OFF` (`CMakeLists.txt`). `GMAO_hermes` shrinks to a small set of modules (`GMAO_hermes/CMakeLists.txt`). [Verified]
- Likely reason: *Inferred:* a minimal build for a configuration that needs only a few dynamics utilities, without the assimilation libraries.
- Question for maintainers: which configurations use `HERMES_LIGHT`?
- On this line (`release/MAPL-v3`) the option is forced ON: `CMakeLists.txt` does `set (HERMES_LIGHT ON)` and adds `GMAO_mpeu`, `GMAO_hermes`, `GEOS_Shared` unconditionally; the `GMAO_hermes` light library gains `m_topo_remap.F90`, `m_dyn.f90`, `m_const.f90` and depends on `GMAO_mpeu` (`GMAO_hermes/CMakeLists.txt`). [Verified]

### `GMAO_gfio` built twice
- Pattern: `add_subdirectory (GMAO_gfio GMAO_gfio_r4)` and `... GMAO_gfio_r8`, `GMAO_gfio_r4` depends on `GMAO_gfio_r8` (`CMakeLists.txt`; `GMAO_gfio/CMakeLists.txt`, `precision`). [Verified]
- Likely reason: *Inferred:* consumers need 32-bit and 64-bit real interfaces of the same code (Chemistry `r4`; hermes, stoch `r8`); the ordering avoids a shared-directory race.
- Question for maintainers: is the r4 dependency on r8 an ordering fix or a code dependency?

### Sparse checkout in the GEOSgcm fixture
- Pattern: `config/GMAO_Shared.sparse` excludes nine directories (`GMAO_ncdiag`, `GMAO_ods`, `GMAO_psas`, `GMAO_pyobs`, `GEOS_Pert`, `arpack`, `pnagpack`, `GMAO_pilgrim`, `gmao2ioda`); `CMakeLists.txt` lists them anyway and skips missing ones (`esma_add_subdirectories`; `EXISTS` guard on `gmao2ioda`). [Verified]
- Likely reason: *Inferred:* one source repository serves several fixtures; the GCM does not need the assimilation code, and the build system tolerates absence.
- Question for maintainers: is `CMakeLists.txt` the shared list for all fixtures?

### CICE4 kept alongside CICE6
- Pattern: `LANL_Shared/CICE4` is built as target `CICE4` and linked by several components; `CICE_GEOSPlug` builds its own `cice6` (`LANL_Shared/CMakeLists.txt`; `CICE_GEOSPlug/CMakeLists.txt` in GEOSgcm_GridComp, `DEPENDENCIES MAPL cice6`). [Verified]
- Likely reason: *Inferred:* both sea-ice models are supported in different configurations during a transition.
- Question for maintainers: can both be linked into one executable?

### Two mpeu targets
- Pattern: `GMAO_mpeu` (real MPI) and `GMAO_eu` (stub `mpi0`) from the same sources (`GMAO_mpeu/CMakeLists.txt`). [Verified]
- Likely reason: *Inferred:* serial tools that must not need MPI.
- Question for maintainers: which executables are meant to use `GMAO_eu`?

### Restricted build while porting to MAPL3 (this line only)
- Pattern: `CMakeLists.txt` comment: "not yet ported to MAPL3 - only build GMAO_mpeu, GMAO_hermes (light), and GEOS_Shared for FVdycoreCubed_GridComp (needs m_set_eta, m_topo_remap, tropovars)"; "restore full build (and set HERMES_LIGHT OFF) when all components are ported". [Verified]
- Question for maintainers: which libraries will be ported first, and does `GMAO_hermes` keep a light mode afterwards?
