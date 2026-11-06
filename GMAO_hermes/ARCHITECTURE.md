# GMAO_hermes

**Audience:** coding agents and developers working on the GMAO_hermes library and tools in GMAO_Shared.
**Describes:** the MAPL2 line, branch `main`.
**Verified against:** `GMAO_Shared` @ `b0ebdca` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`; this directory is `GMAO_hermes/`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Library-level document: purpose, contents and users. Read the source for routine detail.

## 1. Purpose

Library and tool set for the "dyn-vector" (FVGCM dynamics state) of the GEOS
data-assimilation heritage: reading and writing dyn-vector files, vertical
and horizontal remapping, eta-level (ak/bk) definition, interpolation to
pressure levels, and ensemble/analysis utilities. Most modules carry
1990s-2000s DAO/GEOS-DAS headers. In the GCM the live parts are small
(`m_set_eta`, `m_topo_remap`); the rest serves analysis and post-processing.

## 2. Build (`CMakeLists.txt`)

- Library `GMAO_hermes`, deps `GMAO_mpeu`, `GMAO_gfio_r8`. Defines `-DHERMES
  -DDEBUG_SHAVE -Dunix -D__unix__`; compiled with the `FREAL8` flags (real*8).
- Option `HERMES_LIGHT` (default OFF, declared in `GMAO_Shared/CMakeLists.txt`):
  when ON the library is only `m_set_eta.F90` + `shared_topo_remap.F90`, with no
  deps and no executables except `write_eta.x`.
- Normal build: 29 sources (the module list below plus `progsig.f90`) and
  42 `ecbuild_add_executable` targets `*.x` (41 in the full build plus `write_eta.x`;
  one per program, linked to the library; `lcv2prs.x` also to `MAPL`, `GEOS_Shared`;
  `dyndiff.x` also to `GEOS_Shared`; `GFIO_remap.x`, `GFIO_Shave.x`, `lcv2prs.x` to
  `MAPL_cfio_r8`). `write_eta.x` is built in both modes.
- Installed: `cnv2prs.pl`, `ana5sfc`, `echorc.pl` to `bin`, inside the full build only (skipped under `HERMES_LIGHT`); every
  `*.rc` (glob) to `etc`, unconditionally.
- Present, not compiled: `fv2prs.F90`, `extract_stations.F90`, `blendq.f90`,
  `rout2prs.f`, `geos2fv.f`/`.h`, `ut_*.f90` unit programs
  (with `.ctl` files), `alias.h`, `unused_dummy.H`, `GMAO_hermes_arch.mk`, `AnIntro`
  (BOI documentation text for the PSAS-FVGCM interface).

## 3. Contents

### Modules (library)

| Module / file | Role (from header) |
|---|---|
| `m_dyn` | FVGCM dynamics state class (dyn-vector type, get/put) |
| `m_ecdyn` | Variant of the dyn-state class (same header title) |
| `m_dynp` | Dynamics state perturbation class |
| `m_dyn_util` | Utilities on dyn-vectors |
| `m_dyn2dyn` | Dyn-vector to dyn-vector vertical interpolation (header text names `m_mapz`) |
| `m_set_eta`, `m_ec_set_eta` | Define eta-level `ak`/`bk` for a layer count (EC variant) |
| `m_mapz`, `m_mapz_pert` | Vertical interpolation of rst-like fields (and perturbations) |
| `m_maph`, `m_maph_pert` | Horizontal interpolation of rst fields (and perturbations) |
| `m_remap` | Remapping of dyn-vectors |
| `m_topo_remap`, `shared_topo_remap` | Remap dyn-vector after topography change |
| `m_interp` | Interpolation class for FVCCM |
| `m_interpack`, `m_interpack_ad` | Interpolation package (R. Todling) and its adjoint |
| `m_insitu` | In-situ vector simulators (grid to pressure levels / stations) |
| `m_ana2dyn` | PSAS analysis to FVGCM coupler |
| `m_lp` | Large-Pond surface-layer parameterizations |
| `m_vdc2vdc` | Map vertical diffusion coefficients between resolutions |
| `m_qsat`, `m_const`, `m_spline`, `m_tick` | Saturation humidity; constants; splines; date/time increment |
| `m_nc_akbk`, `m_nc_JEDIinc` | NetCDF ak/bk and JEDI-increment I/O |

### Programs (grouped by role, from header comments)

| Role | Programs |
|---|---|
| Dyn-vector conversion | `dyn2dyn`, `dyn52dyn`, `rs52dyn`, `dyn2rs5`, `drs2dyn`, `dyn2drs`, `dyn_ncf2dyn`, `diag2dyn`, `dyn_boot` |
| Pressure-level / eta output | `dyn2prs` (uses `m_insitu`), `lcv2prs` (eta to pressure, subsetting), `fvsubset`, `dyn2real_eta`, `cnv2prs.pl` |
| Remapping | `mapz`, `maph`, `maph_pert`, `dyn_pert_remap`, `GFIO_remap` (55-level mass fields to given levels), `GFIO_Shave` |
| Analysis / increments | `ana5sfc`, `gcmbkg2ana`, `ibc_upd`, `dyn_iupd`, `dyn_jediupd`, `recalcsfc` |
| Ensemble / sensitivity | `dyn_recenter`, `dyn_inflate`, `dyn_cov`, `dyn_rankhist`, `dyn_efsens`, `dyn_fsens_conv`, `dynp` (add/subtract perturbation) |
| Inspection / small tools | `dyndiff`, `dyndims`, `dyndot`, `dyn_hydro`, `dyn_blob`, `getgfiodim`, `reset_time`, `echorc`, `eta_echo`, `write_eta` |

### Resource files (`*.rc`, 85 files)

Read by the standalone post-processing tools. `fv2prs.rc` states "resource
file for the fv2prs utility"; `fv2prs.F90` defaults `rcfile` to `fv2prs.rc` and
accepts `-rc` (e.g. `fv2prs.lcv.rc`), and `fvsubset` accepts class names such as
`@tavg3d_met_e`. The files are HISTORY-style: `tavg*`/`inst*` per-collection
(e.g. `tavg3d_dyn_p`, `tavg3d_met_e`), `const_2d_*`, `dyn2prs.rc`, `dyn2dyn.rc`,
`ana.prs.rc`, `fcst_fv2prs*.rc`. This `CMakeLists.txt` installs all of them
to `etc` by glob; the readers are the tools above. *Inferred:* which `.rc` each
tool reads is set by its own `-rc` option or default, not by CMake.

## 4. Who uses it

- `GMAO_Shared/@GEOS_Util/post/CMakeLists.txt` links `GMAO_hermes` for
  `convert_eta.x`, `rs_vinterp.x`, `rsg3_vinterp.x`, `eta2rst.x`; sources there
  `use m_set_eta` (`ec_prs2eta.F`, `ncep_prs2fv.F`, `rs_vinterp.F90`); `m_tick` in
  `pre/prepare_ocnExtData/*_SST_FRACI_eight.F90`.
- `GMAO_transf` uses `m_dyn`, `m_mapz`, `m_remap`, `m_const`, `m_set_eta`
  (`sfc_regrid.f90`, `ec2fv.F90`, `g4tog5.f90`, `m_xform.f90`).
- `GMAO_stoch` uses `m_const` and includes `include_GMAO_hermes`
  (`GMAO_stoch/CMakeLists.txt`).
- `FVdycoreCubed_GridComp`: `m_set_eta` in `DynCore_GridCompMod.F90`,
  `AdvCore_GridCompMod.F90`, `interp_restarts.F90`, `interp_restarts_bin.F90`;
  `m_topo_remap` (`dyn_topo_remap`) in `DynCore_GridCompMod.F90`.
- `GEOSchem_GridComp`: `m_set_eta` in `Shared/Chem_Shared/Lightning_mod.F90`.
- `GEOSgcm_GridComp`: `m_set_eta` in `GEOSmkiau_GridComp/GEOS_mkiauGridComp.F90`.
- No `use` of `m_interp`, `m_maph`, `m_qsat`, `m_spline`, `m_insitu`, `m_lp`
  found outside this directory in the GEOSgcm fixture (*Unverified:* GEOSadas).

## 5. Known traps

- Real*8 is forced (`FREAL8`); `.F`/`.F90` files are preprocessed; with
  `EXTENDED_SOURCE`, eight files get extended-line flags (`CMakeLists.txt` l.24-33).
- `HERMES_LIGHT` changes the library to two modules; code needing others fails there.
- Module names `m_dyn`, `m_const`, `m_tick` etc. are common to the GMAO_Shared
  libraries; check which library defines the one in use.

## 6. Observations (noted, not fixed)

- `fv2prs.F90`, which the `fv2prs*.rc` files describe, has no executable target here.
- The `CMakeLists.txt` sets flags for `ibc_upd.F`, `gcmbkg2ana.F`, `diag2dyn.F` that are executables, not library sources.
- `m_dyn2dyn.f90` header names `m_mapz`; `m_ecdyn.f90` header repeats `m_dyn`'s title.
- `unused_dummy.H` and the `ut_*` programs are present but not built.
