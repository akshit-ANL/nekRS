# dgComp in nekRS: implementation notes and Fortran map

This document describes how the CompNek5000 DGSEM solver is laid out in nekRS. For every Fortran routine it gives the nekRS function or kernel that replaces it.
- Plan and design: [ROADMAP.md](ROADMAP.md).
- Status of the work items (`T-nn`): [TRACKER.md](TRACKER.md).

**Rule:** when you port, move or stub a Fortran routine, update its row here in the same commit.

Last updated: 2026-11-25

---

## 1. Why there are only three plugin files

`CompNek5000/src` has 14 Fortran files plus the `COMPRESS` include, about 4.3k lines. The plugin has 1.1k lines in three files.

| file | lines | role |
|---|---|---|
| `src/app/nrs/plugins/dgComp.hpp` | 145 | public API (`dgComp::` namespace), primitive field indices `dgComp::prim` |
| `src/app/nrs/plugins/dgComp.cpp` | 703 | host side: options, thermo input, memory, face exchange setup, RK driver, checkpoint fields |
| `src/app/nrs/plugins/kernels/dgComp.okl` | 299 | device kernels, one per pointwise or per-element operation |

The code is short for three reasons.

1. **Most of the physics is not ported yet.** The plugin is at milestone M0/M1. The volume flux, the interface flux, the boundary conditions, the CFL, the viscous terms, entropy viscosity, the limiter and the chemistry are stubs or absent. That accounts for about 3k of the 4.3k Fortran lines; see the *Status* column in §3. The plugin grows as those land. The kernels stay in `dgComp.okl`, which can be split per topic later (for example `dgFlux.okl`, `dgViscous.okl`) without changing the API.
2. **nekRS supplies what the Fortran built by hand.** These include `full2face`, `gs_op` on `dg_hndlx`, `bm1`/`jacmi`, `rxm1..tzm1`, `area`/`unx`, `dg_setup` and the `.fld` I/O. In nekRS they are `mesh->vmapM`, `oogs` on `mesh->globalFaceIds`, `mesh->o_vgeo`/`o_sgeo`, `mesh->o_LMM` and the checkpoint writer.
3. **Small Fortran helpers become one line inside a kernel.** For example `eval_EOS_D` is `rho = p / (Ru T sum(Y/W))` in `dgConserved`, and `eval_dgeh` is `h = (rhoE + p)/rho` in `dgPrimitives`.

Case-specific code from `base.usr` (initial condition, statistics) lives in the example's `.udf`, not in the plugin. See §4.

## 2. Layout and call flow

### State (replaces the `COMPRESS` common blocks)

| Fortran (`COMPRESS`) | nekRS |
|---|---|
| `U_EQ, V_EQ, W_EQ, E_EQ, N_EQ = 4 + nspec` (1-based; ρ is not a separate equation) | conserved fields 0, 1, 2 = ρu, ρv, ρw; 3 = ρE; 4.. = ρY_k (0-based, stride `nrs->fieldOffset`), `p_Neq`. ρ = Σ ρY_k, as in Fortran |
| `R_u = 8.3144598` | `Ru` in `dgComp.cpp`, `p_Ru` in kernels (same value) |
| `rho, ux/uy/uz, pr, temp1, gamma, aspd, dgeh, wmean, Yi` (`comp_vars*`) | primitives `o_primv`, indices `dgComp::prim::RHO, U, V, W, P, T, GAMMA, ASPD, H, WMEAN, Y0+k`, kernel defines `p_PRHO`.. |
| `bi, aa, h0, s0, mw` (`comp_vars11`, `comp_vars2`) | `o_bi`, `o_aa`, `o_mw` (device), `dgComp::thermo_t` (host) |
| `*f` (own trace) and `*pi` (self + neighbour) face arrays (`comp_flux2/3`, `comp_temp1`) | `o_faceM`, `o_faceP`, `p_Ntrace = Neq + Nprim` fields of stride `faceOffset` |
| `ifviscous, ifavisc, ifsplitchem, ifclip, ifpplim, ifIntLLFfull, ifVolEC` (`comp_flags*`) | `opt` struct in `dgComp.cpp`, filled from the `[DGCOMP]` par section (`viscous`, `artificialViscosity`, `splitChemistry`, `clip`, `positivityLimiter`, `fullInterfaceLLF`, `volumeFlux = kep / ec`) |
| `uparam(1..5)` | `[DGCOMP] dtChem, avCstar, avCmax, avP, riemannSolver` |
| `lmf, enmu, enhe, rhoS1/2, vis1, cond1, diff1, wdot, ss` | not allocated yet (M2 to M4) |

Primitives with no nekRS field yet: `spd`, `mach`, `ss`/`si` (entropy), `Xi` (mole fractions), transport, `wdot`.

### Call flow per time step

```
nekRS time loop (solver = none, so nekRS only advances time and writes checkpoints)
└─ UDF_ExecuteStep(time, tstep)                 [case .udf, replaces userchk]
   ├─ dgComp::step(time - dt, dt)               [runcomp.f::runcompdg + rungeKutta.f::SSPRK3]
   │   ├─ copy U -> U0
   │   └─ 3 stages:
   │       ├─ computeRhs(t)                     [dgrhs.f::getdgrhs, advect.f::advect]
   │       │   ├─ dgTraceGather + oogs add      [fluxLF.f::sum_neighbor_vars]
   │       │   ├─ dgVolumeRhsHex3D   (STUB)     [advect.f::advect_volInt]
   │       │   └─ dgSurfaceRhsHex3D  (STUB)     [advect.f::advect_interface, riemann.f]
   │       ├─ dgSSPRKStage                      [SSPRK3 stage update]
   │       └─ dgPrimitives                      [compvar.f::eval_dependent_vars]
   │   └─ copyToFluid (u, p into nrs->fluid for checkpoints)
   └─ case diagnostics (stat.dat, cons.dat)     [base.usr calc_stat, check_cons]
```

**Sign convention:** the Fortran `rkRHS` is −dU/dt. `o_rhs` in nekRS is +dU/dt. Keep this in mind when you compare dumps.

**While the two RHS kernels are stubs, `o_rhs` = 0 and the solution never changes.** u, p and T stay at the initial condition for the whole run. In particular p stays at p0 = 1 bar, and the velocity keeps its initial urms (120.5 m/s for compHIT) while the Nek5000 run decays to 58.7 m/s by t = 4.5 τ0.

## 3. Fortran routine → nekRS map

Status values:
- **ported**: the nekRS code does the same computation.
- **partial**: part of the routine is ported; the note says what is missing.
- **stub**: a kernel or hook exists but does nothing yet.
- **todo**: nothing in nekRS yet.
- **nekRS**: replaced by an existing nekRS facility.
- **n/a**: not needed in nekRS.

Names without a file are in `dgComp.cpp` (host) or `kernels/dgComp.okl` (kernels starting with `dg`).

### runcomp.f — driver

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `runcompdg` (setup branch: `dg_setup`, `dg_setup2`, `eval_enhe`) | `dgComp::setup()` → `setupImpl` | partial | no `eval_enhe` (T-38) |
| `runcompdg` (unsplit branch: `compress_cfl` + `SSPRK3(dt)`) | `dgComp::step(time, dt)`, called from `UDF_ExecuteStep` | partial | dt is fixed from the par file, no CFL (T-20) |
| `runcompdg` (Strang branch: flow dt/2, chem dt, flow dt/2) | — | todo | T-39; setup aborts if `reacting = yes` |
| `check_mass_int` | `dgCompHIT.udf::checkConservation` (cons.dat) | partial | case level, not in the plugin |
| `check_mass_min` | — | todo | commented out in the Fortran driver too |
| `printFlag` | `printf` block at the end of `setupImpl` | ported | |

### rungeKutta.f — time integration, chemistry, limiter

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `SSPRK3` | `dgComp::step` + `dgSSPRKStage` | ported | same Shu-Osher α, β; stage times `cstage = {0, 1, 1/2}` |
| `integrate_chem`, `chemRHS` | — | todo | T-39 (host Cantera, M4), T-42 (GPU, M5) |
| `pp_limiter`, `eval_rU_shift`, `eval_rS_shift` | hook comment in `dgComp::step` | todo | T-37 (M3); setup aborts if `positivityLimiter = yes` |

### compvar.f — primitives and thermodynamics

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `eval_dependent_vars` | `dgComp::computePrimitives` → `dgPrimitives` | partial | ported: ρ, Y (with clip), W̄, u, T (Newton), γ, a, h, p. Not yet: μ, λ, D_k (T-32), ω̇ (T-40), entropy `ss`/`si` (T-37/T-38), Mach |
| (Newton for T inside `eval_dependent_vars`) | `dgPrimitives`, `#if !p_perfectGas` | ported | same start (500 K), tolerance (1e-8), cap (1000). The species sum is folded into mixture coefficients before the loop, which is algebraically identical |
| — | `dgPrimitives`, `#if p_perfectGas` | new | `thermo = perfectGas`: T = e/cv, constant γ |
| `eval_EOS_D` | inline in `dgConserved` | ported | ρ = p W̄ / (Ru T) |
| `eval_EOS_P` | inline in `dgPrimitives` | ported | p = ρ Ru T / W̄ |
| `eval_EOS_T` | — | n/a | not used on the current path |
| `eval_dgeh` | inline in `dgPrimitives` (`prim::H`) | ported | |
| `get_bi` | `computeBi` | ported | |
| (`useric` conserved-from-primitives, case code) | `dgComp::conservedFromPrimitives` → `dgConserved` | ported | generic, so every case can use it |

### fluxLF.f — face traces and LLF pieces

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `sum_neighbor_vars`, `sum_neighbor` (`full2face` + `gs_op add`) | `dgTraceGather` + `oogs::startFinish(o_faceP, ...)` in `computeRhs` | ported | gives M + P like the `*pi` arrays; the neighbour is `faceP − faceM`. Periodic pairing check: T-13 |
| `eval_lmf` (λ_max on faces, `gs_op max`) | `dgWaveSpeed` (pointwise \|u\|+a only) | partial | no face max yet; part of T-19. See B-01 |
| `compute_flux` | — | todo | inside the T-18/T-19 kernels |

### advect.f — inviscid volume and interface terms

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `advect` | `computeRhs` | ported | orchestration only |
| `advect_volInt` | `dgVolumeRhsHex3D` | **stub** (rhs = 0) | T-18, main performance kernel |
| `getFlux_KEP` | — | todo | T-18 (`volumeFlux = kep`) |
| `getFlux_EC`, `ecAvg`, `ecLog` | — | todo | T-18 (`volumeFlux = ec`) |
| `advect_interface` (LLF) | `dgSurfaceRhsHex3D` | **stub** (adds nothing) | T-19 |

### riemann.f — HLL / HLLC

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `advect_interface_riemann` | `dgSurfaceRhsHex3D` (HLLC branch) | todo | T-19; setup aborts unless `riemannSolver = llf` |
| `riemann_flux`, `phys_flux`, `cons_state`, `flux_hll`, `flux_hllc` | — | todo | T-19, as device functions in `dgComp.okl` |

### dgrhs.f, compf.f — RHS assembly and source term

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `getdgrhs` | `computeRhs` | partial | inviscid stubs only; viscous, AV and source hooks are TODO comments |
| `compf` (ω̇ in the species RHS, unsplit) | — | todo | T-40 (M4) |

### compbc.f — boundary conditions

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `compbc`, `compbc_setup`, `compbc_flux` | `EToB > 0` branch in `dgSurfaceRhsHex3D` | todo | T-36 (M3). Not needed for periodic HIT |
| `compbc_visc` | — | todo | T-36 with T-31 |

### compcfl.f — time-step limits

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `compress_cfl`, `cmp_acfl` | `dgComp::maxWaveSpeed` (max \|u\|+a, reported only) | partial | no dt control yet (T-20) |
| `cmp_dcfl` | — | todo | M2, with the viscous terms |
| `mindr`, `dist3`, `dist2` | — | todo | T-20 |

### dgdiv.f, viscous.f — BR1 viscous terms

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `dgGrad_strong`, `dudxyz_strong`, `dgGrad_weak` | — | todo | T-31 |
| `viscous`, `eval_dudx`, `eval_visc_U`, `eval_visc_EY` | hook comment in `computeRhs` | todo | T-31; `viscous = yes` currently prints a warning and runs inviscid |

### entrovisc.f — entropy viscosity

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `entrovisc` | — | todo | T-38 |
| `liftGrad`, `divQ_strong`, `divQ_weak` | — | todo | shared with BR1 (T-31) |
| `eval_enmu`, `eval_enhe`, `eval_lmf_visc_reg` | — | todo | T-38; setup aborts if `artificialViscosity = yes` |

### cantera_stubs.f and core/ftnlib.cpp — Cantera

| Fortran / C++ | nekRS | status | notes / task |
|---|---|---|---|
| `newIdealGasMix`, `getmolecularweights` (in `usrdat`) | `dgComp::readThermo` reads `mw.dat` | partial | temporary until T-32/T-34 |
| `viscosity`, `thermalconductivity`, `getmixdiffcoeffs` | — | todo | T-32 (host, M2), T-41 (GPU, M5) |
| `getnetproductionrates`, `setstate_uv`, `setmassfractions` | — | todo | T-39/T-42 |
| `cantera_stub_abort` and the other stubs | — | todo | T-33 (stubs when built without Cantera) |

### Nek5000 core routines that src/*.f calls

| Nek5000 | nekRS |
|---|---|
| `dg_setup`, `dg_setup2`, `dg_hndlx` | `oogs::setup` on `mesh->globalFaceIds` (`setupImpl`) |
| `full2face` | `mesh->o_vmapM` (`dgTraceGather`) |
| `gs_op(..., add)` / `gs_op(..., max)` | `oogs::startFinish(..., ogsAdd)`; max not used yet |
| `rxm1..tzm1`, `jacm1`, `jacmi`, `bm1` | `mesh->o_vgeo` (`p_RXID..p_TZID`, `p_JID`, `p_IJWID`, `p_JWID`), `mesh->o_LMM` |
| `area`, `unx/uny/unz` | `mesh->o_sgeo` (`p_NXID..`, `p_WSJID`) |
| `dxm1` | `mesh->o_D` |
| `glsc2(f, bm1, n)`, `volvm1` | `platform->linAlg->innerProd(N, mesh->o_LMM, f)`, `mesh->volume` |
| `outpost`, `.fld` writer | nekRS checkpoint; `dgComp::addCheckpointFields` registers ρ and Y_k (ADIOS names or `scalarNN`) |

## 4. base.usr (compHIT) → dgCompHIT.udf

| base.usr / hit_fourier.f | dgCompHIT.udf | status |
|---|---|---|
| `usrdat` (Cantera init, `aa/h0/s0.dat`, `get_bi`) | `dgComp::setup()` (plugin) | ported (`mw.dat` replaces Cantera) |
| `usrdat2` (`rescale_x` to [0, L]) | `linAlg->rescale` + `mesh->update()` in `UDF_Setup` | ported |
| `usrdat3` (feature flags) | `[DGCOMP]` par section | ported |
| `usrdat3` → `init_hit_direct` (read `hit_modes.dat`) | `readModes` | ported |
| `hit_eval` (hit_fourier.f) | `fourierVelocity` kernel | ported (the analytic `diag` derivatives are not computed) |
| `useric` (uniform T, p, Y; ρ, ρu, ρE, ρY) | `UDF_Setup` IC block + `dgComp::conservedFromPrimitives` | ported; T0, p0, Y0 come from `[CASEDATA]` |
| `hit_ic_report` | — | todo |
| `userchk` | `UDF_ExecuteStep` | ported |
| `calc_stat` (stat.dat) | `hitStatistics` | ported |
| `hit_ref_stat` (stat_ref.dat) | `hitStatistics` | ported |
| `check_cons` (cons.dat) | `checkConservation` | ported |
| `post_inst`, `get_xyzpts` | — | todo (post-processing) |
| `uservp`, `userf`, `userq`, `userbc` | — | n/a (the DG path does not call them) |
| `userbc_comp` | — | todo with T-36 |

## 5. Comparing with a Nek5000 run

- **Step 0 should match exactly.**
  - The nekRS startup prints urms0, lambda0 and tau0. The first `stat_ref.dat` row has the full set.
  - For compHIT the expected values are urms0 = 120.50, lambda0 = 1.3367e-5 and tau0 = 1.109e-7.
  - p = p0 and T = T0 are uniform in both codes at t = 0.
- **Later times cannot match until T-18 and T-19 land.**
  - The nekRS field stays frozen at the initial condition.
  - In Nek5000, urms decays (to 58.7 m/s at 4.5 τ0), and prms and Trms grow from 0.
- **Once the fluxes are in, use the harness (T-22).**
  - Dump one step of CompNek5000 and compare it field by field.
  - Mind the rkRHS sign flip (§2).
