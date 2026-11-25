# dgComp in nekRS: implementation notes and Fortran map

This document describes how the CompNek5000 DGSEM solver is laid out in nekRS. For every Fortran routine it gives the nekRS function or kernel that replaces it.
- Plan and design: [ROADMAP.md](ROADMAP.md).
- Status of the work items (`T-nn`): [TRACKER.md](TRACKER.md).

**Rule:** when you port, move or stub a Fortran routine, update its row here in the same commit.

Last updated: 2026-11-25

---

## 1. Why there are only three plugin files

`CompNek5000/src` has 14 Fortran files plus the `COMPRESS` include, about 4.3k lines. The plugin has 1.6k lines in three files.

| file | lines | role |
|---|---|---|
| `src/app/nrs/plugins/dgComp.hpp` | 148 | public API (`dgComp::` namespace), primitive field indices `dgComp::prim` |
| `src/app/nrs/plugins/dgComp.cpp` | 765 | host side: options, thermo input, memory, face exchange setup, RK driver, checkpoint fields |
| `src/app/nrs/plugins/kernels/dgComp.okl` | 698 | device kernels, one per pointwise or per-element operation |

The code is short for three reasons.

1. **Most of the physics is not ported yet.** The plugin is at milestone M1. The inviscid fluxes are in (KEP volume flux; LLF, HLL, HLLC and blended interface fluxes), and so are the BR1 viscous terms with mixture-averaged transport. The EC volume flux, the boundary conditions, the CFL, entropy viscosity, the limiter and the chemistry are absent. That accounts for about 2.5k of the 4.3k Fortran lines; see the *Status* column in §3. The plugin grows as those land. The kernels stay in `dgComp.okl`, which can be split per topic later (for example `dgFlux.okl`, `dgViscous.okl`) without changing the API.
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
| `uparam(1..5)` | `[DGCOMP] dtChem, avCstar, avCmax, avP, riemannSolver` (`llf, hll, hllc, hllcBlend` = 0..3) |
| `uparam(6..8)` (blend switch, `riemann.f`) | `[DGCOMP] blendSw0, blendSw1, blendSwRho`, same fallbacks (0.01, 10·sw0, 0.1) |
| `vis1, cond1, diff1` (μ, λ, ρD_k) | `o_trans`, `p_Nspec + 2` fields: μ, λ, ρD_k (viscous only) |
| BR1 gradients (`gradU`, `gradT`, ...), viscous fluxes, their face traces | `o_grad` (`3·Nvisc` fields, `[(3g+d)·offset]`), `o_vflux` (same layout), `o_vfaceM` / `o_vfaceP` (`Nvisc` normal-flux fields of stride `faceOffset`) |
| Cantera transport fits (inside Cantera) | `o_viscFit`, `o_condFit`, `o_diffFit` from `visc/cond/diff.dat` |
| `lmf, enmu, enhe, rhoS1/2, wdot, ss` | not allocated yet (M3 to M4) |

`Nvisc = 4 + (nspec > 1 ? nspec : 0)`: the gradient variables are u, v, w, T and, for a mixture, X_k. The viscous fluxes reuse the same count: ρu, ρv, ρw, ρE and ρY_k.

Primitives with no nekRS field yet: `spd`, `mach`, `ss`/`si` (entropy), `Xi` (mole fractions; computed on the fly in the viscous kernels), `wdot`.

### Call flow per time step

```
nekRS time loop (solver = none, so nekRS only advances time and writes checkpoints)
└─ UDF_ExecuteStep(time, tstep)                 [case .udf, replaces userchk]
   ├─ dgComp::step(time - dt, dt)               [runcomp.f::runcompdg + rungeKutta.f::SSPRK3]
   │   ├─ copy U -> U0
   │   └─ 3 stages:
   │       ├─ computeRhs(t)                     [dgrhs.f::getdgrhs, advect.f::advect]
   │       │   ├─ dgTraceGather + oogs add      [fluxLF.f::sum_neighbor_vars]
   │       │   ├─ dgVolumeRhsHex3D  (KEP)       [advect.f::advect_volInt, getFlux_KEP]
   │       │   ├─ dgSurfaceFluxHex3D            [advect.f::advect_interface, fluxLF.f::eval_lmf, riemann.f]
   │       │   ├─ dgSurfaceLiftHex3D            [the area/bm1 lift inside advect_interface]
   │       │   └─ viscous = yes:                [dgrhs.f -> viscous.f::viscous]
   │       │       ├─ dgTransport               [compvar.f vis1/cond1/diff1, Cantera MixTransport]
   │       │       ├─ dgViscGradHex3D × Nvisc   [dgdiv.f::dgGrad_weak + entrovisc.f::liftGrad]
   │       │       ├─ dgViscFlux                [viscous.f::eval_visc_U, eval_visc_EY]
   │       │       ├─ dgViscTraceGather + oogs  [full2face of F·n + gs_op add]
   │       │       └─ dgViscDivHex3D × Nvisc    [entrovisc.f::divQ_weak]
   │       ├─ dgSSPRKStage                      [SSPRK3 stage update]
   │       └─ dgPrimitives                      [compvar.f::eval_dependent_vars]
   │   └─ copyToFluid (u, p into nrs->fluid for checkpoints)
   └─ case diagnostics (stat.dat, cons.dat)     [base.usr calc_stat, check_cons]
       └─ dgComp::computeTransport()            [rungeKutta.f eval_dependent_vars after the step: vis1 for Re_lambda]
```

**Sign convention:** the Fortran `rkRHS` is −dU/dt. `o_rhs` in nekRS is +dU/dt. Keep this in mind when you compare dumps.

**Status.** With `viscous = no` a nekRS run matches a Nek5000 run with `ifviscous = .false.`; the inviscid HLLC compHIT matches to 12-13 digits (TRACKER B-05). With `viscous = yes` (BR1, T-31; fit-based transport, T-50) the viscous HLLC compHIT matches CompNek5000 to about 1e-7 relative over 4.5 τ0. The transport inside `computeRhs` is evaluated at stage states; `dgComp::computeTransport()` evaluates it at the end-of-step state, which is what the Fortran `vis1` holds when `userchk` runs. Setup aborts on any boundary face (fully periodic meshes only, until T-36).

**BR1 in nekRS form.**
- Gradient: strong derivative in each element plus the central lift ½(q_P − q_M)·n·WSJ·IJW on each face node. The exchanged q are the M+P primitives already in `o_faceP`, so there is no second exchange for the gradients.
- Divergence: weak form, −IJW·Σ_l (D_li F̂r + D_lj F̂s + D_lk F̂t) with F̂r = JW·(rx Fx + ry Fy + rz Fz), plus the central face flux ½(F_M + F_P)·n. Only F·n_M is exchanged, one field per conserved equation.
- Full metrics at every node. The Fortran `dgGrad_weak` uses node-i metrics and a single weight ratio, which is the same on the affine HIT box but not on curved meshes (TRACKER B-04).

**Surface term in two kernels.** `dgSurfaceFluxHex3D` runs once per face node and writes the numerical-flux correction `−(F* − F_M)·WSJ·IJW` (= Fortran `(F* − F_L)·area/bm1` with the sign flip) to `o_faceFlux`. `dgSurfaceLiftHex3D` runs per element; thread (i,j) owns column k, so it adds the faces touching each node without atomics or a race.

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
| `eval_dependent_vars` | `dgComp::computePrimitives` → `dgPrimitives` | partial | ported: ρ, Y (with clip), W̄, u, T (Newton), γ, a, h, p. μ, λ, ρD_k are in the separate `dgTransport` kernel (T-50). Not yet: ω̇ (T-40), entropy `ss`/`si` (T-37/T-38), Mach |
| (transport part of `eval_dependent_vars`: `viscosity`, `thermalconductivity`, `getmixdiffcoeffs`) | `dgTransport` | ported | T-50. Mode 1 (`transport = fits`): Cantera `MixTransport` formulas from the fits in `visc/cond/diff.dat`. Mode 0 (`constant`): μ, λ = μ cp/Pr, ρD = μ/(Pr Le). Writes ρD_k, as `diff1` holds |
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
| `eval_lmf` (λ_max on faces, `gs_op max`) | inside `dgSurfaceFluxHex3D` | ported | λ = max(a + \|u·n\|) over M and P, computed pointwise from both traces, so no `gs_op max` is needed. `dgWaveSpeed` is the separate \|u\|+a reduction for the CFL |
| `compute_flux` | `physFlux` (full LLF) / `kepFluxAdd` (partial LLF) in `dgComp.okl` | ported | `fullInterfaceLLF` picks ½(F_M + F_P) instead of the KEP average |

### advect.f — inviscid volume and interface terms

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `advect` | `computeRhs` | ported | orchestration only |
| `advect_volInt` | `dgVolumeRhsHex3D` | ported (KEP) | T-18. Primitives in shared memory, thread (i,j) loops over k, Σ_l 2 D_il F#(i,l)·(metrics at i) in r, s, t; writes rhs = −rkRHS |
| `getFlux_KEP` | `kepFluxAdd` | ported | T-18 |
| `getFlux_EC`, `ecAvg`, `ecLog` | — | todo | T-49; setup aborts on `volumeFlux = ec` |
| `advect_interface` (LLF) | `dgSurfaceFluxHex3D` (`riemannSolver = llf`) + `dgSurfaceLiftHex3D` | ported | T-19. Dissipation −½λ(U_P − U_M) |

### riemann.f — HLL / HLLC

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `advect_interface_riemann` | `dgSurfaceFluxHex3D` (`riemannSolver = hll, hllc, hllcBlend`) + `dgSurfaceLiftHex3D` | ported | T-19. L = own trace, R = `faceP − faceM`, Davis wave speeds |
| `riemann_flux` (incl. the HLLC/HLL pressure-jump blend) | `dgSurfaceFluxHex3D` | ported | blend weights from `blendSw0/blendSw1/blendSwRho` |
| `phys_flux`, `flux_hll`, `flux_hllc` | `physFlux`, `fluxHLL`, `fluxHLLC` | ported | |
| `cons_state` | — | n/a | the conserved states are read directly from `faceM`/`faceP` |

### dgrhs.f, compf.f — RHS assembly and source term

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `getdgrhs` | `computeRhs` | partial | inviscid and viscous (BR1) parts ported; AV and source hooks are TODO comments |
| `compf` (ω̇ in the species RHS, unsplit) | — | todo | T-40 (M4) |

### compbc.f — boundary conditions

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `compbc`, `compbc_setup`, `compbc_flux` | (planned: `EToB > 0` branch in `dgSurfaceFluxHex3D`) | todo | T-36 (M3). Not needed for periodic HIT; setup aborts if the mesh has any boundary face |
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
| `dgGrad_weak` + `liftGrad` (BR1 gradient) | `dgViscGradHex3D` (`br1GradLift` for the face term) | ported | T-31. Strong form with central lift, equivalent to weak form + lift on any mesh; full metrics (B-04) |
| `dgGrad_strong`, `dudxyz_strong` | — | n/a | not on the BR1 path used by compHIT |
| `viscous` | `computeRhs` (`opt.viscous` block) | ported | T-31 |
| `eval_dudx` | `dgViscGradHex3D` on u, v, w (`viscVar`) | ported | |
| `eval_visc_U` (τ, u·τ, λ∇T) | `dgViscFlux` | ported | τ = μ(∇u + ∇uᵀ) − ⅔μ(∇·u)I |
| `eval_visc_EY` (species diffusion, Σ h_k J_k) | `dgViscFlux`, `#if p_Nspec > 1` | ported | viscous species flux ρD_k (W_k/W̄)∇X_k minus Y_k times its sum over k (zero net mass flux); the energy flux adds Σ h_k of it, h_k from `bi` plus Ru T / W_k |

### entrovisc.f — entropy viscosity

| Fortran | nekRS | status | notes / task |
|---|---|---|---|
| `entrovisc` | — | todo | T-38 |
| `liftGrad` | `br1GradLift` in `dgViscGradHex3D` | ported | T-31; reused by T-38 |
| `divQ_weak` | `dgViscDivHex3D` (`br1DivFace` for the face term) + `dgViscTraceGather` + oogs | ported | T-31; reused by T-38 |
| `divQ_strong` | — | n/a | not on the BR1 path |
| `eval_enmu`, `eval_enhe`, `eval_lmf_visc_reg` | — | todo | T-38; setup aborts if `artificialViscosity = yes` |

### cantera_stubs.f and core/ftnlib.cpp — Cantera

| Fortran / C++ | nekRS | status | notes / task |
|---|---|---|---|
| `newIdealGasMix`, `getmolecularweights` (in `usrdat`) | `dgComp::readThermo` reads `mw.dat` | partial | `exportTransport.py` writes `mw.dat` from Cantera; T-34 for a run-time bridge |
| `viscosity`, `thermalconductivity`, `getmixdiffcoeffs` | `dgTransport` (GPU) from the fits `exportTransport.py` writes | ported | T-50; replaces T-32/T-41 for transport. Cantera default (non-CK) mixture-averaged model |
| `getnetproductionrates`, `setstate_uv`, `setmassfractions` | — | todo | T-39/T-42 |
| `cantera_stub_abort` and the other stubs | — | todo | T-33 (stubs when built without Cantera) |

### Nek5000 core routines that src/*.f calls

| Nek5000 | nekRS |
|---|---|
| `dg_setup`, `dg_setup2`, `dg_hndlx` | `oogs::setup` on `mesh->globalFaceIds` (`setupImpl`) |
| `full2face` | `mesh->o_vmapM` (`dgTraceGather`) |
| `gs_op(..., add)` / `gs_op(..., max)` | `oogs::startFinish(..., ogsAdd)`; max is not needed (λ is taken pointwise from M and P) |
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
- **Inviscid runs (`viscous = no`).**
  - Compare against a Nek5000 run with `ifviscous = .false.`, the same `ifVolEC = .false.` (KEP), the same `uparam(5)` / `riemannSolver`, `ifIntLLFfull` / `fullInterfaceLLF`, and the same fixed dt.
  - HLLC at dt = 1e-10 agrees to 12-13 digits until CompNek5000 goes NaN at 3.99e-8 s; nekRS goes NaN at 4.51e-8 s. Inviscid KEP is not entropy stable, so the exact crash time is set by round-off (TRACKER B-05).
- **Viscous runs (`viscous = yes`).** The reference is urms 58.7 m/s at 4.5 τ0. The transport must be the same Cantera mixture-averaged model on both sides (`transport = fits`).
- **Step-0 RHS is the sharpest check (T-22).**
  - Dump `rkRHS` from CompNek5000 after the first `getdgrhs` and `o_rhs` from nekRS after the first `computeRhs`; they should agree to round-off with the sign flipped (§2).
  - Then compare one full SSPRK3 step field by field.
