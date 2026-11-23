# dgComp port: task tracker

The plan and design live in [ROADMAP.md](ROADMAP.md). This file tracks the status of the work.

## How to edit

- **Adding a row:** copy an existing row and give it the next free ID in that table. IDs never get reused.
  - `T-nn` work item
  - `C-nn` example or validation case
  - `B-nn` bug, or a build or run failure
- **When something finishes:** move its row from *Pending* to *Implemented*. Keep the ID and fill in the commit.
- **Status values:**
  - `todo`: not started
  - `wip`: in progress
  - `blocked`: waiting on something; say what in Notes
  - `done`: code is written
  - `verified`: done and checked on a machine; say which in *Verified on*
- **Milestones:** `M0` to `M6` are the ROADMAP §4 milestones. Use `-` for items outside them.
- **Machines:** `polaris` (CUDA), `aurora` (DPCPP), `local`.
- **Bug reports:** use the template at the end of this file.

Last updated: 2026-11-23

---

## 1. Implemented

None of this code has been compiled yet.

| ID | Item | M | Status | Verified on | Commit | Notes |
|---|---|---|---|---|---|---|
| T-01 | Porting roadmap (`doc/dg/ROADMAP.md`) | - | done | - | 56afba11 | |
| T-02 | Plugin skeleton: `dgComp.{hpp,cpp}`, `dgComp.okl`, `cmake/nekrs.cmake`, `pluginTable` autoload | M1 | done | - | 7b2cc73d | Including `dgComp.hpp` in a .udf calls `dgComp::buildKernel` |
| T-03 | `[DGCOMP]` options: nspec, npa, flux choices, feature flags, species names | M1 | done | - | 7b2cc73d | Unimplemented features abort at setup |
| T-04 | `[DGCOMP]` userparam01..04: `dtChem`, `avCstar`, `avCmax`, `avP` | M1 | done | - | uncommitted | Parsed and printed only; used from M3/M4 on |
| T-05 | Thermo input: `aa/h0/s0/mw.dat` reader, `bi` polynomials (`get_bi`) | M1 | done | - | 7b2cc73d | `mw.dat` is temporary until the Cantera bridge |
| T-06 | Primitives kernel (EOS, Newton for T) and conserved-from-primitives kernel | M1 | done | - | 7b2cc73d | |
| T-07 | Trace gather and face exchange (oogs on `globalFaceIds`, gives M+P) | M1 | done | - | 7b2cc73d | Periodic pairing not checked yet (T-13) |
| T-08 | SSPRK3 driver and stage kernel | M1 | done | - | 7b2cc73d | |
| T-09 | Max wave speed reduction | M1 | done | - | 7b2cc73d | Not used for dt yet (T-20) |
| T-10 | Checkpoint fields: named with ADIOS (`density`, `Y_O2`, ...), numbered with .fld | M1 | done | - | 7b2cc73d | |
| T-11 | `dgCompHIT` case: par, udf, Fourier IC kernel, mesh rescale, `input.box`, README | M2 | done | - | 7b2cc73d | |
| T-12 | `dgCompHIT` `[CASEDATA]`: boxLength, T0, p0, Y0 | M2 | done | - | uncommitted | Read in `UDF_Setup0`. Replaces the constants that were hard-coded in the udf |

## 2. Pending (current: M0 and M1)

| ID | Item | M | Status | Notes |
|---|---|---|---|---|
| T-13 | Check that `globalFaceIds` pairs periodic faces | M1 | todo | Needed by the HIT box |
| T-14 | Separate `dg_dev` install on Polaris, run `examples/ethier` | M0 | todo | Prefix `$HOME/.local/nekrs-dg-polaris` |
| T-15 | Separate `dg_dev` install on Aurora, run `examples/ethier` | M0 | todo | Prefix `$HOME/.local/nekrs-dg-aurora` |
| T-16 | Module and env scripts per machine under `doc/dg/` | M0 | todo | |
| T-17 | Compile `dgCompHIT`, run the M0 stub; the field should stay unchanged | M0 | todo | Also needs `hit_modes.dat`, `aa/h0/s0/mw.dat` copied in |
| T-18 | KEP and EC flux-differencing volume kernel (`advect_volInt`) | M1 | todo | Main performance kernel |
| T-19 | LLF and HLLC surface flux, lift (`advect_interface`, `riemann.f`) | M1 | todo | |
| T-20 | CFL and adaptive dt (`compcfl.f`) | M1 | todo | |
| T-21 | Restart from checkpoint (read rho, u, p, T, Y_k by name) | M1 | todo | The udf aborts on restart for now |
| T-22 | Verification harness: CompNek5000 one-step dump plus compare script | M1 | todo | |
| T-23 | Test: EC flux conserves total entropy to round-off | M1 | todo | |
| T-24 | Test: isentropic vortex convergence (C-03) | M1 | todo | |
| T-25 | Cantera source build on Polaris | M1 | todo | Must be done by the end of M1 |
| T-26 | Cantera source build on Aurora | M1 | todo | Must be done by the end of M1 |

## 3. Future pipeline

| ID | Item | M | Status | Notes |
|---|---|---|---|---|
| T-30 | HIT statistics (`calc_stat` to `stat.dat`: urms, lambda, Re_lambda, Mt, Trms, prms) | M2 | todo | Needed for validation against `stat_ref.dat` |
| T-31 | BR1 viscous kernels (gradient lift, weak divergence) | M2 | todo | |
| T-32 | `dgCantera` host class (mu, lambda, D_k), Phase A | M2 | todo | |
| T-33 | `ENABLE_CANTERA` / `CANTERA_DIR` in CMake, with stubs when Cantera is absent | M2 | todo | |
| T-34 | Molar masses and thermo fits from Cantera instead of the `.dat` files | M2 | todo | |
| T-35 | Validate compHIT, LLF and HLLC, against `stat_ref.dat` and Bielawski Fig. A.15 | M2 | todo | |
| T-36 | Boundary-condition kernel (`compbc.f`, keyed on EToB) | M3 | todo | |
| T-37 | Zhang-Shu positivity limiter | M3 | todo | |
| T-38 | Entropy viscosity (`eval_enmu`, `eval_enhe`) | M3 | todo | Uses avCstar, avCmax, avP |
| T-39 | Strang-split chemistry (`integrate_chem`), host Cantera | M4 | todo | Uses dtChem |
| T-40 | Unsplit chemistry source (`compf.f`) | M4 | todo | |
| T-41 | GPU transport (Phase B) | M5 | todo | |
| T-42 | GPU kinetics, generated OKL (Phase C) | M5 | todo | Ask the nekRS developers about existing code first |
| T-43 | Performance: shared-memory volume flux, overlapped exchange, fused kernels | M6 | todo | |
| T-44 | Scaling runs on Polaris and Aurora | M6 | todo | |

## 4. Examples and validation cases

| ID | Case | Purpose | Source | Needs | Status |
|---|---|---|---|---|---|
| C-01 | `dgCompHIT` (Re_lambda=100, Mt=0.6) | Main 3D validation | `compHIT_direct_Re_lamda_100/LLF` | M1 (inviscid), M2 (full) | wip |
| C-02 | `tgvComp` | Optional bring-up, analytic IC | CompNek5000 `tgvComp` | M1 | todo |
| C-03 | Isentropic vortex | Convergence order | new | M1 | todo |
| C-04 | `dmr2d`, extruded one element in z | Shocks, BCs, limiter | CompNek5000 `dmr2d` | M3 | todo |
| C-05 | `det2d`, extruded (H2, 10 species) | Reacting | CompNek5000 `det2d` | M4 | todo |
| C-06 | HIT at 64^3 and 128^3 | Weak scaling | C-01 | M5 | todo |
| C-07 | 3D detonation | End goal | new | M6 | todo |

## 5. Bugs and build or run failures

Stage is one of `configure`, `build`, `udf-compile`, `kernel-jit` or `run`.

| ID | Date | Machine | Stage | Summary | Status | Fix / commit |
|---|---|---|---|---|---|---|
| B-01 | 2026-11-23 | - | - | Possible bug in the Fortran `fluxLF.f::eval_lmf`: `k` is incremented but never set to 0 first. The port initialises it. If results differ from Nek5000, check this first. | watch | - |

### Bug report template

Copy this into an issue or under a `### B-nn` heading here, and add a row to the table above.

```
ID:        B-nn
Date:
Machine:   polaris | aurora | local
Stage:     configure | build | udf-compile | kernel-jit | run
Commit:    (git rev-parse --short HEAD)
Modules:   (module list)
Command:
Error (first lines that matter):

Notes / suspected cause:
```
