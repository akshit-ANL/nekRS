# DGSEM compressible multi-species plugin for nekRS — porting roadmap

Branch: `dg_dev`
Source being ported: `CompNek5000/src/*.f` (about 4.3k lines of Fortran) plus `core/ftnlib.cpp` (the Cantera bridge).

---

## 1. What the Nek5000 code actually does

| Fortran (CompNek5000/src) | Role | GPU porting notes |
|---|---|---|
| `runcomp.f::runcompdg` | Driver, called from `userchk`. Does Strang splitting (flow dt/2 → chem dt → flow dt/2) or unsplit. | Host-side driver in the plugin. |
| `rungeKutta.f::SSPRK3` | 3-stage SSP-RK update of all `N_EQ = 4 + nspec` conserved fields. | One axpby kernel (or `platform->linAlg`). |
| `rungeKutta.f::integrate_chem / chemRHS` | Explicit RK2 chemistry sub-cycling with `dt_chem = uparam(1)`. **Calls Cantera per point.** | Phase A: on the host. Phase C: on the GPU. |
| `rungeKutta.f::pp_limiter` | Zhang–Shu style per-element positivity limiter (ρ, ρY, ρu_shift, ρs_shift). Needs a Newton solve for T. | One work-group per element, with element reductions. |
| `compvar.f::eval_dependent_vars` | Primitive recovery: ρ = Σρ_k Y_k; Newton solve for T from the refit-NASA `bi` polynomials; p, h, γ, c, s. **Calls Cantera for μ, λ, D_k (and ω̇ if unsplit).** | Thermo is pure polynomial work and maps to one thread per node. Transport goes to Cantera. |
| `fluxLF.f::sum_neighbor_vars` | `full2face` + `gs_op(dg_hndlx, add)` gives *self+neighbour* face traces (the `*pi` arrays). `eval_lmf` uses `gs_op max`. | Maps directly onto nekRS face-oogs (see §2.3). |
| `advect.f::advect_volInt` | Flux-differencing volume term with KEP (Pirozzoli) or EC (Gouasmi 2020) two-point fluxes. O(N⁴) per element. | **Main performance kernel.** Use shared memory per element. |
| `advect.f::advect_interface`, `riemann.f` | Surface term: LLF (default), or HLL/HLLC when `uparam(5)>0`. | One face-node kernel, then lift into the volume using 1/bm1. |
| `compbc.f` | Ghost states for `W`, `SYM`, `V`/`v`, `O`/`o`, `ON` boundaries; `userbc_comp` sets user data. | Kernel keyed on `EToB`. The cbc strings need mapping to nekRS boundary IDs. |
| `viscous.f`, `entrovisc.f::liftGrad, divQ_weak` | BR1-type gradient lifting and weak divergence of viscous / AV fluxes (τ, q, species diffusion). | Gradient, flux and divergence kernels, reusing the face exchange. |
| `entrovisc.f::eval_enmu, eval_enhe` | Entropy-viscosity coefficient from the residual of ρs over time (`rhoS1`, `rhoS2`). | Pointwise work plus element max. |
| `compcfl.f` | Acoustic and diffusive CFL. | Global max reduction. |
| `compf.f` | Adds ω̇ to the species RHS (unsplit mode). | Pointwise. |

**Key finding: Cantera is *not* in the EOS hot path.**
- EOS, T(e) inversion, cp, h_k, s_k and γ already use the code's own NASA-polynomial refit:
  - `aa.dat`, `h0.dat` and `s0.dat` are read in `usrdat`;
  - `bi` is built by `get_bi`.
  - All of this is plain arithmetic, so it is fully GPU-portable.
- Cantera is called in only three places, all per GLL point, through `setmassfractions` + `setstate_uv`:
  1. `getnetproductionrates`: ω̇_k (`compvar.f:146`, `rungeKutta.f:174`);
  2. `viscosity`, `thermalconductivity`: μ, λ (`compvar.f:161`);
  3. `getmixdiffcoeffs`: D_k (`compvar.f:163`).
- These three are what need a GPU strategy (§3).

**Possible latent bug worth checking in the Fortran:**
- In `fluxLF.f::eval_lmf`, `k` is incremented (`k = k+1`, line 24) but never set to 0 first.
- It probably works by luck of stack contents or compiler flags. The port initialises it, so if results differ, check this first.

---

## 2. Mapping onto nekRS

### 2.1 Where the code lives
- New plugin files:
  - `src/app/nrs/plugins/dgComp.{hpp,cpp}`, following the pattern of `lowMach.cpp`;
  - kernels in `src/app/nrs/plugins/kernels/dgComp/*.okl`.
- `cmake/install_kernels.cmake` installs the kernels automatically.
- Build edits:
  - add one line to `cmake/nekrs.cmake` (the plugin source list);
  - add one entry to the `pluginTable` in `src/core/udf/udfMake.hpp:140`, so that including `dgComp.hpp` in a `.udf` auto-registers its kernels.
- Later, if it outgrows the plugin model, promote it to its own `app_t`, selected through `APPLICATION` (`src/lib/nekrs.cpp:177`).

### 2.2 Time stepping
- Leave `[FLUID VELOCITY]` and all scalars out of the `.par` file. `fluid` is then never constructed and the nekRS step is effectively a no-op.
- Drive SSPRK3 and the chemistry from `UDF_ExecuteStep`. This mirrors how `userchk → runcompdg` works today.

### 2.3 DG trace exchange
- nekRS has `vmapM` but no filled `vmapP`.
- It does have the exact primitive the Fortran uses: a face-trace oogs handle built on `mesh->globalFaceIds` with `ogsAdd`. See `src/core/gjp.hpp:36-49`.
- That handle returns self+neighbour, which is identical to `gs_op(dg_hndlx, …, add)`, so `*pi` arrays port 1:1.
- Use `ogsMax` for `lmf`.
- Use split `oogs::start/finish` to overlap the volume kernel with communication.
- Batch all fields into one exchange (`N_EQ + primitives`, stride = face-array size) rather than calling the exchange once per variable as the Fortran does.
- **Verify early** that `globalFaceIds` pairs periodic faces correctly. The TGV case needs it.

### 2.4 Data layout and geometry
- Conserved state: `o_U[N_EQ * fieldOffset]`. Primitives: `o_rho, o_vel(3), o_p, o_T, o_h, o_gamma, o_c, o_Y(nspec), o_s …`.
- Kernel defines: `p_Nspec`, `p_Neq`, `p_Npa`, plus the mesh defines `p_Nq`, `p_Np`, `p_Nfp`.
- Thermo model, a compile-time branch (`p_perfectGas`), chosen with `[DGCOMP] thermo`:
  - `perfectGas` (default): one component, constant `gamma` and `gasConstant`, closed-form T. nspec = 1, so field 4 is rho. This makes the plugin a general single-component compressible DG solver.
  - `nasa`: the CompNek5000 multi-species mixture with NASA cp fits and Newton for T.
  - Both share the same state layout, so the flux kernels are written once.
- Geometry:
  - volume metrics from `o_vgeo` (RX..TZ, JID, JWID);
  - face normals and area weights from `o_sgeo` (NXID.., WSJID);
  - derivative matrix `o_D`, GLL weights `o_gllw`;
  - boundary IDs from `o_EToB`.
- Parameters:
  - read them from a par section or `userParam` instead of `uparam(1..8)`;
  - replace the `ifviscous`, `ifavisc`, … flags with options.
- **Use FP64 only.** Do not build or run the `nekrs-lib-fp32` variant for this.

### 2.5 nekRS is 3D only
- `dmr2d` and `det2d` must be extruded to one element thick in z, with periodic or symmetry z-faces.
- Compare them against the 2D Nek5000 runs.
- The 3D validation case is compHIT (§2.6).

### 2.6 Validation cases
**Primary 3D case: `Simulations/CompNek5000/compHIT/compHIT_direct_Re_lamda_100`.** It is decaying compressible HIT:
- Re_λ0 = 100, Mt0 = 0.6.
- Air (O₂/N₂ = 0.23/0.77) at 300 K and 1 bar.
- 32³ periodic elements on [0,1]³, rescaled to L = 1.6467e-4 m in `usrdat2`.
- dt = 1.5e-10, endTime = 5e-7 (about 4.5 τ0).

Why it is preferred over `tgvComp`:
- It has published reference curves: Bielawski Fig. A.15 in `ref/*.csv`, covering K/K0, enstrophy, dilatation and temperature variance.
- It already has Nek5000 runs with both LLF (`userparam05=0`) and HLLC (`userparam05=2`).
- At Mt0 = 0.6 it develops shocklets, so it exercises dilatation and Riemann-solver robustness. It already runs without AV or the positivity limiter.
- Each physical statistic is volume-weighted with `bm1` (`stat_ref.dat`), so it maps directly to nekRS mass-matrix reductions.
- 32³ elements is already a reasonable single-node GPU problem. It is also the natural weak-scaling base (64³, 128³).

What porting it requires:
- **IC:** read `hit_modes.dat` on the host (3576 mode pairs) and evaluate the Fourier sum at the GLL nodes in a small kernel or on the host. No `.re2` changes are needed: nekRS supports `usrdat2`/`rescale_x` in the `.usr` file, as `examples/ktauChannel` shows.
- **Statistics:** port `calc_stat` (the 12-column `stat_ref.dat`) to a UDF routine, using `o_LMM` weights and SEM gradients.
- **Fluxes:** HLLC (`riemann.f`) moves into the inviscid milestone so that both flux paths can be compared.

**Optional bring-up case: `tgvComp`.**
- It has an analytic IC and needs no input files.
- Use it only if HIT is hard to debug.

**Reacting cases:**
- `det2d`, extruded to one element thick, for comparison with Nek5000.
- Then a genuinely 3D detonation, which is the end goal.

---

## 3. Cantera strategy (three phases)

Cantera is a host-only C++ library and cannot run inside OCCA kernels. Calling it per point on every RK stage would also need a device↔host round trip each time.

**Phase A — host Cantera (correctness first).**
- Port `ftnlib.cpp` to a small C++ class (`dgCantera`) that holds `Solution`, `ThermoPhase`, `Kinetics` and `Transport`, with no Fortran name mangling.
- For μ, λ, D_k and ω̇: copy Y, u and ρ to the host, loop with Cantera, then copy back.
- This is slow but bit-comparable to Nek5000. It is the reference for Phases B and C.
- Also use Cantera at setup to *generate* `aa`, `h0`, `s0` and `mw`, replacing the `.dat` files, and to set the initial condition (`setState_TPY`, as `tgvComp/base.usr:112` does).

**Phase B — GPU transport.**
- Cantera's mixture-averaged model uses per-species polynomial fits in ln T (viscosity, conductivity, binary diffusion) combined with Wilke and mixture-averaged mixing rules.
- At setup, extract the fit coefficients from the Cantera `GasTransport` object.
  - Check the accessor names for your Cantera version (e.g. `getViscosityPolynomial`, `getConductivityPolynomial`, `getBinDiffusivityPolynomial`).
  - Upload the coefficients to the device and evaluate the mixing rules in an OKL kernel.
- Validate against Phase A pointwise.

**Phase C — GPU kinetics.**
- At setup, walk `Kinetics` (elementary/Arrhenius, three-body, falloff/Troe, reversibility through equilibrium constants from the NASA polynomials).
- From it, **generate OKL source** that is JIT-compiled with the mechanism baked in. OCCA's JIT is a natural fit for this.
- Keep the explicit RK2 sub-cycling from `integrate_chem` first, which ports trivially. Add point-implicit or Rosenbrock later if stiffness demands it.
- **Do not depend on nekRS CVODE.** `cmake/cvode.cmake` rejects HIP and DPCPP, so it would not work on Aurora.
- Before writing a generator, ask the nekRS developers whether there is existing GPU chemistry code from the nekRS reacting-flow work that can be reused.

**Linking Cantera into the build.**
1. *Quick path, no CMake changes:*
   - set `NEKRS_UDF_INCLUDES` and `NEKRS_UDF_LIBS` to point at the Cantera install, or use a case-local `udf.cmake`;
   - only the UDF/host code includes Cantera headers.
2. *Proper path:*
   - add `ENABLE_CANTERA` and `CANTERA_DIR`, modelled on the ADIOS pattern (`cmake/adios.cmake`, `CMakeLists.txt:291-296`);
   - `target_link_libraries(nekrs-lib … cantera)` and `-DNEKRS_ENABLE_CANTERA`;
   - when Cantera is absent, provide stub implementations that abort, the same idea as `cantera_stubs.f`.
3. *Building Cantera itself on the HPC machines:*
   - build it from source with SCons (shared library, C++17) using the same compiler and MPI as nekRS:
     - Polaris: `PrgEnv-gnu`;
     - Aurora: oneAPI `icpx`, or GNU with a compatible libstdc++.
   - Dependencies are yaml-cpp, fmt, Eigen, and Boost headers depending on the Cantera version. Check that version's compile docs.
   - Avoid relying on a conda env on compute nodes. The current `Cantera.mak` points to `~/.conda/envs/ct-dev`.

---

## 4. Milestones

The status of each item, plus case ideas and bug reports, is tracked in [TRACKER.md](TRACKER.md).

| # | Goal | Done when |
|---|---|---|
| 0 | Baseline | A separate `dg_dev` install (§5) builds and runs `examples/ethier` on Polaris (CUDA) and Aurora (DPCPP). Module and env scripts are committed under `doc/dg/`. |
| 1 | Inviscid, non-reacting skeleton, no Cantera | `dgComp` plugin: primitives kernel (Newton T), face exchange, KEP/EC volume flux, LLF **and HLLC** surface fluxes, SSPRK3, CFL; periodic box. `aa`/`h0`/`s0`/`mw` are read from files. Tests: EC flux with no dissipation conserves total entropy to round-off; isentropic vortex convergence; one-step diff against CompNek5000 on the compHIT mesh, N and IC with viscous terms off. |
| 2 | Viscous + host Cantera transport → **compHIT** | BR1 viscous kernels, the `dgCantera` class (μ, λ, D_k), the HIT IC reader and `stat_ref.dat` statistics. Test: `compHIT_direct_Re_lamda_100`, LLF and HLLC, against the Nek5000 `stat_ref.dat` and Bielawski `ref/*.csv`. |
| 3 | Shocks and BCs | `compbc` kernel, positivity limiter, entropy viscosity. Test: extruded `dmr2d` against Nek5000. |
| 4 | Reacting, host Cantera chemistry | Strang splitting with `integrate_chem`. Test: extruded `det2d` (H₂, 10 species). |
| 5 | GPU transport and kinetics | Phases B and C. Kernel-level comparison with Phase A. HIT at 64³ and 128³ as a weak-scaling check. |
| 6 | 3D detonation and performance | First genuinely 3D detonation run. Fused kernels, shared-memory volume flux, overlapped exchange, scaling on Aurora and Polaris. |

Milestone 2 needs Cantera on Polaris and Aurora, so the Cantera source builds (§3, item 3) must be done by the end of Milestone 1.

---

## 5. Separate nekRS install for this work

The other project's nekRS stays untouched. The default `CMAKE_INSTALL_PREFIX` is `$HOME/.local/nekrs` (`CMakeLists.txt:10-11`). If the existing build used that default, a plain `./build.sh` here would **overwrite it**. Always pass a prefix of its own, one per machine:

```bash
# Polaris
CC=cc CXX=CC FC=ftn ./build.sh -DCMAKE_INSTALL_PREFIX=$HOME/.local/nekrs-dg-polaris
# Aurora
CC=mpicc CXX=mpicxx FC=mpif77 ./build.sh -DCMAKE_INSTALL_PREFIX=$HOME/.local/nekrs-dg-aurora
```

- Set `NEKRS_HOME` to the matching prefix in a per-machine `env` script.
  - Never put both installs on `PATH` at once.
  - Check with `which nekrs` before running.
- The `build/` directory inside this clone is separate from the other project's clone, so it does not conflict.
- The compiler wrappers above are typical defaults. The exact modules for each machine still need to be confirmed and recorded in `doc/dg/`.

### Verification harness (build it in milestone 1)
- Add a debug dump to CompNek5000 (`userchk`) that writes the state, primitives and RHS after a single `getdgrhs` call.
- Write the same quantities from nekRS.
- Compare them with a small script.
- Use this for every kernel added, not only full runs.
