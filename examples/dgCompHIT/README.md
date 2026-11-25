# dgCompHIT: compressible HIT with the dgComp plugin

This is the nekRS port of the CompNek5000 case
`compHIT/compHIT_direct_Re_lamda_100/LLF`:
- Initial Re_lambda = 100 and Mt = 0.6.
- Air (O2/N2 = 0.23/0.77) at 300 K and 1e5 Pa.
- Box length L = 1.6467320829664873e-4 m, 32^3 elements, N = 7.
- dt = 1.5e-10 s, endTime = 5e-7 s.

It is the main validation case for the DG port (see `doc/dg/ROADMAP.md`).

**Status: M1, inviscid.**
- The inviscid RHS is in: KEP volume flux, and the LLF, HLL, HLLC or blended HLLC/HLL interface flux (`[DGCOMP] riemannSolver`). Not compiled or run yet.
- Viscous terms (BR1 and the Cantera transport properties, M2) are not implemented, so `viscous = no` for now.
  - The run is therefore inviscid. It does **not** reproduce the viscous Nek5000 reference (urms decays to 58.7 m/s at 4.5 tau0); only numerical dissipation removes energy here.
  - To compare now, run CompNek5000 with `ifviscous = .false.` and the same flux settings (`ifVolEC = .false.`, `uparam(5)` = `riemannSolver`, `ifIntLLFfull` = `fullInterfaceLLF`, same dt). Or compare the step-0 RHS, see `doc/dg/IMPLEMENTATION.md` §5.
- No CFL control yet: dt is fixed from the par file.
- See `doc/dg/IMPLEMENTATION.md` for which Fortran routines are ported.

## Parameters

The par file has two user sections:
- `[DGCOMP]` is read by the plugin: thermo model, species, flux choices, feature switches, and CompNek5000 `userparam01..08` (`dtChem`, `avCstar`, `avCmax`, `avP`, `riemannSolver`, `blendSw0`, `blendSw1`, `blendSwRho`).
- `[CASEDATA]` is read by `dgCompHIT.udf` in `UDF_Setup0`, with plain `platform->par->extract` calls as in the other examples. These used to be hard-coded in `base.usr`. `checkCaseData()` checks them in `UDF_Setup`, once the plugin knows nspec.

| key | meaning | default |
|---|---|---|
| `boxLength` | box length L [m] (`base.usr` L0) | required |
| `T0`, `p0` | initial temperature [K] and pressure [Pa] | required |
| `Y0` | initial mass fractions, same order as `species` | required for nspec > 1 |
| `mu0` | dynamic viscosity [Pa s], only used for Re_lambda in `stat.dat` | 0 (Re_lambda written as 0, with a warning) |
| `statInterval` | steps between `stat.dat` / `stat_ref.dat` rows (`base.usr` sfreq) | 1 |
| `consInterval` | steps between `cons.dat` rows (`base.usr` cfreq) | 10 |

## Statistics output

The udf ports `calc_stat`, `hit_ref_stat` and `check_cons` from `base.usr`. The files have the same columns and number formats as the CompNek5000 run, so `plot_hit_stats.py` and the Bielawski Fig. A.15 comparison read them unchanged.

- Averages are volume averages, sum(Jw f) / volume, the same as `glsc2(f,bm1,n)/volvm1`.
- Velocity gradients are element-local, as with `gradm1`.
- Variances are computed in two passes, about the mean.
- Everything is computed at tstep 0 from the initial condition, which also sets the reference scales. The scales are printed at startup.

| file | written | columns |
|---|---|---|
| `stat.dat` | every `statInterval` steps, from step 1 | `time/tau0  <\|u\|^2>  urms  lambda  Re_lambda  prms  Trms` (E13.5, 6 x E15.7) |
| `stat_ref.dat` | header and a row at step 0, then every `statInterval` steps | `time  t/tau0  K/K0  omega2*tau0^2  Tvar/Tscale^2  theta2*tau0^2  K  omega2  Tvar  theta2  theta_mean  theta_var` (12 x ES24.16E3). The header line `# scales:` holds K0, tau0, Tscale, Mt0, gamma0, c0 |
| `cons.dat` | every `consInterval` steps, from step 1 | `time  (m0-m)/m0  (e0-e)/e0` (E13.5, 2 x E15.7) |

Definitions, as in `base.usr`:
- urms = sqrt(<\|u - <u>\|^2> / 3); lambda = sqrt(urms^2 / <(du/dx)^2>); Re_lambda = urms lambda <rho> / mu0.
- K = 1.5 urms^2; omega2 = <\|curl u\|^2> (no factor 1/2); theta = div u, theta2 = <theta^2>; Tvar = Trms^2.
- tau0 = lambda/urms, K0, c0 = <a>, gamma0 = <gamma>, Mt0 = sqrt(3) urms / c0 and Tscale = (gamma0 - 1) <T> Mt0^2, all at step 0.
- m and e are the integrals of rho and rho E.

Until the transport properties arrive (M2) the viscosity is not computed, so Re_lambda uses the constant `mu0`. The par file sets the Cantera value for air at 300 K and 1e5 Pa, which is what the CompNek5000 run reports at step 0. The other columns do not depend on mu.

### Thermo model

| `thermo` | gas | needs | Y0 |
|---|---|---|---|
| `perfectGas` (plugin default) | one component, constant `gamma` and `gasConstant` | nothing extra | ignored |
| `nasa` (this case) | mixture, NASA cp fits per species | `nspec`, `aa/h0/s0/mw.dat` | one mass fraction per species, summing to 1 |

The par file uses `nasa` to match the CompNek5000 reference. Switching to `perfectGas` with gamma = 1.4 and R = 287 gives single-component air without any `.dat` files.

## Files to copy in from the CompNek5000 case

The `.dat` files are only needed with `thermo = nasa`.

| file | source |
|---|---|
| `hit_modes.dat` | `compHIT_direct_Re_lamda_100/LLF/hit_modes.dat` (same file for LLF and HLLC) |
| `aa.dat`, `h0.dat`, `s0.dat` | thermo refit files the Fortran `usrdat` reads (air1.yaml, O2 then N2) |
| mesh | `dgCompHIT.re2`, see [Mesh](#mesh) |
`mw.dat` is already in this directory: molar masses in g/mol, one per line, O2 then N2 (31.998, 28.014, the values Cantera's `getmolecularweights` returns). Its order must match `aa.dat`, `h0.dat` and `s0.dat`. To check the order, or to write it for another mechanism:

```bash
python -c "import cantera as ct; g = ct.Solution('air1.yaml'); print(g.species_names); print(*g.molecular_weights, sep='\n')"
```

The CompNek5000 run also prints the order at startup (`species       mol wt`). `mw.dat` is temporary. It goes away once the Cantera host bridge (M2) supplies the molar masses. No `.yaml` file is read at run time.

## Mesh

```bash
genbox <<< input.box
```

```bash
gencon dgCompHIT.re2    # optional; nekRS can also build connectivity with parRSB
```

`genbox` writes `box.re2`. Rename it to `dgCompHIT.re2`. The udf then rescales [0,1]^3 to [0,L]^3 and recomputes the geometric factors.

## Run

Use the separate dg_dev install, not the other project's nekRS:

```bash
export NEKRS_HOME=$HOME/.local/nekrs-dg
```

```bash
nekrs --setup dgCompHIT --build-only <nranks>
```

```bash
mpiexec -n <nranks> nekrs --setup dgCompHIT
```

## Checkpoint contents

The par file sets `checkpointEngine = adios`. nekRS then writes one ADIOS2 BP directory instead of `.f0000x` files, and the fields carry real names:

| field | content |
|---|---|
| `velocity` | u, v, w |
| `pressure` | p |
| `temperature` | T |
| `density` | rho |
| `Y_O2`, `Y_N2` | mass fractions, named from `[DGCOMP] species` |

- ADIOS2 is built into nekRS by default (`ENABLE_ADIOS=ON`, bundled in `3rd_party/adios`). With `-DENABLE_ADIOS=OFF` the par file is rejected.
- Open the `.bp` directory in ParaView with the **ADIOS2VTXReader**.
- The setup log prints the field list.

With `checkpointEngine = nek` (the nekRS default) the run writes `.fld` files instead. That format only allows the names velocity, pressure, temperature and `scalarNN`, so density and the mass fractions are numbered. The setup log prints the mapping:

| .fld field | ParaView name | content |
|---|---|---|
| `scalar00` | S01 | density |
| `scalar01` | S02 | Y_O2 |
| `scalar02` | S03 | Y_N2 |

If `species` is left out, the mass fractions are named `Y_00`, `Y_01`, ...

With `thermo = perfectGas` there is one component, so no mass fraction is written; the scalars are density only.
