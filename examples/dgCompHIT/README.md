# dgCompHIT: compressible HIT with the dgComp plugin

This is the nekRS port of the CompNek5000 case
`compHIT/compHIT_direct_Re_lamda_100/LLF`:
- Initial Re_lambda = 100 and Mt = 0.6.
- Air (O2/N2 = 0.23/0.77) at 300 K and 1e5 Pa.
- Box length L = 1.6467320829664873e-4 m, 32^3 elements, N = 7.
- dt = 1.5e-10 s, endTime = 5e-7 s.

It is the main validation case for the DG port (see `doc/dg/ROADMAP.md`).

**Status: M0 scaffold.**
- The volume and surface RHS kernels are still stubs that return zero.
- A run should therefore leave the initial field unchanged.
- This exercises setup, the EOS/Newton primitives, the face exchange, SSPRK3 and checkpointing, not the physics.
- Viscous terms (M2) are not implemented, so `viscous = no` for now.

## Parameters

The par file has two user sections:
- `[DGCOMP]` is read by the plugin: species, flux choices, feature switches, and CompNek5000 `userparam01..05` (`dtChem`, `avCstar`, `avCmax`, `avP`, `riemannSolver`).
- `[CASEDATA]` is read by `dgCompHIT.udf`: `boxLength`, `T0`, `p0` and `Y0` (one mass fraction per species, summing to 1). These used to be hard-coded in `base.usr`.

## Files to copy in from the CompNek5000 case

| file | source |
|---|---|
| `hit_modes.dat` | `compHIT_direct_Re_lamda_100/LLF/hit_modes.dat` (same file for LLF and HLLC) |
| `aa.dat`, `h0.dat`, `s0.dat` | thermo refit files the Fortran `usrdat` reads (air1.yaml, O2 then N2) |
| `mw.dat` | molar masses in g/mol, one per line, same species order (Cantera `getmolecularweights`) |

`mw.dat` is temporary. It goes away once the Cantera host bridge (M2) supplies the molar masses.

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
