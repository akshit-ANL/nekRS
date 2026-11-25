#!/usr/bin/env python3
"""Export Cantera's mixture-averaged transport fits for the dgComp plugin.

The plugin ([DGCOMP] transport = fits) evaluates Cantera's MixTransport model
on the GPU from the same polynomial fits Cantera uses internally, so no
Cantera library is needed at run time. This script writes, in the species
order of the mechanism:

  visc.dat  nspec rows of 5: sqrt(mu_k) / T^(1/4)        as a polynomial in ln T
  cond.dat  nspec rows of 5: lambda_k / sqrt(T)          as a polynomial in ln T
  diff.dat  nspec*nspec rows of 5 (full matrix, row j*nspec+k):
            p D_jk / T^(3/2)                             as a polynomial in ln T
  mw.dat    molar masses [g/mol]

Coefficients are in ascending powers of ln T, SI units. The species order must
match aa.dat, h0.dat, s0.dat and [DGCOMP] species.

Usage:
  python exportTransport.py air1.yaml [--phase NAME] [--out DIR]
  python exportTransport.py air1.yaml --check --T 300 --p 1e5 --Y O2:0.23,N2:0.77

--check evaluates the formulas of dgComp.okl::dgTransport in Python and
compares them with Cantera's own viscosity, thermal_conductivity and
mix_diff_coeffs at the given state.
"""

import argparse
import math
import os
import sys

import numpy as np
import cantera as ct


def load(mech, phase):
    gas = ct.Solution(mech, phase) if phase else ct.Solution(mech)
    for model in ("mixture-averaged", "Mix"):
        try:
            gas.transport_model = model
            break
        except Exception:
            continue
    else:
        sys.exit("cannot set the mixture-averaged transport model")
    if getattr(gas, "CK_mode", False):
        sys.exit("CK-mode transport fits are not supported (the plugin uses the Cantera default, non-CK)")
    return gas


def fits(gas):
    ns = gas.n_species
    visc = np.array([gas.get_viscosity_polynomial(k) for k in range(ns)])
    cond = np.array([gas.get_thermal_conductivity_polynomial(k) for k in range(ns)])
    diff = np.zeros((ns, ns, 5))
    for j in range(ns):
        for k in range(ns):
            diff[j, k] = gas.get_binary_diff_coeffs_polynomial(min(j, k), max(j, k))
    for name, a in (("viscosity", visc), ("conductivity", cond), ("diffusion", diff)):
        if a.shape[-1] != 5:
            sys.exit(f"{name} fits have {a.shape[-1]} coefficients, the plugin expects 5")
    return visc, cond, diff


def write(path, rows):
    with open(path, "w") as f:
        for r in rows:
            f.write(" ".join("%.17e" % c for c in r) + "\n")


def kernel_eval(gas, visc, cond, diff, T, p, Y):
    """Python copy of dgComp.okl::dgTransport, mode 1."""
    mw = gas.molecular_weights
    ns = gas.n_species
    W = 1.0 / np.sum(Y / mw)
    X = np.maximum(1e-20, Y * W / mw)
    lt = np.array([1.0, math.log(T), math.log(T) ** 2, math.log(T) ** 3, math.log(T) ** 4])

    sqmu = T ** 0.25 * visc @ lt
    muk = sqmu ** 2
    lamk = math.sqrt(T) * cond @ lt
    Djk = T ** 1.5 * diff @ lt / p

    mu = 0.0
    for k in range(ns):
        den = 0.0
        for j in range(ns):
            f1 = 1 + (sqmu[k] / sqmu[j]) * (mw[j] / mw[k]) ** 0.25
            den += X[j] * f1 * f1 / math.sqrt(8 * (1 + mw[k] / mw[j]))
        mu += X[k] * muk[k] / den
    lam = 0.5 * (np.sum(X * lamk) + 1 / np.sum(X / lamk))

    D = np.zeros(ns)
    for k in range(ns):
        D[k] = Djk[k, k]
        if ns > 1:
            s = sum(X[j] / Djk[j, k] for j in range(ns) if j != k)
            if s > 0:
                D[k] = (W - X[k] * mw[k]) / (W * s)
    return mu, lam, D


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mech", help="Cantera mechanism file (.yaml)")
    ap.add_argument("--phase", default=None, help="phase name in the mechanism")
    ap.add_argument("--out", default=".", help="output directory (default .)")
    ap.add_argument("--check", action="store_true", help="compare the kernel formulas with Cantera")
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--p", type=float, default=1e5)
    ap.add_argument("--Y", default=None, help="mass fractions, e.g. O2:0.23,N2:0.77")
    args = ap.parse_args()

    gas = load(args.mech, args.phase)
    visc, cond, diff = fits(gas)
    ns = gas.n_species
    print("species (order of visc/cond/diff/mw.dat):", ", ".join(gas.species_names))

    if args.check:
        Y = args.Y if args.Y else {gas.species_names[0]: 1.0}
        gas.TPY = args.T, args.p, Y
        mu, lam, D = kernel_eval(gas, visc, cond, diff, args.T, args.p, gas.Y)
        rel = lambda a, b: abs(a - b) / abs(b)
        print(f"T = {args.T} K, p = {args.p} Pa, Y = {gas.Y}")
        print(f"mu     kernel {mu:.10e}  cantera {gas.viscosity:.10e}  rel {rel(mu, gas.viscosity):.2e}")
        print(f"lambda kernel {lam:.10e}  cantera {gas.thermal_conductivity:.10e}  "
              f"rel {rel(lam, gas.thermal_conductivity):.2e}")
        for k, name in enumerate(gas.species_names):
            Dc = gas.mix_diff_coeffs[k]
            print(f"D_{name:<6} kernel {D[k]:.10e}  cantera {Dc:.10e}  rel {rel(D[k], Dc):.2e}")
        return

    os.makedirs(args.out, exist_ok=True)
    write(os.path.join(args.out, "visc.dat"), visc)
    write(os.path.join(args.out, "cond.dat"), cond)
    write(os.path.join(args.out, "diff.dat"), diff.reshape(ns * ns, 5))
    with open(os.path.join(args.out, "mw.dat"), "w") as f:
        for w in gas.molecular_weights:
            f.write("%.17g\n" % w)
    print(f"wrote visc.dat, cond.dat, diff.dat, mw.dat to {args.out}")


if __name__ == "__main__":
    main()
