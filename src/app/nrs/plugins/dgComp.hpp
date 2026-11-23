#if !defined(nekrs_dgComp_hpp_)
#define nekrs_dgComp_hpp_

// Compressible multi-species DGSEM plugin (port of CompNek5000).
// See doc/dg/ROADMAP.md for the design and milestones.
//
// The plugin owns its own state on meshV. The nekRS fluid solver should be
// disabled in the .par file:
//
//   userSections = DGCOMP
//   [FLUID VELOCITY]
//   solver = none
//
//   [DGCOMP]
//   nspec = 2                 # number of species
//   npa = 5                   # NASA polynomial coefficients per species
//   riemannSolver = llf       # llf | hllc      (CompNek5000 userparam05 = 0 | 2)
//   volumeFlux = kep          # kep | ec        (ifVolEC)
//   viscous = yes             # ifviscous
//   reacting = no             # ifforces
//   splitChemistry = no       # ifsplitchem
//   clip = no                 # ifclip
//   positivityLimiter = no    # ifpplim
//   artificialViscosity = no  # ifavisc
//   fullInterfaceLLF = no     # ifIntLLFfull
//   species = O2, N2          # optional, names the Y_k output fields
//   dtChem = 2e-9             # chemistry sub-step      (userparam01), reacting only
//   avCstar = 1.0             # entropy viscosity c*    (userparam02)
//   avCmax = 0.3              # entropy viscosity cmax  (userparam03)
//   avP = 0.7                 # entropy viscosity P     (userparam04)
//
// Case-specific inputs (box size, initial state, ...) belong in the case's
// own [CASEDATA] section and are read by the .udf, not by the plugin.
//
// The kernels are compiled with the species count baked in (p_Nspec), so
// buildKernel reads [DGCOMP] nspec/npa before setup() is called.
//
// Conserved-variable layout, field stride = nrs->fieldOffset:
//   0,1,2      rho*u, rho*v, rho*w   (CompNek5000 U_EQ..W_EQ)
//   3          rho*E                 (E_EQ)
//   4 .. 4+ns  rho*Y_k               (E_EQ+k)
//
// Primitive layout (see dgComp::prim):
//   rho, u, v, w, p, T, gamma, a, h(=(rhoE+p)/rho), W(mean molar mass), Y_k ...

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

namespace dgComp
{

namespace prim
{
constexpr int RHO = 0;
constexpr int U = 1;
constexpr int V = 2;
constexpr int W = 3;
constexpr int P = 4;
constexpr int T = 5;
constexpr int GAMMA = 6;
constexpr int ASPD = 7;
constexpr int H = 8;
constexpr int WMEAN = 9;
constexpr int Y0 = 10; // first species mass fraction
} // namespace prim

// Thermodynamic data per species, same layout as CompNek5000 (COMPRESS):
//   aa[k*npa + c]  cp refit, cp_k(T) = sum_{c=0}^{npa-1} aa[k*npa+c] T^c  [J/(kg K)]
//   h0[k]          enthalpy at 298.15 K                                   [J/kg]
//   s0[k]          entropy at 298.15 K, 101325 Pa                         [J/(kg K)]
//   mw[k]          molar mass                                             [kg/mol]
struct thermo_t {
  int npa = 5;
  std::vector<dfloat> aa;
  std::vector<dfloat> h0;
  std::vector<dfloat> s0;
  std::vector<dfloat> mw;
};

// Read the CompNek5000 case files from caseDir: aa.dat, h0.dat, s0.dat
// (list-directed, same files the Fortran usrdat reads) and mw.dat (one molar
// mass per species in g/mol, as Cantera's getmolecularweights returns).
// TODO(M2): take mw from the Cantera host bridge instead of mw.dat.
thermo_t readThermo(const std::string &caseDir, int nspec, int npa = 5);

void buildKernel(occa::properties kernelInfo);

void setup(const thermo_t &thermo);

// Host conserved state, size Neq()*nrs->fieldOffset, layout as above.
void setSolution(const std::vector<dfloat> &U);

// Build conserved variables from host primitives on meshV nodes:
// velocity (3*fieldOffset), p, T (fieldOffset each), Y (nspec*fieldOffset).
std::vector<dfloat> conservedFromPrimitives(const std::vector<dfloat> &vel,
                                            const std::vector<dfloat> &p,
                                            const std::vector<dfloat> &T,
                                            const std::vector<dfloat> &Y);

// Advance o_U by one SSPRK3 step of size dt.
void step(double time, dfloat dt);

// Recompute primitives from o_U (Newton iteration for T).
void computePrimitives();

// Largest |u| + a on this mesh, reduced over all ranks.
dfloat maxWaveSpeed();

// Register the primitives as checkpoint fields. velocity, pressure and
// temperature always keep those names. Density and the mass fractions are
// named only with checkpointEngine = adios, because the nek (.fld) writer
// accepts nothing but scalarNN:
//   adios:     density, Y_<species> ...
//   nek .fld:  scalar00 = density, scalar01.. = Y_k (ParaView: S01, S02, ...)
// speciesNames overrides [DGCOMP] species; with neither, Y_00, Y_01, ...
// The mapping is printed at setup.
void addCheckpointFields(const std::vector<std::string> &speciesNames = {});

int Nspec();
int Neq();
int Nprim();
occa::memory o_U();
occa::memory o_prim();

} // namespace dgComp

#endif
