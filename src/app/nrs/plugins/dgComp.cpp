#include <algorithm>
#include <cctype>
#include <fstream>
#include <iterator>
#include <sstream>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "linAlg.hpp"

#include "dgComp.hpp"

// Status (M2): primitives (perfect gas, or NASA mixture with Newton for T),
// trace gather and face exchange, KEP flux-differencing volume term, LLF /
// HLL / HLLC / blended interface flux, BR1 viscous terms with constant or
// Cantera-fit (mixture-averaged) transport, SSPRK3 and the wave-speed
// reduction. Not yet: EC volume flux, boundary conditions, chemistry,
// positivity limiter, artificial viscosity.
// See doc/dg/ROADMAP.md.

namespace
{

constexpr dfloat Ru = 8.3144598; // J/(mol K), same value as CompNek5000

nrs_t *_nrs = nullptr;
mesh_t *mesh = nullptr;

// [DGCOMP] thermo = perfectGas (single component, default) | nasa (mixture)
bool isPerfectGas = true;
dfloat gammaGas = 1.4; // perfectGas only
dfloat Rgas = 287.0;   // perfectGas only, J/(kg K)

int nspec = 0;
int npa = 5;
int nEq = 0;
int nPrim = 0;
int nTrace = 0;
// BR1: gradient variables u, v, w, T (+ X_k for nspec > 1) and viscous-flux
// rows (3 momentum, energy, + rho Y_k for nspec > 1); same count
int nVisc = 0;

dlong fieldOffset = 0;
dlong faceOffset = 0;

struct options_t {
  // 0 LLF, 1 HLL, 2 HLLC, 3 HLLC blended into HLL (CompNek5000 userparam05)
  int riemannSolver = 0;
  int volumeFlux = 0; // 0 KEP (Pirozzoli), 1 EC (Gouasmi)
  // riemann.f blend switch: alpha ramps 0 -> 1 as |dp|/(pL+pR) goes from
  // blendSw0 to blendSw1, damped by |drho|/(rhoL+rhoR) / blendSwRho at
  // subsonic faces
  dfloat blendSw0 = 0.01;
  dfloat blendSw1 = 0; // <= blendSw0 means 10 blendSw0
  dfloat blendSwRho = 0.1;
  bool viscous = false;
  bool reacting = false;
  bool splitChemistry = false;
  bool clip = false;
  bool positivityLimiter = false;
  bool artificialViscosity = false;
  bool fullInterfaceLLF = false;

  // viscous = yes: transport = constant (mu, prandtl, lewis) or fits
  // (Cantera mixture-averaged fits from visc.dat, cond.dat, diff.dat)
  int transport = 0; // 0 constant, 1 fits
  dfloat mu = 0;
  dfloat prandtl = 0.72;
  dfloat lewis = 1.0;

  // CompNek5000 userparam01..04 (base.par of compHIT)
  dfloat dtChem = 0;    // uparam(1), chemistry sub-step (M4)
  dfloat avCstar = 1.0; // uparam(2), entropy-viscosity c*  (M3)
  dfloat avCmax = 0.3;  // uparam(3), entropy-viscosity cmax (M3)
  dfloat avP = 0.7;     // uparam(4), entropy-viscosity P   (M3)
} opt;

occa::memory o_Ucons; // conserved state
occa::memory o_U0;    // state at the start of the step
occa::memory o_rhs;   // dU/dt
occa::memory o_primv; // primitives
occa::memory o_faceM; // own trace, nTrace fields of faceOffset
occa::memory o_faceP; // M + P after the face exchange
occa::memory o_faceFlux; // lifted interface flux, nEq fields of faceOffset
occa::memory o_lambda;

occa::memory o_bi; // internal-energy polynomial, nspec*(npa+1)
occa::memory o_aa; // cp polynomial, nspec*npa
occa::memory o_mw;

oogs_t *gshFace = nullptr;

// viscous terms (allocated only with viscous = yes)
occa::memory o_trans;  // mu, lambda, rho D_k: (2 + nspec) fields
occa::memory o_grad;   // BR1 gradients, 3*nVisc fields, row 3g+d
occa::memory o_vflux;  // viscous fluxes, 3*nVisc fields, row 3q+d
occa::memory o_vfaceM; // F_q.n_M on the trace, nVisc fields of faceOffset
occa::memory o_vfaceP; // F_M.n_M + F_P.n_P after the exchange
occa::memory o_viscFit; // transport fits, 5 per species (pair)
occa::memory o_condFit;
occa::memory o_diffFit;
oogs_t *gshVisc = nullptr;

occa::kernel primitivesKernel;
occa::kernel conservedKernel;
occa::kernel traceGatherKernel;
occa::kernel volumeRhsKernel;
occa::kernel surfaceFluxKernel;
occa::kernel surfaceLiftKernel;
occa::kernel ssprkStageKernel;
occa::kernel waveSpeedKernel;
occa::kernel transportKernel;
occa::kernel viscGradKernel;
occa::kernel viscFluxKernel;
occa::kernel viscTraceGatherKernel;
occa::kernel viscDivKernel;

bool buildKernelCalled = false;
bool setupCalled = false;

std::string toLower(std::string s)
{
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
  return s;
}

bool parFlag(const std::string &key, bool defaultValue)
{
  std::string val;
  if (!platform->par->extract("dgcomp", key, val)) {
    return defaultValue;
  }
  val = toLower(val);
  return val == "yes" || val == "true" || val == "on" || val == "1";
}

std::string parString(const std::string &key, const std::string &defaultValue)
{
  std::string val;
  if (!platform->par->extract("dgcomp", key, val)) {
    return defaultValue;
  }
  return toLower(val);
}

void readOptions()
{
  const auto riemann = parString("riemannsolver", "llf");
  const std::vector<std::string> riemannNames = {"llf", "hll", "hllc", "hllcblend"};
  const auto it = std::find(riemannNames.begin(), riemannNames.end(), riemann);
  nekrsCheck(it == riemannNames.end(),
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "[DGCOMP] riemannSolver = %s is not supported (llf | hll | hllc | hllcBlend)\n",
             riemann.c_str());
  opt.riemannSolver = static_cast<int>(it - riemannNames.begin());

  // CompNek5000 userparam06..08, with the same fallbacks as riemann.f
  opt.blendSw1 = 0;
  platform->par->extract("dgcomp", "blendsw0", opt.blendSw0);
  platform->par->extract("dgcomp", "blendsw1", opt.blendSw1);
  platform->par->extract("dgcomp", "blendswrho", opt.blendSwRho);
  if (!(opt.blendSw0 > 0)) {
    opt.blendSw0 = 0.01;
  }
  if (!(opt.blendSw1 > opt.blendSw0)) {
    opt.blendSw1 = 10 * opt.blendSw0;
  }
  if (!(opt.blendSwRho > 0)) {
    opt.blendSwRho = 0.1;
  }

  const auto vflux = parString("volumeflux", "kep");
  nekrsCheck(vflux != "kep" && vflux != "ec",
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "[DGCOMP] volumeFlux = %s is not supported (kep | ec)\n",
             vflux.c_str());
  opt.volumeFlux = (vflux == "ec") ? 1 : 0;

  opt.viscous = parFlag("viscous", false);
  opt.reacting = parFlag("reacting", false);
  opt.splitChemistry = parFlag("splitchemistry", false);
  opt.clip = parFlag("clip", false);
  opt.positivityLimiter = parFlag("positivitylimiter", false);
  opt.artificialViscosity = parFlag("artificialviscosity", false);
  opt.fullInterfaceLLF = parFlag("fullinterfacellf", false);

  const auto transport = parString("transport", isPerfectGas ? "constant" : "fits");
  nekrsCheck(transport != "constant" && transport != "fits",
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "[DGCOMP] transport = %s is not supported (constant | fits)\n",
             transport.c_str());
  opt.transport = (transport == "fits") ? 1 : 0;
  platform->par->extract("dgcomp", "mu", opt.mu);
  platform->par->extract("dgcomp", "prandtl", opt.prandtl);
  platform->par->extract("dgcomp", "lewis", opt.lewis);
  nekrsCheck(opt.viscous && opt.transport == 0 && !(opt.mu > 0),
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "%s\n",
             "[DGCOMP] viscous = yes with transport = constant needs mu > 0");
  nekrsCheck(opt.viscous && opt.transport == 0 && (!(opt.prandtl > 0) || !(opt.lewis > 0)),
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "%s\n",
             "[DGCOMP] prandtl and lewis must be > 0");

  // keys missing from the par file keep the defaults above
  platform->par->extract("dgcomp", "dtchem", opt.dtChem);
  platform->par->extract("dgcomp", "avcstar", opt.avCstar);
  platform->par->extract("dgcomp", "avcmax", opt.avCmax);
  platform->par->extract("dgcomp", "avp", opt.avP);

  nekrsCheck(opt.reacting && opt.dtChem <= 0,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "%s\n",
             "[DGCOMP] reacting = yes needs dtChem > 0 (CompNek5000 userparam01)");
}

// [DGCOMP] species = O2, N2 (comma or space separated). The par reader
// lower-cases values, so the names are upper-cased back (O2, N2, H2O, OH).
std::vector<std::string> speciesNamesFromPar()
{
  std::vector<std::string> names;
  std::string val;
  if (!platform->par->extract("dgcomp", "species", val)) {
    return names;
  }
  std::replace(val.begin(), val.end(), ',', ' ');
  std::istringstream ss(val);
  for (std::string s; ss >> s;) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::toupper(c); });
    names.push_back(s);
  }
  return names;
}

std::vector<dfloat> readList(const std::string &fileName, size_t n)
{
  std::ifstream f(fileName);
  nekrsCheck(!f.is_open(), platform->comm.mpiComm(), EXIT_FAILURE, "cannot open %s\n", fileName.c_str());

  // Fortran list-directed files may use D exponents.
  std::string txt((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  std::replace(txt.begin(), txt.end(), 'D', 'E');
  std::replace(txt.begin(), txt.end(), 'd', 'e');
  std::replace(txt.begin(), txt.end(), ',', ' ');

  std::istringstream is(txt);
  std::vector<dfloat> v;
  double x;
  while (v.size() < n && is >> x) {
    v.push_back(x);
  }
  nekrsCheck(v.size() != n,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "%s: expected %zu values, found %zu\n",
             fileName.c_str(),
             n,
             v.size());
  return v;
}

// Port of compvar.f::get_bi. Internal energy of species k:
//   u_k(T) = sum_{c=0}^{npa} bi[k*(npa+1)+c] T^c
std::vector<dfloat> computeBi(const dgComp::thermo_t &th)
{
  std::vector<dfloat> bi(nspec * (npa + 1));
  for (int k = 0; k < nspec; k++) {
    const dfloat *a = th.aa.data() + k * npa;
    dfloat sum1 = 0;
    for (int c = 0; c < npa; c++) {
      sum1 += a[c] * std::pow(298.15, c + 1) / (c + 1);
    }
    dfloat *b = bi.data() + k * (npa + 1);
    b[0] = th.h0[k] - sum1;
    b[1] = a[0] - Ru / th.mw[k];
    for (int c = 2; c <= npa; c++) {
      b[c] = a[c - 1] / c;
    }
  }
  return bi;
}

// Single species with constant cp = gamma R/(gamma-1) and W = Ru/R, so the
// shared EOS code gives p = rho R T. bi is set directly in setupImpl.
dgComp::thermo_t perfectGasThermo()
{
  dgComp::thermo_t th;
  th.npa = 1;
  th.aa = {gammaGas * Rgas / (gammaGas - 1)};
  th.h0 = {0};
  th.s0 = {0};
  th.mw = {Ru / Rgas};
  return th;
}

void setupImpl(const dgComp::thermo_t &th); // defined after dgComp::setup

void checkBuilt()
{
  nekrsCheck(!buildKernelCalled,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "%s\n",
             "dgComp::buildKernel was not called (include dgComp.hpp in the .udf)");
}

bool checkpointFieldsAdded = false;

bool writeThroughFluid()
{
  return _nrs->fluid && platform->options.compareArgs(upperCase(_nrs->fluid->name) + " CHECKPOINTING", "TRUE");
}

void copyToFluid()
{
  if (!checkpointFieldsAdded || !writeThroughFluid()) {
    return;
  }
  using namespace dgComp;
  auto &fluid = _nrs->fluid;
  for (int d = 0; d < 3; d++) {
    fluid->o_U.copyFrom(o_primv, mesh->Nlocal, d * fluid->fieldOffset, (prim::U + d) * fieldOffset);
  }
  fluid->o_P.copyFrom(o_primv, mesh->Nlocal, 0, prim::P * fieldOffset);
}

void checkReady(const char *caller)
{
  nekrsCheck(!setupCalled || !buildKernelCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "dgComp::%s called prior to dgComp::setup()!\n",
             caller);
}

// dU/dt = volume + surface terms (advect.f::advect, inviscid part).
// Assumes o_primv is current for o_Ucons.
void computeRhs(double time)
{
  traceGatherKernel(mesh->Nelements, fieldOffset, faceOffset, mesh->o_vmapM, o_Ucons, o_primv, o_faceM);
  o_faceP.copyFrom(o_faceM, nTrace * faceOffset);
  oogs::startFinish(o_faceP, nTrace, faceOffset, ogsDfloat, ogsAdd, gshFace);

  // KEP flux differencing (advect_volInt); sets o_rhs
  volumeRhsKernel(mesh->Nelements, fieldOffset, mesh->o_vgeo, mesh->o_D, o_primv, o_rhs);

  // interface flux (advect_interface / advect_interface_riemann), interior
  // faces only, then lift into o_rhs
  const dlong Ntr = mesh->Nelements * mesh->Nfaces * mesh->Nfp;
  surfaceFluxKernel(Ntr,
                    faceOffset,
                    opt.riemannSolver,
                    static_cast<int>(opt.fullInterfaceLLF),
                    opt.blendSw0,
                    opt.blendSw1,
                    opt.blendSwRho,
                    mesh->o_sgeo,
                    o_faceM,
                    o_faceP,
                    o_faceFlux);
  surfaceLiftKernel(mesh->Nelements, fieldOffset, faceOffset, o_faceFlux, o_rhs);

  if (opt.viscous) {
    // transport properties (compvar.f vis1, cond1, diff1)
    transportKernel(mesh->Nlocal,
                    fieldOffset,
                    opt.transport,
                    opt.mu,
                    opt.prandtl,
                    opt.lewis,
                    o_aa,
                    o_mw,
                    o_viscFit,
                    o_condFit,
                    o_diffFit,
                    o_primv,
                    o_trans);

    // BR1 gradients (liftGrad), from the primitive trace exchanged above
    for (int g = 0; g < nVisc; g++) {
      viscGradKernel(mesh->Nelements,
                     fieldOffset,
                     faceOffset,
                     g,
                     mesh->o_vgeo,
                     mesh->o_sgeo,
                     mesh->o_D,
                     o_mw,
                     o_primv,
                     o_faceM,
                     o_faceP,
                     o_grad);
    }

    // viscous fluxes (viscous.f) and their central face average
    viscFluxKernel(mesh->Nlocal, fieldOffset, o_bi, o_mw, o_primv, o_trans, o_grad, o_vflux);
    viscTraceGatherKernel(mesh->Nelements, fieldOffset, faceOffset, mesh->o_vmapM, mesh->o_sgeo, o_vflux, o_vfaceM);
    o_vfaceP.copyFrom(o_vfaceM, nVisc * faceOffset);
    oogs::startFinish(o_vfaceP, nVisc, faceOffset, ogsDfloat, ogsAdd, gshVisc);

    // weak divergence (divQ_weak), added to o_rhs
    for (int q = 0; q < nVisc; q++) {
      viscDivKernel(mesh->Nelements,
                    fieldOffset,
                    faceOffset,
                    q,
                    mesh->o_vgeo,
                    mesh->o_sgeo,
                    mesh->o_D,
                    o_vflux,
                    o_vfaceM,
                    o_vfaceP,
                    o_rhs);
    }
  }

  // TODO(M4): unsplit chemistry source when reacting && !splitChemistry
}

} // namespace

dgComp::thermo_t dgComp::readThermo(const std::string &caseDir, int ns, int na)
{
  thermo_t th;
  th.npa = na;
  th.aa = readList(caseDir + "/aa.dat", ns * na);
  th.h0 = readList(caseDir + "/h0.dat", ns);
  th.s0 = readList(caseDir + "/s0.dat", ns);
  th.mw = readList(caseDir + "/mw.dat", ns);
  for (auto &w : th.mw) {
    w *= 1e-3; // g/mol -> kg/mol, as CompNek5000 usrdat does
  }
  return th;
}

void dgComp::buildKernel(occa::properties kernelInfo)
{
  const auto thermo = parString("thermo", "perfectgas");
  nekrsCheck(thermo != "perfectgas" && thermo != "nasa",
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "[DGCOMP] thermo = %s is not supported (perfectGas | nasa)\n",
             thermo.c_str());
  isPerfectGas = (thermo == "perfectgas");

  if (isPerfectGas) {
    nspec = 1;
    platform->par->extract("dgcomp", "nspec", nspec);
    nekrsCheck(nspec != 1,
               platform->comm.mpiComm(),
               EXIT_FAILURE,
               "%s\n",
               "[DGCOMP] thermo = perfectGas is single-component, nspec must be 1 (use thermo = nasa for "
               "mixtures)");
    npa = 1; // cp constant
    platform->par->extract("dgcomp", "gamma", gammaGas);
    platform->par->extract("dgcomp", "gasconstant", Rgas);
    nekrsCheck(!(gammaGas > 1) || !(Rgas > 0),
               platform->comm.mpiComm(),
               EXIT_FAILURE,
               "[DGCOMP] perfectGas needs gamma > 1 and gasConstant > 0 (got %g, %g)\n",
               gammaGas,
               Rgas);
  } else {
    nekrsCheck(!platform->par->extract("dgcomp", "nspec", nspec) || nspec < 1,
               platform->comm.mpiComm(),
               EXIT_FAILURE,
               "%s\n",
               "[DGCOMP] thermo = nasa needs nspec (and userSections = DGCOMP in the par file)");
    platform->par->extract("dgcomp", "npa", npa);
  }

  nEq = 4 + nspec;
  nPrim = prim::Y0 + nspec;
  nTrace = nEq + nPrim;
  nVisc = 4 + (nspec > 1 ? nspec : 0);

  kernelInfo["defines/p_Nspec"] = nspec;
  kernelInfo["defines/p_Nvisc"] = nVisc;
  kernelInfo["defines/p_Npa"] = npa;
  kernelInfo["defines/p_Neq"] = nEq;
  kernelInfo["defines/p_Nprim"] = nPrim;
  kernelInfo["defines/p_Ntrace"] = nTrace;
  kernelInfo["defines/p_Ru"] = Ru;
  kernelInfo["defines/p_perfectGas"] = isPerfectGas ? 1 : 0;
  kernelInfo["defines/p_gamma"] = gammaGas;
  kernelInfo["defines/p_cv"] = Rgas / (gammaGas - 1);

  kernelInfo["defines/p_PRHO"] = prim::RHO;
  kernelInfo["defines/p_PU"] = prim::U;
  kernelInfo["defines/p_PV"] = prim::V;
  kernelInfo["defines/p_PW"] = prim::W;
  kernelInfo["defines/p_PP"] = prim::P;
  kernelInfo["defines/p_PT"] = prim::T;
  kernelInfo["defines/p_PGAMMA"] = prim::GAMMA;
  kernelInfo["defines/p_PASPD"] = prim::ASPD;
  kernelInfo["defines/p_PH"] = prim::H;
  kernelInfo["defines/p_PWMEAN"] = prim::WMEAN;
  kernelInfo["defines/p_PY0"] = prim::Y0;

  auto buildKernel = [&kernelInfo](const std::string &kernelName) {
    const auto path = getenv("NEKRS_KERNEL_DIR") + std::string("/app/nrs/plugins/");
    const auto fileName = path + "dgComp.okl";
    const auto reqName = "dgComp::";
    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      platform->kernelRequests.add(reqName, fileName, kernelInfo);
      return occa::kernel();
    } else {
      buildKernelCalled = true;
      return platform->kernelRequests.load(reqName, kernelName);
    }
  };

  primitivesKernel = buildKernel("dgPrimitives");
  conservedKernel = buildKernel("dgConserved");
  traceGatherKernel = buildKernel("dgTraceGather");
  volumeRhsKernel = buildKernel("dgVolumeRhsHex3D");
  surfaceFluxKernel = buildKernel("dgSurfaceFluxHex3D");
  surfaceLiftKernel = buildKernel("dgSurfaceLiftHex3D");
  ssprkStageKernel = buildKernel("dgSSPRKStage");
  waveSpeedKernel = buildKernel("dgWaveSpeed");
  transportKernel = buildKernel("dgTransport");
  viscGradKernel = buildKernel("dgViscGradHex3D");
  viscFluxKernel = buildKernel("dgViscFlux");
  viscTraceGatherKernel = buildKernel("dgViscTraceGather");
  viscDivKernel = buildKernel("dgViscDivHex3D");
}

void dgComp::setup()
{
  if (setupCalled) {
    return;
  }
  checkBuilt();
  setupImpl(isPerfectGas ? perfectGasThermo() : readThermo(".", nspec, npa));
}

void dgComp::setup(const thermo_t &th)
{
  if (setupCalled) {
    return;
  }
  checkBuilt();
  nekrsCheck(isPerfectGas,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "%s\n",
             "dgComp::setup(thermo) needs [DGCOMP] thermo = nasa; use dgComp::setup() for a perfect gas");
  setupImpl(th);
}

bool dgComp::perfectGas()
{
  return isPerfectGas;
}

namespace
{
void setupImpl(const dgComp::thermo_t &th)
{
  _nrs = dynamic_cast<nrs_t *>(platform->app);
  mesh = _nrs->meshV;
  fieldOffset = _nrs->fieldOffset;

  nekrsCheck(th.npa != npa || th.aa.size() != static_cast<size_t>(nspec * npa) ||
                 th.h0.size() != static_cast<size_t>(nspec) || th.mw.size() != static_cast<size_t>(nspec),
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "%s\n",
             "dgComp::setup: thermo data does not match [DGCOMP] nspec/npa");

  readOptions();

  nekrsCheck(opt.reacting || opt.positivityLimiter || opt.artificialViscosity,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "%s\n",
             "dgComp: chemistry, positivity limiter and AV are not implemented yet; see doc/dg/ROADMAP.md");
  nekrsCheck(opt.volumeFlux != 0,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "%s\n",
             "dgComp: volumeFlux = ec is not implemented yet, use kep");
  // No boundary states yet: the interface kernels (inviscid and BR1) treat
  // every face as an interior (or periodic) one.
  {
    int nBoundaryFaces = 0;
    for (dlong f = 0; f < mesh->Nelements * mesh->Nfaces; f++) {
      if (mesh->EToB[f] > 0) {
        nBoundaryFaces++;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &nBoundaryFaces, 1, MPI_INT, MPI_SUM, platform->comm.mpiComm());
    nekrsCheck(nBoundaryFaces > 0,
               platform->comm.mpiComm(),
               EXIT_FAILURE,
               "dgComp: mesh has %d boundary faces, but boundary conditions are not implemented yet "
               "(fully periodic meshes only)\n",
               nBoundaryFaces);
  }

  // trace layout: e*Nfaces*Nfp + f*Nfp + n, nTrace fields
  const dlong Ntr = mesh->Nelements * mesh->Nfaces * mesh->Nfp;
  faceOffset = alignStride<dfloat>(Ntr);

  o_Ucons = platform->device.malloc<dfloat>(nEq * fieldOffset);
  o_U0 = platform->device.malloc<dfloat>(nEq * fieldOffset);
  o_rhs = platform->device.malloc<dfloat>(nEq * fieldOffset);
  o_primv = platform->device.malloc<dfloat>(nPrim * fieldOffset);
  o_faceM = platform->device.malloc<dfloat>(nTrace * faceOffset);
  o_faceP = platform->device.malloc<dfloat>(nTrace * faceOffset);
  o_faceFlux = platform->device.malloc<dfloat>(nEq * faceOffset);
  o_lambda = platform->device.malloc<dfloat>(fieldOffset);

  platform->linAlg->fill(nEq * fieldOffset, 0.0, o_Ucons);
  platform->linAlg->fill(nEq * fieldOffset, 0.0, o_rhs);
  platform->linAlg->fill(nPrim * fieldOffset, 0.0, o_primv); // T = 0 -> Newton starts at 500 K
  platform->linAlg->fill(nTrace * faceOffset, 0.0, o_faceM);

  // perfect gas: u(T) = cv T exactly, matching p_cv in dgPrimitives
  const auto bi = isPerfectGas ? std::vector<dfloat>{0, Rgas / (gammaGas - 1)} : computeBi(th);
  o_bi = platform->device.malloc<dfloat>(bi.size());
  o_bi.copyFrom(bi.data());
  o_aa = platform->device.malloc<dfloat>(th.aa.size());
  o_aa.copyFrom(th.aa.data());
  o_mw = platform->device.malloc<dfloat>(th.mw.size());
  o_mw.copyFrom(th.mw.data());

  // Face exchange: globalFaceIds pairs each face node with its neighbour's
  // (periodicity included), so a gs-add gives M + P.
  gshFace = oogs::setup(Ntr,
                        mesh->globalFaceIds,
                        nTrace,
                        faceOffset,
                        ogsDfloat,
                        platform->comm.mpiComm(),
                        1,
                        platform->device.occaDevice(),
                        NULL,
                        OOGS_AUTO);

  // BR1 viscous terms. The gradient lift reuses the primitive trace in
  // o_faceM/o_faceP; the central viscous flux needs its own exchange of
  // F_q.n, one field per row.
  {
    // fit arrays are kernel arguments in either transport mode
    std::vector<dfloat> viscFit(5 * nspec, 0.0), condFit(5 * nspec, 0.0), diffFit(5 * nspec * nspec, 0.0);
    if (opt.viscous && opt.transport == 1) {
      // exportTransport.py: Cantera MixTransport fits, ascending powers of
      // ln T, species in the same order as mw.dat (diff.dat: full matrix)
      viscFit = readList("./visc.dat", 5 * nspec);
      condFit = readList("./cond.dat", 5 * nspec);
      diffFit = readList("./diff.dat", 5 * nspec * nspec);
    }
    o_viscFit = platform->device.malloc<dfloat>(viscFit.size());
    o_viscFit.copyFrom(viscFit.data());
    o_condFit = platform->device.malloc<dfloat>(condFit.size());
    o_condFit.copyFrom(condFit.data());
    o_diffFit = platform->device.malloc<dfloat>(diffFit.size());
    o_diffFit.copyFrom(diffFit.data());
  }
  if (opt.viscous) {
    o_trans = platform->device.malloc<dfloat>((2 + nspec) * fieldOffset);
    o_grad = platform->device.malloc<dfloat>(3 * nVisc * fieldOffset);
    o_vflux = platform->device.malloc<dfloat>(3 * nVisc * fieldOffset);
    o_vfaceM = platform->device.malloc<dfloat>(nVisc * faceOffset);
    o_vfaceP = platform->device.malloc<dfloat>(nVisc * faceOffset);
    platform->linAlg->fill((2 + nspec) * fieldOffset, 0.0, o_trans);
    platform->linAlg->fill(3 * nVisc * fieldOffset, 0.0, o_grad);
    platform->linAlg->fill(3 * nVisc * fieldOffset, 0.0, o_vflux);
    platform->linAlg->fill(nVisc * faceOffset, 0.0, o_vfaceM);

    gshVisc = oogs::setup(Ntr,
                          mesh->globalFaceIds,
                          nVisc,
                          faceOffset,
                          ogsDfloat,
                          platform->comm.mpiComm(),
                          1,
                          platform->device.occaDevice(),
                          NULL,
                          OOGS_AUTO);
  }

  if (platform->comm.mpiRank() == 0) {
    if (isPerfectGas) {
      printf("dgComp: thermo=perfectGas gamma=%g R=%g J/(kg K)\n", gammaGas, Rgas);
    } else {
      printf("dgComp: thermo=nasa\n");
    }
    const char *riemannNames[] = {"llf", "hll", "hllc", "hllcBlend"};
    printf("dgComp: nspec=%d npa=%d nEq=%d riemann=%s volumeFlux=%s clip=%d\n",
           nspec,
           npa,
           nEq,
           riemannNames[opt.riemannSolver],
           opt.volumeFlux ? "ec" : "kep",
           opt.clip);
    if (opt.riemannSolver == 0) {
      printf("dgComp: fullInterfaceLLF=%d\n", opt.fullInterfaceLLF);
    } else if (opt.riemannSolver == 3) {
      printf("dgComp: blendSw0=%g blendSw1=%g blendSwRho=%g\n", opt.blendSw0, opt.blendSw1, opt.blendSwRho);
    }
    printf("dgComp: viscous=%d reacting=%d splitChemistry=%d dtChem=%g\n",
           opt.viscous,
           opt.reacting,
           opt.splitChemistry,
           opt.dtChem);
    if (opt.viscous) {
      if (opt.transport == 0) {
        printf("dgComp: BR1 viscous, transport=constant mu=%g prandtl=%g lewis=%g\n",
               opt.mu,
               opt.prandtl,
               opt.lewis);
      } else {
        printf("dgComp: BR1 viscous, transport=fits (Cantera mixture-averaged; visc.dat, cond.dat, "
               "diff.dat)\n");
      }
    }
    printf("dgComp: positivityLimiter=%d artificialViscosity=%d avCstar=%g avCmax=%g avP=%g\n",
           opt.positivityLimiter,
           opt.artificialViscosity,
           opt.avCstar,
           opt.avCmax,
           opt.avP);
  }

  setupCalled = true;
}
} // namespace

void dgComp::setSolution(const std::vector<dfloat> &U)
{
  checkReady("setSolution");
  nekrsCheck(U.size() != static_cast<size_t>(nEq * fieldOffset),
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "%s\n",
             "dgComp::setSolution: size must be nEq*fieldOffset");
  o_Ucons.copyFrom(U.data(), U.size());
  computePrimitives();
  copyToFluid();
}

std::vector<dfloat> dgComp::conservedFromPrimitives(const std::vector<dfloat> &vel,
                                                    const std::vector<dfloat> &p,
                                                    const std::vector<dfloat> &T,
                                                    const std::vector<dfloat> &Y)
{
  checkReady("conservedFromPrimitives");

  nekrsCheck(Y.empty() && nspec != 1,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "dgComp::conservedFromPrimitives: Y is required with nspec = %d\n",
             nspec);
  std::vector<dfloat> ones;
  if (Y.empty()) {
    ones.assign(fieldOffset, 1.0);
  }
  const auto &Yk = Y.empty() ? ones : Y;
  nekrsCheck(Yk.size() < static_cast<size_t>(nspec * fieldOffset),
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "%s\n",
             "dgComp::conservedFromPrimitives: Y must hold nspec*fieldOffset values");

  auto o_vel = platform->deviceMemoryPool.reserve<dfloat>(3 * fieldOffset);
  auto o_p = platform->deviceMemoryPool.reserve<dfloat>(fieldOffset);
  auto o_T = platform->deviceMemoryPool.reserve<dfloat>(fieldOffset);
  auto o_Y = platform->deviceMemoryPool.reserve<dfloat>(nspec * fieldOffset);
  o_vel.copyFrom(vel.data(), 3 * fieldOffset);
  o_p.copyFrom(p.data(), fieldOffset);
  o_T.copyFrom(T.data(), fieldOffset);
  o_Y.copyFrom(Yk.data(), nspec * fieldOffset);

  auto o_Utmp = platform->deviceMemoryPool.reserve<dfloat>(nEq * fieldOffset);
  conservedKernel(mesh->Nlocal, fieldOffset, o_mw, o_bi, o_vel, o_p, o_T, o_Y, o_Utmp);

  std::vector<dfloat> U(nEq * fieldOffset, 0.0);
  o_Utmp.copyTo(U.data(), U.size());
  return U;
}

void dgComp::computePrimitives()
{
  checkReady("computePrimitives");
  primitivesKernel(mesh->Nlocal, fieldOffset, static_cast<int>(opt.clip), o_mw, o_bi, o_aa, o_Ucons, o_primv);
}

void dgComp::step(double time, dfloat dt)
{
  checkReady("step");

  // SSPRK3 in Shu-Osher form, U <- alpha*U0 + beta*(U + dt*rhs).
  // Same coefficients as CompNek5000 rungeKutta.f; note the Fortran rkRHS is
  // -dU/dt, here o_rhs is +dU/dt.
  constexpr dfloat alpha[3] = {0.0, 0.75, 1.0 / 3.0};
  constexpr dfloat beta[3] = {1.0, 0.25, 2.0 / 3.0};
  constexpr dfloat cstage[3] = {0.0, 1.0, 0.5};

  o_U0.copyFrom(o_Ucons, nEq * fieldOffset);

  for (int s = 0; s < 3; s++) {
    computeRhs(time + cstage[s] * dt);
    ssprkStageKernel(mesh->Nlocal, fieldOffset, alpha[s], beta[s], dt, o_U0, o_rhs, o_Ucons);
    // TODO(M3): Zhang-Shu positivity limiter when opt.positivityLimiter
    computePrimitives();
  }

  // TODO(M4): Strang-split chemistry when reacting && splitChemistry

  copyToFluid();
}

dfloat dgComp::maxWaveSpeed()
{
  checkReady("maxWaveSpeed");
  waveSpeedKernel(mesh->Nlocal, fieldOffset, o_primv, o_lambda);
  return platform->linAlg->max(mesh->Nlocal, o_lambda, platform->comm.mpiComm());
}

void dgComp::addCheckpointFields(const std::vector<std::string> &speciesNames)
{
  checkReady("addCheckpointFields");
  if (checkpointFieldsAdded) {
    return;
  }

  // single component: Y = 1 everywhere, not worth writing
  const int nY = (nspec == 1) ? 0 : nspec;

  auto names = speciesNames.empty() ? speciesNamesFromPar() : speciesNames;
  if (names.empty()) {
    for (int k = 0; k < nY; k++) {
      names.push_back(scalarDigitStr(k));
    }
  }
  nekrsCheck(nY > 0 && names.size() != static_cast<size_t>(nspec),
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "dgComp: %zu species names given, nspec = %d\n",
             names.size(),
             nspec);

  // The nek (.fld) writer only accepts velocity/pressure/temperature and
  // scalarNN, so named fields need checkpointEngine = adios.
  const bool namedFields = platform->options.compareArgs("CHECKPOINT ENGINE", "ADIOS");
  nekrsCheck(!namedFields && nY + 1 > NSCALAR_MAX,
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "dgComp: density + %d mass fractions exceed the %d scalars of the .fld format, use "
             "checkpointEngine = adios\n",
             nY,
             NSCALAR_MAX);

  auto slice = [](int f) {
    return deviceMemory<dfloat>(o_primv.slice(f * fieldOffset, mesh->Nlocal));
  };

  // With [FLUID VELOCITY] solver = none nekRS still writes fluid->o_U/o_P, so
  // copy into them (copyToFluid) instead of adding a second "velocity" field.
  checkpointFieldsAdded = true;
  if (writeThroughFluid()) {
    copyToFluid();
  } else {
    _nrs->addUserCheckpointField("velocity", {slice(prim::U), slice(prim::V), slice(prim::W)});
    _nrs->addUserCheckpointField("pressure", {slice(prim::P)});
  }
  _nrs->addUserCheckpointField("temperature", {slice(prim::T)});

  // scalar fields: density, then one mass fraction per species
  std::vector<std::pair<std::string, int>> scalars = {{"density", prim::RHO}};
  for (int k = 0; k < nY; k++) {
    scalars.push_back({"Y_" + names[k], prim::Y0 + k});
  }

  std::ostringstream map;
  for (int i = 0; i < static_cast<int>(scalars.size()); i++) {
    const auto &[label, field] = scalars[i];
    const auto fieldName = namedFields ? label : "scalar" + scalarDigitStr(i);
    _nrs->addUserCheckpointField(fieldName, {slice(field)});
    // ParaView's Nek5000 reader shows scalarNN as S<NN+1>
    map << "\n  " << (namedFields ? label : "S" + scalarDigitStr(i + 1) + " = " + label);
  }

  if (platform->comm.mpiRank() == 0) {
    std::cout << "dgComp: checkpoint fields (" << (namedFields ? "adios" : "nek .fld")
              << "): velocity, pressure, temperature" << map.str() << std::endl;
  }
}

int dgComp::Nspec()
{
  return nspec;
}

int dgComp::Neq()
{
  return nEq;
}

int dgComp::Nprim()
{
  return nPrim;
}

occa::memory dgComp::o_U()
{
  return o_Ucons;
}

occa::memory dgComp::o_prim()
{
  return o_primv;
}
