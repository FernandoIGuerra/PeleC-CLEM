#ifndef CLEMSTIRRING_H
#define CLEMSTIRRING_H

#include <cmath>
#include <limits>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include "ClemIndex.H"
#include "ClemScratch.H"

namespace clem {

// Fernando-Clem: the STIRRING CLOSURE of the LEM line. It answers three
// questions and nothing else:
//
//   1. WHEN does the next triplet map fire?   computeDtStirring()
//   2. HOW BIG is its eddy?                   sampleEddySize()
//   3. WHERE on the line does it sit?         sampleEddyCenter()
//
// It does NOT apply the map (that is the rearrangement operator) and it does
// NOT sub-cycle (that is the caller's loop, which already exists in
// ClemDiffusion).
//
// THE CHAIN (1). Everything is taken from the RESOLVED LES field - nothing is
// prescribed a priori. nu_t comes from the LES model PeleC is actually running
// (computeTurbViscosity), |S| from the resolved velocity gradients:
//
//   D_T     = nu_t / Sc_t = nu_t                          Sc_t = 1
//   epsilon = nu_t |S|^2                                  equilibrium
//   eta     = ( nu^3 / epsilon )^(1/4)                    Kolmogorov scale
//   Lambda  = 27 D_T / ( 2 I(eta, Delta) )                [CM 2] inverted
//   dt_stir = 1 / ( Lambda Delta )                        one line, one Delta
//
// with the eddy-size PDF over the inertial range eta <= l <= Delta [CM 5]
//
//   f(l) = (5/3) l^(-8/3) / ( eta^(-5/3) - Delta^(-5/3) )
//
// and its third moment, in closed form,
//
//   I(eta, Delta) = Int_eta^Delta l^3 f(l) dl
//                 = (5/4) ( Delta^(4/3) - eta^(4/3) )
//                       / ( eta^(-5/3) - Delta^(-5/3) )
//
// epsilon is written through nu_t rather than as the Smagorinsky-specific
// (C_s Delta)^2 |S|^3 so that it stays correct for the dynamic model and WALE.
//
// THE SAMPLERS (2, 3) invert that same f(l): sampleEddySize() is its inverse
// CDF, sampleEddyCenter() places the segment uniformly on the (periodic) line.
//
// DIMENSIONS. CGS throughout, as everywhere in PeleC: [nu] = [nu_t] = cm^2/s,
// [Delta] = [eta] = cm, [|S|] = 1/s, [Lambda] = 1/(cm s), [dt_stir] = s.
//
// [CM n] = Chakravarthy & Menon, Phys. Fluids 13(2), 2001, Eqs. (2) and (5) -
// used in preference to Maxwell [M 3.29-3.33], who states the same chain
// non-dimensionally (lengths over the laminar flame thickness), so his
// constants cannot be lifted into a dimensional code.

// Fernando-Clem: runtime configuration ("clem.stir_*"), read once by
// ClemManager::readParams. There is only a master switch: every number in the
// chain above is either a resolved-field quantity or a fixed constant of the
// [CM 5] PDF, so there is nothing else to tune
struct StirParams
{
  bool do_stir = false; // clem.do_stir

  static void readParams();
  static const StirParams& get() { return m_params; }

private:
  static StirParams m_params;
};

// Fernando-Clem: a read-only view of PeleC's ACTIVE LES model, handed in the
// same way the transport-parameter block is (ClemDiffusion.H), so Clem keeps
// its one-way dependency and never includes PeleC.H. Filled from PeleC::do_les,
// PeleC::les_model, PeleC::Cs, PeleC::Cw and PeleC::LES_Coeffs
struct LesView
{
  bool do_les = false;
  int model = -1;       // PeleC::les_model: 0 Smag, 1 dynamic, 2 WALE, 3 Vreman
  amrex::Real Cs = 0.0; // PeleC::Cs, model 0
  amrex::Real Cw = 0.0; // PeleC::Cw, model 2
  // Fernando-Clem: PeleC::LES_Coeffs, model 1 only - the per-cell, already
  // filtered C_s^2 (NOT C_s) in component cs2_comp. nullptr otherwise
  const amrex::MultiFab* coeffs = nullptr;
  int cs2_comp = 0;
};

// Fernando-Clem: layout of the per-cell field built by computeTurbViscosity
struct NutComp
{
  static constexpr int nut = 0;    // turbulent viscosity nu_t [cm^2/s]
  static constexpr int strain = 1; // |S| [1/s]
  static constexpr int ncomp = 2;
};

// Fernando-Clem: Delta = (dx dy dz)^(1/3), the geometric mean cell size
// (dimension-generic: (dx dy)^(1/2) in 2D, dx in 1D)
amrex::Real filterWidth(const amrex::Geometry& geom);

// Fernando-Clem: fill nut_strain with (nu_t, |S|) per cell using the LES model
// PeleC is actually running - NOT a private Smagorinsky coefficient, so the
// subgrid stirring and the resolved subgrid stress cannot disagree about how
// much turbulence there is. state MUST carry >= 1 ghost cell (central
// differences of the resolved velocity); nut_strain is on the same BoxArray
// with NutComp::ncomp components and 0 ghosts. Built once per step, then read
// cell by cell by computeDtStirring
void computeTurbViscosity(
  const amrex::MultiFab& state,
  int rho_indx,
  int mom_indx,
  const amrex::Geometry& geom,
  const LesView& les,
  amrex::MultiFab& nut_strain);

// Fernando-Clem: setup-time check that the LES settings can actually feed the
// chain - i.e. that there is a non-zero eddy-viscosity COEFFICIENT and a model
// this code implements. Note it does NOT require pelec.do_les: nu_t is formed
// here, from that coefficient and the resolved velocity, and PeleC's own LES
// term aborts for NUM_SPECIES > 2 - see the note in checkLesCompatibility.
// Aborts rather than letting a run proceed with nu_t identically zero
void checkLesCompatibility(const LesView& les);

// Fernando-Clem: the MOLECULAR kinematic viscosity of a gathered line,
// nu = mu/rho mass-averaged over its elements [cm^2/s]. This is the only
// molecular quantity the chain needs (it is what sets eta), and it is the third
// input to computeDtStirring, so it belongs here and not in the diffusion
// operator: DiffusionOperator::updateProperties deliberately does NOT ask
// transport for mu, because the line equations never use it. Evaluate ONCE per
// line, before the sub-cycle - dt_stir is latched once per LES step
amrex::Real
lineViscosity(const LineScratch& lscr, const TransParmType* tparm);

// Fernando-Clem: eta = (nu^3/epsilon)^(1/4) with epsilon = nu_t |S|^2 [cm].
// Returns 0 for a cell with no subgrid turbulence
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
kolmogorovScale(
  const amrex::Real nu, const amrex::Real nu_t, const amrex::Real strain_mag)
{
  const amrex::Real eps = nu_t * strain_mag * strain_mag;
  if (nu <= 0.0 || eps <= 0.0) {
    return 0.0;
  }
  return std::pow(nu * nu * nu / eps, 0.25);
}

// Fernando-Clem: THE STIRRING TIMESCALE [s] - the whole chain in one place.
//   nu         MOLECULAR kinematic viscosity of the cell [cm^2/s] (mu/rho)
//   delta      filter width Delta [cm], from filterWidth()
//   nu_t       turbulent viscosity [cm^2/s]   } both per-cell, from
//   strain_mag |S| [1/s]                      } computeTurbViscosity
//
// Returns +max when the cell must not be stirred - no subgrid turbulence, or
// eta has reached the filter width so there is no inertial range left to draw
// eddies from (the cell is resolved). Scheduling the next event at +max is the
// correct answer for both: it never arrives
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
computeDtStirring(
  const amrex::Real nu,
  const amrex::Real delta,
  const amrex::Real nu_t,
  const amrex::Real strain_mag)
{
  constexpr amrex::Real never = std::numeric_limits<amrex::Real>::max();

  // (1) the Kolmogorov scale - the SMALLEST eddy
  const amrex::Real eta = kolmogorovScale(nu, nu_t, strain_mag);
  if (delta <= 0.0 || eta <= 0.0 || eta >= delta) {
    return never;
  }

  // (2) the event rate. I = Int l^3 f(l) dl in closed form, then Lambda from
  // D_T = (2/27) Lambda I [CM 2] with D_T = nu_t (Sc_t = 1)
  const amrex::Real integral =
    1.25 * (std::pow(delta, 4.0 / 3.0) - std::pow(eta, 4.0 / 3.0)) /
    (std::pow(eta, -5.0 / 3.0) - std::pow(delta, -5.0 / 3.0));
  if (integral <= 0.0 || nu_t <= 0.0) {
    return never;
  }
  const amrex::Real lambda = 27.0 * nu_t / (2.0 * integral);

  // (3) Lambda is a rate PER UNIT LENGTH, so Lambda*Delta is the rate on a
  // line of one filter width [CM IIB]
  return 1.0 / (lambda * delta);
}

// Fernando-Clem: draw an eddy size l [cm] from f(l) [CM 5] by inverting its
// CDF, F(l) = (eta^-5/3 - l^-5/3)/(eta^-5/3 - Delta^-5/3), with u uniform on
// [0,1). f goes as l^(-8/3), so the draws pile up just above eta
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
sampleEddySize(
  const amrex::Real eta, const amrex::Real delta, const amrex::Real u)
{
  if (eta <= 0.0 || delta <= 0.0 || eta >= delta) {
    return 0.0;
  }
  const amrex::Real a = std::pow(eta, -5.0 / 3.0);
  const amrex::Real b = std::pow(delta, -5.0 / 3.0);
  return std::pow(a - u * (a - b), -0.6); // -3/5
}

// Fernando-Clem: how many LEM elements an eddy of size l spans. Rounded to a
// MULTIPLE OF 3 - the triplet map lays down three compressed copies, so the
// segment has to split into three equal thirds - and capped at the line.
//
// Returns 0 when the eddy is smaller than the 6 elements a triplet map needs
// [CM IIB]: the event is then REJECTED, not stretched up to 6. That rejection
// is not a detail - f(l) piles up at the small end, so at high Delta/eta most
// draws land there and the realised D_T falls below the target unless n_lem is
// large enough (n_lem ~ 6 Delta/eta resolves the whole range)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
eddyElements(
  const amrex::Real l, const amrex::Real delta, const int n_lem = config::n_lem)
{
  if (l <= 0.0 || delta <= 0.0 || n_lem < 6) {
    return 0;
  }
  // Fernando-Clem: l/Delta is the fraction of the line the eddy covers; round
  // that element count to the nearest multiple of 3
  const amrex::Real n_exact = l / delta * static_cast<amrex::Real>(n_lem);
  int n_seg = 3 * static_cast<int>(std::lround(n_exact / 3.0));
  const int n_max = 3 * (n_lem / 3); // largest multiple of 3 that fits
  if (n_seg > n_max) {
    n_seg = n_max;
  }
  return (n_seg >= 6) ? n_seg : 0;
}

//Here this must be a uniform pdf that span from 0 to the number of lem -1 wiht the same probability.
//that is any particle can be the center of the vortex 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
sampleEddyCenter(const amrex::Real u, const int n_lem = config::n_lem)
{
  const int c = static_cast<int>(u * static_cast<amrex::Real>(n_lem));
  return (c < 0) ? 0 : ((c >= n_lem) ? n_lem - 1 : c);
}

} // namespace clem
#endif
