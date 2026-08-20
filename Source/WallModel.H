#ifndef WALLMODEL_H
#define WALLMODEL_H

#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_REAL.H>

#include "IndexDefines.H"
#include "PelePhysics.H"
#include "PeleCUtilities.H"
#include "TransCoeff.H"

// Fernando-Clem: EQUILIBRIUM WALL-STRESS MODEL (WMLES) for no-slip domain
// boundaries.
//
// WHY. On this grid the first cell centre sits at y+ ~ 30 and the cells are
// ~60 wall units in every direction, so the viscous sublayer is not resolved
// and the two-point molecular gradient mu * du/dy|wall UNDER-predicts the true
// wall stress by roughly a factor of two in a turbulent layer. Undamped
// constant Smagorinsky "fixes" that by dumping nu_t ~ 50 nu into the first
// cell, which thickens the boundary layer by an order of magnitude and (in the
// CrossJet case) opens a subsonic channel along the wall that lets a flame
// walk upstream. The right answer is to model the wall stress directly and
// leave the eddy viscosity alone.
//
// THE LAW. Spalding's single blended profile, valid from the wall through the
// log region:
//
//   y+ = u+ + exp(-kappa B) [ exp(kappa u+) - 1 - kappa u+
//                             - (kappa u+)^2/2 - (kappa u+)^3/6 ]
//
// Spalding rather than a bare log law ON PURPOSE: this boundary layer starts
// as a slug profile at the leading edge and is LAMINAR and sub-cell thick for
// the first centimetres. Spalding degenerates to u+ = y+ there, i.e. exactly
// the molecular stress, so the model is inert until y+ actually grows. A bare
// log law would impose turbulent drag on a laminar layer from x = 0 and
// manufacture a fresh version of the very artefact this is meant to remove.
//
// HOW IT IS APPLIED. Let y1 = dx[dir]/2 be the first cell-centre distance,
// u_t the resolved tangential speed there and nu_w the kinematic viscosity at
// the wall. With
//
//   R = y1 u_t / nu_w        (known),        y+ = R / u+
//
// Spalding becomes one monotone equation in u+ alone, solved by bisection
// (bulletproof: G(0+) = -inf, G(sqrt(R)) >= 0, so the root is bracketed by
// [0, sqrt(R)] always). Then u_tau = u_t / u+ and tau_w = rho u_tau^2.
//
// Rather than overwrite the face flux - which would require reproducing
// PeleC's sign and area conventions - we SCALE the tangential momentum flux
// the diffusion operator already computed. With the standard reflect_odd ghost
// the operator produced tau_2pt = mu_w u_t / y1, so the required factor is
//
//   scale = tau_w / tau_2pt = R / (u+)^2
//
// which is identically 1 in the laminar limit (u+ = sqrt(R)) and grows to ~2
// once the layer is a genuine log layer. Sign-agnostic, and it degrades to a
// no-op exactly where the flow is resolved.
//
// ENERGY. Nothing is done to UEDEN. For a stationary adiabatic wall the
// correct total-energy flux is zero, and it already is: the reflect_even
// ghost gives dT/dn = 0 (no Fourier flux) and the reflect_odd ghost gives
// u_face = 0 (no viscous work). Because the equations are advanced in
// conservative TOTAL energy, the kinetic energy that tau_w removes from the
// first cell reappears as internal energy automatically - the adiabatic-wall
// viscous heating comes out right without a separate dissipation term.
namespace pele::wallmodel {

// Spalding's y+(u+).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
spalding_yplus(
  const amrex::Real uplus, const amrex::Real kappa, const amrex::Real B)
{
  const amrex::Real ku = kappa * uplus;
  return uplus +
         std::exp(-kappa * B) *
           (std::exp(ku) - 1.0 - ku - 0.5 * ku * ku - ku * ku * ku / 6.0);
}

// Solve Spalding for u+ given R = y1 u_t / nu_w, and return the stress
// amplification factor R/(u+)^2 (1 = molecular/laminar, >1 = log layer).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
stress_amplification(
  const amrex::Real R, const amrex::Real kappa, const amrex::Real B)
{
  // R <= 0 means no tangential motion: nothing to model.
  if (!(R > 0.0)) {
    return 1.0;
  }

  // G(u+) = spalding_yplus(u+) - R/u+ is monotonically increasing, negative
  // as u+ -> 0 and >= 0 at u+ = sqrt(R) (the laminar root), so [lo, hi]
  // brackets the root for every R.
  amrex::Real lo = 0.0;
  amrex::Real hi = std::sqrt(R);

  // The exponential overflows for large u+; Spalding's u+ never exceeds ~35
  // in any practical wall layer, so clip the bracket rather than risk inf.
  constexpr amrex::Real uplus_max = 40.0;
  hi = amrex::min<amrex::Real>(hi, uplus_max);
  if (spalding_yplus(hi, kappa, B) - R / hi < 0.0) {
    // R so large that even u+ = uplus_max undershoots: saturate there.
    return R / (hi * hi);
  }

  for (int it = 0; it < 60; ++it) {
    const amrex::Real mid = 0.5 * (lo + hi);
    if (mid <= 0.0) {
      break;
    }
    const amrex::Real g = spalding_yplus(mid, kappa, B) - R / mid;
    if (g > 0.0) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  const amrex::Real uplus = 0.5 * (lo + hi);
  if (!(uplus > 0.0)) {
    return 1.0;
  }

  // Never REDUCE the stress below the molecular value: on a resolved/laminar
  // patch the model must be inert, not a drag reducer.
  return amrex::max<amrex::Real>(1.0, R / (uplus * uplus));
}

} // namespace pele::wallmodel

// Rescale the tangential momentum flux on one no-slip domain-boundary face so
// that it carries the modelled wall stress. Mirrors the calling convention of
// pc_isothermal_wall_fluxes (Diffterm.H).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_wall_model_fluxes(
  const int i,
  const int j,
  const int k,
  const int dir,
  const int normal,
  const amrex::Real kappa,
  const amrex::Real Bconst,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<amrex::EBCellFlag const>& flags,
  const amrex::Array4<amrex::Real>& flx,
  const amrex::GeometryData& geomdata,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::PhysicsType::transport_type> const* tparm,
  ProbParmDevice const& prob_parm)
{
  // iv: face-centred index where the flux lives
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  // ivm: cell-centred index of the first interior cell
  amrex::IntVect ivm{AMREX_D_DECL(i, j, k)};
  ivm[dir] += ((normal > 0) ? -1 : 0);
  // ivr: boundary direction relative to that interior cell
  amrex::IntVect ivr{0};
  ivr[dir] = normal;

  if (!flags(ivm).isConnected(ivr)) {
    return;
  }

  // The case decides where this is a wall at all - e.g. the CrossJet injector
  // hole punched through the ylo wall must be left alone.
  if (!ProblemSpecificFunctions::wall_model_active(
        ivm, dir, normal, geomdata, prob_parm)) {
    return;
  }

  // Resolved tangential velocity in the first interior cell
  const amrex::Real uvec[AMREX_SPACEDIM] = {
    AMREX_D_DECL(q(ivm, QU), q(ivm, QV), q(ivm, QW))};
  amrex::Real ut2 = 0.0;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (d != dir) {
      ut2 += uvec[d] * uvec[d];
    }
  }
  const amrex::Real ut = std::sqrt(ut2);
  if (ut <= 0.0) {
    return;
  }

  // Wall properties. The wall is adiabatic, so the wall state is the first
  // interior cell state.
  const amrex::Real rho_w = q(ivm, QRHO);
  const amrex::Real T_w = q(ivm, QTEMP);
  amrex::Real Y_w[NUM_SPECIES] = {0.0};
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    Y_w[ns] = q(ivm, ns + QFS);
  }

  const bool wtr_get_xi = false;
  const bool wtr_get_mu = true;
  const bool wtr_get_lam = false;
  const bool wtr_get_Ddiag = false;
  const bool wtr_get_chi = false;
  amrex::Real* dummy_Ddiag = nullptr;
  amrex::Real* dummy_chi_mix = nullptr;
  amrex::Real mu_w, dummy_xi, dummy_lam;
  const amrex::RealVect x = pc_cmp_loc(ivm, geomdata);
  pc_transcoeff(
    wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi, T_w,
    rho_w, Y_w, dummy_Ddiag, dummy_chi_mix, mu_w, dummy_xi, dummy_lam, tparm,
    prob_parm, x);

  if (!(mu_w > 0.0) || !(rho_w > 0.0)) {
    return;
  }

  const amrex::Real y1 = 0.5 * geomdata.CellSize(dir);
  const amrex::Real nu_w = mu_w / rho_w;
  const amrex::Real R = y1 * ut / nu_w;

  const amrex::Real scale =
    pele::wallmodel::stress_amplification(R, kappa, Bconst);

  // Scale only the TANGENTIAL momentum components; the wall-normal component
  // carries the normal stress, which no wall model touches.
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (d != dir) {
      flx(i, j, k, UMX + d) *= scale;
    }
  }
}

#endif
