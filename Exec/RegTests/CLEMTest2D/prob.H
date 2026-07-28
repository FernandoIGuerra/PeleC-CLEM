#ifndef PROB_H
#define PROB_H

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"

// H2_ID, O2_ID, N2_ID come from mechanism.H (LiDryer: H2, O2, H2O, H, O, OH,
// HO2, H2O2, N2)

// Fill the conserved state at (i,j) from a (T, u, v, Y) primitive set at the
// ambient pressure prob_parm.pamb - shared by the IC and both inflow BCs so
// the crossflow/jet states are defined in exactly one place.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fill_state(
  const amrex::Real T,
  const amrex::Real u,
  const amrex::Real v,
  const amrex::Real massfrac[NUM_SPECIES],
  amrex::Real* s,
  ProbParmDevice const& prob_parm)
{
  amrex::Real rho, energy;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(prob_parm.pamb, massfrac, T, rho, energy);

  s[URHO] = rho;
  s[UMX] = rho * u;
  s[UMY] = rho * v;
  s[UMZ] = 0.0;
  s[UEINT] = rho * energy;
  s[UEDEN] = rho * (energy + 0.5 * (u * u + v * v));
  s[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s[UFS + n] = rho * massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
crossflow_massfrac(amrex::Real Y[NUM_SPECIES], ProbParmDevice const& prob_parm)
{
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y[n] = 0.0;
  }
  Y[O2_ID] = prob_parm.Y_O2_cross;
  Y[N2_ID] = prob_parm.Y_N2_cross;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
jet_massfrac(amrex::Real Y[NUM_SPECIES])
{
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y[n] = 0.0;
  }
  Y[H2_ID] = 1.0;
}

// t=0: the whole domain starts filled with the (quiescent-jet) crossflow -
// the jet only enters through the y-lo slot once the run starts.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  amrex::Real Y[NUM_SPECIES];
  crossflow_massfrac(Y, prob_parm);

  amrex::Real s[NVAR] = {0.0};
  fill_state(prob_parm.T_cross, prob_parm.u_cross, 0.0, Y, s, prob_parm);
  for (int n = 0; n < NVAR; n++) {
    state(i, j, k, n) = s[n];
  }
}

// Boundary conditions (only "Hard" faces reach here; bcnormal is not called
// for "FOExtrap"):
//   x-lo (idir=0, sgn=+1): crossflow inflow, uniform in y.
//   y-lo (idir=1, sgn=+1): the reactor wall, with a jet slot cut into it.
//     s_ext arrives PRE-FILLED as an odd-velocity/even-scalar reflection of
//     the interior state (see BCfill.cpp) - i.e. it is ALREADY an adiabatic
//     no-slip wall. Outside the slot we do nothing and keep that wall;
//     inside the slot [jet_xlo, jet_xhi] we overwrite it with the jet inflow.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  if (idir == 0 && sgn == 1) {
    // x-lo: crossflow inflow
    amrex::Real Y[NUM_SPECIES];
    crossflow_massfrac(Y, prob_parm);
    fill_state(prob_parm.T_cross, prob_parm.u_cross, 0.0, Y, s_ext, prob_parm);
  } else if (idir == 1 && sgn == 1) {
    // y-lo: wall, with a jet slot in [jet_xlo, jet_xhi]
    if (x[0] >= prob_parm.jet_xlo && x[0] <= prob_parm.jet_xhi) {
      amrex::Real Y[NUM_SPECIES];
      jet_massfrac(Y);
      fill_state(prob_parm.T_jet, 0.0, prob_parm.v_jet, Y, s_ext, prob_parm);
    }
    // else: leave s_ext as the pre-filled no-slip wall reflection.
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
